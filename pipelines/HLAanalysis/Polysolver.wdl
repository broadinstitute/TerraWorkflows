version 1.0
workflow HLAAnalysis {
    
    input{
      File normalBam
      File normalBamIndex
      File tumorBam
      File tumorBamIndex
      String race
      String build
      String format
      String indiv
      Int includeFreq
      Int insertCalc
      String cloud_provider
    }
      
    # docker images
    String polysolver_docker = "polysolver:v4"

    String gcr_docker_prefix = "sachet/"
    String acr_docker_prefix = "terraworkflows.azurecr.io/"

    # choose docker prefix based on cloud provider
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

    call PolysolverType {
        input:  
            normalBam = normalBam,
            normalBamIndex = normalBamIndex,
            race = race,
            build = build,
            format = format,
            includeFreq = includeFreq,
            insertCalc = insertCalc,
            polysolver_docker_path = docker_prefix + polysolver_docker
    }
   
   call PolysolverMut {
        input:
            normalBam = normalBam,
            normalBamIndex = normalBamIndex,
            tumorBam = tumorBam,
            tumorBamIndex = tumorBamIndex,
            winners = PolysolverType.hla_winner,
            build = build,
            polysolver_docker_path = docker_prefix + polysolver_docker
    }
    
     call PolysolverAnnot {
        input:
            tarZipDir = PolysolverMut.hlaMut,
            indiv = indiv,
            polysolver_docker_path = docker_prefix + polysolver_docker
    }

    output {
        File winners = PolysolverType.hla_winner
        File hlaMut = PolysolverMut.hlaMut
        File hlaTypeBam = PolysolverMut.hlaTypeBam
        Array[File] hla_annot = PolysolverAnnot.hla_annot_out
    }
}

task PolysolverType {
      
    input {
      File normalBam
      File normalBamIndex
      String race
      String build
      String format
      Int includeFreq
      Int insertCalc
      String polysolver_docker_path

      Int disk_size = 100
      Int mem_size = 16
      Int preemptible_tries = 3
      Int cpu = 4
    }

    command <<<
        #!/bin/sh
        set -euo pipefail
        set -x

        # Step 0. Make output directory and set variables
        outDir="$(pwd)/hla_out/"
        mkdir -pv $outDir

        data_ids=/home/polysolver/data/ids
        tag_file=/home/polysolver/data/abc_v14.uniq
        samtools=/home/polysolver/binaries/samtools

        var=`cat $data_ids`

        # Step 1. Check if bam is paired end
        echo "Step 1: Check if bam is paired end"
        $samtools view ~{normalBam} | head -10000 | cut -f2   > $outDir/temp.checkpairs
        count=`wc -l $outDir/temp.checkpairs`
        echo "temp.checkpairs # lines = $count"

        /home/polysolver/scripts/check_bam_flag_pairs_at_least_one.pl $outDir/temp.checkpairs 0 1 $outDir
        status=`tail -1 $outDir/check.status.out.txt | cut -f2`
        echo "check_bam_flag_pairs status = $status"
        if [ $status == 0 ]; then
                echo "bam=~{normalBam} file is not paired"
                exit 1
        fi

        # Step 2. Calculate insert size distribution
        echo "Step 2: Calculate insert size distribution"
        if [ ~{insertCalc} == 1 ]; then
            echo "Calculating insert size distribution"
            /usr/bin/java -Xmx12g -Xms5g -jar /home/polysolver/binaries/CollectInsertSizeMetrics.jar I=~{normalBam} O=$outDir/insertsize.txt H=$outDir/insertsize.hist.pdf VALIDATION_STRINGENCY=SILENT TMP_DIR=$outDir
        else
            echo "Skipping insert size distribution"
            iFile=0
        fi

        # Step 3. Get matching tag sequences
        echo "Step 3: Get matching tag sequences"
        # Extract the header from the BAM file and save it to tag.sam
        $samtools view -H ~{normalBam} > $outDir/tag.sam
        # Filter the BAM file using the tag file and append the results to tag.sam
        $samtools view ~{normalBam} | grep -F -f $tag_file >> $outDir/tag.sam
        # Convert the SAM file to BAM format
        $samtools view -bS -o $outDir/tag.bam $outDir/tag.sam
        # Convert BAM to fastq
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/tag.bam F=$outDir/tag.1.fastq F2=$outDir/tag.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.2.fastq

        # Step 4. Get chr6 region
        echo "Step 4: Get chr6 region"
        # Append header to chr6region.sam file
        $samtools view -H ~{normalBam} > $outDir/chr6region.sam

        # Identify the regions based on the build
        declare -A regions
        if [ ~{build} == "hg38" ]; then
            echo "build=hg38"
            regions=(
                ["region1"]="6:29941260-29945884"
                ["region2"]="6:31353872-31357187"
                ["region3"]="6:31268749-31272105"
            )
        elif [ ~{build} == "hg19" ]; then
            echo "build=hg19"
            regions=(
                ["region1"]="6:29909037-29913661"
                ["region2"]="6:31321649-31324964"
                ["region3"]="6:31236526-31239869"
            )
        else
            echo "build=hg18"
            regions=(
                ["region1"]="chr6:30016016-30022640"
                ["region2"]="chr6:31428628-31433943"
                ["region3"]="chr6:31343505-31348848"
            )
        fi
        
        # Extract regions from bam files and append to chr6region.sam
        for region in ${regions[@]}; 
        do
            $samtools view ~{normalBam} $region >> $outDir/chr6region.sam
        done

        # Convert SAM to BAM
        $samtools view -bS -o $outDir/chr6region.bam $outDir/chr6region.sam
        # Convert BAM to fastq
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/chr6region.bam F=$outDir/chr6region.1.fastq F2=$outDir/chr6region.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.2.fastq

        # Step 5. Merge the two sets of fastqs
        echo "Step 5: Merge the two sets of fastqs"
        /home/polysolver/scripts/merge_fastq.pl $outDir/tag $outDir/chr6region $outDir/merged

        # Step 6. Alignment
        echo "Step 6: Alignment to HLA library"
        cat /home/polysolver/data/novoalign_complete_header.sam > $outDir/nv.complete.chr6region.R0k6.sam
        /home/polysolver/scripts/align_fork_fh.pl $outDir/merged.1.fastq $outDir/merged.2.fastq 8 ~{format} /home/polysolver/data/abc_complete.nix $outDir/nv.complete.chr6region.R0k6.sam 0 /home/polysolver/binaries
        $samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.bam $outDir/nv.complete.chr6region.R0k6.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.sam
        echo "Sort bam file"
        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.R0k6.bam O=$outDir/nv.complete.chr6region.R0k6.csorted.bam SORT_ORDER=coordinate TMP_DIR=/home/polysolver
        t=`ls -lh $outDir/nv.complete.chr6region.R0k6.csorted.bam`
        echo -n -e "Size of bam = $t\n"
        # Index bam file
        $samtools index $outDir/nv.complete.chr6region.R0k6.csorted.bam

        # Step 7. Get first allele calculations
        echo "Step 7. Compute first allele calculations."
        echo "Calculating lik1"
        /home/polysolver/scripts/first_allele_calculations_fork.pl $data_ids /home/polysolver $samtools 8 ~{race} $iFile $outDir
        echo "Get the first winners"
        # Loop through the variants and get the best allele
        for i in $var; do
                a=`tail -1 $outDir/$i.lik1 | cut -f2`
                echo -n -e "$i\t$a\n" >> $outDir/counts1.R0k6
        done

        # Assign the winners to individual variables
        winner1_a=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_a | awk 'NR==1{print $1}' | cut -f1`
        winner1_b=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_b | awk 'NR==1{print $1}' | cut -f1`
        winner1_c=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_c | awk 'NR==1{print $1}' | cut -f1`
   
        # Step 8. Get second allele calculations
        echo "Step 8. Compute second allele calculations."
        echo "Calculating lik2"
        /home/polysolver/scripts/second_allele_calculations.pl ~{race} $outDir/counts1.R0k6 $data_ids 1 /home/polysolver $outDir
        echo "Get the second winners"
        for i in $var; do
                a=`tail -1 $outDir/$i.lik2 | cut -f2`
                echo -n -e "$i\t$a\n" >> $outDir/counts2.R0k6
        done
             
        # Assign the winners to individual variables
        winner2_a=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_a | awk 'NR==1{print $1}' | cut -f1`
        winner2_b=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_b | awk 'NR==1{print $1}' | cut -f1`
        winner2_c=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_c | awk 'NR==1{print $1}' | cut -f1`

        # Echo winners1 and winnners2
        echo -n -e "winners1\t$winner1_a\t$winner1_b\t$winner1_c\n"
        echo -n -e "winners2\t$winner2_a\t$winner2_b\t$winner2_c\n"

        # Step 9. Write the winners to a file
        echo "Step 9. Write the winners to a file"
        # Define the output filename
        output_file="$outDir/winners.hla.txt"

        # Prepare the winners text
        winners_text=$(printf "HLA-A\t%s\t%s\nHLA-B\t%s\t%s\nHLA-C\t%s\t%s\n" "$winner1_a" "$winner2_a" "$winner1_b" "$winner2_b" "$winner1_c" "$winner2_c")

        # Write the winners text to the appropriate file
        echo -e "$winners_text" > $output_file

        # Optionally create the no frequency file -- this is never used so why produce?
        if [ ~{includeFreq} -ne 1 ]; then
            cp $output_file $outDir/winners.hla.nofreq.txt
        fi
    >>>
    
    output {
        File hla_winner = "./hla_out/winners.hla.txt"
        File status = "./hla_out/check.status.out.txt"

    }

    runtime {
        docker: polysolver_docker_path
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
}

task PolysolverMut {
    input {
      File normalBam
      File normalBamIndex
      File tumorBam
      File tumorBamIndex
      File winners
      String build
      String polysolver_docker_path

      Int disk_size = 100
      Int mem_size = 16
      Int preemptible_tries = 3
      Int cpu = 4
    }

    command <<<
        #!/bin/sh
        set -e
        set -x
        #/home/polysolver/scripts/shell_call_hla_mutations_from_type ${normalBam} ${tumorBam} ${winners} ${build} ${format} $(pwd)/hla_mut_out ${indiv}
       
        # Step 0. Make output directory and set variables
        outDir="$(pwd)/hla_mut_out"
        mkdir -p "$outDir"

        tag_file=/home/polysolver/data/abc_38_both_pm_update.uniq
        samtools=/home/polysolver/binaries/samtools
        
        # Step 1. Check if normal bam is paired end
        echo "Step 1: Check if bam is paired end"
        $samtools view ~{normalBam} | head -10000 | cut -f2 > temp.normal.checkpairs
        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.normal. checkpairs 0 1
        status=$(tail -1 check.status.out.txt | cut -f2)
        echo "check_bam_flag_pairs: normal status = $status"
        if [ $status == 0 ]; then
                echo "normal_bam_hla=~{normalBam} file is not paired"
                exit 1
        fi

        # Step 2. Check if tumor bam is paired end
        echo "Step 2: Check if tumor bam is paired end"
        $samtools view ~{tumorBam} | head -10000 | cut -f2 > temp.tumor.checkpairs
        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.tumor.checkpairs 0 1
        status=$(tail -1 check.status.out.txt | cut -f2)
        echo "check_bam_flag_pairs: tumor status = $status"
        if [ $status == 0 ]; then
                echo "tumor_bam_hla=~{tumorBam} file is not paired"
                exit 1
        fi

        ######################################## Align normalBam and tumorBam #####################################################

        # Function to process bam file
        process_bam() {
          local bam_file=$1
          local prefix=$2
          local out_dir=$3

          echo "Processing $bam_file with prefix $prefix"

          # Step a: Get matching tag sequences
          echo "(a) Get matching tag sequences"
          $samtools view -H "$bam_file" > "$out_dir/tag.$prefix.sam"
          $samtools view "$bam_file" | grep -F -f $tag_file >> "$out_dir/tag.$prefix.sam"
          $samtools view -bS -o "$out_dir/tag.$prefix.bam" "$out_dir/tag.$prefix.sam"
          /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I="$out_dir/tag.$prefix.bam" F="$out_dir/tag.$prefix.1.fastq" F2="$out_dir/tag.$prefix.2.fastq" VALIDATION_STRINGENCY=SILENT
          /home/polysolver/scripts/clean_unpaired_fastq.pl "$out_dir/tag.$prefix.1.fastq"
          /home/polysolver/scripts/clean_unpaired_fastq.pl "$out_dir/tag.$prefix.2.fastq"

          # Step b: Get chr6 region
          echo "(b) Get chr6 region"
          $samtools view -H "$bam_file" > "$out_dir/chr6region.$prefix.sam"
          for region in "${regions[@]}"; do
            $samtools view "$bam_file" "$region" >> "$out_dir/chr6region.$prefix.sam"
          done

          # Step c: Convert SAM to BAM
          echo "(c) Convert SAM to BAM"
          $samtools view -bS -o "$out_dir/chr6region.$prefix.bam" "$out_dir/chr6region.$prefix.sam"

          # Step d: Convert BAM to fastq
          echo "(d) Convert BAM to fastq"
          /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I="$out_dir/chr6region.$prefix.bam" F="$out_dir/chr6region.$prefix.1.fastq" F2="$out_dir/chr6region.$prefix.2.fastq" VALIDATION_STRINGENCY=SILENT
          /home/polysolver/scripts/clean_unpaired_fastq.pl "$out_dir/chr6region.$prefix.1.fastq"
          /home/polysolver/scripts/clean_unpaired_fastq.pl "$out_dir/chr6region.$prefix.2.fastq"

          # Step e: Merge fastqs
          echo "(e) Merge fastqs"
          /home/polysolver/scripts/merge_fastq.pl "$out_dir/tag.$prefix" "$out_dir/chr6region.$prefix" "$out_dir/merged.$prefix"
          rm -f "$out_dir/*sam"

          # Step f: Alignment to HLA library
          echo "(f) Alignment to HLA library"
          cat /home/polysolver/data/novoalign_complete_header.sam > "$out_dir/nv.complete.chr6region.$prefix.R0k6.sam"
          /home/polysolver/scripts/align_fork_fh.pl "$out_dir/merged.$prefix.1.fastq" "$out_dir/merged.$prefix.2.fastq" 8 STDFQ /home/polysolver/data/abc_complete.nix "$out_dir/nv.complete.chr6region.$prefix.R0k6.sam" 0 /home/polysolver/binaries
          $samtools view -bS -o "$out_dir/nv.complete.chr6region.$prefix.R0k6.bam" "$out_dir/nv.complete.chr6region.$prefix.R0k6.sam"
          rm -f "$out_dir/nv.complete.chr6region.$prefix.R0k6.sam"

          # Step g: Remove duplicates
          echo "(g) Remove duplicates"
          /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I="$out_dir/nv.complete.chr6region.$prefix.R0k6.bam" O="$out_dir/nv.complete.chr6region.$prefix.R0k6.csorted.bam" SORT_ORDER=coordinate
          $samtools rmdup "$out_dir/nv.complete.chr6region.$prefix.R0k6.csorted.bam" "$out_dir/nv.complete.chr6region.$prefix.R0k6.csorted.nodup.bam"
          local t=`ls -lh "$out_dir/nv.complete.chr6region.$prefix.R0k6.csorted.nodup.bam"`
          echo -n -e "Size of bam = $t\n"

          # Step h: Index aligned file
          echo "(h) Index aligned file"
          $samtools index "$out_dir/nv.complete.chr6region.$prefix.R0k6.csorted.nodup.bam"
       }

        # Define regions based on build
        declare -A regions
        if [ ~{build} == "hg38" ]; then
            echo "build=hg38"
            regions=(
                ["region1"]="6:29941260-29945884"
                ["region2"]="6:31353872-31357187"
                ["region3"]="6:31268749-31272105"
            )
        elif [ ~{build} == "hg19" ]; then
            echo "build=hg19"
            regions=(
                ["region1"]="6:29909037-29913661"
                ["region2"]="6:31321649-31324964"
                ["region3"]="6:31236526-31239869"
            )
        else
            echo "build=hg18"
            regions=(
                ["region1"]="chr6:30016016-30022640"
                ["region2"]="chr6:31428628-31433943"
                ["region3"]="chr6:31343505-31348848"
            )
        fi

        # Align normal and tumor BAMs
        echo "Step 3. Align normal BAM"
        process_bam "~{normalBam}" "normal" "$outDir"
        
        echo "Step 4. Align tumor BAM"
        process_bam "~{tumorBam}" "tumor" "$outDir"


        ################################################# Call HLA mutations #########################################################
        echo "Step 5. Call HLA mutations"

        normalBamChr6=$outDir/nv.complete.chr6region.normal.R0k6.csorted.nodup.bam
        tumorBamChr6=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam
        
        # Step 5.1 Call AddOrReplaceReadGroups.jar 
        echo "Step 5.1 Call AddOrReplaceReadGroups.jar"
        /usr/bin/java -jar /home/polysolver/binaries/AddOrReplaceReadGroups.jar I=$normalBamChr6 O=$outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam RGID=foo RGLB=foo RGPL=foo RGPU=foo RGSM=foo
        /usr/bin/java -jar /home/polysolver/binaries/AddOrReplaceReadGroups.jar I=$tumorBamChr6 O=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam RGID=foo RGLB=foo RGPL=foo RGPU=foo RGSM=foo
        
        echo "Changing flags and mapping quality" #? where is this done?
        
        # Step 5.2 Filter reads
        echo "Step 5.2 Filter reads from normal and tumor bam"
        /home/polysolver/scripts/filterReads.pl $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam 1
        $samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam

        /home/polysolver/scripts/filterReads.pl $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam 1
        $samtools view -bS -o $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam.temp.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam

        # Step 5.3 Index the filtered reads
        echo "Step 5.3 Index the filtered reads"
        $samtools index $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam
        $samtools index $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam
        /home/polysolver/scripts/hlaOneline.pl ~{winners} $outDir/hla.intervals

        # Step 5.4 Get type bams
        echo "Step 5.4 Get type bams"
        normalBam=$outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam
        tumorBam=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam
        
        # Get the HLA intervals - Allele format looks like this : "hla_b_39_01_01_02l"
        var=`cat $outDir/hla.intervals`

        for allele in $var; do
                echo -n -e "creating type bams for $allele\n"
                cat /home/polysolver/data/novoalign_complete_header.sam | grep "\@HD" > $outDir/temp.$allele.tumor.sam
                cat /home/polysolver/data/novoalign_complete_header.sam | grep "\@PG" >> $outDir/temp.$allele.tumor.sam
                cat /home/polysolver/data/novoalign_complete_header.sam | grep -P "$allele\t" >> $outDir/temp.$allele.tumor.sam
                echo -n -e "@RG\tID:foo\tPL:foo\tPU:foo\tLB:foo\tSM:foo\n" >> $outDir/temp.$allele.tumor.sam
                /home/samtools/samtools view $tumorBam | grep -P "\t$allele\t" | sort -k1 > $outDir/temp.$allele.tumor
                /home/polysolver/scripts/keep_only_read_pairs.pl $outDir/temp.$allele.tumor $outDir/temp.$allele.tumor2
                cat $outDir/temp.$allele.tumor2 >> $outDir/temp.$allele.tumor.sam
                /home/samtools/samtools view -bS -o $outDir/temp.$allele.tumor.bam $outDir/temp.$allele.tumor.sam
                /home/samtools/samtools sort $outDir/temp.$allele.tumor.bam $outDir/type.$allele.tumor
                /home/samtools/samtools index $outDir/type.$allele.tumor.bam
                rm -rf $outDir/temp.$allele.tumor.*

                cat /home/polysolver/data/novoalign_complete_header.sam | grep "\@HD" > $outDir/temp.$allele.sam
                cat /home/polysolver/data/novoalign_complete_header.sam | grep "\@PG" >> $outDir/temp.$allele.sam
                cat /home/polysolver/data/novoalign_complete_header.sam | grep -P "$allele\t" >> $outDir/temp.$allele.sam
                echo -n -e "@RG\tID:foo\tPL:foo\tPU:foo\tLB:foo\tSM:foo\n" >> $outDir/temp.$allele.sam
                /home/samtools/samtools view $normalBam | grep -P "\t$allele\t" | sort -k1 > $outDir/temp.$allele
                /home/polysolver/scripts/keep_only_read_pairs.pl $outDir/temp.$allele $outDir/temp.$allele.2
                cat $outDir/temp.$allele.2 >> $outDir/temp.$allele.sam
                /home/samtools/samtools view -bS -o $outDir/temp.$allele.bam $outDir/temp.$allele.sam
                /home/samtools/samtools sort $outDir/temp.$allele.bam $outDir/type.$allele
                /home/samtools/samtools index $outDir/type.$allele.bam
                rm -rf $outDir/temp.$allele.*

                # run mutect
                echo -n -e "mutect: $allele\n"
                /usr/bin/java -Xmx2g -jar /home/polysolver/binaries/muTect-1.1.6.jar \
                --analysis_type MuTect \
                --reference_sequence /home/polysolver/data/complete/$allele.fasta \
                --input_file:normal $outDir/type.$allele.bam \
                --input_file:tumor $outDir/type.$allele.tumor.bam \
                --out $outDir/call_stats.$allele.out \
                --coverage_file $outDir/coverage.$allele.txt \
                --vcf $outDir/hla.output.$allele.vcf
                #--normal_panel /cga/wu/sachet/hla/hla_caller/capture/hla_abc_hapmap133_normal_panel.vcf \

                # run strelka
                echo -n -e "strelka: $allele\n"
                /home/polysolver/binaries/bin/configureStrelkaWorkflow.pl --normal=$outDir/type.$allele.bam --tumor=$outDir/type.$allele.tumor.bam --ref=/home/polysolver/data/abc_complete.fasta --config=/home/strelka_workflow-1.0.10/etc/strelka_config_eland_default.ini --output-dir=$outDir/strelka.$allele
                make -C $outDir/strelka.$allele
                cp $outDir/strelka.$allele/results/all.somatic.indels.vcf $outDir/$allele.all.somatic.indels.vcf
                rm -rf $outDir/strelka.$allele
        done

         ########################################  Make directories and tar the output files ############################################
         echo "Step 6. Make directories and tar the output files"
         mkdir -p hla_type_bams
         cp $outDir/type* hla_type_bams
         tar -czvf hla_type_bams.tar.gz hla_type_bams

         mkdir -p hla_mut
         cp $outDir/* hla_mut/
         tar -czvf hla_mut.tar.gz hla_mut

         mkdir -p hla_type_bams
         cp $outDir/type* hla_type_bams
         tar -czvf hla_type_bams.tar.gz hla_type_bams
    >>>

    output {
        File hlaMut = "./hla_mut.tar.gz"
        File hlaTypeBam = "./hla_type_bams.tar.gz"
        Array[File] hlaMutout= glob("./hla_mut_out/*")
    }

    runtime {
        docker: polysolver_docker_path
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
}

task PolysolverAnnot {
    
    input {
      File tarZipDir
      String indiv
      String polysolver_docker_path

      Int disk_size = 100
      Int mem_size = 16
      Int preemptible_tries = 3
      Int cpu = 4
    }
    
    command <<<
      set -euo pipefail
      set -x
        
      # /home/polysolver/scripts/shell_annotate_hla_mutations ${indiv} ${tarZipDir} $outDir

      # set PERL5LIB
      export PERL5LIB=$PERL5LIB:/home/polysolver/scripts
        
      # make output directory
      outDir=$(pwd)/hla_annot_out
      mkdir -p "$outDir"

      # untar the input directory
      tar -xvfz ~{tarZipDir}

      # set dirPrefix
      dirPrefix=`/home/polysolver/scripts/get_prefix.pl ~{tarZipDir}`
      echo "dirPrefix=$dirPrefix"
        
      # set complete fasta files
      a_complete_fasta=/home/polysolver/data/a_complete.3100.new.eb.fasta
      b_complete_fasta=/home/polysolver/data/b_complete.3100.new.eb.fasta
      c_complete_fasta=/home/polysolver/data/c_complete.3100.new.eb.fasta

      # call annotate_hla_mutect
      /home/polysolver/scripts/annotate_hla_mutect.pl ~{indiv} $dirPrefix $a_complete_fasta $b_complete_fasta $c_complete_fasta /home/polysolver $outDir
        
      # call filter_hla_mutect
      /home/polysolver/scripts/filter_hla_mutect.pl "$outDir/~{indiv}.mutect.unfiltered.annotated" ~{indiv} $outDir 0 /home/polysolver
        
      # call remove_synonymous_hla_mutect
      /home/polysolver/scripts/remove_synonymous_hla_mutect.pl "$outDir/~{indiv}.mutect.filtered.annotated" ~{indiv} $outDir
        
      # call annotate_hla_strelka
      /home/polysolver/scripts/annotate_hla_strelka_indels.pl $~{indiv} $dirPrefix $a_complete_fasta $b_complete_fasta $c_complete_fasta /home/polysolver $outDir

      # call filter_hla_strelka_indels
      /home/polysolver/scripts/filter_hla_strelka_indels.pl "$outDir/~{indiv}.strelka_indels.unfiltered.annotated" ~{indiv} $outDir /home/polysolver

    >>>

    output {
        Array[File] hla_annot_out= glob("./hla_annot_out/*")
     }

    runtime {
        docker: polysolver_docker_path
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
}
