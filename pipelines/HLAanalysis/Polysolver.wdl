workflow HLAAnalysis {
    input {
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
            format = format,
            indiv = indiv,
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

    command <<<
        #!/bin/sh
        set -e 

        bam=${normalBam}
        race=${race}
        includeFreq=${includeFreq}
        build=${build}
        format=${format}
        insertCalc=${insertCalc}
        outDir="$(pwd)/hla_out/"

        ids=/home/polysolver/data/ids
        tag_file=/home/polysolver/data/abc_v14.uniq

        mkdir -pv $outDir

        var=`cat $ids`

        # check if the bam is paired end

        /home/polysolver/binaries/samtools view $bam | head -10000 | cut -f2   > $outDir/temp.checkpairs
        count=`wc -l $outDir/temp.checkpairs`
        echo "temp.checkpairs # lines = $count"

        /home/polysolver/scripts/check_bam_flag_pairs_at_least_one.pl $outDir/temp.checkpairs 0 1 $outDir
        status=`tail -1 $outDir/check.status.out.txt | cut -f2`
        echo "check_bam_flag_pairs status = $status"
        if [ $status == 0 ]; then
                echo "bam=$bam file is not paired"
                exit 1
        fi

        # calculate insert size distribution
        if [ $insertCalc == 1 ]; then
            echo -n -e "calculating insert size distribution\n"
            /usr/bin/java -Xmx12g -Xms5g -jar /home/polysolver/binaries/CollectInsertSizeMetrics.jar I=$bam O=$outDir/insertsize.txt H=$outDir/insertsize.hist.pdf VALIDATION_STRINGENCY=SILENT TMP_DIR=$outDir
        else
            echo -n -e "skipping insert size distribution\n"
            iFile=0
        fi

        # getting matching tag sequences
        echo -n -e "getting matching tags\n"
        /home/polysolver/binaries/samtools view -H $bam > $outDir/tag.sam
        /home/polysolver/binaries/samtools view $bam | grep -F -f $tag_file >> $outDir/tag.sam
        /home/polysolver/binaries/samtools view -bS -o $outDir/tag.bam $outDir/tag.sam
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/tag.bam F=$outDir/tag.1.fastq F2=$outDir/tag.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.2.fastq

        #getting chr6 region
        echo -n -e "getting chr6 region\n"
        /home/polysolver/binaries/samtools view -H $bam > $outDir/chr6region.sam
        if [ $build == "hg38" ]; then
                echo -n -e "build=hg38\n"
                /home/polysolver/binaries/samtools view $bam 6:29941260-29945884 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31353872-31357187 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31268749-31272105 >> $outDir/chr6region.sam
        elif [ $build == "hg19" ]; then
                echo -n -e "build=hg19\n"
                /home/polysolver/binaries/samtools view $bam 6:29909037-29913661 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31321649-31324964 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31236526-31239869 >> $outDir/chr6region.sam
        else
                echo -n -e "build=hg18\n"
                /home/polysolver/binaries/samtools view $bam chr6:30016016-30022640 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam chr6:31428628-31433943 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam chr6:31343505-31348848 >> $outDir/chr6region.sam
        fi

        /home/polysolver/binaries/samtools view -bS -o $outDir/chr6region.bam $outDir/chr6region.sam
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/chr6region.bam F=$outDir/chr6region.1.fastq F2=$outDir/chr6region.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.2.fastq

        # merge the two sets of fastqs
        echo -n -e "merge the two sets of fastqs\n"
        /home/polysolver/scripts/merge_fastq.pl $outDir/tag $outDir/chr6region $outDir/merged


        #rm -f $outDir/*sam

        # alignment
        echo -n -e "aligning to HLA library\n"
        cat /home/polysolver/data/novoalign_complete_header.sam > $outDir/nv.complete.chr6region.R0k6.sam
        /home/polysolver/scripts/align_fork_fh.pl $outDir/merged.1.fastq $outDir/merged.2.fastq 8 $format /home/polysolver/data/abc_complete.nix $outDir/nv.complete.chr6region.R0k6.sam 0 /home/polysolver/binaries
        /home/polysolver/binaries/samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.bam $outDir/nv.complete.chr6region.R0k6.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.sam
        echo -n -e "sorting\n"
        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.R0k6.bam O=$outDir/nv.complete.chr6region.R0k6.csorted.bam SORT_ORDER=coordinate TMP_DIR=/home/polysolver
        t=`ls -lh $outDir/nv.complete.chr6region.R0k6.csorted.bam`
        echo -n -e "size of bam = $t\n"
        #rm -f $outDir/nv.complete.chr6region.R0k6.bam
        /home/polysolver/binaries/samtools index $outDir/nv.complete.chr6region.R0k6.csorted.bam

        # first allele calculations

        echo -n -e "calculating lik1\n"
        date
        /home/polysolver/scripts/first_allele_calculations_fork.pl /home/polysolver/data/ids /home/polysolver /home/polysolver/binaries/samtools 8 $race $iFile $outDir
        date
        echo -n -e "get first winners\n"
        date
        #rm $outDir/counts1.R0k6
        for i in $var; do
                a=`tail -1 $outDir/$i.lik1 | cut -f2`
                echo -n -e "$i\t$a\n" >> $outDir/counts1.R0k6
        done

        winner1_a=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_a | head -1 | cut -f1`
        winner1_b=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_b | head -1 | cut -f1`
        winner1_c=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_c | head -1 | cut -f1`
        date


        # second allele calculations
        echo -n -e "calculating lik2\n"
        /home/polysolver/scripts/second_allele_calculations.pl $race $outDir/counts1.R0k6 /home/polysolver/data/ids 1 /home/polysolver $outDir
        date
        echo -n -e "get second winners\n"
        #rm $outDir/counts2.R0k6
        for i in $var; do
                a=`tail -1 $outDir/$i.lik2 | cut -f2`
                echo -n -e "$i\t$a\n" >> $outDir/counts2.R0k6
        done

        winner2_a=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_a | head -1 | cut -f1`
        winner2_b=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_b | head -1 | cut -f1`
        winner2_c=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_c | head -1 | cut -f1`

        echo -n -e "winners1\t$winner1_a\t$winner1_b\t$winner1_c\n"
        echo -n -e "winners2\t$winner2_a\t$winner2_b\t$winner2_c\n"

        touch $outDir/winners.hla.nofreq.txt
        touch $outDir/winners.hla.txt
        if [ $includeFreq -eq 1 ]; then
            echo -n -e "HLA-A\t$winner1_a\t$winner2_a\nHLA-B\t$winner1_b\t$winner2_b\nHLA-C\t$winner1_c\t$winner2_c\n" > $outDir/winners.hla.txt
        else
            echo -n -e "HLA-A\t$winner1_a\t$winner2_a\nHLA-B\t$winner1_b\t$winner2_b\nHLA-C\t$winner1_c\t$winner2_c\n" > $outDir/winners.hla.nofreq.txt
            cp $outDir/winners.hla.nofreq.txt $outDir/winners.hla.txt
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
    
    File normalBam
    File normalBamIndex
    File tumorBam
    File tumorBamIndex
    File winners
    String build
    String format
    String indiv
    String polysolver_docker_path

    Int disk_size = 100
    Int mem_size = 16
    Int preemptible_tries = 3
    Int cpu = 4

    command <<<
        set -x
        #/home/polysolver/scripts/shell_call_hla_mutations_from_type ${normalBam} ${tumorBam} ${winners} ${build} ${format} $(pwd)/hla_mut_out ${indiv}

        
        #### check if an appropriate number of arguments were passed ####
        normal_bam_hla=${normalBam}
        tumor_bam_hla=${tumorBam}
        hla=${winners}
        build=${build}
        format=${format}
        outDir="$(pwd)/hla_mut_out"
        indiv=${indiv}

        mkdir -p "$outDir"
        tag_file=/home/polysolver/data/abc_38_both_pm_update.uniq



        # check if the normal bam is paired end
        /home/polysolver/binaries/samtools view $normal_bam_hla | head -10000 | cut -f2 > temp.normal.checkpairs
        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.normal. checkpairs 0 1

        status=$(tail -1 check.status.out.txt | cut -f2)
        echo "check_bam_flag_pairs: normal status = $status"
        if [ $status == 0 ]; then
                echo "normal_bam_hla=$normal_bam_hla file is not paired"
                exit 1
        fi


        # check if the tumor bam is paired end
        /home/polysolver/binaries/samtools view $tumor_bam_hla | head -10000 | cut -f2 > temp.tumor.checkpairs
        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.tumor.checkpairs 0 1
        status=$(tail -1 check.status.out.txt | cut -f2)
        echo "check_bam_flag_pairs: tumor status = $status"
        if [ $status == 0 ]; then
                echo "tumor_bam_hla=$tumor_bam_hla file is not paired"
                exit 1
        fi

        ######################################### align normal_bam_hla #####################################################

        # getting matching tag sequences
        echo -n -e "getting matching tags\n"
        /home/samtools/samtools view -H $normal_bam_hla > "$outDir/tag.sam"
        /home/samtools/samtools view $normal_bam_hla | grep -F -f $tag_file >> "$outDir/tag.sam"
        /home/samtools/samtools view -bS -o "$outDir/tag.bam" "$outDir/tag.sam"
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I="$outDir/tag.bam" F="$outDir/tag.1.fastq" F2="$outDir/tag.2.fastq" VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl "$outDir/tag.1.fastq"
        /home/polysolver/scripts/clean_unpaired_fastq.pl "$outDir/tag.2.fastq"

        #getting chr6 region
        echo -n -e "getting chr6 region\n"
        /home/samtools/samtools view -H $normal_bam_hla > $outDir/chr6region.sam

        if [ $build == "hg38" ]; then
                echo -n -e "build=hg38\n"
                /home/polysolver/binaries/samtools view $bam 6:29941260-29945884 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam 6:31353872-31357187 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam 6:31268749-31272105 >> "$outDir/chr6region.sam"
        elif [ $build == "hg19" ]; then
                echo -n -e "build=hg19\n"
                /home/polysolver/binaries/samtools view $bam 6:29909037-29913661 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam 6:31321649-31324964 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam 6:31236526-31239869 >> "$outDir/chr6region.sam"
        else
                echo -n -e "build=hg18\n"
                /home/polysolver/binaries/samtools view $bam chr6:30016016-30022640 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam chr6:31428628-31433943 >> "$outDir/chr6region.sam"
                /home/polysolver/binaries/samtools view $bam chr6:31343505-31348848 >> "$outDir/chr6region.sam"
        fi

        /home/samtools/samtools view -bS -o $outDir/chr6region.bam $outDir/chr6region.sam
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/chr6region.bam F=$outDir/chr6region.1.fastq F2=$outDir/chr6region.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.2.fastq

        # merge the two sets of fastqs
        echo -n -e "merge the two sets of fastqs\n"
        /home/polysolver/scripts/merge_fastq.pl $outDir/tag $outDir/chr6region $outDir/merged


        #alignment

        echo -n -e "aligning to HLA library\n"
        cat /home/polysolver/data/novoalign_complete_header.sam > $outDir/nv.complete.chr6region.R0k6.sam
        /home/polysolver/scripts/align_fork_fh.pl $outDir/merged.1.fastq $outDir/merged.2.fastq 8 STDFQ /home/polysolver/data/abc_complete.nix $outDir/nv.complete.chr6region.R0k6.sam 0 /home/polysolver/binaries
        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.bam $outDir/nv.complete.chr6region.R0k6.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.sam


        echo -n -e "removing duplicates\n"
        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.R0k6.bam O=$outDir/nv.complete.chr6region.R0k6.csorted.bam SORT_ORDER=coordinate
        /home/samtools/samtools rmdup $outDir/nv.complete.chr6region.R0k6.csorted.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam
        t=`ls -lh $outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam`
        echo -n -e "size of bam = $t\n"


        echo -n -e "indexing aligned file\n"
        date
        /home/samtools/samtools index $outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam


        ########################################## align tumor_bam_hla #####################################################
        # getting matching tag sequences
        echo -n -e "getting matching tags\n"
        /home/samtools/samtools view -H $tumor_bam_hla > $outDir/tag.tumor.sam
        /home/samtools/samtools view $tumor_bam_hla | grep -F -f $tag_file >> $outDir/tag.tumor.sam
        /home/samtools/samtools view -bS -o $outDir/tag.tumor.bam $outDir/tag.tumor.sam
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/tag.tumor.bam F=$outDir/tag.tumor.1.fastq F2=$outDir/tag.tumor.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.tumor.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.tumor.2.fastq

        #getting chr6 region
        echo -n -e "getting chr6 region\n"
        /home/samtools/samtools view -H $tumor_bam_hla > $outDir/chr6region.sam

        if [ $build == "hg38" ]; then
                echo -n -e "build=hg38\n"
                /home/polysolver/binaries/samtools view $bam 6:29941260-29945884 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31353872-31357187 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31268749-31272105 >> $outDir/chr6region.sam
        elif [ $build == "hg19" ]; then
                echo -n -e "build=hg19\n"
                /home/polysolver/binaries/samtools view $bam 6:29909037-29913661 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31321649-31324964 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam 6:31236526-31239869 >> $outDir/chr6region.sam
        else
                echo -n -e "build=hg18\n"
                /home/polysolver/binaries/samtools view $bam chr6:30016016-30022640 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam chr6:31428628-31433943 >> $outDir/chr6region.sam
                /home/polysolver/binaries/samtools view $bam chr6:31343505-31348848 >> $outDir/chr6region.sam
        fi

        /home/samtools/samtools view -bS -o $outDir/chr6region.bam $outDir/chr6region.sam
        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/chr6region.bam F=$outDir/chr6region.tumor.1.fastq F2=$outDir/chr6region.tumor.2.fastq VALIDATION_STRINGENCY=SILENT
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.tumor.1.fastq
        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/chr6region.tumor.2.fastq

        # merge the two sets of fastqs
        echo -n -e "merge the two sets of fastqs\n"
        /home/polysolver/scripts/merge_fastq.pl $outDir/tag.tumor $outDir/chr6region.tumor $outDir/merged.tumor
        rm -f $outDir/*sam

        #alignment
        echo -n -e "aligning to HLA library\n"
        cat /home/polysolver/data/novoalign_complete_header.sam > $outDir/nv.complete.chr6region.tumor.R0k6.sam
        /home/polysolver/scripts/align_fork_fh.pl $outDir/merged.tumor.1.fastq $outDir/merged.tumor.2.fastq 8 STDFQ /home/polysolver/data/abc_complete.nix $outDir/nv.complete.chr6region.tumor.R0k6.sam 0 /home/polysolver/binaries
        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.tumor.R0k6.bam $outDir/nv.complete.chr6region.tumor.R0k6.sam
        rm -f $outDir/nv.complete.chr6region.tumor.R0k6.sam

        echo -n -e "removing duplicates\n"
        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.tumor.R0k6.bam O=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.bam SORT_ORDER=coordinate
        /home/samtools/samtools rmdup $outDir/nv.complete.chr6region.tumor.R0k6.csorted.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam
        t=`ls -lh $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam`
        echo -n -e "size of bam = $t\n"


        echo -n -e "indexing aligned file\n"
        date
        /home/samtools/samtools index $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam

        ################################## call hla mutations ########################################

        normalBam=$outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam
        tumorBam=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam
        hlaFile=$hla
        /usr/bin/java -jar /home/polysolver/binaries/AddOrReplaceReadGroups.jar I=$normalBam O=$outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam RGID=foo RGLB=foo RGPL=foo RGPU=foo RGSM=foo
        /usr/bin/java -jar /home/polysolver/binaries/AddOrReplaceReadGroups.jar I=$tumorBam O=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam RGID=foo RGLB=foo RGPL=foo RGPU=foo RGSM=foo
        echo -n -e "changing flags and mapping quality\n"

        /home/polysolver/scripts/filterReads.pl $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam 1
        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam

        /home/polysolver/scripts/filterReads.pl $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam 1
        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam.temp.sam
        rm -f $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam.temp.sam

        echo -n -e "indexing\n"
        /home/samtools/samtools index $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam
        /home/samtools/samtools index $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam
        /home/polysolver/scripts/hlaOneline.pl $hlaFile $outDir/hla.intervals

        echo -n -e "get type bams\n"
        normalBam=$outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam
        tumorBam=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam

        var=`cat $outDir/hla.intervals`


        #allele="hla_b_39_01_01_02l"

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
    
    File tarZipDir
    String indiv
    String polysolver_docker_path

    Int disk_size = 100
    Int mem_size = 16
    Int preemptible_tries = 3
    Int cpu = 4

    command <<<
        set -x
        outDir=$(pwd)/hla_annot_out

        mkdir -p "$outDir"
        #/home/polysolver/scripts/shell_annotate_hla_mutations ${indiv} ${tarZipDir} $outDir
        ls -lah

        #### check if an appropriate number of arguments were passed ####
        indiv=${indiv}
        tarZipDir=${tarZipDir}

        export PERL5LIB=$PERL5LIB:/home/polysolver/scripts
        tar xvfz $tarZipDir 
        dirPrefix=`/home/polysolver/scripts/get_prefix.pl $tarZipDir`
        echo "dirPrefix=$dirPrefix"
        /home/polysolver/scripts/annotate_hla_mutect.pl $indiv $dirPrefix /home/polysolver/data/a_complete.3100.new.eb.fasta /home/polysolver/data/b_complete.3100.new.eb.fasta /home/polysolver/data/c_complete.3100.new.eb.fasta /home/polysolver $outDir
        /home/polysolver/scripts/filter_hla_mutect.pl "$outDir/$indiv.mutect.unfiltered.annotated" $indiv $outDir 0 /home/polysolver
        /home/polysolver/scripts/remove_synonymous_hla_mutect.pl "$outDir/$indiv.mutect.filtered.annotated" $indiv $outDir
        /home/polysolver/scripts/annotate_hla_strelka_indels.pl $indiv $dirPrefix /home/polysolver/data/a_complete.3100.new.eb.fasta /home/polysolver/data/b_complete.3100.new.eb.fasta /home/polysolver/data/c_complete.3100.new.eb.fasta /home/polysolver $outDir
        /home/polysolver/scripts/filter_hla_strelka_indels.pl "$outDir/$indiv.strelka_indels.unfiltered.annotated" $indiv $outDir /home/polysolver

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