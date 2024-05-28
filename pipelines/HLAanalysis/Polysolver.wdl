workflow HLAAnalysis {
    
    File normalBam
    File normalBamIndex
    File tumorBam
    File tumorBamIndex
    String race
    String build
    String format
    String indiv
    String container


    call PolysolverType {
        input:
            normalBam = normalBam,
            normalBamIndex = normalBamIndex,
            race = race,
            build = build,
            format = format,
            container = container
    }

    call PolysolverMut {
        input:
            normalBam = normalBam,
            normalBamIndex = normalBamIndex,
            tumorBam = tumorBam,
            tumorBamIndex = tumorBamIndex,
            winners = PolysolverType.winners,
            build = build,
            format = format,
            indiv = indiv,
            container = container
    }

    call PolysolverAnnot {
        input:
            tarZipDir = PolysolverMut.hlaMut,
            indiv = indiv,
            build = build,
            format = format,
            container = container
    }

    output {
        File winners = PolysolverType.winners
        File hlaMut = PolysolverMut.hlaMut
        File hlaTypeBam = PolysolverMut.hlaTypeBam
        Array[File] mutUnfiltAnnot = PolysolverAnnot.mutUnfiltAnnot
        Array[File] mutFiltNonSynAnnot = PolysolverAnnot.mutFiltNonSynAnnot
        Array[File] mutAmbAnnot = PolysolverAnnot.mutAmbAnnot
        Array[File] mutFiltSynAnnot = PolysolverAnnot.mutFiltSynAnnot
        Array[File] strelkaIndelUnfiltAnnot = PolysolverAnnot.strelkaIndelUnfiltAnnot
        Array[File] strelkaIndelAmbAnnot = PolysolverAnnot.strelkaIndelAmbAnnot
        Array[File] strelkaIndelFiltAnnot = PolysolverAnnot.strelkaIndelFiltAnnot
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
    String container

    Int disk_size = 100
    Int mem_size = 2
    Int preemptible_tries = 3
    Int cpu = 2

    command <<<
        #!/bin/bash
        set -x

        # Use optional race to 'compute' effective race
        RACE_WC_C_COUNT=$(echo -ne "${race}" | wc -c | grep -Po '\\d+')
        if [ "$RACE_WC_C_COUNT" -eq "0" ]; then
            # If character count is 0, then use unknown
            EFFECTIVE_POLYSOLVER_RACE="Unknown"
        else
            EFFECTIVE_POLYSOLVER_RACE="${race}"
        fi
        echo "EFFECTIVE_POLYSOLVER_RACE IS $EFFECTIVE_POLYSOLVER_RACE"

        # Create output dir
        mkdir -pv hla_out

        # Call the shell_call_hla_type script with the required arguments, ignoring errors
        /home/polysolver/scripts/shell_call_hla_type_test ${normalBam} "$EFFECTIVE_POLYSOLVER_RACE" ${includeFreq} ${build} ${format} ${insertCalc} $(pwd)/hla_out 
    
    >>>
    
    output {
        File winners = "hla_out/winners.hla.txt"
    }

    runtime {
        docker: container
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
    String container
    Int disk_size = 100
    Int mem_size = 16
    Int preemptible_tries = 3
    Int cpu = 4

    command <<<
        set -x
        mkdir -p hla_mut_out
        #/home/polysolver/scripts/shell_call_hla_mutations_from_type ${normalBam} ${tumorBam} ${winners} ${build} ${format} $(pwd)/hla_mut_out ${indiv}

        
        #### check if an appropriate number of arguments were passed ####
        normal_bam_hla=${normalBam}
        tumor_bam_hla=${tumorBam}
        hla="winners.hla.nofreq.txt"
        build=${build}
        format=${format}
        outDir=$(pwd)/hla_mut_out
        indiv=${indiv}

        tag_file=/home/polysolver/data/abc_38_both_pm_update.uniq


        echo "Environment variables"
        echo "  -PSHOME: POLYSOLVER home directory = /home/polysolver"
        echo "  -SAMTOOLS_DIR: directory containing the samtools executable = /home/samtools"
        echo "  -JAVA_DIR: directory containing the JAVA_DIR executable"
        echo "  -NOVOALIGN_DIR: directory containing the Novoalign executables = /home/polysolver/binaries"
        echo "  -GATK_DIR: directory containing the GATK jar files = /home/polysolver/binaries"
        echo "  -MUTECT_DIR: directory containing the MuTect executable = /home/polysolver/binaries"
        echo "  -STRELKA_DIR: directory containing the Strelka file = /home/polysolver/binaries"
        echo "  -NUM_THREADS: number of threads = 8"



        # check if the normal bam is paired end

        /home/polysolver/binaries/samtools view $normal_bam_hla | head -10000 | cut -f2 > temp.normal.checkpairs

        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.normal. checkpairs 0 1

        status=`tail -1 check.status.out.txt | cut -f2`

        echo "check_bam_flag_pairs: normal status = $status"

        if [ $status == 0 ]; then
                echo "normal_bam_hla=$normal_bam_hla file is not paired"
                exit 1
        fi

        rm -f temp.normal.checkpairs
        rm -f check.status.out.txt


        # check if the tumor bam is paired end

        /home/polysolver/binaries/samtools view $tumor_bam_hla | head -10000 | cut -f2 > temp.tumor.checkpairs

        /home/polysolver/scripts/check_bam_flag_pairs.pl temp.tumor.checkpairs 0 1

        status=`tail -1 check.status.out.txt | cut -f2`

        echo "check_bam_flag_pairs: tumor status = $status"

        if [ $status == 0 ]; then
                echo "tumor_bam_hla=$tumor_bam_hla file is not paired"
                exit 1
        fi

        rm -f temp.tumor.checkpairs
        #rm -f check.status.out.txt


        ######################################### align normal_bam_hla #####################################################

        # getting matching tag sequences

        echo -n -e "getting matching tags\n"

        /home/samtools/samtools view -H $normal_bam_hla > $outDir/tag.sam

        /home/samtools/samtools view $normal_bam_hla | grep -F -f $tag_file >> $outDir/tag.sam

        /home/samtools/samtools view -bS -o $outDir/tag.bam $outDir/tag.sam

        /usr/bin/java -jar /home/polysolver/binaries/SamToFastq.jar I=$outDir/tag.bam F=$outDir/tag.1.fastq F2=$outDir/tag.2.fastq VALIDATION_STRINGENCY=SILENT

        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.1.fastq

        /home/polysolver/scripts/clean_unpaired_fastq.pl $outDir/tag.2.fastq

        #getting chr6 region

        echo -n -e "getting chr6 region\n"

        /home/samtools/samtools view -H $normal_bam_hla > $outDir/chr6region.sam

        if [ $build == "hg38" ]; then

                echo -n -e "build=hg38\n"

                $SAMTOOLS_DIR/samtools view $normal_bam_hla 6:29941260-29945884 >> $outDir/chr6region.sam

                $SAMTOOLS_DIR/samtools view $normal_bam_hla 6:31353872-31357187 >> $outDir/chr6region.sam

                $SAMTOOLS_DIR/samtools view $normal_bam_hla 6:31268749-31272105 >> $outDir/chr6region.sam

        elif [ $build == "hg19" ]; then

                echo -n -e "build=hg19\n"

                /home/samtools/samtools view $normal_bam_hla 6:29909037-29913661 >> $outDir/chr6region.sam

                /home/samtools/samtools view $normal_bam_hla 6:31321649-31324964 >> $outDir/chr6region.sam

                /home/samtools/samtools view $normal_bam_hla 6:31236526-31239869 >> $outDir/chr6region.sam

        else

                echo -n -e "build=hg18\n"

                /home/samtools/samtools view $normal_bam_hla chr6:30016016-30022640 >> $outDir/chr6region.sam

                /home/samtools/samtools view $normal_bam_hla chr6:31428628-31433943 >> $outDir/chr6region.sam

                /home/samtools/samtools view $normal_bam_hla chr6:31343505-31348848 >> $outDir/chr6region.sam

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

        #/home/polysolver/binaries/novoalign -d /home/polysolver/data/abc_complete.nix -f merged.1.fastq merged.2.fastq -F $format -R 0 -r all -o SAM -o FullNW | grep -P '\thla' >> nv.complete.chr6region.R0k6.sam

        ##time /broad/software/free/Linux/redhat_5_x86_64/pkgs/novocraft-2.07.18/novoalign -d /home/polysolver/data/abc_complete.nix -f merged.1.fastq merged.2.fastq -F $format -R 0 -r all -o SAM -o FullNW | grep -P '\thla' >> nv.complete.chr6region.R0k6.sam

        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.bam $outDir/nv.complete.chr6region.R0k6.sam

        rm -f $outDir/nv.complete.chr6region.R0k6.sam

        echo -n -e "removing duplicates\n"

        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.R0k6.bam O=$outDir/nv.complete.chr6region.R0k6.csorted.bam SORT_ORDER=coordinate

        #/usr/bin//usr/bin/java -Xmx12g -Xms5g -jar /home/polysolver/binaries/MarkDuplicates.jar I=nv.complete.chr6region.R0k6.csorted.bam O=nv.complete.chr6region.R0k6.csorted.nodup.bam M=dup.txt REMOVE_DUPLICATES=TRUE TMP_DIR=/broad/hptmp/firehose VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

        /home/samtools/samtools rmdup $outDir/nv.complete.chr6region.R0k6.csorted.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam

        t=`ls -lh $outDir/nv.complete.chr6region.R0k6.csorted.nodup.bam`
        echo -n -e "size of bam = $t\n"

        #rm -f nv.complete.chr6region.R0k6.csorted.bam
        #rm -f nv.complete.chr6region.R0k6.bam

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

                $SAMTOOLS_DIR/samtools view $tumor_bam_hla 6:29941260-29945884 >> $outDir/chr6region.sam

                $SAMTOOLS_DIR/samtools view $tumor_bam_hla 6:31353872-31357187 >> $outDir/chr6region.sam

                $SAMTOOLS_DIR/samtools view $tumor_bam_hla 6:31268749-31272105 >> $outDir/chr6region.sam

        elif [ $build == "hg19" ]; then

                echo -n -e "build=hg19\n"

                /home/samtools/samtools view $tumor_bam_hla 6:29909037-29913661 >> $outDir/chr6region.sam

                /home/samtools/samtools view $tumor_bam_hla 6:31321649-31324964 >> $outDir/chr6region.sam

                /home/samtools/samtools view $tumor_bam_hla 6:31236526-31239869 >> $outDir/chr6region.sam

        else

                echo -n -e "build=hg18\n"

                /home/samtools/samtools view $tumor_bam_hla chr6:30016016-30022640 >> $outDir/chr6region.sam

                /home/samtools/samtools view $tumor_bam_hla chr6:31428628-31433943 >> $outDir/chr6region.sam

                /home/samtools/samtools view $tumor_bam_hla chr6:31343505-31348848 >> $outDir/chr6region.sam

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

        #/home/polysolver/binaries/novoalign -d /home/polysolver/data/abc_complete.nix -f merged.tumor.1.fastq merged.tumor.2.fastq -F $format -R 0 -r all -o SAM -o FullNW | grep -P '\thla' >> nv.complete.chr6region.tumor.R0k6.sam

        #time /broad/software/free/Linux/redhat_5_x86_64/pkgs/novocraft-2.07.18/novoalign -d /home/polysolver/data/abc_complete.nix -f merged.tumor.1.fastq merged.tumor.2.fastq -F $format -R 0 -r all -o SAM -o FullNW | grep -P '\thla' >> nv.complete.chr6region.tumor.R0k6.sam

        /home/samtools/samtools view -bS -o $outDir/nv.complete.chr6region.tumor.R0k6.bam $outDir/nv.complete.chr6region.tumor.R0k6.sam

        rm -f $outDir/nv.complete.chr6region.tumor.R0k6.sam

        echo -n -e "removing duplicates\n"

        /usr/bin/java -jar /home/polysolver/binaries/SortSam.jar I=$outDir/nv.complete.chr6region.tumor.R0k6.bam O=$outDir/nv.complete.chr6region.tumor.R0k6.csorted.bam SORT_ORDER=coordinate

        #/usr/bin//usr/bin/java -Xmx12g -Xms5g -jar /home/polysolver/binaries/MarkDuplicates.jar I=nv.complete.chr6region.tumor.R0k6.csorted.bam O=nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam M=dup.txt REMOVE_DUPLICATES=TRUE TMP_DIR=/broad/hptmp/firehose VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

        /home/samtools/samtools rmdup $outDir/nv.complete.chr6region.tumor.R0k6.csorted.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam

        t=`ls -lh $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.bam`
        echo -n -e "size of bam = $t\n"

        #rm -f nv.complete.chr6region.tumor.R0k6.csorted.bam
        #rm -f nv.complete.chr6region.tumor.R0k6.bam

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


        #mv $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.R0k6.csorted.nodup.RG.rm.tag.bam 
        #mv $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.bam $outDir/nv.complete.chr6region.tumor.R0k6.csorted.nodup.RG.rm.tag.bam


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
        File hlaMut = "hla_mut_out/hla_mut.tar.gz"
        File hlaTypeBam = "hla_mut_out/hla_type_bams.tar.gz"
        Array[File] hlaMutout= glob("./hla_mut_out/*")
    }

    runtime {
        docker: container
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
}

task PolysolverAnnot {
    
    File tarZipDir
    String indiv
    String build
    String format
    String container

    Int disk_size = 100
    Int mem_size = 16
    Int preemptible_tries = 3
    Int cpu = 4

    command <<<
        set -x
        mkdir -p hla_annot_out
        /home/polysolver/scripts/shell_annotate_hla_mutations ${indiv} ${tarZipDir} $(pwd)/hla_annot_out
        ls -lah
    >>>

    output {
        Array[File] mutUnfiltAnnot = glob("hla_annot_out/*.mutect.unfiltered.annotated")
        Array[File] mutFiltNonSynAnnot = glob("hla_annot_out/*.mutect.filtered.nonsyn.annotated")
        Array[File] mutAmbAnnot = glob("hla_annot_out/*.mutect.ambiguous.annotated")
        Array[File] mutFiltSynAnnot = glob("hla_annot_out/*.mutect.filtered.syn.annotated")
        Array[File] strelkaIndelUnfiltAnnot = glob("hla_annot_out/*.strelka_indels.unfiltered.annotated")
        Array[File] strelkaIndelAmbAnnot = glob("hla_annot_out/*.strelka_indels.ambiguous.annotated")
        Array[File] strelkaIndelFiltAnnot = glob("hla_annot_out/*.strelka_indels.filtered.annotated")
    }

    runtime {
        docker: container
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
}

