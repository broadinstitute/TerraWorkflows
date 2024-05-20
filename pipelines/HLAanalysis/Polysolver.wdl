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

        python /home/process_monitor.py /bin/bash /home/polysolver/scripts/shell_call_hla_type ${normalBam} "$EFFECTIVE_POLYSOLVER_RACE" ${includeFreq} ${build} ${format} ${insertCalc} $(pwd)/hla_out
    >>>
    
    output {
        File winners = "hla_out/winners.hla.txt"
    }

    runtime {
        docker: container
        memory: "2G"
        cpu: 2
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


    command <<<
        set -x
        mkdir -p hla_mut_out
        /home/polysolver/scripts/shell_call_hla_mutations_from_type ${normalBam} ${tumorBam} ${winners} ${build} ${format} $(pwd)/hla_mut_out ${indiv}
    >>>

    output {
        File hlaMut = "hla_mut_out/hla_mut.tar.gz"
        File hlaTypeBam = "hla_mut_out/hla_type_bams.tar.gz"
    }

    runtime {
        docker: container
        memory: "2G"
        cpu: 2
    }
}

task PolysolverAnnot {
    
    File tarZipDir
    String indiv
    String build
    String format
    String container


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
        memory: "2G"
        cpu: 2
    }
}

