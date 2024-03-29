version 1.0

# assume hg38

workflow CyriusStarAlleleCalling {
    input {
        File input_cram
        File input_cram_index
        String sample_name
        File ref_fasta
    }
    
    call RunCyrius {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            sample_name = sample_name,
            ref_fasta = ref_fasta
    }
    
    output {
        File cyrius_output = RunCyrius.cyrius_output
        File cyrius_details = RunCyrius.cyrius_details
    }
}

task RunCyrius {
    input {
        File input_cram
        File input_cram_index
        String sample_name
        File ref_fasta
    }
    
    command {
        /usr/bin/seq_cache_populate.pl -root ./ref/cache ~{ref_fasta} > /dev/null #silence ~3500 lines of contigs
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s
        
        echo ~{input_cram} > manifest.txt
        
        python /Cyrius/star_caller.py --manifest manifest.txt \
            --genome 38 --prefix ~{sample_name} --outDir ~{sample_name} \
            --threads 1 --reference ~{ref_fasta}

        mv ~{sample_name}/~{sample_name}.tsv ~{sample_name}_cyp2d6_genotype-calls.tsv
    }
    
    output {
        File cyrius_output = "~{sample_name}_cyp2d6_genotype-calls.tsv"
        File cyrius_details = "~{sample_name}/~{sample_name}.json"
    }
    
    runtime {
        memory: "3750 MiB"
        disks: "local-disk 100 HDD"
        disk: "100 GB"
        bootDiskSizeGb: 15
        preemptible: 3
        docker: "terraworkflows.azurecr.io/cyrius@sha256:5b634a57ec667046a6b71223653a12fa466e75f88df561af4988b2811cbfe15f"
    }    
}
