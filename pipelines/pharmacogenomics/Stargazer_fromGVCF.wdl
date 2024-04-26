version 1.0

#GRCh38 only

workflow StargazerStarAlleleCalling {
    input {
        File input_gvcf
        File input_gvcf_index
        File? panel_vcf_override
        File intervals_file
        String sample_name
        Array[String] gene_names
        String control_gene_name = "vdr"
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String? stargazer_docker
        Int? stargazer_mem_in_mb
    }

    call GenotypeGVCF {
        input:
            input_gvcf = input_gvcf,
            input_gvcf_index = input_gvcf_index,
            intervals_file = intervals_file,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
            }


   call RunStargazer {
        input:
            input_vcf = GenotypeGVCF.output_vcf,
            input_vcf_index = GenotypeGVCF.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample_name = sample_name,
            gene_names = gene_names,
            control_gene_name = control_gene_name,
            panel_vcf_override = panel_vcf_override,
            stargazer_docker = stargazer_docker
    }



    call RunStargazer as Stargazer_DPYD {
        input:
            input_vcf = GenotypeGVCF.output_vcf,
            input_vcf_index = GenotypeGVCF.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample_name = sample_name,
            gene_names = ["dpyd"],
            control_gene_name = control_gene_name,
             panel_vcf_override = panel_vcf_override,
            stargazer_docker = stargazer_docker,
            memory_in_mb = 7500
    }

    Array[File] star_allele_output = flatten([RunStargazer.stargazer_output, Stargazer_DPYD.stargazer_output])
    Array[File] star_allele_details = flatten([RunStargazer.stargazer_details, Stargazer_DPYD.stargazer_details])


    output {
        Array[File] stargazer_output = star_allele_output
        Array[File] stargazer_details = star_allele_details
        File genotyped_vcf = GenotypeGVCF.output_vcf
        File genotyped_vcf_index = GenotypeGVCF.output_vcf_index
    }
}

task GenotypeGVCF {
    input {
        File input_gvcf
        File input_gvcf_index
        File intervals_file
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int disk_size_gb = 100
        Int machine_mem_mb = 3750
        String gatk_docker = "terraworkflows.azurecr.io/gatk:4.5.0.0"
    }

    #TODO: try streaming with GATK 4.5
    parameter_meta {
    input_gvcf: {
      localization_optional: true
    }
    input_gvcf_index: {
      localization_optional: true
    }
    }

    command {
     /gatk/gatk GenotypeGVCFs -V ~{input_gvcf} -R ~{ref_fasta} \
        --force-output-intervals ~{intervals_file} -L ~{intervals_file} \
        --allow-old-rms-mapping-quality-annotation-data \
        -O  ~{sub(basename(input_gvcf), ".g.vcf.gz", ".vcf.gz")}
    }

    output {
        File output_vcf = sub(basename(input_gvcf), ".g.vcf.gz", ".vcf.gz")
        File output_vcf_index = sub(basename(input_gvcf_index), ".g.vcf.gz.tbi", ".vcf.gz.tbi")
    }

    runtime {
        memory: "~{machine_mem_mb} MiB"
        preemptible: 1
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
    disk: disk_size_gb + " GB"
        docker: gatk_docker
     }
}

task RunStargazer {
    input {
           File input_vcf
        File input_vcf_index
        File ref_fasta
        File ref_fasta_index
        String sample_name
        Array[String] gene_names
        String control_gene_name
        File? panel_vcf_override
        String? stargazer_docker
        Int? memory_in_mb
    }

    command <<<
        ~{if defined(panel_vcf_override) then "mv " + panel_vcf_override + " /stargazer-grc38-v.2.0.2/stargazer/1kgp_vcf/grc38/" else ""}

        # Create REF_CACHE. Used when indexing a CRAM
        /usr/bin/seq_cache_populate.pl -root ./ref/cache ~{ref_fasta} > /dev/null #don't clog the log with 3600 contig names
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s

        gene_names_array=(~{sep=" " gene_names})

        control_name = ~{control_gene_name}

        for gene_name in ${gene_names_array[@]}; do
          if [[ $gene_name -eq vdr ]]
            then
            control_name=egfr
          fi

          python /stargazer-grc38-v.2.0.2/stargazer  -t ${gene_name} -o ~{sample_name}_${gene_name} -i ~{input_vcf} -a grc38

          mv ~{sample_name}_${gene_name}/report.tsv ~{sample_name}_${gene_name}_report.tsv
          mv ~{sample_name}_${gene_name}/genotype-calls.tsv ~{sample_name}_${gene_name}_genotype-calls.tsv

        done
    >>>

    output {
        Array[File] stargazer_output = glob("~{sample_name}_*_report.tsv")
        Array[File] stargazer_details = glob("~{sample_name}_*_genotype-calls.tsv")
    }

    runtime {
        memory: select_first([memory_in_mb, 3750]) + " MiB"
        disks: "local-disk 100 HDD"
        disk: "100 GB"
        bootDiskSizeGb: 15
        preemptible: 3
        docker: select_first([stargazer_docker, "terraworkflows.azurecr.io/stargazer@sha256:2b43c8753a0fe4c2e26f007dfad13050e5e581a7348bbe562fcfe98af2a6752f"])
    }
}

