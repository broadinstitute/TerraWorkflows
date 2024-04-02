version 1.0

workflow regenie_bgen {
    input {
        File input_bgen
        File input_samples

        Int? prefilter_mac_cutoff
        Int? regenie_mac_threshold

        File pheno_csv
        String pheno_col

        File? covariate_csv

        String output_prefix

        Int? step1_block_size
        Int? step2_block_size
    }

    call filter_variants_for_gwas {
      input:
        bgen = input_bgen,
        samples = input_samples,
        mac_cutoff = prefilter_mac_cutoff,
        output_prefix = output_prefix
    }

        output {
           String step2_output = "Done"
        }
}

task filter_variants_for_gwas {
    input {
        File bgen
        File samples

        Int mac_cutoff = 100

        String output_prefix
    }

    command <<<
        set -e
        set -x
        mkdir -p plink
        ls
        echo "Processing bgen/samples: ~{bgen}"
        ls
   
    >>>

    output {
        File ldpruned_data_bgen = "plink/~{output_prefix}_ldpruned_data.bgen"
        File ldpruned_data_samples = "plink/~{output_prefix}_ldpruned_data.sample"
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "31 GB"
        cpu: "4"
        disks: "local-disk 800 HDD"
    }
}
