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

    call regenie_steps {
        input:
            input_bgen=filter_variants_for_gwas.ldpruned_data_bgen,
            input_samples=filter_variants_for_gwas.ldpruned_data_samples,
            pheno_csv=pheno_csv,
            pheno_col = pheno_col,
            covariate_csv = covariate_csv,
            mac_threshold = regenie_mac_threshold,
            output_prefix = output_prefix,
            step1_block_size = step1_block_size,
            step2_block_size = step2_block_size
    }

        output {
           File step2_output = regenie_steps.step2_output
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
        mkdir -p plink
        echo "Processing bgen/samples: ~{bgen}"
        #######

        # modify filter on missingness -- AC/AF filter yields ~25% missing GTs, probably due to multiallelic issues
        # one weird almost all het site is causing problem, but it's call rate is abyssmal, so just take it out
        #“–geno 0.1” tells PLINK to throw out every variant where more than 10% of the genotype calls are “NA”s.
       plink2 --bgen ~{bgen} --sample ~{samples} \
        --geno 0.9 \
        --make-bed \
        --out plink/missingness_filtered_data

       # skip remove duplicates variants and related samples
       # skip allele frequency table no one uses

        # filter on minor allele frequency
        # --maf 0.05 removes all variants with a minor allele frequency less than 0.05
        plink2 --bfile plink/missingness_filtered_data \
        --maf 0.05 \
        --make-bed \
        --out plink/maf_filtered_data

        #regenie still reports some low variance SNPs, so use AC
        plink2 --bfile plink/maf_filtered_data \
        --mac ~{mac_cutoff} \
        --make-bed \
        --out plink/max_filtered_data

        # hardy weinberg filtering
        # --hwe 1e-25 removes all variants with a Hardy-Weinberg p-value greater than 1e-25
        plink2 --bfile plink/max_filtered_data \
        --hwe 1e-25 keep-fewhet \
        --make-bed \
        --out plink/hwe_filtered_data

        plink2 --bfile plink/hwe_filtered_data \
        --set-missing-var-ids @:#\$1,\$2 \
        --make-bed --out plink/hwe_filtered_data.newIDs \
        --new-id-max-allele-len 1000

        # linkage disequilibrium
        # --ld-window-r2 0.5 sets the window size to 0.5
        plink2 --bfile plink/hwe_filtered_data.newIDs \
        --indep-pairwise 200kb 1 0.5 \
        --out plink/ldpruned_snplist \
        --rm-dup force-first

        # prune the data
        plink2 --bfile plink/hwe_filtered_data.newIDs \
        --extract plink/ldpruned_snplist.prune.in \
        --export bgen-1.2 \
        --out plink/~{output_prefix}_ldpruned_data
    >>>

    output {
        File ldpruned_data_bgen = "plink/~{output_prefix}_ldpruned_data.bgen"
        File ldpruned_data_samples = "plink/~{output_prefix}_ldpruned_data.sample"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:latest"
        memory: "31 GB"
        cpu: "4"
        disks: "local-disk 800 HDD"
    }
}

task regenie_steps {

    input {
        File input_bgen
        File input_samples

        File pheno_csv
        String pheno_col

        File? covariate_csv

        Int step1_block_size = 100
        Int step2_block_size = 200

        Int mac_threshold = 100
        String output_prefix
    }

    command <<<
        # LTL Added the next lines
        set -e

        TEMP_IN_LD=~{input_bgen}
        echo $TEMP_IN_LD

        mkdir -p plink

        # make regenie output directory
        mkdir -p regenie

        plink2 \
          --bgen ~{input_bgen} \
          --sample ~{input_samples} \
          --mac ~{mac_threshold} \
          --make-bed \
          --out high_mac_variants

       regenie \
        --step 1 \
        --bed high_mac_variants \
        --phenoFile ~{pheno_csv} \
        --phenoCol ~{pheno_col} \
        ~{"--covarFile " + covariate_csv} \
        --bsize ~{step1_block_size} \
        --lowmem \
        --out regenie/~{output_prefix}

        echo "STEP 1 COMPLETE..."
        ls -la regenie/*

        regenie \
        --step 2 \
        --bgen ~{input_bgen} \
        --sample ~{input_samples} \
        --phenoFile ~{pheno_csv} \
        --phenoCol ~{pheno_col} \
        ~{"--covarFile " + covariate_csv} \
        --bsize ~{step2_block_size} \
        --firth --approx \
        --pThresh 0.01 \
        --pred regenie/~{output_prefix}_pred.list \
        --out regenie/~{output_prefix}_firth

        ls -la regenie/*
    >>>

    output {
        File step2_output = "regenie/~{output_prefix}_firth_~{pheno_col}.regenie"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:latest"
        memory: "14 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
    }
}