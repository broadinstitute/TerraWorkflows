version 1.0

import "regenie_bgen-1.2.wdl" as RegenieSingleWorkflow

#Genome-wide association study with Regenie
workflow regenie_multiple_contigs {
  input {
    File per_contig_sample_files
    File pheno_csv
    String pheno_col
      File? covariate_csv
    String prefix_in
    Int regenie_mac_threshold=200
    Int? step1_block_size=1000
    Int? step2_block_size=2000
  }

  Array[String] input_contig_files = read_lines(per_contig_sample_files)

  scatter(contig_file in input_contig_files){
    call RegenieSingleWorkflow.regenie_bgen as regenie_steps{
      input: input_bgen = sub(contig_file, ".sample", ".bgen"),
              input_samples = contig_file,
              pheno_csv = pheno_csv,
              pheno_col = pheno_col,
              covariate_csv = covariate_csv,
              output_prefix = basename(contig_file, ".sample"),
              step1_block_size = step1_block_size,
              step2_block_size = step2_block_size
              }
  }

    call gather_per_contig_files {
        input: per_contig_files = regenie_steps.step2_output,
            output_name = prefix_in
    }

  output {
      File workflow_output = gather_per_contig_files.gathered_output
  }
}

task gather_per_contig_files {
    input {
        Array[File] per_contig_files
        String output_name
    }

    command <<<
        set -euo pipefail
        head -n1 ~{per_contig_files[0]} >> ~{output_name}
        while read per_contig_file;
            do tail -n+2 $per_contig_file >> ~{output_name};
        done < ~{write_lines(per_contig_files)}
    >>>

    output {
        File gathered_output = output_name
    }

    runtime {
        docker: "debian:bookworm-slim"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 800 HDD"
    }
}

