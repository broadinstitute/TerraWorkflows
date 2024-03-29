version 1.0

workflow MakePGxReport {
    input {
        Array[File] pharmacogenomics_details
        String sample_name
        Array[String] gene_names
        String? docker_override
    }

    call MakePGxReport {
        input:
            pharmacogenomics_details = pharmacogenomics_details,
            sample_name = sample_name,
            gene_names = gene_names,
            docker_override = docker_override
    }

    output {
        File pgx_report = MakePGxReport.pgx_report
    }
}

task MakePGxReport {
    input {
        Array[File] pharmacogenomics_details
        String sample_name
        Array[String] gene_names
        String? docker_override
    }

    command <<<
        python /merge_calls.py ~{sample_name} ~{write_lines(gene_names)} ~{write_lines(pharmacogenomics_details)}
    >>>

    output {
        File pgx_report = sample_name + ".pgx_report.tsv"
    }

    runtime {
        memory: "3750 MiB"
        disks: "local-disk 100 HDD"
        disk: "100 GB"
        bootDiskSizeGb: 15
        preemptible: 3
        docker: select_first([docker_override, "terraworkflows.azurecr.io/make_pgx_report@sha256:d66d1c0fef3e3d46a008ff2033771eb5d245cdbfc70bb79e7612b6f2fe0892b3"])
    }
}