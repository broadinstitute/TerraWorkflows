version 1.0

workflow PrintReads {
    input {
        String input_bam
        String input_bam_index
        String intervals
        String output_name
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        Int preemptible_tries = 1
    }

    call RunPrintReads {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            intervals = intervals,
            output_name = output_name,
            ref_fasta = ref_fasta,
            ref_fasta_index =ref_fasta_index,
            ref_dict = ref_dict
    }

    output {
        File output_bam = RunPrintReads.output_bam
        File output_bam_index = RunPrintReads.output_bam_index
    }

}


task RunPrintReads {
    input {
        String input_bam
        String input_bam_index
        String intervals
        String output_name
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        Int preemptible_tries = 1
    }

    Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Int disk_size = ceil(((size(input_bam, "GiB") + 30)) + ref_size) + 20

    # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
    # the interval because the assembly uses them.
    command <<<
    set -euo pipefail

        /gatk/gatk --java-options "-Xms2000m -Xmx9000m"\
        PrintReads \
        --input ~{input_bam}?$AZURE_STORAGE_SAS_TOKEN \
        --interval-padding 500 \
        -L ~{intervals} \
        -O ~{output_name}.bam \
        -R ~{ref_fasta}
    >>>
    runtime {
        docker: "terraworkflows.azurecr.io/gatk:4.5.0.0"
        preemptible: preemptible_tries
        memory: "10000 MiB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        azureSasEnvironmentVariable: "AZURE_STORAGE_SAS_TOKEN"
    }
    output {
        File output_bam = "~{output_name}.bam"
        File output_bam_index = "~{output_name}.bai"
    }
}