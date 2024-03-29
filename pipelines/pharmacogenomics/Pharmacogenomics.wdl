version 1.0

import "Cyrius_localize.wdl" as CyriusWorkflow
import "Stargazer_fromGVCF.wdl" as StargazerWorkflow
import "MakePGxReport.wdl" as MakePGxReport

workflow Pharmacogenomics {
    input {
        File input_cram
        File input_cram_index
        String sample_name
        File intervals_file
        File input_gvcf
        File input_gvcf_index
        File? panel_vcf_override
        String sample_name
        Array[String] gene_names
        String control_gene_name = "vdr"
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String? stargazer_docker
        String? report_docker_override
    }

    call CyriusWorkflow.CyriusStarAlleleCalling {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            sample_name = sample_name,
            ref_fasta = ref_fasta
    }

    call StargazerWorkflow.StargazerStarAlleleCalling{
        input:
            input_gvcf = input_gvcf,
            input_gvcf_index = input_gvcf_index,
            intervals_file = intervals_file,
            panel_vcf_override = panel_vcf_override,
            sample_name = sample_name,
            gene_names = gene_names,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            stargazer_docker = stargazer_docker
    }

    call StargazerWorkflow.StargazerStarAlleleCalling as Stargazer_DPYD {
        input:
            input_gvcf = input_gvcf,
            input_gvcf_index = input_gvcf_index,
            intervals_file = intervals_file,
            panel_vcf_override = panel_vcf_override,
            sample_name = sample_name,
            gene_names = ["dpyd"],
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            stargazer_docker = stargazer_docker,
            stargazer_mem_in_mb = 7500
    }

    Array[File] star_allele_calls = flatten([StargazerStarAlleleCalling.stargazer_details, Stargazer_DPYD.stargazer_details, [CyriusStarAlleleCalling.cyrius_output]])

    call MakePGxReport.MakePGxReport {
        input:
            pharmacogenomics_details = star_allele_calls,
            sample_name = sample_name,
            gene_names = gene_names,
            docker_override = report_docker_override
    }

    output {
        File pgx_report = MakePGxReport.pgx_report
    }
}