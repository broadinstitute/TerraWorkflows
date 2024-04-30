version 1.0

workflow AnnotateVCFWorkflow {
    input {
        File input_vcf
        File bed_file
        String cloud_provider
        File omim_annotations
        Int batch_size
    }

    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "terraworkflows.azurecr.io/"


    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

    #TODO which gatk docker image to use for azure?
    String gatk_docker_path = if cloud_provider == "gcp" then gatk_gcr_docker_path else gatk_acr_docker_path
    String gatk_gcr_docker_path= "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    String gatk_acr_docker_path= "dsppipelinedev.azurecr.io/gatk_reduced_layers:latest"

    String ubuntu_docker_path = if cloud_provider == "gcp" then gcr_ubuntu_docker_path else azure_ubuntu_docker_path
    String gcr_ubuntu_docker_path = "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    String azure_ubuntu_docker_path = "dsppipelinedev.azurecr.io/ubuntu_16_0_4:latest"

    # Define docker images
    String nirvana_docker_image = "nirvana:np_add_nirvana_docker"
    String variantreport_docker_image = "variantreport:latest"



    call AnnotateVCF as AnnotateVCF {
            input:
                input_filtered_vcf_tars = input_vcf,
                cloud_provider = cloud_provider,
                omim_annotations = omim_annotations,
                docker_path = docker_prefix + nirvana_docker_image
        }


    call VariantReport {
        input:
            positions_annotation_json = AnnotateVCF.positions_annotation_json,
            docker_path = docker_prefix + variantreport_docker_image
    }

    output {
        File positions_annotation_json = AnnotateVCF.positions_annotation_json
        File genes_annotation_json = AnnotateVCF.genes_annotation_json
        Array[File] variant_report_pdf = VariantReport.pdf_report
        Array[File] variant_table_tsv = VariantReport.tsv_file
    }
}



task BatchVCFs {
    input {
        File input_vcfs
        Int batch_size
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
        String docker_path
    }

    Int mem_size = ceil(size(input_vcfs, "MiB")) + 2000 + additional_memory_mb
    Int disk_size = ceil(size(input_vcfs, "GiB")) + 200 + additional_disk_gb

    command <<<
        set -euo pipefail

        # Calculate the number of batches based on the batch size
        batch_size=~{batch_size}
        num_vcfs=1
        num_batches=$(( ($num_vcfs + $batch_size - 1) / $batch_size ))

        echo "Total VCF files: $num_vcfs"
        echo "Batch size (the number of files we will process in a batch): ~{batch_size}"
        echo "Number of batches: $num_batches"

        # Process each batch. Move the VCF files to a batch directory and tar the directory

        mkdir batch_1

        mv ~{input_vcfs} batch_1

        tar_name="batch_1.tar.gz"
        tar czf $tar_name batch_1
        echo "Created tar: $tar_name"
        cp $tar_name /cromwell-executions/

    >>>
    runtime {
        docker: docker_path
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " MiB"
    }

    output {
        File batch_tar = "batch_1.tar.gz"
    }
}

task FilterVCF {
    input {
        File batch_tars
        File bed_file
        Int batch_size

        Int disk_size_gb = ceil(2*size(batch_tars, "GiB")) + 50
        Int cpu = 1
        Int memory_mb = 8000
        String docker_path
    }

    Int command_mem = memory_mb - 1000
    Int max_heap = memory_mb - 500
    command <<<
        set -euo pipefail

        # Untar the tarred inputs
        tar -xf $batch_tars --strip-components=1

        vcf_file=$(ls | grep ".vcf$")

        # Perform indexing
        echo "Indexing VCF file: $vcf_file"
        gatk \
        IndexFeatureFile \
        -I "$vcf_file"

        # Perform filtering
        echo "Filtering VCF file: $vcf_file"
        gatk \
        SelectVariants \
        -V  "$vcf_file" \
        -L ~{bed_file} \
        -O "diabetes_pathogenic_variant.vcf"


        # Create a directory for each batch and move the filtered VCF files to the corresponding directory
        batch_size=~{batch_size}
        mkdir batch1

        filtered_vcf_files=($(ls | grep ".filtered.vcf$"))
        echo "filtered_vcf_files: ${filtered_vcf_files[@]}"

        mv $filtered_vcf_files batch1/$filtered_vcf_files
        tar -czf "1.filtered_vcf_files.tar.gz" "batch1"/*.filtered.vcf

    >>>
    runtime {
        docker: docker_path
        memory: "${memory_mb} MiB"
        cpu: cpu
        disks: "local-disk ${disk_size_gb} HDD"
    }

    output {
       File tarred_filtered_vcfs = "1.filtered_vcf_files.tar.gz"
    }
}

task AnnotateVCF {
    input {
        File input_filtered_vcf_tars
        File omim_annotations
        String cloud_provider
        String docker_path
    }

    String nirvana_location = "/Nirvana/Nirvana.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"


    command <<<
        set -euo pipefail

        # Declare array of filtered vcf tars
        filtered_vcf_files=(~{input_filtered_vcf_tars})

        # Choose where to create the Nirvana DATA_SOURCES_FOLDER based on cloud_provider
        if [[ "~{cloud_provider}" == "azure" ]]; then
            DATA_SOURCES_FOLDER=/cromwell-executions/nirvana_references
        elif [[ "~{cloud_provider}" == "gcp" ]]; then
            DATA_SOURCES_FOLDER=/cromwell_root/nirvana_references
        else
            >&2 echo "Invalid cloud_provider value. Please specify either 'azure' or 'gcp'."
            exit 1
        fi

        mkdir ${DATA_SOURCES_FOLDER}

        # Download the references
        dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

        # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
        # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
        # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
        ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/

        # Bash function to perform the annotation task on a single VCF file

            # Create Nirvana annotations:
            dotnet ~{nirvana_location} \
                -i diabetes_pathogenic_variant.filtered.vcf \
                -c $DATA_SOURCES_FOLDER~{path} \
                --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
                -r $DATA_SOURCES_FOLDER~{path_reference} \
                -o diabetes_pathogenic_variant.vcf

              # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
              # Parse out the Genes section into a separate annotated json
              dotnet  ~{jasix_location} \
                  --in diabetes_pathogenic_variant.json.gz \
                  --section genes \
                  --out diabetes_pathogenic_variant.genes.json.gz

          # Parse out the Positions section into a separate annotated json
          dotnet  ~{jasix_location} \
              --in $diabetes_pathogenic_variant.json.gz \
              --section positions \
              --out diabetes_pathogenic_variant.positions.json.gz
        }

         # Define lists of vcf files
          vcf_files=($(ls | grep ".vcf$"))
          ls -lh
          echo "vcf_files: ${vcf_files[@]}"


        # Tar up the genes.json and positions.json files
        tar -czf genes_annotation_json.tar.gz *.genes.json.gz
        tar -czf positions_annotation_json.tar.gz *.positions.json.gz
    >>>

    runtime {
        docker: docker_path
        memory: "64 GB"
        cpu: "4"
        preemptible: 3
        maxRetries: 2
        disks: "local-disk 2000 HDD"
    }

    output {
        File genes_annotation_json = "genes_annotation_json.tar.gz"
        File positions_annotation_json = "positions_annotation_json.tar.gz"
    }
}

task VariantReport {

    input {
        File positions_annotation_json

        String docker_path
        Int memory_mb = 4000
        Int disk_size_gb = 15
        Int cpu = 1
    }

    command <<<

        set -euo pipefail

        # Loop through the json files and unpack them
        declare -a input_jsons=(~{positions_annotation_json})

        for json in "${input_jsons[@]}"; do
            tar -xzf $json
        done

        # Loop through the json.gz files and get the sample_id by using basename
        for json in $(ls *.json.gz); do
            sample_id=$(basename "$json" ".positions.json.gz")
            echo "Creating the report for $sample_id"
            python3 /src/variants_report.py \
            --positions_json $json \
            --sample_identifier $sample_id
        done

    >>>

    runtime {
        docker: docker_path
        memory: "${memory_mb} MiB"
        cpu: cpu
        disks: 'local-disk ${disk_size_gb} HDD'
        maxRetries: 2
    }
    output {
        Array[File] pdf_report = glob("*.pdf")
        Array[File] tsv_file = glob("*.tsv")
    }
}
