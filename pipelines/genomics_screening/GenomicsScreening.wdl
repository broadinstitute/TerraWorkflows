version 1.0

workflow GenomicsScreening {
    input {
        Array[File] input_vcfs
        File bed_file
        String cloud_provider
        Int batch_size
        String output_prefix
    }

    String pipeline_version = "1.0"

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
    String variantreport_docker_image = "genomics_variant_report:0e50556"


    call BatchVCFs as batch_vcfs {
        input:
            input_vcfs = input_vcfs,
            batch_size = batch_size,
            docker_path = ubuntu_docker_path
    }

    scatter(tar in batch_vcfs.batch_tars) {
        call FilterVCF as filter_vcf {
                input:
                    batch_tars = tar,
                    bed_file = bed_file,
                    batch_size = batch_size,
                    docker_path = gatk_docker_path
            }
    }


    scatter(tar in filter_vcf.tarred_filtered_vcfs) {
        call AnnotateVCF as AnnotateVCF {
                input:
                    input_filtered_vcf_tars = tar,
                    cloud_provider = cloud_provider,
                    docker_path = docker_prefix + nirvana_docker_image,
                    output_prefix = output_prefix
            }
        }

    call VariantReport {
        input:
            positions_annotation_json = AnnotateVCF.positions_annotation_json,
            docker_path = docker_prefix + variantreport_docker_image,
            output_prefix = output_prefix,
            bed_file = bed_file
    }

    call PipelineMetadata {
        input:
            pipeline_version = pipeline_version,
            output_prefix = output_prefix,
            input_vcfs = input_vcfs,
            pdf_report = VariantReport.pdf_report,
            bed_file = bed_file
    }

    output {
        Array[File] positions_annotation_json = AnnotateVCF.positions_annotation_json
        Array[File] genes_annotation_json = AnnotateVCF.genes_annotation_json
        Array[File] variant_report_pdf = VariantReport.pdf_report
        Array[File] variant_table_tsv = VariantReport.tsv_file
        File pipeline_metadata = PipelineMetadata.pipeline_metadata
    }
}



task BatchVCFs {
    input {
        Array[File] input_vcfs
        Int batch_size
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
        String docker_path
    }

    Int mem_size = ceil(size(input_vcfs, "MiB")) + 2000 + additional_memory_mb
    Int disk_size = ceil(size(input_vcfs, "GiB")) + 200 + additional_disk_gb

    command <<<
        set -euo pipefail

        # Get the size of the input_vcfs array
        declare -a input_vcfs=(~{sep=' ' input_vcfs})
        num_vcfs=${#input_vcfs[@]}

        # Calculate the number of batches based on the batch size
        batch_size=~{batch_size}
        num_batches=$(( ($num_vcfs + $batch_size - 1) / $batch_size ))

        echo "Total VCF files: $num_vcfs"
        echo "Batch size (the number of files we will process in a batch): ~{batch_size}"
        echo "Number of batches: $num_batches"

        # Process each batch. Move the VCF files to a batch directory and tar the directory
        for ((i=0; i<num_batches; i++)); do
            start_idx=$((i * batch_size))
            end_idx=$((start_idx + batch_size - 1))

            # Create a directory for the current batch
            batch_dir="batch_${i}"
            mkdir $batch_dir

            echo "Processing batch $i..."
            for ((j=start_idx; j<=end_idx && j<num_vcfs; j++)); do
                vcf="${input_vcfs[j]}"
                mv $vcf $batch_dir
            done

            # Tar the batch directory
            tar_name="batch_${i}.tar.gz"
            tar czf $tar_name $batch_dir
            echo "Created tar: $tar_name"
        done
    >>>
    runtime {
        docker: docker_path
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " MiB"
    }

    output {
        Array[File] batch_tars = glob("*.tar.gz")
    }
}

task FilterVCF {
    input {
        File batch_tars
        File bed_file
        Int batch_size

        Int disk_size_gb = ceil(2*size(batch_tars, "GiB")) + 50
        Int cpu = 1
        Int memory_mb = 64000
        String docker_path
    }

    Int command_mem = memory_mb - 1000
    Int max_heap = memory_mb - 500
    command <<<
        set -euo pipefail

        # Bash function to perform the filtering/indexing task on a single VCF file

        task() {
            local vcf_file=$1

            echo "Starting task for $vcf_file.."

            # Determine the file extension and extract the sample ID
            if [[ "$vcf_file" == *.vcf.gz ]]; then
                sample_id=$(basename "$vcf_file" ".vcf.gz")
            elif [[ "$vcf_file" == *.vcf ]]; then
                sample_id=$(basename "$vcf_file" ".vcf")
            else
                echo "Unsupported file extension"
                exit 1
            fi

            # Perform sorting
            echo "Sorting VCF file: $vcf_file"
            gatk \
            SortVcf \
            -I "$vcf_file" \
            -O "$sample_id.sorted.vcf"
            echo "Sorting done."

            # Perform indexing
            echo "Indexing VCF file: $sample_id.sorted.vcf"
            gatk \
            IndexFeatureFile \
            -I "$sample_id.sorted.vcf"
            echo "Indexing done."

            # Perform filtering
            echo "Filtering VCF file: $sample_id.sorted.vcf"
            gatk \
            SelectVariants \
            -V  "$sample_id.sorted.vcf" \
            -L ~{bed_file} \
            -O "$sample_id.filtered.vcf"
            echo "Filtering done."
        }

        # Declare array of input batched tars
        declare -a batch_tars=(~{sep=' ' batch_tars})
        echo "batch_tars: ${batch_tars[@]}"

        # Untar the tarred inputs
        tar -xf $batch_tars --strip-components=1

        # Declare array of vcf or .vcf.gz files
        declare -a input_vcfs=($(ls | grep -E "\.vcf$|\.vcf\.gz$"))
        echo "input_vcfs: ${input_vcfs[@]}"

        # Launch tasks in parallel for each VCF file
        for ((i = 0; i < ${#input_vcfs[@]}; i += 2)); do

            task "${input_vcfs[i]}" &

            # Check if there's another file available in the batch
            if [[ $((i + 1)) -lt ${#input_vcfs[@]} ]]; then
                echo "Next file to process: ${input_vcfs[i + 1]}"
                task "${input_vcfs[i + 1]}" &
            fi

            # Limit the number of concurrent tasks to 2 (adjust as needed)
            if [[ $(jobs -p | wc -l) -ge 2 ]]; then
                wait -n # Wait for any background job to finish
            fi
        done

        wait
        echo "Tasks all done."

        # Create a directory for each batch and move the filtered VCF files to the corresponding directory
        batch_size=~{batch_size}
        for i in $(seq 1 "${batch_size}"); do
            mkdir -p "batch${i}"
        done

        folder_index=1

        filtered_vcf_files=($(ls | grep ".filtered.vcf$"))
        echo "filtered_vcf_files: ${filtered_vcf_files[@]}"

        for file in "${filtered_vcf_files[@]}"; do
            mv $file batch$((folder_index))/$file
            folder_index=$(( (folder_index % $batch_size) + 1 ))
        done

        # Create a tar for each batch directory
        for i in $(seq 1 "${batch_size}"); do
        # Check if files exist in batch directory before creating a tar
            if [ -n "$(find "batch${i}" -maxdepth 1 -name '*.filtered.vcf')" ]; then
                tar -czf "${i}.filtered_vcf_files.tar.gz" "batch${i}"/*.filtered.vcf
            else
                echo "No files found in batch${i}. Skipping tar creation."
            fi
        done

    >>>
    runtime {
        docker: docker_path
        memory: "${memory_mb} MiB"
        cpu: cpu
        disks: "local-disk ${disk_size_gb} HDD"
    }

    output {
        Array[File] tarred_filtered_vcfs = glob("*.tar.gz")
    }
}

task AnnotateVCF {
    input {
        Array[File] input_filtered_vcf_tars
        String cloud_provider
        String docker_path
        String output_prefix
    }

    String nirvana_location = "/Nirvana/Nirvana.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"


    command <<<
        set -euo pipefail

        # Declare array of filtered vcf tars
        filtered_vcf_files=(~{sep=' ' input_filtered_vcf_tars})

        # Loop through the filtered_vcf_files and untar them
        for file in "${filtered_vcf_files[@]}"; do
        (
            echo "extracting $file"
            tar -xf $file --strip-components=1
            rm $file
        )
        done

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


        # Bash function to perform the annotation task on a single VCF file
        task() {
            local file=$1
            sample_id=$(basename "$file" ".filtered.vcf")
            echo $sample_id

            # Create Nirvana annotations:
            dotnet ~{nirvana_location} \
                -i ${sample_id}.filtered.vcf \
                -c $DATA_SOURCES_FOLDER~{path} \
                --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
                -r $DATA_SOURCES_FOLDER~{path_reference} \
                -o ${sample_id}

              # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
              # Parse out the Genes section into a separate annotated json
              dotnet  ~{jasix_location} \
                  --in ${sample_id}.json.gz \
                  --section genes \
                  --out ~{output_prefix}_${sample_id}.genes.json.gz

          # Parse out the Positions section into a separate annotated json
          dotnet  ~{jasix_location} \
              --in ${sample_id}.json.gz \
              --section positions \
              --out ~{output_prefix}_${sample_id}.positions.json.gz
        }

         # Define lists of vcf files
          vcf_files=($(ls | grep "filtered.vcf$"))
          echo "vcf_files: ${vcf_files[@]}"

         # Run 2 instances of the task in parallel
         for file in "${vcf_files[@]}"; do
           (
             echo "starting task $file.."
             task "$file"
             sleep $(( (RANDOM % 3) + 1))
           ) &
           # allow to execute up to 2 jobs in parallel
           if [[ $(jobs -r -p | wc -l) -ge 2 ]]; then
             wait -n
           fi
         done

         wait
         echo "Tasks all done."

        # Tar up the genes.json and positions.json files
        tar -czf ~{output_prefix}.genes_annotation_json.tar.gz *.genes.json.gz
        tar -czf ~{output_prefix}.positions_annotation_json.tar.gz *.positions.json.gz
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
        File genes_annotation_json = "~{output_prefix}.genes_annotation_json.tar.gz"
        File positions_annotation_json = "~{output_prefix}.positions_annotation_json.tar.gz"
    }
}

task VariantReport {

    input {
        Array[File] positions_annotation_json
        File bed_file
        String output_prefix


        String docker_path
        Int memory_mb = 8000
        Int disk_size_gb = 64
        Int cpu = 1
    }

    command <<<

        set -euo pipefail

        # Loop through the json files and unpack them
        declare -a input_jsons=(~{sep=' ' positions_annotation_json})

        for json in "${input_jsons[@]}"; do
            tar -xzf $json
        done

        # Loop through the json.gz files and get the sample_id by using basename
        for json in $(ls *.json.gz); do
            sample_id=$(basename "$json" ".positions.json.gz")
            echo "Creating the report for $sample_id"
            python3 /src/variants_report.py \
            --positions_json $json \
            --sample_identifier $sample_id \
            --output_prefix ~{output_prefix} \
            --bed_file ~{bed_file}
        done

        ls -lh

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

task PipelineMetadata {
    input {
        String pipeline_version
        String output_prefix
        Array[File] input_vcfs
        Array[File] pdf_report
        File bed_file

        String memory_mb = 4000
        String disk_size_gb = 10
        String cpu = 1
    }
    command <<<

        # Variables
        pipeline_version=v~{pipeline_version}
        declare -a input_vcfs=(~{sep=' ' input_vcfs})
        declare -a pdf_reports=(~{sep=' ' pdf_report})

        echo -e "input_vcf\tvariant_report_pdf" >> ~{output_prefix}_pipeline_metadata.txt

        for ((i=0; i<${#pdf_reports[@]}; i++)); do
            basename_pdf=$(basename "${pdf_reports[i]}")
            basename_vcf=$(basename "${input_vcfs[i]}")
            echo -e "$basename_vcf\t$basename_pdf" >> ~{output_prefix}_pipeline_metadata.txt
        done

        echo "" >> ~{output_prefix}_pipeline_metadata.txt

        # Write the pipeline version to the file
        date=$(date)
        bed_file=$(basename ~{bed_file})
        echo "The variant report PDFs were generated on $date using the $bed_file file and $pipeline_version version of this pipeline" >> ~{output_prefix}_pipeline_metadata.txt

    >>>
    runtime {
        docker: "ubuntu:latest"
        memory: "${memory_mb} MiB"
        cpu: cpu
        disks: 'local-disk ${disk_size_gb} HDD'
    }
    output {
        File pipeline_metadata = "~{output_prefix}_pipeline_metadata.txt"
    }
}

