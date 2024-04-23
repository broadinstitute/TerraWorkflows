version 1.0

workflow AnnotateVCFWorkflow {
    input {
        Array[File] input_vcf
        File bed_file
        String output_annotated_file_name
        Boolean use_reference_disk
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


call BatchVCFs as batch_vcfs {
    input:
        input_vcfs = input_vcf,
        batch_size = batch_size,
        docker_path = ubuntu_docker_path
}

scatter(tar in batch_vcfs.batch_tars) {
    call FilterVCF as filter_vcf {
            input:
                batch_tars = tar,
                bed_file = bed_file,
                output_annotated_file_name = output_annotated_file_name,
                batch_size = batch_size,
                docker_path = gatk_docker_path
        }
}


scatter(tar in filter_vcf.tarred_filtered_vcfs) {
    call AnnotateVCF as AnnotateVCF {
            input:
                input_filtered_vcf_tars = tar,
                output_annotated_file_name = output_annotated_file_name,
                use_reference_disk = use_reference_disk,
                cloud_provider = cloud_provider,
                omim_annotations = omim_annotations,
                docker_path = docker_prefix + nirvana_docker_image
        }
    }

    output {
        Array[File] positions_annotation_json = AnnotateVCF.positions_annotation_json
        Array[File] genes_annotation_json = AnnotateVCF.genes_annotation_json
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
    Int disk_size = ceil(size(input_vcfs, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        declare -a input_vcfs=(~{sep=' ' input_vcfs})

        # Get the size of the input_vcfs array
        num_vcfs=${#input_vcfs[@]}

        # Split input_vcfs into batches based on batch_sizes
        # Calculate the number of batches based on the batch size
        batch_size=~{batch_size}
        num_batches=$(( ($num_vcfs + $batch_size - 1) / $batch_size ))


        echo "Total VCF files: $num_vcfs"
        echo "Batch size (the number of files we will process in a batch): ~{batch_size}"
        echo "Number of batches: $num_batches"

        # Process each batch
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
        String output_annotated_file_name
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

        # Function to perform the task on a single VCF file
        task() {
            local vcf_file=$1

            echo "Starting task for $vcf_file.."
            sample_id=$(basename "$vcf_file" ".rb.g.vcf")

            # perform indexing
            echo "Indexing VCF file: $vcf_file"
            gatk \
            IndexFeatureFile \
            -I "$vcf_file"

            # perform filtering
            echo "Filtering VCF file: $vcf_file"
            gatk \
            SelectVariants \
            -V  "$vcf_file" \
            -L ~{bed_file} \
            -O "$sample_id.filtered.vcf"
        }

        # declare array of vcf tars
        declare -a batch_tars=(~{sep=' ' batch_tars})
        echo "batch_tars: ${batch_tars[@]}"

        # untar the tarred inputs
        tar -xf $batch_tars --strip-components=1

        declare -a input_vcfs=($(ls | grep ".rb.g.vcf$"))
        echo "input_vcfs: ${input_vcfs[@]}"


        for ((i = 0; i < ${#input_vcfs[@]}; i += 2)); do
            # Launch tasks in parallel for the current batch of VCF files
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
        ls -lR

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


        for i in $(seq 1 "${batch_size}"); do
        # Check if files exist in batch directory before creating tar archive
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
        String output_annotated_file_name
        Boolean use_reference_disk
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

        #declare array of filtered vcf tars

        filtered_vcf_files=(~{sep=' ' input_filtered_vcf_tars})

        # untar the tarred inputs
        tar -xf $filtered_vcf_files --strip-components=1
        rm $filtered_vcf_files

      if [[ "~{use_reference_disk}" == "true" ]]
      then
          # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
          # disk image was built: while all the reference files do exist on the image they are not at the expected
          # locations. The following code works around this issue and should continue to work even after a corrected
          # version of the Nirvana reference image is deployed into Terra.

          # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
          # mount point" because we do not actually confirm this is a reference disk until the following code block.
          CANDIDATE_MOUNT_POINT=$(lsblk | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
          if [[ -z ${CANDIDATE_MOUNT_POINT} ]]; then
              >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
              exit 1
          fi

          # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
          # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
          REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
          REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
          if [[ -z ${REFERENCE_PATH} ]]; then
              >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
              exit 1
          fi

          # Take the parent of the parent directory of this file as root of the locally mounted  references:
          DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
      else
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

      fi

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
                  --out ${sample_id}.genes.json.gz

          # Parse out the Positions section into a separate annotated json
          dotnet  ~{jasix_location} \
              --in ${sample_id}.json.gz \
              --section positions \
              --out ${sample_id}.positions.json.gz
        }

         # define lists of vcf files
          vcf_files=($(ls | grep "filtered.vcf$"))

         # run 2 instances of task in parallel
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

        #tar up the genes.json and positions.json files
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