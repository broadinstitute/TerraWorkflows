version 1.0

workflow variantsreport{
	input {
    	File positions_annotation_json
		String sample_id
    }

    call parsejson {
    	input:
    		positions_annotation_json = positions_annotation_json,
		    sample_id = sample_id
    }

    output {
    	File report = parsejson.pdf_report
    	File tsv = parsejson.tsv_file
    }
}

task parsejson {

	input {
		File positions_annotation_json
        String sample_id

        String docker = 'terraworkflows.azurecr.io/variantreport:testing'
        Int mem_gb = 4
        Int disk_gb = 15
	}

	command {

		python3 /src/variants_report.py \
			--positions_json ~{positions_annotation_json} \
			--sample_identifier ~{sample_id}

	}

	runtime {
		docker: docker
		memory: '${mem_gb} GB'
		disks: 'local-disk ${disk_gb} HDD'
	}

	output {
		File pdf_report = "~{sample_id}_mody_variants_report.pdf"
		File tsv_file = "~{sample_id}_mody_variants_table.tsv"
	}
}