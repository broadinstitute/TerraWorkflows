version 1.0

workflow variantsreport{
	input {
    	File positions_annotation_json
		String sampleid
    }

    call parsejson {
    	input:
    		positions_annotation_json = positions_annotation_json,
		    sampleid = sampleid,
    }

    output {
    	File report = parsejson.pdf_report
    	File tsv = parsejson.tsv_file
    }
}

task parsejson {

	input {
		File positions_annotation_json
        String output_file_name = 'filtered_positions.json'
        String sampleid

        String docker = 'samclairehv/variantsreport:v5'
        Int mem_gb = 4
        Int disk_gb = 15
	}

	command {

		python3 /src/variants_report.py \
			--positions_json ~{positions_annotation_json} \
			--output_file_name ~{output_file_name} \
			--sample_identifier ~{sampleid}

	}

	runtime {
		docker: docker
		memory: '${mem_gb} GB'
		disks: 'local-disk ${disk_gb} HDD'
	}

	output {
		File pdf_report = "{sampleid}_mody_variants_report.pdf"
		File tsv_file = "{sampleid}_mody_variants.tsv"
	}
}