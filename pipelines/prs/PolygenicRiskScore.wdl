version 1.0

workflow PolygenicRiskScore {
	input { 
	File linear_weights

	File vcf  # VCF for scoring
	Int scoring_mem = 16
	Int vcf_to_plink_mem = 8

	String basename # for naming the output of array scoring and the array projection files

    Boolean assume_missing_hom_ref = false # this can be used when using a whole genome vcf, where any uncalled sites can be assumed to be hom-ref.
	#  In this case, you must also provide ref_fasta and ref_fai
	File? ref_fasta
    File? ref_fai
  }
	call DetermineChromosomeEncoding {
		input:
			weights = linear_weights
	}

	call ScoreVcf {
		input:
			vcf = vcf,
			basename = basename,
			weights = linear_weights,
			base_mem = scoring_mem,
			chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding,
			assume_missing_hom_ref = assume_missing_hom_ref,
			ref_fasta = ref_fasta,
			ref_fai = ref_fai
	}

  output {
	File raw_scores = ScoreVcf.score
	File sites_scored = ScoreVcf.sites_scored
  }
}

# score with plink2
task ScoreVcf {
  input {
    File vcf
    String basename
    File weights
    Int base_mem = 8
    String? extra_args
    File? sites
    String? chromosome_encoding
    Boolean use_ref_alt_for_ids = false
    File? ref_fasta
    File? ref_fai
    Boolean assume_missing_hom_ref = false
  }

  Int runtime_mem = base_mem + 2
  Int plink_mem = ceil(base_mem * 0.75 * 1000)
  Int disk_space =  3*ceil(size(vcf, "GB")) + 20
  String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

  command <<<
    /plink2 --score ~{weights} header ignore-dup-ids list-variants no-mean-imputation \
    cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --set-all-var-ids ~{var_ids_string} --allow-extra-chr ~{extra_args} -vcf ~{vcf} dosage=DS \
    --new-id-max-allele-len 1000 missing ~{"--extract " + sites} --out ~{basename} --memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}

    python3 << "EOF"
    if ~{if assume_missing_hom_ref then "True" else "False"}:
      if ~{if defined(ref_fasta) then "False" else "True"}:
        raise FileNotFoundError("ref_fasta must be defined when assume_missing_hom_ref is set to true.")
      if ~{if defined(ref_fai) then "False" else "True"}:
        raise FileNotFoundError("ref_fai must be defined when assume_missing_hom_ref is set to true.")
      import csv
      import pysam
      import pandas as pd
      from sortedcontainers import SortedList

      with open("~{basename}.sscore.vars") as sites_file:
        sites = set(sites_file.read().splitlines())
      
      def sites_filter(s):
        if ~{if defined(sites) then "True" else "False"}:
          with open("~{sites}") as sites_filter_file:
            sites_filter = set(sites_filter_file.read().splitlines())
            return s in sites_filter
        else:
            return True

      def sort_fn(s):
        s_split = s.split(":")
        return (s_split[0], int(s_split[1]), s_split[2], s_split[3])
      
      all_sites = SortedList(sites, key=sort_fn)

      with open("~{weights}") as weights:
          reader=csv.reader(weights, delimiter='\t')
          next(reader)
          with pysam.FastaFile("~{ref_fasta}") as ref:
              hom_ref_score_addition=0
              hom_ref_dosage_addition = 0
              for row in reader:
                  if row[0] not in sites and sites_filter(row[0]):
                      all_sites.add(row[0])
                      id_split = row[0].split(":")
                      contig = id_split[0]
                      pos = int(id_split[1])
                      ref_base = ref.fetch(contig, pos-1, pos)
                      if ref_base == row[1]:
                          hom_ref_score_addition+=2*float(row[2])
                          hom_ref_dosage_addition+=2
      
      scores = pd.read_csv("~{basename}.sscore",sep="\t")
      scores['SCORE1_SUM']=scores['SCORE1_SUM'] + hom_ref_score_addition
      scores['NAMED_ALLELE_DOSAGE_SUM']=scores['NAMED_ALLELE_DOSAGE_SUM'] + hom_ref_dosage_addition
      scores['SCORE1_AVG'] = scores['SCORE1_SUM']/len(all_sites)/2

      scores.to_csv("~{basename}.sscore",sep="\t", index=False)
      with open("~{basename}.sscore.vars","w") as sites_out_file:
        for s in all_sites:
          sites_out_file.write(f'{s}\n') 


    EOF
  >>>

  output {
    File score = "~{basename}.sscore"
    File log = "~{basename}.log"
    File sites_scored = "~{basename}.sscore.vars"
  }

  runtime {
    docker: "terraworkflows.azurecr.io/plink2:@sha256:0ced514e651b08a93684ec6b0670e1a21380fecb6e07e725a9f9ec0d2c4a0e72"
    disks: "local-disk " + disk_space + " HDD"
    memory: runtime_mem + " GB"
  }
}


#plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
  input {
    File weights
  }

  command <<<
    python3 << "EOF"
    code = 'MT'
    with open("~{weights}") as weights_file:
      chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
      if any('chr' in c for c in chroms):
          if 'chrM' in chroms:
              code = 'chrM'
          else:
              code = 'chrMT'
      elif 'M' in chroms:
          code = 'M'

    with open('chr_encode_out.txt', 'w') as write_code_file:
        write_code_file.write(f'{code}\n')
    EOF
  >>>

  runtime {
    docker : "terraworkflows.azurecr.io/python:3.9.10"
  }

  output {
    String chromosome_encoding = read_string("chr_encode_out.txt")
  }
}