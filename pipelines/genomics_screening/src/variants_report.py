import gzip
import json
import argparse
import logging
import pandas as pd
import pdfkit
import re
from gtfparse import read_gtf
from tqdm import tqdm

# static file because we are using this for filtering and do not want it to change.
def get_MANE():
    df = read_gtf('/src/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz')
    df = df['transcript_id'].to_pandas()  # polar dataframe to pandas dataframe
    df = pd.DataFrame(df)
    df[['transcript', 'version']] = df['transcript_id'].str.split('.', expand=True)
    return df['transcript']


def isMANE(transcript_id, MANE_list):
    transcript = transcript_id.split('.')[0]  # remove version
    return MANE_list.str.contains(transcript).any()


def get_impact_score(consequence, impact_table):
    return impact_table.loc[impact_table['consequence'] == consequence, 'impact_score'].item()


def filter_transcripts(positions):
    variants_field = 'variants'
    transcripts_field = 'transcripts'

    print('Filtering variant transcripts...')
    total_iterations = len(positions)*7
    progress_bar = tqdm(total=total_iterations, unit="iteration")

    # filter to canonical transcripts only
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                transcripts_to_keep = []
                if transcripts_field in variant_dict:
                    for transcript_dict in variant_dict[transcripts_field]:
                        if 'isCanonical' in transcript_dict and True:
                            transcripts_to_keep.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_to_keep
        progress_bar.set_description('Selecting canonical transcripts...')
        progress_bar.update(1)

    # select transcripts in MANE
    # if none are in MANE, continue with transcripts
    MANE = get_MANE()
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                transcripts_to_keep = []
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field]) > 1:
                    for transcript_dict in variant_dict[transcripts_field]:
                        if isMANE(transcript_dict['transcript'], MANE):
                            transcripts_to_keep.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_to_keep if len(transcripts_to_keep) > 1 \
                        else variant_dict[transcripts_field]
        progress_bar.set_description('Selecting transcripts in MANE...')
        progress_bar.update(1)

    # select transcript with highest impact
    # impact table was manually created from
    # https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
    impact_table = pd.read_csv('/src/impact_table_with_score.csv', index_col=False)
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field]) > 1:
                    transcripts_with_highest_impact = []
                    highest_impact = 0
                    for transcript_dict in variant_dict[transcripts_field]:
                        transcript_impact = max([get_impact_score(consequence, impact_table)
                                                 for consequence in transcript_dict['consequence']])
                        if transcript_impact > highest_impact:
                            highest_impact = transcript_impact
                            transcripts_with_highest_impact = [transcript_dict]
                        elif transcript_impact == highest_impact:
                            transcripts_with_highest_impact.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_with_highest_impact
        progress_bar.set_description('Selecting transcript with highest impact...')
        progress_bar.update(1)

    # choose the ensembl transcript
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                transcripts_to_keep = []
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field]) > 1:
                    for transcript_dict in variant_dict[transcripts_field]:
                        if 'source' in transcript_dict and transcript_dict['source'] == 'Ensembl':
                            transcripts_to_keep.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_to_keep
        progress_bar.set_description('Selecting Ensembl transcript...')
        progress_bar.update(1)

    # if there are any variants with >1 transcript left, sort by transcript id and select the first one
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field]) > 1:
                    sorted_transcripts = sorted(variant_dict[transcripts_field], key=lambda x: x['transcript'])
                    variant_dict[transcripts_field] = [sorted_transcripts[0]]
        progress_bar.set_description('Selecting first transcript when sorted alphabetically...')
        progress_bar.update(1)

    # remove any variants that don't have any transcripts left
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field][:]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field]) == 0:
                    position[variants_field].remove(variant_dict)
        progress_bar.set_description('Remove any variants that are empty...')
        progress_bar.update(1)

    # remove any positions that don't have any variants left
    for position in positions[:]:
        if variants_field in position and len(position[variants_field]) == 0:
            positions.remove(position)
        progress_bar.update(1)

    progress_bar.close()

    return positions


def filter_variants_for_report(filtered_positions):
    variants_field = 'variants'
    clinvar_field = 'clinvar'

    removed_positions = 0

    print('Filtering variants for inclusion in report based on ClinVar significance...')
    # Variant Report Inclusion Criteria 1.a: Include if the variant has a ClinVar classification
    for position in filtered_positions[:]:
        for variant_dict in position[variants_field][:]:
            # remove those variants that do not have a ClinVar classification
            if clinvar_field not in variant_dict:
                logging.info(f'no ClinVar classification found for variant, /'
                             f'removing variant' + str(variant_dict))
                position[variants_field].remove(variant_dict)
            # remove those classifications that have a review status of
            # no assertion criteria provided,
            # no classification provided,
            # or no classification for the individual variant
            elif clinvar_field in variant_dict:
                for clinvar_dict in variant_dict[clinvar_field][:]:
                    review_status_disregard = ['no assertion criteria provided',
                                               'no classification provided',
                                               'no classification for the individual variant']
                    if clinvar_dict['reviewStatus'] in review_status_disregard:
                        logging.info(f'review status does not pass inclusion filters, /'
                                     f'removing clinvar dict' + str(clinvar_dict))
                        variant_dict[clinvar_field].remove(clinvar_dict)
                    # Include if ClinVar significance includes pathogenic or likely pathogenic,
                    # but does not include benign
                    elif 'benign' in clinvar_dict['significance']:
                        logging.info(f'benign in significance, removing clinvar dict {clinvar_dict}')
                        variant_dict[clinvar_field].remove(clinvar_dict)
                    elif 'significance' in clinvar_dict and ('pathogenic' not in clinvar_dict['significance']
                                                             and 'likely pathogenic'
                                                             not in clinvar_dict['significance']):
                        logging.info(f'ClinVar significance does not pass inclusion filters, /'
                                     f'removing clinvar dict' + str(variant_dict))
                        variant_dict[clinvar_field].remove(clinvar_dict)

    # remove any variants that don't have any populated clinvar fields left
    for position in filtered_positions[:]:
        for variant_dict in position[variants_field][:]:
            if len(variant_dict[clinvar_field]) == 0:
                logging.info(f'no ClinVar fields remaining for variant, /'
                             f'removing variant' + str(variant_dict))
                position[variants_field].remove(variant_dict)

    # remove any positions that don't have any variants left
    for position in filtered_positions[:]:
        if variants_field in position and len(position[variants_field]) == 0:
            logging.info(f'no variants found for position, removing position' + str(position))
            filtered_positions.remove(position)
            removed_positions += 1

    # Criteria 2 & 3
    # remove any variants where genotype is 0/0, ./N, N/., or ./.
    print('Filtering variants for inclusion in report based on genotype...')
    pattern = r"\./|./\."
    for position in filtered_positions[:]:
        for samples_dict in position['samples'][:]:
            # if any of the genotypes match the pattern, remove the variant
            if re.match(pattern, samples_dict['genotype']) or samples_dict['genotype'] == '0/0':
                logging.info(f'genotype does not pass inclusion filters, removing sample' + str(samples_dict))
                position['samples'].remove(samples_dict)

    # remove any positions that don't have any samples left
    for position in filtered_positions[:]:
        if 'samples' in position and len(position['samples']) == 0:
            logging.info(f'position has no samples left, removing position' + str(position))
            filtered_positions.remove(position)
            removed_positions += 1

    # check len transcripts vs len variants
    print('Checking number of variants == number of transcripts...')
    filtered_positions_str = str(filtered_positions)
    num_variants = filtered_positions_str.count('variants')
    num_transcripts = filtered_positions_str.count('transcripts')
    print('Number of variants: ' + str(num_variants))
    print('Number of transcripts: ' + str(num_transcripts))

    if not num_variants == num_transcripts:
        raise ValueError('number of variants and transcripts do not match')
    elif num_variants == 0:
        variants_to_include = []
    else:
        print(f'removed_positions = {removed_positions}')
        variants_to_include = filtered_positions

    return variants_to_include


# calculate the zygosity of a genotype
def calculate_zygosity(genotype):
    alleles = genotype.split("/")
    if len(set(alleles)) == 1:
        if alleles[0] == "0":
            return "Homozygous reference"
        else:
            return "Homozygous alternate"
    else:
        if "0" in alleles:
            return "Heterozygous"
        else:
            return "Heterozygous alternate"


# get an index for sorting the classifications
def classification_sort(classification):
    classification_order = ["pathogenic",
                            "likely pathogenic",
                            "risk factor",
                            "uncertain significance",
                            "benign",
                            "likely benign",
                            "affects",
                            "association",
                            "conflicting data from submitters",
                            "drug response",
                            "not provided",
                            "other",
                            "protective",
                            "'-'"]
    return classification_order.index(classification)


def map_variants_to_table(variants_to_include):
    variants = []
    for variant in variants_to_include:

        genotype = variant["samples"][0]["genotype"]
        zygosity = calculate_zygosity(genotype)

        # make list of classifications
        clinvar_classifications = []
        for clinvar_dict in variant.get("variants", [{}])[0].get("clinvar", [{}]):
            clinvar_classifications.append(clinvar_dict.get("significance", "NA"))
        # Convert the list of lists into a list of strings
        classification_strings = [item[0] for item in clinvar_classifications]
        # Remove duplicates from the list
        unique_classifications = list(set(classification_strings))
        # sort the unique list by the classification
        sorted_classifications = sorted(unique_classifications, key=classification_sort)
        # Flatten the list
        classifications_list = ', '.join(sorted_classifications)

        # make list of phenotypes
        phenotypes = []
        for clinvar_dict in variant.get("variants", [{}])[0].get("clinvar", [{}]):
            phenotypes.append(clinvar_dict.get("phenotypes", "NA"))

        print(f'this is phenotypes {phenotypes}')

        # Convert the mixed list of lists and strings into a list of strings
        phenotype_strings = [item[0] if isinstance(item, list) else item for item in phenotypes]
        print(f'this is phenotype_strings {phenotype_strings}')
        # drop NA and not provided
        filtered_phenotypes = [item for item in phenotype_strings if item not in ['NA', 'not provided']]
        # Remove duplicates from the list
        unique_phenotypes = list(set(filtered_phenotypes))
        print(f'this is unique_phenotypes {unique_phenotypes}')
        # Flatten the list
        phenotypes_list = ', '.join(unique_phenotypes)
        print(f'phenotypes_list {phenotypes_list}')

        # map the variant data to a dictionary
        # nullable uses get() to avoid KeyError
        data = {
            "gene": variant["variants"][0]["transcripts"][0]["hgnc"],
            "contig": variant["variants"][0]["chromosome"],
            "position": variant["position"],
            "Reference allele": variant["refAllele"],
            "Alternate Allele": ', '.join(map(str, variant["altAlleles"])),
            "dbSNP": ', '.join(map(
                str, variant.get("variants", [{}])[0].get("dbsnp", ["NA"]))
            ),
            "HGVSC": variant.get("variants", [{}])[0].get("transcripts", [{}])[0].get("hgvsc", "NA"),
            "Zygosity": zygosity,
            "consequence": ', '.join(map(
                str, variant.get("variants", [{}])[0].get("transcripts", [{}])[0].get("consequence", "NA"))
            ),
            "Protein change": variant.get("variants", [{}])[0].get("transcripts", [{}])[0].get("hgvsp", "NA"),
            "gnomAD AF": variant.get("variants", [{}])[0].get("gnomad", {}).get("allAf", "NA"),
            "ClinVar Classification": classifications_list,
            "ClinVar Phenotype": phenotypes_list,
        }
        variants.append(data)

    # Convert the list of variant dictionaries to a DataFrame
    mapped_variants = pd.DataFrame(variants)
    return mapped_variants


# formats the report
def format_report(mapped_variants):
    sample_identifier = args.sample_identifier
    output_prefix = args.output_prefix
    bed_file = args.bed_file

    if mapped_variants.empty:
        logging.info('no variants to report, creating empty report')
        identified_variants_message = ""
        df_var = mapped_variants
        df_var_html = "<h3>No variants marked for further investigation</h3>"
    else:
        identified_variants_message = "We have identified the following variants that warrant further investigation:"
        df_var = mapped_variants

        # Apply styles to the DataFrame
        style = [{'selector': 'th',
                  'props': [('border', '1px solid black'),
                            ('background-color', '#f2f2f2'),
                            ('font-size', '12px'),
                            ('padding', '5px'),
                            ('text-align', 'left')]
                  },
                 {'selector': 'td',
                  'props': [('border', '1px solid black'),
                            ('border-collapse', 'collapse'),
                            ('border-spacing', '0'),
                            ('font-size', '12px'),
                            ('padding', '5px'),
                            ('text-align', 'left')]
                  }]
        # Convert the DataFrame to an HTML table
        df_var_html = df_var.style.hide().set_table_styles(style).to_html()

    # MODY (Maturity-Onset Diabetes of the Young) genes
    # obtained from Naylor 2018, https://www.ncbi.nlm.nih.gov/books/NBK500456/#mody-ov.Genetic_Causes_of_MODY
    # FAQ:
    # Why are HNF1A and KLF11 missing from the MODY list?
    # While these two genes are included in the list of fourteen MODY genes,
    # they are also included in the list of Challenging Medically Relevant Genes (CMRG) (bed file).
    # The CMRG genes are difficult to call in hg38,
    # which is the reference genome for the proposed reports in this document.
    # In order to call variants in these genes effectively,
    # we would need to implement and run another annotation pipeline, which is out of scope for this project.
    mody_genes = ["ABCC8", "APPL1", "BLK", "CEL", "GCK", "HNF1B", "HNF4A", "INS", "KCNJ11", "NEUROD1", "PAX4", "PDX1"]

    genes = []
    # list of genes is 4th column of bed file
    with open(bed_file)as f:
        for line in f:
            genes.append(line.strip().split()[3])
    genes.sort()

    if genes == mody_genes:
        gene_list = mody_genes
        mody = True
    else:
        gene_list = genes
        mody = False

    # Generate HTML report
    html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Genomic Variant Report</title>
        </head>
        <body style="font-family: Arial, sans-serif;">
            <h1 style="font-size: 24px;">Genomic Variant Report</h1>
            <p>Sample identifier: {sample_identifier}</p>
            <p style="font-size: 14px;">Please note:<br>
            <ul style="font-size: 14px;">
                <li>This report is for research purposes only.</li>
                <li>This report is not meant for direct disclosure to patients.</li>
                <li>This report does not mean that the person with these variants has diabetes.</li>
                <li>Clinicians should confirm this result before using it as part of clinical care.</li>
            </ul>
            </p>
            <p>{identified_variants_message}</p>
            {df_var_html}
        """

    if mody:
        html_content += """
                <h3 style="font-size: 16px;">We have examined the following genes, which have been implicated in 
                Maturity-Onset Diabetes of the Young (MODY):</h3>
                <ul>
            """
    else:
        html_content += """
                <h3 style="font-size: 16px;">We have examined the following genes, which have been implicated in 
                the phenotype:</h3>
                <ul>
            """

    for gene in gene_list:
        html_content += f"<li>{gene}</li>"

    if mody:
        html_content += """
                </ul>
                <p style="font-size: 14px;">The analysis to identify these variants has limitations, 
                which include potentially missing pathogenic variants, 
                missing variants that are not yet associated with MODY, and missing variants that are associated 
                with MODY as part of a polygenic effect.  This analysis did not use any familial genomic information 
                to confirm variants.</p>
                <p style="font-size: 14px;">Please note that we do not include HNF1A and KLF11, two MODY genes, 
                in this report, due to difficulties calling these genes on the hg38 reference 
                (<a href="http://dx.doi.org/10.1038/s41587-021-01158-1">10.1038/s41587-021-01158-1</a>).</p>
            </body>
            </html>
            """

    # write out the HTML content to a file
    with open('report.html', 'w') as html_file:
        html_file.write(html_content)
    logging.info("HTML report generated successfully!")

    # convert the HTML report to a PDF file
    pdf_report_name = f'{sample_identifier}_{"mody" if mody else output_prefix}_variants_report.pdf'
    pdfkit.from_file('report.html', pdf_report_name)
    logging.info("PDF report generated successfully!")

    # Export variants dataframe as tsv
    table_name = f'{sample_identifier}_{"mody" if mody else output_prefix}_variants_table.tsv'
    df_var.to_csv(table_name, sep='\t', index=False)
    logging.info("variants.tsv exported successfully!")


def report(args):
    # read in the positions json
    positions_json = gzip.open(args.positions_json, 'rt')
    positions = json.load(positions_json)

    # do the filtering
    filtered_positions = filter_transcripts(positions)

    # output the filtered positions to a json file if specified
    if args.output_json is False:
        pass
    else:
        with open("filtered_positions.json", 'w') as outfile:
            json.dump(filtered_positions, outfile, indent=4)

    variants_to_include = filter_variants_for_report(filtered_positions)
    mapped_variants = map_variants_to_table(variants_to_include)
    format_report(mapped_variants)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse NIRVANA JSON file and create a Variant Report.')
    parser.add_argument('--positions_json', type=str, help='path to nirvana positions json', required=True)
    parser.add_argument('--sample_identifier', type=str, help='sample id', required=True)
    parser.add_argument('--bed_file', type=str, help='path to bed file', required=True)
    parser.add_argument('--output_prefix', type=str, help='output file prefix', required=True)
    parser.add_argument('--output_json', type=bool, default=False,
                        help='Specify to output the initial filtered json', required=False)

    args = parser.parse_args()

    report(args)
