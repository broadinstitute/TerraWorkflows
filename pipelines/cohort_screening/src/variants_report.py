import gzip
import json
import argparse
import logging

import pandas as pd
import pdfkit
from gtfparse import read_gtf
from tqdm import tqdm


def get_MANE():
    df = read_gtf('/src/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz') #  TODO: change back to /src/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz
    df = df['transcript_id'].to_pandas() # polar dataframe to pandas dataframe
    df = pd.DataFrame(df)
    df[['transcript', 'version']] = df['transcript_id'].str.split('.', expand=True)
    return df['transcript']


def isMANE(transcript_id, MANE_list):
    transcript = transcript_id.split('.')[0] # remove version
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
    MANE = get_MANE()
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                transcripts_to_keep = []
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])>1:
                    for transcript_dict in variant_dict[transcripts_field]:
                        if isMANE(transcript_dict['transcript'], MANE):
                            transcripts_to_keep.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_to_keep
        progress_bar.set_description('Selecting transcripts in MANE...')
        progress_bar.update(1)

    # select transcript with highest impact
    # impact table was manually  created from https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
    impact_table = pd.read_csv('/src/impact_table_with_score.csv', index_col=False) # TODO: change back to /src/impact_table_with_score.csv

    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])>1:
                    transcripts_with_highest_impact = []
                    highest_impact = 0
                    for transcript_dict in variant_dict[transcripts_field]:
                        transcript_impact = max([get_impact_score(consequence, impact_table) for consequence in transcript_dict['consequence']])
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
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])>1:
                    for transcript_dict in variant_dict[transcripts_field]:
                        if 'source' in transcript_dict and transcript_dict['source']=='Ensembl':
                            transcripts_to_keep.append(transcript_dict)
                    variant_dict[transcripts_field] = transcripts_to_keep
        progress_bar.set_description('Selecting Ensembl transcript...')
        progress_bar.update(1)

    # if there are any variants with >1 transcript left, sort by transcript id and select the first one
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])>1:
                    sorted_transcripts = sorted(variant_dict[transcripts_field], key=lambda x: x['transcript'])
                    variant_dict[transcripts_field] = sorted_transcripts[0]
        progress_bar.set_description('Selecting first transcript when sorted alphabetically...')
        progress_bar.update(1)

    # remove any variants that don't have any transcripts left
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field][:]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])==0:
                    position[variants_field].remove(variant_dict)
        progress_bar.set_description('Remove any variants that are empty...')
        progress_bar.update(1)

    # remove any positions that don't have any variants left
    for position in positions[:]:
        if variants_field in position and len(position[variants_field])==0:
            positions.remove(position)
        progress_bar.update(1)

    progress_bar.close()

    return positions


# TODO: Bobbie to implement this function
# do this before we map to report table
# selects the variants that meet the inclusion criteria
def select_variants_for_report(filtered_positions):
    selected_variants = []
    for position in filtered_positions:
        # check to see if we have any variants to report
        if 'variants' not in position:
            logging.warning('no variants found for position, creating empty report')
            selected_variants = None
            return selected_variants
        else:
            for variant_dict in position['variants']:
                transcript = variant_dict['transcripts']
                # check to see if we have too few or too many transcripts - should be 1 per variant
                if 'transcripts' in variant_dict and len(transcript) == 0:
                    logging.warning('no transcripts selected for variant, skipping variant')
                    continue
                elif len(transcript) > 1:
                    logging.warning('more than one transcript selected for variant, skipping variant')
                    continue
                else:
                    # Ignore if ClinVar classification is not present
                    if 'clinvar' not in 'variants':
                        logging.warning('no ClinVar classification found for variant, skipping variant')
                        continue
                    # Include if ClinVar significance includes pathogenic or likely pathogenic,
                    # but does not include benign
                    # TODO FIX ME
                    elif 'clinvar' in 'variants' and 'significance' in variants['ClinVar'] and 'pathogenic' in transcript['ClinVar']['classification'] and 'benign' not in transcript['ClinVar']['classification']:
                        selected_variants.append(variant_dict)
                    # Disregard ClinVar classifications that are tagged with a review status of:
                    # - no assertion criteria provided
                    # - no classification provided
                    # - no classification for the individual variant
                    if 'ClinVar' in 'variants' and 'classification' in transcript['ClinVar']:
                        if 'review_status' in transcript['ClinVar'] and transcript['ClinVar']['review_status'] not in ['no assertion criteria provided', 'no classification provided', 'no classification for the individual variant']:
                            selected_variants.append(variant_dict)
    if len(selected_variants) == 0:
        selected_variants = None
    return selected_variants
# returns selected_variants

# TODO: Bobbie to implement this function
# maps the filtered variants to the report table
# def map_to_report(selected_variants)
# returns mapped_variants
# TODO: create a for loop to iterate over the variants and *append* them to the DataFrame
# NB this is part of the mapping
# df_var = df_var.append(mapped_variants, ignore_index=True)


# formats the report
def format_report():
# def format_report(filtered_positions):
    # dummy data
    # sample_identifier = args.sample_identifier)
    # comes from the data table - can get from args
    sample_identifier = args.sample_identifier

    # mapped_variants = map_to_report(select_variants_for_report(filtered_positions))
    # TODO - replace with mapped_variants from map_to_report
    variants = [
        {
            "Gene": "Gene 1",
            "contig": "contig 1",
            "position": "position 1",
            "Ref allele": "Ref allele 1",
            "Alt Allele": "Alt Allele 1",
            "dbSNP": "dbSNP 1",
            "Zygosity": "Zygosity 1",
            "consequence": "consequence 1",
            "Protein change": "Protein change 1",
            "gnomAD AF": "gnomAD AF 1",
            "ClinVar classification": "ClinVar classification 1",
            "ClinVar phenotypes": "ClinVar phenotypes 1",
        },
        {
            "Gene": "Gene 2",
            "contig": "contig 2",
            "position": "position 2",
            "Ref allele": "Ref allele 2",
            "Alt Allele": "Alt Allele 2",
            "dbSNP": "dbSNP 2",
            "Zygosity": "Zygosity 2",
            "consequence": "consequence 2",
            "Protein change": "Protein change 2",
            "gnomAD AF": "gnomAD AF 2",
            "ClinVar classification": "ClinVar classification 2",
            "ClinVar phenotypes": "ClinVar phenotypes 2",
        },
    ]

    # Convert the list of dummy variant dictionaries to a DataFrame
    # TODO replace with mapped_variants
    if mapped_variants is None:
        logging.info('no variants to report, creating empty report')
        identified_variants_message = "No variants marked for further investigation"
    else:
        df_var = pd.DataFrame(variants)
        identified_variants_message = "We have identified the following variants that warrant further investigation:"


    # Apply styles to the DataFrame
    style = [{'selector': 'th',
              'props': [('border', '1px solid black'),
                        ('background-color', '#f2f2f2'),
                        ('font-size', '14px'),
                        ('padding', '5px'),
                        ('text-align', 'left')]
              },
             {'selector': 'td',
              'props': [('border', '1px solid black'),
                        ('border-collapse', 'collapse'),
                        ('border-spacing', '0'),
                        ('font-size', '14px'),
                        ('padding', '5px'),
                        ('text-align', 'left')]
              }]
    # Convert the DataFrame to an HTML table
    df_var_html = df_var.style.hide().set_table_styles(style).to_html()
    # df_var_html = df_var.style.hide().to_html()

    # MOBY (Maturity-Onset Diabetes of the Young) genes
    # obtained from Naylor 2018, https://www.ncbi.nlm.nih.gov/books/NBK500456/#mody-ov.Genetic_Causes_of_MODY
    # TODO: Add FAQ#2 in this comment and to report
    # see spec
    mody_genes = ["ABCC8", "APPL1", "BLK", "CEL", "GCK", "HNF1B", "HNF4A", "INS", "KCNJ11", "NEUROD1", "PAX4", "PDX1"]

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
        <p>Please note:<br>
        <ul style="font-size: 16px;">
            <li>This report is for research purposes only.</li>
            <li>This report is not meant for direct disclosure to patients.</li>
            <li>This report does not mean that the person with these variants has diabetes.</li>
            <li>Clinicians should confirm this result before using it as part of clinical care.</li>
        </ul>
        <p>{identified_variants_message}</p>
        {df_var_html}
    """
    html_content += """
        <h2 style="font-size: 20px;">We have examined the following genes, which have been implicated in 
        Maturity-Onset Diabetes of the Young (MODY):</h2>
        <ul>
    """

    for gene in mody_genes:
        html_content += f"<li>{gene}</li>"

    html_content += """
        </ul>
        <p>The analysis to identify these variants has limitations, which include potentially missing pathogenic 
        variants, missing variants that are not yet associated with MODY, and missing variants that are associated 
        with MODY as part of a polygenic effect.  This analysis did not use any familial genomic information 
        to confirm variants.</p>
    </body>
    </html>
    """
    # write out the HTML content to a file
    with open('report.html', 'w') as html_file:
        html_file.write(html_content)
    logging.info("HTML report generated successfully!")

    # convert the HTML report to a PDF file
    # TODO - match with WDL
    pdf_report_name = args.sample_identifier + '_mody_variants_report.pdf'
    pdfkit.from_file('report.html', pdf_report_name)
    logging.info("PDF report generated successfully!")

    # Export variants dataframe as tsv
    # TODO - match with WDL
    table_name = args.sample_identifier + '_mody_variants_table.tsv'
    df_var.to_csv(table_name, sep='\t', index=False)
    logging.info("variants.tsv exported successfully!")


def report(args):
    # read in the positions json
    positions_json = gzip.open(args.positions_json, 'rt')
    positions = json.load(positions_json)

    # do the filtering
    filtered_positions = filter_transcripts(positions)

    # output the filtered positions to a json file if specified
    if args.output_file_name is None:
        pass
    else:
        with open("filtered_positions_final.json", 'w') as outfile:
            json.dump(filtered_positions, outfile)
            # json.dump(filtered_positions, outfile, indent=4) # add indent for pretty print, remove for py reading

    # for testing
    format_report()
    # format_report(
    #     map_to_report(
    #         select_variants_for_report(filtered_positions)
    #     )
    # )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse NIRVANA JSON file and create a Variant Report.')
    parser.add_argument('--positions_json', type=str, help='path to nirvana positions json', required=True)
    parser.add_argument('--sample_identifier', type=str, help='sample id', required=True)

    args = parser.parse_args()

    report(args)
