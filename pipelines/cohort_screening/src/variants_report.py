import gzip
import json
import argparse
import pandas as pd
from gtfparse import read_gtf

def get_MANE():
    df = read_gtf('/src/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz')
    df = df['transcript_id'].to_pandas() #polar dataframe to pandas dataframe
    df = pd.DataFrame(df)
    df[['transcript','version']] = df['transcript_id'].str.split('.', expand=True)
    return df['transcript']


def isMANE(transcript_id, MANE_list):
    transcript = transcript_id.split('.')[0] #remove version
    return MANE_list.str.contains(transcript).any()

def get_impact_score(consequence, impact_table):
    return impact_table.loc[impact_table['consequence'] == consequence, 'impact_score'].item()

def filter_transcripts(positions):
    variants_field = 'variants'
    transcripts_field = 'transcripts'

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

    # select transcripts in MANE
    # download MANE list or put it in the docker???
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

    # select transcript with highest impact
    # impact table was manually  created from https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html --> download to docker??
    impact_table = pd.read_csv('/src/impact_table_with_score.csv', index_col=False)

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

    # if there are any variants with >1 transcript left, sort by transcript id and select the first one
    for position in positions:
        if variants_field in position:
            for variant_dict in position[variants_field]:
                if transcripts_field in variant_dict and len(variant_dict[transcripts_field])>1:
                    sorted_transcripts = sorted(variant_dict[transcripts_field], key=lambda x: x['transcript'])
                    variant_dict[transcripts_field] = sorted_transcripts[0]

    return positions

def report(args):

    # load input json files
    genes_json = gzip.open(args.genes_json, 'rt')
    genes = json.load(genes_json)

    positions_json = gzip.open(args.positions_json, 'rt')
    positions = json.load(positions_json)

    # do the filtering
    filtered_positions = filter_transcripts(positions)

    #for testing
    with open(args.output_file_name, 'w') as outfile:
        json.dump(filtered_positions, outfile)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse NIRVANA JSON file and create a Variant Report.')
    parser.add_argument('--positions_json', type=str, help='path to nirvana positions json', required=True)
    parser.add_argument('--genes_json', type=str, help='path to nirvana genes json', required=True)
    parser.add_argument('--output_file_name', type=str, help='filename for output', required=True)


    args = parser.parse_args()

    report(args)