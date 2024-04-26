import sys
import pandas as pd

sample_name = sys.argv[1]
gene_list_file = sys.argv[2]
file_paths = sys.argv[3]
#this script is intended to be called from a single-sample workflow
#assume only files from one sample

#---------cyrius------------
exit_code = 0

import cyp2d6_parser

with open(file_paths, 'r') as file:
    file_list = [line.strip() for line in file.readlines()]

with open(gene_list_file, 'r') as file:
    called_genes = [line.strip() for line in file.readlines()]
called_genes.append('cyp2d6')
called_genes.append('dpyd')

print("called genes: ", called_genes, file=sys.stderr)

file_lookup = dict()
#this isn't super efficient, but there should only be ~20 genes
for file in file_list:
    for gene in called_genes:
        if file.endswith(sample_name + '_' + gene + '_genotype-calls.tsv'):
            file_lookup[gene] = file
            continue

import re

call_data_cyrius = []
if 'cyp2d6' in file_lookup:
    cyp2d6_parser.parse_genotype(file=file_lookup['cyp2d6'], caller="cyrius", output_file = 'cyp2d6_results.tsv')
    call_data_cyrius = pd.read_csv('cyp2d6_results.tsv', sep='\t')
    match = re.search(r'Score (.*)', call_data_cyrius.iloc[0]['activity_score'])
    if match:
        cyp2d6_activity_score = match.group(1)
    else:
        cyp2d6_activity_score = 'NaN'
    call_data_cyrius.iloc[0]['activity_score'] = cyp2d6_activity_score
else:
    print('ERROR: Call file for cyp2d6 not found')
    exit_code = 1

#----------stargazer-------------



stargazer_csv_lookup_genes = ['2c_cluster', 'abcg2', 'cacna1s',
                    'cyp2c19', 'cyp2c9', 'cyp3a5', 'cyp4f2', 'dpyd', 'g6pd',
                    'nudt15', 'ryr1', 'slco1b1', 'tpmt', 'ugt1a1', 'vkorc1']



samples = [sample_name]
print("sample name: " + sample_name)

stargazer_parser_genes = set(called_genes).intersection(stargazer_csv_lookup_genes)

#subset file_lookup dict to genes in the stargazer parser table
subset_dict = {key: file_lookup[key] for key in stargazer_parser_genes if key in file_lookup}

import stargazer_parser
stargazer_parser.run_me(subset_dict) #stargazer phasing csv doesn't have cftr, cyp2b6, cyp3a4, ifnl3
stargazer_file = "call_data_stargazer.pkl"
call_data_stargazer = pd.read_pickle(stargazer_file)

call_data_stargazer = call_data_stargazer.rename(columns = {'person_id':'sample_id'})

dfs = []
for gene_name, file_path in file_lookup.items():
    if gene_name not in stargazer_csv_lookup_genes and gene_name != 'cyp2d6': #did these already
        df = pd.read_csv(file_path, sep='\t')  # Assuming files are tab-separated
        df['gene'] = gene_name.upper()
        df = df.rename(columns = {'dip_score' : 'activity_score'})
        df = df.rename(columns = {'name' : 'sample_id'})
        df['genotype'] = df['hap1_main'] + '/' + df['hap2_main']
        df['phenotype'] = df['phenotype'].str.replace('_', ' ').str.title()
        if not gene_name.startswith('cyp') or df.iloc[0]['activity_score'] < 0:
            df['activity_score'] = ''
        dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
stargazer_df = pd.concat(dfs, ignore_index=True)

stargazer_df = stargazer_df[['sample_id','gene','genotype','activity_score','phenotype']]

call_data_stargazer = call_data_stargazer[['sample_id','gene','genotype','activity_score','phenotype']]

report = pd.concat([stargazer_df, call_data_stargazer])
report['copy_number'] = '.' #we assume 2, but don't actually check

call_data_cyrius['gene'] = 'CYP2D6'
call_data_cyrius['sample_id'] = sample_name #cyrius seems to get the sample name from the cram path
report = pd.concat([report, call_data_cyrius[['sample_id', 'gene', 'genotype', 'activity_score', 'phenotype', 'copy_number']]])
report.sort_values('gene').to_csv(sample_name + ".pgx_report.tsv", sep = '\t', index=False)

sys.exit(exit_code)