import pandas as pd
import numpy as np
import json
import requests
import os

def send_cpic_request(gene, request_type=None):
    request_types = {
        "allele_freq": f"https://api.cpicpgx.org/v1/population_frequency_view?genesymbol=eq.{gene}",
        "pop_freq": f"https://api.cpicpgx.org/v1/gene_result?genesymbol=eq.{gene}",
        "allele_fxn": f"https://api.cpicpgx.org/v1/allele?genesymbol=eq.{gene}"
    }
    url = request_types[request_type]
    response = requests.get(url)
    if not response.ok:
        raise ValueError(f"Problem accessing data for gene {gene}")
    return response

genes = [
    "2c_cluster",
    "abcg2",
    "cacna1s",
    "cyp2c19",
    "cyp2c9",
    "cyp3a5",
    "cyp4f2",
    'dpyd',
    "g6pd",
    "nudt15",
    'ryr1',
    "slco1b1",
    "tpmt",
    "ugt1a1",
    "vkorc1",
]
genes_with_activity_score = {'cyp2c9', 'dpyd'}
stargazer_converter_path = "/stargazer_converter.csv"
dip_to_phen = pd.read_json("https://api.cpicpgx.org/v1/diplotype").loc[
    :, ["genesymbol", "diplotype", "generesult", "totalactivityscore"]
]
dip_to_phen.columns = ["gene", "diplotype", "phenotype", "activity_score"]
dip_to_phen = dip_to_phen.sort_values(["gene", "diplotype"], ignore_index=True)
dip_to_phen_file =  "cpic_dip_to_phen.pkl"
dip_to_phen.to_pickle(dip_to_phen_file)
dip_to_phen_data = pd.read_pickle(dip_to_phen_file)
dip_to_phen_data_raw = dip_to_phen_data.drop(columns="activity_score")
dip_to_phen = {
    "VKORC1": {
        "rs9923231 reference (C)/rs9923231 reference (C)": "Variant Absent",
        "rs9923231 reference (C)/rs9923231 variant (T)": "Variant Present",
        "rs9923231 variant (T)/rs9923231 variant (T)": "Variant Present",
        "Indeterminate/Indeterminate": "Indeterminate",
    },
    "2C_CLUSTER": {
        "rs12777823 reference (G)/rs12777823 reference (G)": "Variant Absent",
        "rs12777823 reference (G)/rs12777823 variant (A)": "Variant Present",
        "rs12777823 variant (A)/rs12777823 variant (A)": "Variant Present",
        "not called": "Indeterminate",
        "Indeterminate/Indeterminate": "Indeterminate",
    },
    "CYP4F2":{} # Placeholder so later steps aren't overly complicated
}

for _, (gene, dip, phen) in dip_to_phen_data_raw.iterrows():
    try:
        dip_to_phen[gene].update({dip: phen})
    except KeyError:
        dip_to_phen[gene] = {dip: phen}
dip_to_activity_data = dip_to_phen_data.loc[dip_to_phen_data['gene'].isin({"CYP2C9", 'DPYD'})].drop(columns='phenotype')
dip_to_activity_score = {
    "DPYD": {
        "Indeterminate/Indeterminate": "n/a",
    },
    "CYP2C9": {
        "Indeterminate/Indeterminate": "n/a",
    },
}
for _, (gene, dip, activity) in dip_to_activity_data.iterrows():
    dip_to_activity_score[gene].update({dip: activity})
response = send_cpic_request("DPYD", request_type="allele_fxn")
json_data = json.loads(response.text)
dpyd_allele_function = {"Reference": 1.0}
for allele_data in json_data:
    allele = allele_data["name"]
    activity = float(allele_data["activityvalue"])
    dpyd_allele_function[allele] = activity

class StargazerTable:
    def __init__(self, file_path):
        self.file_path = file_path
        self.df = cat_cols(pd.read_csv(file_path, comment='#'))
        self._set_gene_data(attr="gene_data", col="genotype")
        self._set_gene_data(attr="gene_data_possible", col="possible_genotype")

    def _set_gene_data(self, attr, col):
        gene_data = {}
        for gene, group in self.df.groupby("gene"):
            gene_data_individual = {}
            for _, (dip, combo) in group.loc[:, [col, "cols_combined"]].iterrows():
                gene_data_individual[combo] = dip
            gene_data[gene] = gene_data_individual
        setattr(self, attr, gene_data)

def cat_cols(df):
    cols = ['hap1_main', 'hap2_main', 'hap1_cand', 'hap2_cand', 'dip_cand']
    df['cols_combined'] = df[cols].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    return df

def assign_genotype(calls, gene, stargazer_table, attr=None, col=None, ):
    # assigns for both the genotype and possible_genotype columns
    return calls.assign(
        **{col:lambda x: x["cols_combined"].map(getattr(stargazer_table, attr)[gene])}
    )

def assign_phenotype(calls, gene, gene_dip_to_phen):
    if gene == "cyp4f2":
        # The variant that defines *3 allele is the only one with a recommendation
        # This variant is also present in the new *4 allele
        return calls.assign(
            phenotype=np.select(
                [
                    calls["genotype"].str.contains("3|4", regex=True, na=False).to_numpy(dtype=bool),
                    (calls["genotype"] == "Indeterminate/Indeterminate").to_numpy(dtype=bool)
                ],
                ["CYP4F2 Variant Present", "CYP4F2 Indeterminate"],
                default="CYP4F2 Variant Absent",
            )
        )
    return calls.assign(
        #phenotype=gene.upper() + " " + calls["genotype"].map(gene_dip_to_phen)
        phenotype=calls["genotype"].map(gene_dip_to_phen)
    )

def assign_activity(calls, gene):
    if gene not in genes_with_activity_score:
        calls["activity_score"] = np.nan
        return calls
    gene_data = dip_to_activity_score[gene.upper()]
    calls["activity_score"] = (
        #gene.upper()
        #+ " Activity Score "
        calls["genotype"].map(gene_data).replace("Indeterminate", "n/a").fillna("n/a")
    )
    return calls

def get_possible_pheno(data):
    possible_phenotypes = data.unique()
    if possible_phenotypes.shape[0] != 1:
        return np.nan
    return possible_phenotypes[0]


def assign_possible_phenotype(calls, gene, gene_dip_to_phen):
    if gene == 'g6pd':
        # Excluding g6pd because it's not possible to assign any possible phenotypes
        # Inlcuding it would cause XY samples with 2 alleles to be translated
        return calls.assign(possible_phenotype=np.nan)
    calls = (
        calls.astype({"possible_genotype": "object"})
        .assign(possible_geno_split=lambda x: x["possible_genotype"].str.split(";"))
        .explode("possible_geno_split")
        .assign(
            possible_phenotype=lambda x: x["possible_geno_split"].map(gene_dip_to_phen)
        )
        .groupby(
            ["person_id", "genotype", "possible_genotype", "phenotype"],
            as_index=False,
            dropna=False,
        )
        .agg({"possible_phenotype": get_possible_pheno})
    )
    calls.loc[calls["possible_phenotype"].notna(), "possible_phenotype"] = (
        gene.upper() + " " + calls["possible_phenotype"].astype("object")
    )
    return calls

def adjust_dtypes(data):
    return data.astype(
        {
            #"person_id": np.int64,
            "person_id": "string", #TODO: change this back after testing
            "gene": "category",
            "genotype": "string",
            "phenotype": "string",
            "activity_score": "object",
        }
    )

def adjust_g6pd_xy_genotype(calls):
    # need to adjust diplotype to a single allele call for XY samples
    # NA and XX samples will be left as a diplotype
    #calls = calls.merge(imputed_sex, on="person_id", how="outer", validate="m:1")
    #xy_cond = calls["num_x"] == 1
    #TODO FIXME
    xy_cond = calls["person_id"] == "me"
    if sum(xy_cond) == 0:
        return calls #.drop(columns=["imputed_sex", "num_x"])
    num_diff = (
        calls.loc[xy_cond]
        .assign(allele=lambda x: x["genotype"].str.split("/"))
        .explode("allele")
        .groupby("person_id")
        .agg({"allele": "nunique"})
        .query("`allele` != 1")
        .shape[0]
    )
    # need to verify that only one unique allele exists per XY sample
    # Should always be the case since everything is through the DRAGEN pipeline
    # Known issue exists in DRAGEN 3.4.12 that reports hets for males
    # This has been fixed in the VCF before hand
    assert num_diff == 0
    calls.loc[xy_cond, "genotype"] = calls.loc[xy_cond, "genotype"].str.split(
        "/", expand=True
    )[0]
    calls = (
        calls.loc[xy_cond]
        .assign(allele=lambda x: x["possible_genotype"].str.split("/"))
        .explode("allele")
        .groupby("person_id", dropna=False, observed=True, as_index=False, sort=False)
        .agg(
            possible_genotype=(
                "allele",
                lambda x: pd.Series(x.unique()).str.cat(sep="/"),
            )
        )
        .replace("", np.nan)
        .merge(calls, on="person_id", how="outer")
    )
    calls.loc[xy_cond, "possible_genotype_y"] = np.nan # otherwise ffill doesn't work right
    calls.loc[:, "possible_genotype"] = (
        calls.loc[:, ["possible_genotype_x", "possible_genotype_y"]]
        .ffill(axis=1)["possible_genotype_y"]
        .drop(columns=["possible_genotype_y", "possible_genotype_x"])
        .rename({"possible_genotype_x": "possible_genotype"})
    )
    return calls.drop(columns=["imputed_sex", "num_x"])

def get_dpyd_activity_score(alleles):
    # In scenarios where a sample carries >2 alleles for DPYD,
    # CPIC guideline recommends using the 2 alleles with the lowest activity
    # to assign activity score for the gene
    # Therefore a genotype can contain >2 alleles
    # https://pubmed.ncbi.nlm.nih.gov/29152729/
    allele_fxns_sorted = sorted(alleles.map(dpyd_allele_function))
    return sum(allele_fxns_sorted[:2])

def assign_dpyd_cols(calls):
    # DPYD is a special case where samples have >2 alleles
    # See get_dpyd_activity_score for more detail
    dpyd_activity_score = (
        calls.assign(allele=lambda x: x["genotype"].str.split("/"))
        .explode("allele")
        .groupby("person_id")
        .agg(activity_score=("allele", get_dpyd_activity_score))
    )
    calls.loc[:, 'activity_score'] = "DPYD Activity Score " + dpyd_activity_score.reset_index()['activity_score'].astype(str)
    calls.loc[:, 'possible_genotype'] = np.nan
    calls.loc[:, 'possible_phenotype'] = np.nan
    
    return calls

def assign_dpyd_phenotype(calls):
    return calls.assign(
        phenotype=lambda x: x["activity_score"].map(
            {
                "DPYD Activity Score 2.0": "DPYD Normal Metabolizer",
                "DPYD Activity Score 1.5": "DPYD Intermediate Metabolizer",
                "DPYD Activity Score 1.0": "DPYD Intermediate Metabolizer",
                "DPYD Activity Score 0.5": "DPYD Poor Metabolizer",
                "DPYD Activity Score 0.0": "DPYD Poor Metabolizer",
            }
        )
    )


#this is intended for the single-sample case
def run_me(gene_file_dict):
    stargazer_table = StargazerTable(stargazer_converter_path)
    call_data = []
    for gene in gene_file_dict.keys():
        calls = pd.read_table(
            gene_file_dict[gene]
        ).rename(columns={"name": "person_id"})
        calls = cat_cols(calls)
        calls = assign_genotype(calls, gene, stargazer_table, attr="gene_data", col="genotype")
        calls = assign_genotype(
                calls, gene, stargazer_table, attr="gene_data_possible", col="possible_genotype"
            )
        gene_dip_to_phen = dip_to_phen[gene.upper()]
        calls = assign_phenotype(calls, gene, gene_dip_to_phen)
        if gene == 'dpyd':
            calls = assign_dpyd_cols(calls)
            calls = assign_dpyd_phenotype(calls)
        else:
            if gene == "g6pd":
                # This diplotype doesn't exist in CPIC
                gene_dip_to_phen.update({"Farroupilha/B (reference)": "Variable"})
                calls = adjust_g6pd_xy_genotype(calls)
                calls = adjust_g6pd_xy_genotype(calls)
            calls = assign_possible_phenotype(calls, gene, gene_dip_to_phen)
            calls = assign_activity(calls, gene)
        calls["gene"] = gene.upper()
        calls.loc[calls["phenotype"].isna(), "phenotype"] = f"{gene.upper()} Indeterminate"
        call_data.append(
            calls.loc[
                :,
                [
                    "person_id",
                    "gene",
                    "genotype",
                    "phenotype",
                    "activity_score",
                    "possible_genotype",
                    "possible_phenotype",
                ],
            ]
        )
        call_data_stargazer = pd.concat(call_data, axis=0, ignore_index=True)
        call_data_stargazer = adjust_dtypes(call_data_stargazer)
    stargazer_file = "call_data_stargazer.pkl"
    call_data_stargazer.to_pickle(stargazer_file)
    call_data_stargazer = pd.read_pickle(stargazer_file)
