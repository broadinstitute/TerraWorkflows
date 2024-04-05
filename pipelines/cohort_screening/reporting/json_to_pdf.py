
# This script generates an HTML report from a JSON file containing genomic variant data.
import pandas as pd

# dummy data
sample_identifier = "Sample 1"

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
    # Add more variants as needed
]
# Convert the list of variant dictionaries to a DataFrame
# TODO: create a for loop to iterate over the variants and *append* them to the DataFrame
# df_var = df_var.append(data, ignore_index=True)
df_var = pd.DataFrame(variants)

# Apply styles to the DataFrame
borders = [{'selector': 'td, th, table',
            'props': [('border', '1px solid lightgrey'),
                      ('border-collapse', 'collapse')]
            }]
# Convert the DataFrame to an HTML table
df_var_html = df_var.style.hide().set_table_styles(borders).to_html()

# MOBY (Maturity-Onset Diabetes of the Young) genes
# obtained from Naylor 2018, https://www.ncbi.nlm.nih.gov/books/NBK500456/#mody-ov.Genetic_Causes_of_MODY
# TODO: Add FAQ#2 - unclear if this is part of the report or just a comment
# see #m42-pipelines slack
mody_genes = ["ABCC8", "APPL1", "BLK", "CEL", "GCK", "HNF1B", "HNF4A", "INS", "KCNJ11", "NEUROD1", "PAX4", "PDX1"]

# Generate HTML report
html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Genomic Variant Report</title>
</head>
<body>
    <h1>Genomic Variant Report</h1>
    <p>Sample identifier: {sample_identifier}</p>
    <p>Please note:<br>
    <ul>
        <li>This report is for research purposes only.</li>
        <li>This report is not meant for direct disclosure to patients.</li>
        <li>This report does not mean that the person with these variants has diabetes.</li>
        <li>Clinicians should confirm this result before using it as part of clinical care.</li>
    </ul>
    <p>We have identified the following variants that warrant further investigation:</p>
    {df_var_html}
"""
html_content += """
    <h2>We have examined the following genes, which have been implicated in Maturity-Onset Diabetes of the Young (MODY):</h2>
    <ul>
"""

for gene in mody_genes:
    html_content += f"<li>{gene}</li>"

html_content += """
    </ul>
    <p>The analysis to identify these variants has limitations, which include potentially missing pathogenic variants, missing variants that are not yet associated with MODY, and missing variants that are associated with MODY as part of a polygenic effect.  This analysis did not use any familial genomic information to confirm variants.</p>
</body>
</html>
"""

# write out the HTML content to a file
with open('report.html', 'w') as html_file:
    html_file.write(html_content)
print("HTML report generated successfully!")

# TODO: Add code to convert the HTML report to a PDF file


# Export variants dataframe as tsv
df_var.to_csv('variants.tsv', sep='\t', index=False)
print("Variants exported successfully!")
