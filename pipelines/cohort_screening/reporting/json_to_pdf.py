# This script generates an HTML report from a JSON file containing genomic variant data.
import os
import pandas as pd
import logging
import pdfkit

sample_identifier = "TEST"


mapped_variants = [
    {
        "Gene": "BLK",
        "contig": "chr8",
        "position": "11556728",
        "Ref allele": "T",
        "Alt Allele": "C, <NON_REF>",
        "dbSNP": "rs2306234",
        "HGSVG": "NC_000008.11:g.11556728T>C",
        "Zygosity": "homozygous alternate",
        "consequence": "synonymous_variant",
        "Protein change": "ENST00000259089.8:c.843T>C(p.(Phe281=))",
        "gnomAD AF": "0.814359",
        "ClinVar classification": "association",
        "ClinVar phenotypes": "Systemic lupus erythematosus",
    },
    {
        "Gene": "CEL",
        "contig": "chr9",
        "position": "133071212",
        "Ref allele": "C",
        "Alt Allele": "T, <NON_REF>",
        "dbSNP": "rs488087",
        "Zygosity": "heterozygous",
        "consequence": "synonymous_variant",
        "Protein change": "ENST00000372080.6:c.1719C>T(p.(Pro573=))",
        "gnomAD AF": "0.254305",
        "ClinVar classification": "likely benign",
        "ClinVar phenotypes": "not specified",
    },
]

if not mapped_variants:
    logging.info('no variants to report, creating empty report')
    identified_variants_message = ""
    df_var_html = "<h3>No variants marked for further investigation</h3>"
else:
    identified_variants_message = "We have identified the following variants that warrant further investigation:"
    df_var = pd.DataFrame(mapped_variants)

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

# MOBY (Maturity-Onset Diabetes of the Young) genes
# obtained from Naylor 2018, https://www.ncbi.nlm.nih.gov/books/NBK500456/#mody-ov.Genetic_Causes_of_MODY
# FAQ:
# Why are HNF1A and KLF11 missing from the MODY list?
# While these two genes are included in the list of fourteen MODY genes,
# they are also included in the list of Challenging Medically Relevant Genes (CMRG) (bed file).
# The CMRG genes are difficult to call in hg38, which is the reference genome for the proposed reports in this document.
# In order to call variants in these genes effectively, we would need to implement and run another annotation pipeline,
# which is out of scope for this project.
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
html_content += """
        <h3 style="font-size: 16px;">We have examined the following genes*, which have been implicated in 
        Maturity-Onset Diabetes of the Young (MODY):</h3>
        <ul>
    """

for gene in mody_genes:
    html_content += f"<li>{gene}</li>"

html_content += """
        </ul>
        <p style="font-size: 14px;">The analysis to identify these variants has limitations, which include potentially missing pathogenic 
        variants, missing variants that are not yet associated with MODY, and missing variants that are associated 
        with MODY as part of a polygenic effect.  This analysis did not use any familial genomic information 
        to confirm variants.</p>
        <p style="font-size: 14px;">Please note that we do not include HNF1A and KLF11, two MODY genes, in this report, 
        due to difficulties calling these genes on the hg38 reference 
        (<a href="http://dx.doi.org/10.1038/s41587-021-01158-1">10.1038/s41587-021-01158-1</a>).</p>
    </body>
    </html>
    """
# write out the HTML content to a file
with open('report.html', 'w') as html_file:
    html_file.write(html_content)
logging.info("HTML report generated successfully!")

# convert the HTML report to a PDF file
pdf_report_name = sample_identifier + '_mody_variants_report.pdf'
pdfkit.from_file('report.html', pdf_report_name)
logging.info("PDF report generated successfully!")

# Export variants dataframe as tsv
if not mapped_variants:
    pass
else:
    table_name = sample_identifier + '_mody_variants_table.tsv'
    df_var.to_csv(table_name, sep='\t', index=False)
    logging.info("variants.tsv exported successfully!")