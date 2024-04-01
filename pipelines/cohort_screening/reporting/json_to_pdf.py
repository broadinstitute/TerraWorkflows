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
    # Add more variants as needed
]

genes = ["ABCC8", "APPL1", "BLK", "CEL", "GCK", "HNF1B", "HNF4A", "INS", "KCNJ11", "NEUROD1", "PAX4", "PDX1"]

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
    <h2>We have identified the following variants that warrant further investigation:</h2>
    <table>
        <tr>
            <th>Gene</th>
            <th>contig</th>
            <th>position</th>
            <th>Ref allele</th>
            <th>Alt Allele</th>
            <th>dbSNP</th>
            <th>Zygosity</th>
            <th>consequence</th>
            <th>Protein change</th>
            <th>gnomAD AF</th>
            <th>ClinVar classification</th>
            <th>ClinVar phenotypes</th>
        </tr>
"""

for variant in variants:
    html_content += "<tr>"
    for key, value in variant.items():
        html_content += f"<td>{value}</td>"
    html_content += "</tr>"

html_content += """
    </table>
    <h2>We have examined the following genes, which have been implicated in Maturity-Onset Diabetes of the Young (MODY):</h2>
    <ul>
"""

for gene in genes:
    html_content += f"<li>{gene}</li>"

html_content += """
    </ul>
    <p>The analysis to identify these variants has limitations, which include potentially missing pathogenic variants, missing variants that are not yet associated with MODY, and missing variants that are associated with MODY as part of a polygenic effect.  This analysis did not use any familial genomic information to confirm variants.</p>
</body>
</html>
"""

with open('report.html', 'w') as html_file:
    html_file.write(html_content)
