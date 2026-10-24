"""Build Supplementary Table 6 (genome accessions) from the fetched accession list.

Input:  RSSC_142_strain_accessions.tsv (copied unchanged from
        the analysis project folder, MS_REVISION/).
Output: Table_S6_RSSC_genome_accessions.tsv

The script keeps the row order of the input, converts Windows line endings,
and replaces the F1C1 entry "Pending" with the ENA study accession.
Run from this folder: python3 -I make_supp_table6.py
"""
import csv
import re

SOURCE = "RSSC_142_strain_accessions.tsv"
TARGET = "Table_S6_RSSC_genome_accessions.tsv"
F1C1_LABEL = "This study (ENA: PRJEB107717)"

with open(SOURCE, newline="") as handle:
    rows = [(r["Strain"].strip(), r["Accession_Number"].strip())
            for r in csv.DictReader(handle, delimiter="\t")]

genbank = [acc for strain, acc in rows if strain != "F1C1"]
assert len(rows) == 142, len(rows)
assert len(genbank) == 141, len(genbank)
assert all(re.fullmatch(r"GCA_\d{9}\.\d+", acc) for acc in genbank)
assert len(set(genbank)) == 141, "duplicate accession"
assert len({strain for strain, _ in rows}) == 142, "duplicate strain"
assert dict(rows)["F1C1"] == "Pending"

with open(TARGET, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["Strain", "GenBank assembly accession"])
    for strain, acc in rows:
        writer.writerow([strain, F1C1_LABEL if strain == "F1C1" else acc])

print(f"Wrote {TARGET}: {len(rows)} genomes ({len(genbank)} GenBank accessions + F1C1)")
