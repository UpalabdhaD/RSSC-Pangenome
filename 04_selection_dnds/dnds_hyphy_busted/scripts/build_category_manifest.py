#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
from pathlib import Path


source, alignment_dir, manifest_path, genes_path, min_sequences = sys.argv[1:6]
alignment_dir = Path(alignment_dir).resolve()
manifest_path = Path(manifest_path)
genes_path = Path(genes_path)
min_sequences = int(min_sequences)

rows = []
genes = []
seen = set()
with Path(source).open(newline="", encoding="utf-8", errors="replace") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("specific_class") != "Absent in large lineages" or row.get("method") != "NG":
            continue
        gene = row.get("gene", "")
        if not gene or gene in seen:
            continue
        seen.add(gene)
        alignment = alignment_dir / f"{gene}.codon.dedup.fa"
        n_sequences = sum(1 for line in alignment.open() if line.startswith(">")) if alignment.is_file() else 0
        status = "missing_alignment" if not alignment.is_file() else ("eligible" if n_sequences >= min_sequences else "too_few_sequences")
        rows.append({"gene": gene, "n_sequences": n_sequences, "alignment": str(alignment), "status": status})
        if status == "eligible":
            genes.append(gene)

manifest_path.parent.mkdir(parents=True, exist_ok=True)
with manifest_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["gene", "n_sequences", "alignment", "status"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
genes_path.write_text("\n".join(genes) + ("\n" if genes else ""), encoding="utf-8")
print(f"category_rows\t{len(rows)}")
print(f"eligible_min{min_sequences}\t{len(genes)}")
