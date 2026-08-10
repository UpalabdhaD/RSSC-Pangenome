#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
from pathlib import Path


alignment_dir, manifest_path, genes_path, min_sequences = sys.argv[1:5]
alignment_dir = Path(alignment_dir).resolve()
manifest_path = Path(manifest_path)
genes_path = Path(genes_path)
min_sequences = int(min_sequences)

rows = []
eligible = []
for alignment in sorted(alignment_dir.glob("*.codon.dedup.fa")):
    gene = alignment.name.removesuffix(".codon.dedup.fa")
    with alignment.open(encoding="utf-8", errors="replace") as handle:
        n_sequences = sum(1 for line in handle if line.startswith(">"))
    status = "eligible" if n_sequences >= min_sequences else "too_few_sequences"
    rows.append(
        {
            "gene": gene,
            "n_sequences": n_sequences,
            "alignment": str(alignment),
            "input_status": status,
        }
    )
    if status == "eligible":
        eligible.append((n_sequences, gene))

manifest_path.parent.mkdir(parents=True, exist_ok=True)
with manifest_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["gene", "n_sequences", "alignment", "input_status"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

eligible.sort(key=lambda item: (item[0], item[1]))
genes_path.write_text(
    "\n".join(gene for _, gene in eligible) + ("\n" if eligible else ""),
    encoding="utf-8",
)
print(f"all_family_alignments\t{len(rows)}")
print(f"eligible_min{min_sequences}\t{len(eligible)}")
