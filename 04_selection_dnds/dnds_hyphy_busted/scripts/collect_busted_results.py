#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
from pathlib import Path


manifest_path, status_dir, result_dir, output_path = map(Path, sys.argv[1:5])
manifest = {}
with manifest_path.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        manifest[row["gene"]] = row

rows = []
for gene, row in manifest.items():
    summary_path = result_dir / gene / "summary.tsv"
    merged = dict(row)
    merged["input_status"] = merged.pop("status")
    merged["run_status"] = "not_run"
    status_path = status_dir / f"{gene}.status"
    if status_path.is_file():
        fields = status_path.read_text(encoding="utf-8", errors="replace").strip().split("\t")
        if len(fields) > 1:
            merged["run_status"] = fields[1]
    if summary_path.is_file():
        with summary_path.open(newline="", encoding="utf-8") as handle:
            summary_rows = list(csv.DictReader(handle, delimiter="\t"))
        if summary_rows:
            summary_row = dict(summary_rows[0])
            if "status" in summary_row:
                summary_row["busted_status"] = summary_row.pop("status")
            merged.update(summary_row)
    rows.append(merged)

fieldnames = ["gene", "n_sequences", "alignment", "input_status", "run_status"]
for row in rows:
    for key in row:
        if key not in fieldnames:
            fieldnames.append(key)

output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
