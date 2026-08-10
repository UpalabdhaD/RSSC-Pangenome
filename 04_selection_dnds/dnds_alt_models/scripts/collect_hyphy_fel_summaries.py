#!/usr/bin/env python3

import csv
import pathlib
import sys


indir = pathlib.Path(sys.argv[1])
outfile = pathlib.Path(sys.argv[2])

files = sorted(indir.glob("*/summary.tsv"))
if not files:
    raise SystemExit(f"No HyPhy summary.tsv files found under {indir}")

rows = []
fieldnames = None
for path in files:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if fieldnames is None:
            fieldnames = reader.fieldnames
        rows.extend(reader)

outfile.parent.mkdir(parents=True, exist_ok=True)
with outfile.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
