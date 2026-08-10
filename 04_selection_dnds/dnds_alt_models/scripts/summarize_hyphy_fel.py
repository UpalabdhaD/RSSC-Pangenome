#!/usr/bin/env python3

import csv
import json
import math
import sys


gene, json_path, out_path = sys.argv[1:4]

with open(json_path) as handle:
    payload = json.load(handle)

headers = payload.get("MLE", {}).get("headers", [])
content = payload.get("MLE", {}).get("content", {})

header_index = {}
for idx, header in enumerate(headers):
    if isinstance(header, list) and header:
        header_index[str(header[0]).lower()] = idx
    elif isinstance(header, dict):
        name = header.get("name") or header.get("description") or header.get("label")
        if name:
            header_index[str(name).lower()] = idx

def find_col(*patterns):
    for key, idx in header_index.items():
        if all(pattern.lower() in key for pattern in patterns):
            return idx
    return None

alpha_col = find_col("alpha")
beta_col = find_col("beta")
p_col = find_col("p-value") or find_col("p value")

rows = []
for partition_rows in content.values():
    if isinstance(partition_rows, list):
        rows.extend(partition_rows)

tested_sites = len(rows)
positive = 0
negative = 0
omega_values = []

for row in rows:
    try:
        alpha = float(row[alpha_col]) if alpha_col is not None else math.nan
        beta = float(row[beta_col]) if beta_col is not None else math.nan
        pvalue = float(row[p_col]) if p_col is not None else math.nan
    except (TypeError, ValueError, IndexError):
        continue

    if math.isfinite(alpha) and alpha > 0 and math.isfinite(beta):
        omega_values.append(beta / alpha)
    if math.isfinite(pvalue) and pvalue <= 0.1 and math.isfinite(alpha) and math.isfinite(beta):
        if beta > alpha:
            positive += 1
        elif beta < alpha:
            negative += 1

mean_omega = sum(omega_values) / len(omega_values) if omega_values else "NA"

with open(out_path, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "gene",
            "status",
            "tested_sites",
            "mean_site_omega",
            "positive_sites_p_le_0.1",
            "negative_sites_p_le_0.1",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "gene": gene,
            "status": "ok",
            "tested_sites": tested_sites,
            "mean_site_omega": mean_omega,
            "positive_sites_p_le_0.1": positive,
            "negative_sites_p_le_0.1": negative,
        }
    )
