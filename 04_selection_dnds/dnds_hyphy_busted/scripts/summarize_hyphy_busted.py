#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path


gene, json_path, log_path, summary_path = sys.argv[1:5]
json_path = Path(json_path)
log_path = Path(log_path)
summary_path = Path(summary_path)


def walk_numbers(value, path=""):
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            yield from walk_numbers(child, child_path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from walk_numbers(child, f"{path}[{index}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        yield path, float(value)


def find_value(numbers, required_terms, preferred_terms=()):
    candidates = []
    for path, value in numbers:
        lower = path.lower()
        if all(term in lower for term in required_terms):
            preference = sum(term in lower for term in preferred_terms)
            candidates.append((preference, path, value))
    if not candidates:
        return "NA", "not_found"
    _, path, value = sorted(candidates, key=lambda item: (-item[0], item[1]))[0]
    return value, path


with json_path.open(encoding="utf-8") as handle:
    payload = json.load(handle)

numbers = list(walk_numbers(payload))
lrt, lrt_source = find_value(numbers, ("lrt",), ("test", "results"))
pvalue, pvalue_source = find_value(numbers, ("p-value",), ("test", "results"))
if pvalue == "NA":
    pvalue, pvalue_source = find_value(numbers, ("pvalue",), ("test", "results"))

log_text = log_path.read_text(encoding="utf-8", errors="replace")
omega_match = re.findall(r"rate ratio for \*test\*\s*=\s*([0-9.eE+-]+)", log_text)
global_omega = float(omega_match[-1]) if omega_match else "NA"

positive_005 = isinstance(pvalue, (int, float)) and pvalue <= 0.05
positive_010 = isinstance(pvalue, (int, float)) and pvalue <= 0.10

summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "gene",
            "status",
            "busted_lrt",
            "busted_pvalue",
            "busted_positive_p_le_0.05",
            "busted_positive_p_le_0.10",
            "global_omega_from_log",
            "lrt_source",
            "pvalue_source",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "gene": gene,
            "status": "ok",
            "busted_lrt": lrt,
            "busted_pvalue": pvalue,
            "busted_positive_p_le_0.05": positive_005,
            "busted_positive_p_le_0.10": positive_010,
            "global_omega_from_log": global_omega,
            "lrt_source": lrt_source,
            "pvalue_source": pvalue_source,
        }
    )
