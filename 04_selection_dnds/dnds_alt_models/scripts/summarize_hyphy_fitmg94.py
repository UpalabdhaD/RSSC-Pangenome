#!/usr/bin/env python3

import json
import math
import sys


def walk_numbers(obj, path=""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from walk_numbers(value, f"{path}.{key}" if path else str(key))
    elif isinstance(obj, list):
        for idx, value in enumerate(obj):
            yield from walk_numbers(value, f"{path}[{idx}]")
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool) and math.isfinite(obj):
        yield path, float(obj)


gene, json_path, summary_path = sys.argv[1:4]

with open(json_path) as handle:
    payload = json.load(handle)

omega_candidates = []
for path, value in walk_numbers(payload):
    lower = path.lower()
    if "omega" in lower or "dnds" in lower or "dn/ds" in lower:
        omega_candidates.append((path, value))

omega = omega_candidates[0][1] if omega_candidates else "NA"
source = omega_candidates[0][0] if omega_candidates else "not_found"

with open(summary_path, "w") as handle:
    handle.write("gene\tstatus\tomega\tomega_source\n")
    handle.write(f"{gene}\tok\t{omega}\t{source}\n")
