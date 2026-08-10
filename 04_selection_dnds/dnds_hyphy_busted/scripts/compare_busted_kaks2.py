#!/usr/bin/env python3
"""Join all-family HyPhy BUSTED results to the existing KaKs2 estimates."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def parse_float(value: str | None) -> float | None:
    try:
        parsed = float(value) if value not in (None, "", "NA", "nan", "NaN") else math.nan
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def parse_bool(value: str | None) -> bool | None:
    if value == "True":
        return True
    if value == "False":
        return False
    return None


def main() -> None:
    busted_path, kaks_path, output_path = map(Path, sys.argv[1:4])
    summary_path = Path(sys.argv[4]) if len(sys.argv) == 5 else None

    kaks_by_gene: dict[str, float | None] = {}
    with kaks_path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            gene = row.get("gene", "")
            if not gene:
                continue
            if gene in kaks_by_gene:
                raise ValueError(f"Duplicate KaKs2 gene identifier: {gene}")
            kaks_by_gene[gene] = parse_float(row.get("kaks"))

    output_rows: list[dict[str, object]] = []
    with busted_path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            busted_pvalue = parse_float(row.get("busted_pvalue"))
            busted_positive = parse_bool(row.get("busted_positive_p_le_0.05"))
            kaks2 = kaks_by_gene.get(gene)
            kaks2_positive = kaks2 > 1 if kaks2 is not None else None

            if busted_positive is None or kaks2_positive is None:
                comparison = "not_comparable"
            elif busted_positive and kaks2_positive:
                comparison = "both_positive"
            elif busted_positive:
                comparison = "busted_positive_only"
            elif kaks2_positive:
                comparison = "kaks2_positive_only"
            else:
                comparison = "both_not_positive"

            output_rows.append(
                {
                    "gene": gene,
                    "n_sequences": row.get("n_sequences", ""),
                    "hyphy_input_status": row.get("input_status", ""),
                    "hyphy_run_status": row.get("run_status", ""),
                    "hyphy_busted_status": row.get("busted_status", ""),
                    "hyphy_busted_lrt": row.get("busted_lrt", ""),
                    "hyphy_busted_pvalue": busted_pvalue if busted_pvalue is not None else "",
                    "hyphy_busted_positive_p_le_0_05": "" if busted_positive is None else busted_positive,
                    "hyphy_global_omega": row.get("global_omega_from_log", ""),
                    "kaks2_dnds": "" if kaks2 is None else kaks2,
                    "kaks2_positive_gt_1": "" if kaks2_positive is None else kaks2_positive,
                    "busted_vs_kaks2_call": comparison,
                }
            )

    fieldnames = list(output_rows[0]) if output_rows else []
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    if summary_path is not None:
        counts: dict[str, int] = {}
        for row in output_rows:
            category = str(row["busted_vs_kaks2_call"])
            counts[category] = counts.get(category, 0) + 1
        with summary_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["busted_vs_kaks2_call", "families"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                {"busted_vs_kaks2_call": category, "families": count}
                for category, count in sorted(counts.items())
            )


if __name__ == "__main__":
    if len(sys.argv) not in {4, 5}:
        raise SystemExit(
            "usage: compare_busted_kaks2.py BUSTED_SUMMARY KAKS2_CSV OUTPUT_TSV [SUMMARY_TSV]"
        )
    main()
