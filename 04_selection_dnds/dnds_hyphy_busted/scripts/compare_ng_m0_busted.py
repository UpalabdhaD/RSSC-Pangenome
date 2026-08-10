#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = REPO_ROOT / "MS_REVISION/DNDS_HYPHY/results_absent_major_lineages_busted"
CATEGORY_SOURCE = REPO_ROOT / "MS_REVISION/DNDS_SENSITIVITY/results/kaks_sensitivity_joined.tsv"
PAML = REPO_ROOT / "MS_REVISION/DNDS_PAML_VS_KAKS/paml_m0_aggregated.tsv"
HYPHY = OUT_DIR / "hyphy_busted_summary.tsv"


def fraction(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def value_or_none(value):
    if pd.isna(value):
        return None
    return float(value)


def main() -> None:
    category = pd.read_csv(CATEGORY_SOURCE, sep="\t")
    category = category.loc[
        category["specific_class"].eq("Absent in large lineages") & category["method"].eq("NG")
    ].copy()
    category["kaks"] = pd.to_numeric(category["kaks"], errors="coerce")
    category["ng_valid"] = (
        category["valid_numeric"].astype(bool)
        & np.isfinite(category["kaks"])
        & category["kaks"].ge(0)
    )
    category["ng_positive"] = category["ng_valid"] & category["kaks"].gt(1)

    paml = pd.read_csv(PAML, sep="\t")[["gene", "omega", "run_status", "sequence_count"]].rename(
        columns={"omega": "m0_omega", "run_status": "m0_status", "sequence_count": "m0_sequence_count"}
    )
    hyphy = pd.read_csv(HYPHY, sep="\t")[[
        "gene",
        "n_sequences",
        "input_status",
        "run_status",
        "busted_lrt",
        "busted_pvalue",
        "busted_positive_p_le_0.05",
        "busted_positive_p_le_0.10",
        "global_omega_from_log",
    ]].rename(
        columns={
            "input_status": "hyphy_input_status",
            "run_status": "hyphy_run_status",
            "n_sequences": "hyphy_sequence_count",
            "busted_positive_p_le_0.05": "busted_positive_p_le_0_05",
            "busted_positive_p_le_0.10": "busted_positive_p_le_0_10",
        }
    )

    comparison = category.merge(paml, on="gene", how="left", validate="one_to_one")
    comparison = comparison.merge(hyphy, on="gene", how="left", validate="one_to_one")
    comparison["m0_usable"] = comparison["m0_omega"].notna()
    comparison["m0_positive"] = comparison["m0_usable"] & comparison["m0_omega"].gt(1)
    comparison["hyphy_success"] = comparison["hyphy_run_status"].isin(["ok", "skipped_existing_ok"])
    comparison["busted_positive_p_le_0_05"] = comparison["busted_positive_p_le_0_05"].fillna(False).astype(bool)
    comparison["busted_positive_p_le_0_10"] = comparison["busted_positive_p_le_0_10"].fillna(False).astype(bool)
    comparison["ng_m0_hyphy_matched"] = comparison["ng_valid"] & comparison["m0_usable"] & comparison["hyphy_success"]
    comparison.to_csv(OUT_DIR / "ng_m0_busted_comparison.tsv", sep="\t", index=False)

    ng_hyphy = comparison.loc[comparison["ng_valid"] & comparison["hyphy_success"]]
    m0_hyphy = comparison.loc[comparison["m0_usable"] & comparison["hyphy_success"]]
    all_three = comparison.loc[comparison["ng_m0_hyphy_matched"]].copy()
    metrics = {
        "category": "Absent in large lineages",
        "category_rows": int(len(comparison)),
        "category_singletons": int(comparison["n_sequences"].eq(1).sum()),
        "category_valid_ng_n": int(comparison["ng_valid"].sum()),
        "category_ng_positive_n": int(comparison["ng_positive"].sum()),
        "category_ng_positive_fraction": fraction(int(comparison["ng_positive"].sum()), int(comparison["ng_valid"].sum())),
        "hyphy_min4_eligible_n": int(comparison["hyphy_input_status"].eq("eligible").sum()),
        "hyphy_busted_success_n": int(comparison["hyphy_success"].sum()),
        "hyphy_busted_positive_p_le_0.05_n": int(comparison["busted_positive_p_le_0_05"].sum()),
        "hyphy_busted_positive_p_le_0.10_n": int(comparison["busted_positive_p_le_0_10"].sum()),
        "ng_hyphy_matched_n": int(len(ng_hyphy)),
        "ng_positive_fraction_in_hyphy_subset": fraction(int(ng_hyphy["ng_positive"].sum()), int(len(ng_hyphy))),
        "m0_hyphy_matched_n": int(len(m0_hyphy)),
        "m0_positive_fraction_in_hyphy_subset": fraction(int(m0_hyphy["m0_positive"].sum()), int(len(m0_hyphy))),
        "all_three_matched_n": int(len(all_three)),
        "all_three_ng_positive_n": int(all_three["ng_positive"].sum()),
        "all_three_m0_positive_n": int(all_three["m0_positive"].sum()),
        "all_three_busted_positive_p_le_0.05_n": int(all_three["busted_positive_p_le_0_05"].sum()),
    }
    (OUT_DIR / "ng_m0_busted_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    rows = all_three[[
        "gene",
        "n_sequences",
        "kaks",
        "m0_omega",
        "busted_pvalue",
        "global_omega_from_log",
        "ng_positive",
        "m0_positive",
        "busted_positive_p_le_0_05",
    ]].sort_values("gene")
    lines = [
        "# HyPhy BUSTED comparison for “Absent in large lineages”",
        "",
        "## Model and scope",
        "",
        "HyPhy BUSTED was run on the existing codon alignments with at least four deduplicated sequences. BUSTED is a gene-level branch-site test for episodic diversifying selection; its p-value is not another Ka/Ks estimate. A two-class synonymous/nonsynonymous configuration was used for stable overnight execution, with all branches tested and a 0.05 significance threshold reported. Singleton and two-sequence families were excluded because they do not provide a reliable phylogenetic test of episodic selection.",
        "",
        "## Coverage",
        "",
        f"The category contains {metrics['category_rows']} families. The original NG result has {metrics['category_valid_ng_n']} valid numeric values and {metrics['category_ng_positive_n']} values above one ({metrics['category_ng_positive_fraction']:.1%}). {metrics['category_singletons']} families are singletons. HyPhy was eligible for {metrics['hyphy_min4_eligible_n']} families and completed successfully for {metrics['hyphy_busted_success_n']}.",
        "",
        "| Population | n | Positive result |",
        "|---|---:|---:|",
        f"| Original valid NG values | {metrics['category_valid_ng_n']} | {metrics['category_ng_positive_n']} Ka/Ks > 1 ({metrics['category_ng_positive_fraction']:.1%}) |",
        f"| HyPhy BUSTED-eligible families | {metrics['hyphy_min4_eligible_n']} | {metrics['hyphy_busted_positive_p_le_0.05_n']} p ≤ 0.05; {metrics['hyphy_busted_positive_p_le_0.10_n']} p ≤ 0.10 |",
        f"| NG + HyPhy matched subset | {metrics['ng_hyphy_matched_n']} | {int(ng_hyphy['ng_positive'].sum())} NG-positive |",
        f"| M0 + HyPhy matched subset | {metrics['m0_hyphy_matched_n']} | {int(m0_hyphy['m0_positive'].sum())} M0-positive |",
        "",
        "## Direct three-way comparison",
        "",
        "Only families with valid NG, usable M0, and successful BUSTED results enter this table. These are not representative of all 106 valid NG families.",
        "",
        "| Gene | Sequences | NG Ka/Ks | PAML M0 omega | BUSTED p | BUSTED global omega (log) | NG > 1 | M0 > 1 | BUSTED p ≤ 0.05 |",
        "|---|---:|---:|---:|---:|---:|:---:|:---:|:---:|",
    ]
    for row in rows.itertuples(index=False):
        lines.append(
            f"| `{row.gene}` | {int(row.n_sequences)} | {row.kaks:.6f} | {row.m0_omega:.6f} | {row.busted_pvalue:.6g} | {row.global_omega_from_log:.4f} | {'yes' if row.ng_positive else 'no'} | {'yes' if row.m0_positive else 'no'} | {'yes' if row.busted_positive_p_le_0_05 else 'no'} |"
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        f"BUSTED detected episodic diversifying selection in {metrics['hyphy_busted_positive_p_le_0.05_n']} of the {metrics['hyphy_busted_success_n']} eligible families at p ≤ 0.05. The positive BUSTED call is not an NG-positive call: it is `fabF`, whose NG Ka/Ks is invalid, while its M0 omega is 0.04043. The remaining tested families were not BUSTED-positive. This provides no category-wide support for the original 84.9% NG-positive fraction; it is a small, sequence-rich subset test, not a validation of the singleton-heavy category.",
        "",
        "BUSTED p-values answer a different question from the NG and M0 ratios. The BUSTED global omega values are reported for context but should not be treated as interchangeable with PAML M0 or KaKs_Calculator NG. The strongest defensible conclusion is that one multi-sequence family shows evidence compatible with episodic selection under BUSTED, whereas the broad category-level NG claim remains dominated by singleton-derived NG outputs and is not independently confirmed.",
        "",
        "## Outputs",
        "",
        "- [Combined NG/M0/BUSTED table](ng_m0_busted_comparison.tsv)",
        "- [HyPhy summary](hyphy_busted_summary.tsv)",
        "- [Category manifest and coverage](category_manifest.tsv)",
        "- [Runner script](../scripts/run_absent_major_lineages_busted.sh)",
        "- [BUSTED model logs](logs/)",
    ])
    (OUT_DIR / "comparison_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
