#!/usr/bin/env python3
"""Compare prior KaKs_Calculator NG and PAML M0 values for one gene class."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[3]
CATEGORY_SOURCE = REPO_ROOT / "MS_REVISION/DNDS_SENSITIVITY/results/kaks_sensitivity_joined.tsv"
PAML_AGGREGATED = REPO_ROOT / "MS_REVISION/DNDS_PAML_VS_KAKS/paml_m0_aggregated.tsv"
PAML_AUDIT = REPO_ROOT / "MS_REVISION/DNDS_PAML_VS_KAKS/paml_completion_audit.tsv"
OUT_DIR = REPO_ROOT / "MS_REVISION/DNDS_PAML_VS_KAKS/absent_major_lineages"
LOG_PSEUDOCOUNT = 1e-6
CATEGORY = "Absent in large lineages"


def cohen_kappa(table: np.ndarray) -> float:
    total = table.sum()
    if total == 0:
        return math.nan
    observed = np.trace(table) / total
    expected = np.sum(table.sum(axis=0) * table.sum(axis=1)) / (total * total)
    return float((observed - expected) / (1 - expected)) if expected < 1 else math.nan


def correlation(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float]:
    if len(x) < 3:
        return math.nan, math.nan
    result = stats.pearsonr(x, y) if method == "pearson" else stats.spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    category_source = pd.read_csv(CATEGORY_SOURCE, sep="\t")
    category = category_source.loc[
        category_source["specific_class"].eq(CATEGORY) & category_source["method"].eq("NG")
    ].copy()
    if category["gene"].duplicated().any():
        duplicated = category.loc[category["gene"].duplicated(keep=False), "gene"].tolist()
        raise ValueError(f"Category source has duplicate gene rows: {duplicated[:10]}")

    category["kaks"] = pd.to_numeric(category["kaks"], errors="coerce")
    category["ng_valid"] = (
        category["valid_numeric"].astype(bool)
        & np.isfinite(category["kaks"])
        & category["kaks"].ge(0)
    )

    audit = pd.read_csv(PAML_AUDIT, sep="\t")
    audit_subset = audit[["gene", "run_status", "summary_check", "usable_paml"]].rename(
        columns={
            "run_status": "paml_run_status",
            "summary_check": "paml_summary_check",
            "usable_paml": "paml_usable",
        }
    )
    paml = pd.read_csv(PAML_AGGREGATED, sep="\t")
    paml_subset = paml[["gene", "omega", "kappa", "lnL", "np", "run_status", "sequence_count"]].rename(
        columns={"run_status": "paml_aggregated_status"}
    )

    coverage = category.merge(audit_subset, on="gene", how="left", validate="one_to_one")
    coverage = coverage.merge(paml_subset, on="gene", how="left", validate="one_to_one")
    coverage["paml_eligible"] = coverage["paml_run_status"].notna()
    coverage["paml_usable"] = coverage["paml_usable"].fillna(False).astype(bool)
    coverage["paml_status"] = coverage["paml_run_status"].fillna("not_eligible")
    coverage["paml_positive"] = coverage["omega"].gt(1)
    coverage["ng_positive"] = coverage["ng_valid"] & coverage["kaks"].gt(1)
    coverage.to_csv(OUT_DIR / "category_coverage.tsv", sep="\t", index=False)

    matched = coverage.loc[coverage["ng_valid"] & coverage["paml_usable"]].copy()
    matched["paml_positive"] = matched["omega"] > 1
    matched["ng_positive"] = matched["kaks"] > 1
    matched["selection_agreement"] = matched["paml_positive"] == matched["ng_positive"]
    matched["log10_paml"] = np.log10(matched["omega"] + LOG_PSEUDOCOUNT)
    matched["log10_ng"] = np.log10(matched["kaks"] + LOG_PSEUDOCOUNT)
    matched["log10_difference_m0_minus_ng"] = matched["log10_paml"] - matched["log10_ng"]
    matched["absolute_log10_difference"] = matched["log10_difference_m0_minus_ng"].abs()
    matched["omega_minus_kaks"] = matched["omega"] - matched["kaks"]
    matched = matched.sort_values("gene").reset_index(drop=True)
    matched.to_csv(OUT_DIR / "paml_m0_vs_ng_matched.tsv", sep="\t", index=False)

    tn = int((~matched["paml_positive"] & ~matched["ng_positive"]).sum())
    fn = int((~matched["paml_positive"] & matched["ng_positive"]).sum())
    fp = int((matched["paml_positive"] & ~matched["ng_positive"]).sum())
    tp = int((matched["paml_positive"] & matched["ng_positive"]).sum())
    confusion = pd.DataFrame(
        [
            {"m0_class": "omega<=1", "ng_class": "Ka/Ks<=1", "n": tn},
            {"m0_class": "omega<=1", "ng_class": "Ka/Ks>1", "n": fn},
            {"m0_class": "omega>1", "ng_class": "Ka/Ks<=1", "n": fp},
            {"m0_class": "omega>1", "ng_class": "Ka/Ks>1", "n": tp},
        ]
    )
    confusion.to_csv(OUT_DIR / "selection_confusion.tsv", sep="\t", index=False)

    pearson_raw, pearson_raw_p = correlation(matched["kaks"], matched["omega"], "pearson")
    spearman_raw, spearman_raw_p = correlation(matched["kaks"], matched["omega"], "spearman")
    metrics = {
        "category": CATEGORY,
        "category_source_rows": int(len(category)),
        "category_valid_ng_rows": int(category["ng_valid"].sum()),
        "category_invalid_ng_rows": int((~category["ng_valid"]).sum()),
        "category_singleton_rows": int(category["n_sequences"].eq(1).sum()),
        "category_multi_sequence_rows": int(category["n_sequences"].ge(2).sum()),
        "category_paml_eligible_rows": int(coverage["paml_eligible"].sum()),
        "category_paml_completed_rows": int(
            coverage["paml_run_status"].isin({"ok", "skipped_existing_ok"}).sum()
        ),
        "category_paml_usable_rows": int(coverage["paml_usable"].sum()),
        "category_paml_not_run_rows": int(coverage["paml_run_status"].eq("not_run").sum()),
        "category_paml_not_eligible_rows": int((~coverage["paml_eligible"]).sum()),
        "matched_rows": int(len(matched)),
        "matched_share_of_valid_ng": float(len(matched) / category["ng_valid"].sum()),
        "m0_omega_median_matched": float(matched["omega"].median()),
        "ng_kaks_median_matched": float(matched["kaks"].median()),
        "m0_positive_n": int(matched["paml_positive"].sum()),
        "m0_positive_fraction": float(matched["paml_positive"].mean()),
        "ng_positive_n": int(matched["ng_positive"].sum()),
        "ng_positive_fraction": float(matched["ng_positive"].mean()),
        "selection_agreement_n": int(matched["selection_agreement"].sum()),
        "selection_agreement_fraction": float(matched["selection_agreement"].mean()),
        "selection_cohen_kappa": cohen_kappa(np.array([[tn, fn], [fp, tp]])),
        "pearson_raw_r": pearson_raw,
        "pearson_raw_p": pearson_raw_p,
        "spearman_raw_rho": spearman_raw,
        "spearman_raw_p": spearman_raw_p,
        "median_absolute_raw_difference": float(matched["omega_minus_kaks"].abs().median()),
        "median_log10_difference_m0_minus_ng": float(matched["log10_difference_m0_minus_ng"].median()),
        "log10_pseudocount": LOG_PSEUDOCOUNT,
    }
    (OUT_DIR / "comparison_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"metric": key, "value": value} for key, value in metrics.items()]).to_csv(
        OUT_DIR / "comparison_metrics.tsv", sep="\t", index=False
    )
    matched.nlargest(len(matched), "absolute_log10_difference").to_csv(
        OUT_DIR / "largest_method_discrepancies.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
