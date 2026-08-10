#!/usr/bin/env python3
"""Audit PAML M0 completion and compare omega with prior KaKs_Calculator NG values."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[3]
PAML_ROOT = REPO_ROOT / "MS_REVISION/DNDS_ALT_MODELS/results_all_paml_m0"
OLD_KAKS_CSV = REPO_ROOT / "Results/DNDS/RSSC_kaks2_df.csv"
OUT_DIR = REPO_ROOT / "MS_REVISION/DNDS_PAML_VS_KAKS"
PLOTS_DIR = OUT_DIR / "plots"
COMPLETED_STATUSES = {"ok", "skipped_existing_ok"}
LOG_PSEUDOCOUNT = 1e-6


def read_statuses(status_dir: Path) -> dict[str, str]:
    statuses: dict[str, str] = {}
    for filename in sorted(os.listdir(status_dir)):
        if not filename.endswith(".status"):
            continue
        path = status_dir / filename
        fields = path.read_text(encoding="utf-8", errors="replace").strip().split("\t")
        gene = filename.removesuffix(".status")
        statuses[gene] = fields[1] if len(fields) >= 2 else "malformed_status"
    return statuses


def read_one_summary(path: Path, expected_gene: str) -> tuple[dict[str, object] | None, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return None, "missing_summary"
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return None, "unreadable_summary"
    if len(rows) != 1:
        return None, f"summary_row_count_{len(rows)}"
    row = rows[0]
    if row.get("gene") != expected_gene:
        return None, "summary_gene_mismatch"
    if row.get("status") != "ok":
        return None, f"summary_status_{row.get('status', 'missing')}"
    try:
        parsed = {
            "gene": expected_gene,
            "omega": float(row["omega"]),
            "kappa": float(row["kappa"]),
            "lnL": float(row["lnL"]),
            "np": int(float(row["np"])),
        }
    except (KeyError, TypeError, ValueError):
        return None, "nonnumeric_summary"
    numeric_values = [parsed["omega"], parsed["kappa"], parsed["lnL"]]
    if not all(math.isfinite(float(value)) for value in numeric_values):
        return None, "nonfinite_summary"
    if float(parsed["omega"]) < 0 or float(parsed["kappa"]) <= 0 or int(parsed["np"]) <= 0:
        return None, "out_of_domain_summary"
    return parsed, "valid"


def safe_correlation(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if int(valid.sum()) < 3:
        return math.nan, math.nan
    if method == "pearson":
        result = stats.pearsonr(x[valid], y[valid])
    else:
        result = stats.spearmanr(x[valid], y[valid])
    return float(result.statistic), float(result.pvalue)


def cohen_kappa(table: np.ndarray) -> float:
    total = table.sum()
    if total == 0:
        return math.nan
    observed = np.trace(table) / total
    expected = np.sum(table.sum(axis=0) * table.sum(axis=1)) / (total * total)
    return float((observed - expected) / (1 - expected)) if expected < 1 else math.nan


def metric_rows(metrics: dict[str, object]) -> list[dict[str, object]]:
    return [{"metric": key, "value": value} for key, value in metrics.items()]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(PAML_ROOT / "genes_all_min2.sequence_counts.tsv", sep="\t")
    eligible_names = (PAML_ROOT / "genes_all_min2.txt").read_text(encoding="utf-8").split()
    eligible = pd.DataFrame({"gene": eligible_names})
    eligible = eligible.merge(counts, on="gene", how="left", validate="one_to_one")

    statuses = read_statuses(PAML_ROOT / "status")
    eligible["run_status"] = eligible["gene"].map(statuses).fillna("not_run")

    paml_records: list[dict[str, object]] = []
    summary_checks: dict[str, str] = {}
    for gene, run_status in eligible[["gene", "run_status"]].itertuples(index=False):
        if run_status not in COMPLETED_STATUSES:
            summary_checks[gene] = "not_completed"
            continue
        summary_path = PAML_ROOT / "paml_m0" / gene / "summary.tsv"
        record, check = read_one_summary(summary_path, gene)
        summary_checks[gene] = check
        if record is not None:
            record["run_status"] = run_status
            paml_records.append(record)

    eligible["summary_check"] = eligible["gene"].map(summary_checks).fillna("not_checked")
    eligible["usable_paml"] = eligible["summary_check"].eq("valid")
    eligible.to_csv(OUT_DIR / "paml_completion_audit.tsv", sep="\t", index=False)

    paml = pd.DataFrame(paml_records)
    paml = paml.merge(
        counts.rename(columns={"deduplicated_sequence_count": "sequence_count"}),
        on="gene",
        how="left",
        validate="one_to_one",
    )
    paml = paml.sort_values("gene").reset_index(drop=True)
    paml.to_csv(OUT_DIR / "paml_m0_aggregated.tsv", sep="\t", index=False)

    old = pd.read_csv(OLD_KAKS_CSV)
    old["kaks"] = pd.to_numeric(old["kaks"], errors="coerce")
    old_duplicate_genes = int(old["gene"].duplicated(keep=False).sum())
    old_valid_mask = np.isfinite(old["kaks"]) & old["kaks"].ge(0) & old["gene"].notna()
    old_valid = old.loc[old_valid_mask, ["gene", "kaks"]].drop_duplicates("gene", keep="first")

    matched = paml.merge(old_valid, on="gene", how="inner", validate="one_to_one")
    matched["paml_positive"] = matched["omega"] > 1
    matched["kaks_positive"] = matched["kaks"] > 1
    matched["selection_agreement"] = matched["paml_positive"] == matched["kaks_positive"]
    matched["log10_paml"] = np.log10(matched["omega"] + LOG_PSEUDOCOUNT)
    matched["log10_kaks"] = np.log10(matched["kaks"] + LOG_PSEUDOCOUNT)
    matched["log10_difference"] = matched["log10_paml"] - matched["log10_kaks"]
    matched["absolute_log10_difference"] = matched["log10_difference"].abs()
    matched["omega_minus_kaks"] = matched["omega"] - matched["kaks"]
    matched = matched.sort_values("gene").reset_index(drop=True)
    matched.to_csv(OUT_DIR / "paml_vs_kaks_matched.tsv", sep="\t", index=False)

    tn = int((~matched["paml_positive"] & ~matched["kaks_positive"]).sum())
    fp = int((matched["paml_positive"] & ~matched["kaks_positive"]).sum())
    fn = int((~matched["paml_positive"] & matched["kaks_positive"]).sum())
    tp = int((matched["paml_positive"] & matched["kaks_positive"]).sum())
    confusion_array = np.array([[tn, fn], [fp, tp]])
    confusion = pd.DataFrame(
        [
            {"paml_class": "omega<=1", "kaks_class": "Ka/Ks<=1", "n": tn},
            {"paml_class": "omega<=1", "kaks_class": "Ka/Ks>1", "n": fn},
            {"paml_class": "omega>1", "kaks_class": "Ka/Ks<=1", "n": fp},
            {"paml_class": "omega>1", "kaks_class": "Ka/Ks>1", "n": tp},
        ]
    )
    confusion.to_csv(OUT_DIR / "selection_confusion.tsv", sep="\t", index=False)

    pearson_raw, pearson_raw_p = safe_correlation(matched["kaks"], matched["omega"], "pearson")
    spearman_raw, spearman_raw_p = safe_correlation(matched["kaks"], matched["omega"], "spearman")
    positive_both = matched[(matched["kaks"] > 0) & (matched["omega"] > 0)]
    pearson_log, pearson_log_p = safe_correlation(
        np.log10(positive_both["kaks"]), np.log10(positive_both["omega"]), "pearson"
    )

    status_counts = eligible["run_status"].value_counts().to_dict()
    summary_check_counts = eligible["summary_check"].value_counts().to_dict()
    completed_status_count = int(eligible["run_status"].isin(COMPLETED_STATUSES).sum())
    usable_paml_count = int(eligible["usable_paml"].sum())
    metrics: dict[str, object] = {
        "paml_eligible_gene_families": int(len(eligible)),
        "paml_completed_statuses": completed_status_count,
        "paml_usable_summaries": usable_paml_count,
        "paml_completion_rate": completed_status_count / len(eligible),
        "paml_prepare_failed": int(status_counts.get("prepare_failed", 0)),
        "paml_not_run": int(status_counts.get("not_run", 0)),
        "paml_completed_but_invalid_summary": completed_status_count - usable_paml_count,
        "old_kaks_rows": int(len(old)),
        "old_kaks_duplicate_gene_rows": old_duplicate_genes,
        "old_kaks_invalid_rows": int((~old_valid_mask).sum()),
        "old_kaks_valid_unique_genes": int(len(old_valid)),
        "matched_gene_families": int(len(matched)),
        "matched_share_of_usable_paml": len(matched) / usable_paml_count if usable_paml_count else math.nan,
        "matched_share_of_valid_kaks": len(matched) / len(old_valid) if len(old_valid) else math.nan,
        "paml_omega_median": float(matched["omega"].median()),
        "old_kaks_median": float(matched["kaks"].median()),
        "paml_positive_n": int(matched["paml_positive"].sum()),
        "paml_positive_fraction": float(matched["paml_positive"].mean()),
        "paml_omega_at_lower_boundary_n": int((matched["omega"] <= 0.0001).sum()),
        "paml_omega_at_upper_boundary_n": int((matched["omega"] >= 999).sum()),
        "old_kaks_positive_n": int(matched["kaks_positive"].sum()),
        "old_kaks_positive_fraction": float(matched["kaks_positive"].mean()),
        "old_kaks_at_zero_or_below_n": int((matched["kaks"] <= 0).sum()),
        "old_kaks_at_upper_boundary_n": int((matched["kaks"] >= 999).sum()),
        "selection_class_agreement_n": int(matched["selection_agreement"].sum()),
        "selection_class_agreement_fraction": float(matched["selection_agreement"].mean()),
        "selection_cohen_kappa": cohen_kappa(confusion_array),
        "pearson_raw_r": pearson_raw,
        "pearson_raw_p": pearson_raw_p,
        "spearman_raw_rho": spearman_raw,
        "spearman_raw_p": spearman_raw_p,
        "pearson_log10_positive_pairs_r": pearson_log,
        "pearson_log10_positive_pairs_p": pearson_log_p,
        "median_absolute_raw_difference": float(matched["omega_minus_kaks"].abs().median()),
        "median_log10_difference_paml_minus_kaks": float(matched["log10_difference"].median()),
        "log10_pseudocount_for_plots": LOG_PSEUDOCOUNT,
    }

    pd.DataFrame(metric_rows(metrics)).to_csv(OUT_DIR / "comparison_metrics.tsv", sep="\t", index=False)
    matched.nlargest(100, "absolute_log10_difference").to_csv(
        OUT_DIR / "largest_method_discrepancies.tsv", sep="\t", index=False
    )

    rng = np.random.default_rng(20260802)
    if len(matched) > 2500:
        discordant = matched.loc[~matched["selection_agreement"]]
        retained = discordant if len(discordant) <= 1000 else discordant.sample(1000, random_state=20260802)
        remaining = matched.drop(index=retained.index)
        sample_n = max(0, 2500 - len(retained))
        sampled = remaining.iloc[rng.choice(len(remaining), size=sample_n, replace=False)]
        chart_sample = pd.concat([retained, sampled], ignore_index=True)
    else:
        chart_sample = matched.copy()
    chart_sample = chart_sample.sort_values("gene")
    chart_sample.to_csv(OUT_DIR / "chart_scatter_sample.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    hexbin = axes[0].hexbin(
        matched["log10_kaks"], matched["log10_paml"], gridsize=55, mincnt=1, cmap="viridis", bins="log"
    )
    limits = [
        float(min(matched["log10_kaks"].quantile(0.005), matched["log10_paml"].quantile(0.005))),
        float(max(matched["log10_kaks"].quantile(0.995), matched["log10_paml"].quantile(0.995))),
    ]
    axes[0].plot(limits, limits, linestyle="--", color="black", linewidth=1, label="1:1")
    axes[0].set_xlim(limits)
    axes[0].set_ylim(limits)
    axes[0].set_xlabel("KaKs_Calculator NG log10(Ka/Ks + 1e-6)")
    axes[0].set_ylabel("PAML M0 log10(omega + 1e-6)")
    axes[0].set_title("Gene-family dN/dS estimates")
    axes[0].legend(frameon=False)
    fig.colorbar(hexbin, ax=axes[0], label="log10 gene-family count")

    matrix = np.array([[tn, fn], [fp, tp]])
    image = axes[1].imshow(matrix, cmap="Blues")
    for i in range(2):
        for j in range(2):
            axes[1].text(j, i, f"{matrix[i, j]:,}", ha="center", va="center", color="black", fontsize=12)
    axes[1].set_xticks([0, 1], ["Ka/Ks <= 1", "Ka/Ks > 1"])
    axes[1].set_yticks([0, 1], ["PAML <= 1", "PAML > 1"])
    axes[1].set_title("Selection-threshold agreement")
    fig.colorbar(image, ax=axes[1], label="gene-family count")
    fig.suptitle("PAML M0 versus prior KaKs_Calculator NG results", fontsize=15)
    fig.savefig(PLOTS_DIR / "paml_vs_kaks_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(PLOTS_DIR / "paml_vs_kaks_comparison.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "metrics": metrics,
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "summary_check_counts": {str(k): int(v) for k, v in summary_check_counts.items()},
        "inputs": {
            "paml_root": str(PAML_ROOT.relative_to(REPO_ROOT)),
            "old_kaks_csv": str(OLD_KAKS_CSV.relative_to(REPO_ROOT)),
        },
        "method_note": (
            "PAML codeml M0 is a phylogeny-aware one-ratio family estimate. The prior KaKs_Calculator "
            "result is the NG Ka/Ks value from the first row of each family TSV, generated from a multi-sequence AXT record."
        ),
    }
    (OUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
