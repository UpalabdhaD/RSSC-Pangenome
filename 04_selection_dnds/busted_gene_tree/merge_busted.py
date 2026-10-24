"""Merge the original BUSTED results with the gene-tree rerun, and recompute q-values.

Inputs (read only):
  $PROJECT_ROOT/MS_REVISION/DNDS_HYPHY/results_all_busted/hyphy_busted_summary.tsv
  results/busted/<gene>/summary.tsv for the 607 families rerun on gene trees
Output:
  busted_complete.tsv: one row per eligible family (6,811) with tree_source, p-value, q-value
  (Benjamini-Hochberg across all families with a p-value), global omega, and, for the 607 rerun
  families, the original backup-tree p-value and omega.

herA: HyPhy wrote "AIC-c":inf, which is not valid JSON. We replaced inf with Infinity in a copy
(busted.inf_fixed.json) and ran the original summarize_hyphy_busted.py on that copy.

Run from this folder: python3 -I merge_busted.py
"""
import csv
import os
from pathlib import Path

ORIGINAL = Path(os.environ.get("PROJECT_ROOT", "."), "MS_REVISION/DNDS_HYPHY/"
                "results_all_busted/hyphy_busted_summary.tsv")
RERUN = Path("results/busted")
TARGET = Path("busted_complete.tsv")


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


with open(ORIGINAL, newline="") as handle:
    original = [{k.strip(): (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(handle, delimiter="\t")]
eligible = [r for r in original if r["input_status"] == "eligible"]
rerun_genes = [line.strip() for line in open("genes.txt") if line.strip()]
assert len(eligible) == 6811 and len(rerun_genes) == 607

rows = []
for r in eligible:
    row = {"gene": r["gene"], "n_sequences": r["n_sequences"], "tree_source": "pruned_species_tree",
           "busted_pvalue": r["busted_pvalue"] if r["busted_status"] == "ok" else "",
           "global_omega": r["global_omega_from_log"] if r["busted_status"] == "ok" else "",
           "pvalue_backup_tree": "", "omega_backup_tree": ""}
    if r["gene"] in rerun_genes:
        row["tree_source"] = "gene_tree"
        row["pvalue_backup_tree"], row["omega_backup_tree"] = row["busted_pvalue"], row["global_omega"]
        row["busted_pvalue"] = row["global_omega"] = ""
        summary = RERUN / r["gene"] / "summary.tsv"
        if summary.exists():
            with open(summary, newline="") as handle:
                s = next(csv.DictReader(handle, delimiter="\t"))
            if s["status"] == "ok" and number(s["busted_pvalue"]) is not None:
                row["busted_pvalue"], row["global_omega"] = s["busted_pvalue"], s["global_omega_from_log"]
    rows.append(row)

# Benjamini-Hochberg across all families with a p-value
tested = [(number(r["busted_pvalue"]), i) for i, r in enumerate(rows) if number(r["busted_pvalue"]) is not None]
tested.sort()
m = len(tested)
q = [None] * len(rows)
running = 1.0
for rank in range(m, 0, -1):
    p, i = tested[rank - 1]
    running = min(running, p * m / rank)
    q[i] = running
for r, value in zip(rows, q):
    r["busted_qvalue"] = "" if value is None else f"{value:.6g}"

columns = ["gene", "n_sequences", "tree_source", "busted_pvalue", "busted_qvalue", "global_omega",
           "pvalue_backup_tree", "omega_backup_tree"]
with open(TARGET, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

sig = [r for r in rows if r["busted_qvalue"] and float(r["busted_qvalue"]) <= 0.05]
omega_below_1 = [r for r in rows if number(r["global_omega"]) is not None and number(r["global_omega"]) < 1]
print(f"Wrote {TARGET}: {len(rows)} eligible families, {m} with a p-value")
print(f"  p <= 0.05: {sum(number(r['busted_pvalue']) <= 0.05 for r in rows if number(r['busted_pvalue']) is not None)}")
print(f"  q <= 0.05: {len(sig)}")
print(f"  global omega < 1: {len(omega_below_1)} of {m}")

# Backup tree against gene tree, for the rerun families with a p-value on both trees
both = [r for r in rows if r["tree_source"] == "gene_tree"
        and number(r["busted_pvalue"]) is not None and number(r["pvalue_backup_tree"]) is not None]
old_pos = sum(number(r["pvalue_backup_tree"]) <= 0.05 for r in both)
new_pos = sum(number(r["busted_pvalue"]) <= 0.05 for r in both)
lost = sum(number(r["pvalue_backup_tree"]) <= 0.05 < number(r["busted_pvalue"]) for r in both)
gained = sum(number(r["busted_pvalue"]) <= 0.05 < number(r["pvalue_backup_tree"]) for r in both)
print(f"Rerun families with both p-values: {len(both)}")
print(f"  p <= 0.05 on backup tree: {old_pos}; on gene tree: {new_pos}; lost: {lost}; gained: {gained}")
