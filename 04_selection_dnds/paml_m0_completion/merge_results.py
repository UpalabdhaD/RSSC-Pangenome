"""Merge the original PAML M0 table with the final run into one table of all eligible families.

Inputs (read only):
  $PROJECT_ROOT/MS_REVISION/DNDS_PAML_VS_KAKS/paml_m0_aggregated.tsv
  $PROJECT_ROOT/MS_REVISION/DNDS_PAML_VS_KAKS/paml_completion_audit.tsv
  queue.tsv and results_final/ from run_final.sh
Output:
  paml_m0_complete.tsv: the columns of paml_m0_aggregated.tsv plus
    tree_source        pruned_species_tree or gene_tree
    run_batch          original or final_2026-10-14
    omega_backup_tree  the original omega on the backup ladder tree (255 families only)

Rules:
  - Original families fitted on the pruned species tree keep their original values.
  - The 255 original families fitted on the backup tree take their value from the final run
    (gene tree); the original value moves to omega_backup_tree.
  - The 437 families the original run did not finish take their value from the final run.

Run from this folder after run_final.sh finishes:
  PROJECT_ROOT=/path/to/analysis/project python3 -I merge_results.py
"""
import csv
import os
from pathlib import Path

SOURCE = Path(os.environ.get("PROJECT_ROOT", ".")) / "MS_REVISION/DNDS_PAML_VS_KAKS"
RESULTS = Path("results_final")
TARGET = Path("paml_m0_complete.tsv")
COLUMNS = ["gene", "omega", "kappa", "lnL", "np", "run_status", "sequence_count",
           "tree_source", "run_batch", "omega_backup_tree"]
COMPLETED = {"ok", "skipped_existing_ok"}
FINAL_BATCH = "final_2026-10-14"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


original = {r["gene"]: r for r in read_tsv(SOURCE / "paml_m0_aggregated.tsv")}
audit = {r["gene"]: r for r in read_tsv(SOURCE / "paml_completion_audit.tsv")}
queue = {}
with open("queue.tsv") as handle:
    for line in handle:
        n, gene, tree, origin = line.rstrip("\n").split("\t")
        # Five families (accC, aNKYR, purB, sbnB, thiF) appear twice in queue.tsv: the original
        # run prepared them on the backup tree but never finished them. They count as remaining.
        if gene not in queue or origin == "remaining":
            queue[gene] = {"tree": tree, "origin": origin}

assert len(original) == 9187 and len(audit) == 9624 and len(queue) == 692
redo = {g for g, q in queue.items() if q["origin"] == "original_backup_tree"}
remaining = {g for g, q in queue.items() if q["origin"] == "remaining"}
assert len(redo) == 255 and redo <= set(original)
assert len(remaining) == 437 and not remaining & set(original)


def final_row(gene):
    status_file = RESULTS / "status" / f"{gene}.status"
    status = status_file.read_text().split("\t")[1] if status_file.exists() else "not_run"
    row = {"gene": gene, "omega": "NA", "kappa": "NA", "lnL": "NA", "np": "NA",
           "run_status": status, "sequence_count": audit[gene]["deduplicated_sequence_count"],
           "tree_source": "gene_tree" if queue[gene]["tree"] == "gene_tree" else "pruned_species_tree",
           "run_batch": FINAL_BATCH, "omega_backup_tree": ""}
    if status in COMPLETED:
        summary = read_tsv(RESULTS / "paml_m0" / gene / "summary.tsv")[0]
        assert summary["gene"] == gene and summary["status"] == "ok", gene
        row.update({k: summary[k] for k in ("omega", "kappa", "lnL", "np")})
    return row


rows = []
for gene, r in original.items():
    if gene in redo:
        row = final_row(gene)
        row["omega_backup_tree"] = r["omega"]
    else:
        row = dict(r, tree_source="pruned_species_tree", run_batch="original", omega_backup_tree="")
    rows.append(row)
rows += [final_row(gene) for gene in sorted(remaining)]

assert len(rows) == 9624 and len({r["gene"] for r in rows}) == 9624
with open(TARGET, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

counts = {}
for r in rows:
    key = (r["run_status"], r["tree_source"])
    counts[key] = counts.get(key, 0) + 1
print(f"Wrote {TARGET}: {len(rows)} eligible families")
for (status, tree), n in sorted(counts.items()):
    print(f"  {status}\t{tree}\t{n}")
