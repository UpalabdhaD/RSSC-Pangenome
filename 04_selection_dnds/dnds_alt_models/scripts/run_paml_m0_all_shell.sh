#!/usr/bin/env bash
set -euo pipefail

# Run PAML codeml M0 for all eligible deduplicated codon alignments.
# This is shell-only and does not invoke Snakemake or HyPhy.
#
# Defaults:
#   - eligible gene threshold: at least 2 deduplicated sequences
#   - parallel jobs: 8
#   - output directory: MS_REVISION/DNDS_ALT_MODELS/results_all_paml_m0
#
# Examples:
#   bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_paml_m0_all_shell.sh
#   JOBS=12 MIN_SEQUENCES=4 bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_paml_m0_all_shell.sh

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
work_dir="$repo_root/MS_REVISION/DNDS_ALT_MODELS"

out_root="${OUT_ROOT:-$work_dir/results_all_paml_m0}"
alignment_dir="${ALIGNMENT_DIR:-$repo_root/Results/DNDS/PAL2NAL_DEDUP}"
metadata="${METADATA:-$repo_root/Metadata/Full_RSSC_142_Metadata.tsv}"
species_tree="${SPECIES_TREE:-$repo_root/Results/UFbootstrap_iqtree/Midpoint_rooted/ufbootstrap_rssc_coresnpsites_midpointrooted.treefile}"
min_sequences="${MIN_SEQUENCES:-2}"
jobs="${JOBS:-8}"
genes_file="${GENES_FILE:-$out_root/genes_all_min${min_sequences}.txt}"
gene_counts_file="$out_root/genes_all_min${min_sequences}.sequence_counts.tsv"

if [[ -n "${CONDA_ENV:-}" ]]; then
  conda_base="$(conda info --base)"
  # shellcheck disable=SC1091
  source "$conda_base/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
fi

codeml_bin="${CODEML_BIN:-$repo_root/.snakemake/conda/283dac3d4c48743ab09cd5eb6bb0ab9e_/bin/codeml}"
if [[ ! -x "$codeml_bin" ]]; then
  codeml_bin="$(command -v codeml || true)"
fi
if [[ -z "$codeml_bin" || ! -x "$codeml_bin" ]]; then
  echo "ERROR: codeml not found. Set CODEML_BIN=/path/to/codeml or CONDA_ENV=<env_with_paml>." >&2
  exit 1
fi

mkdir -p "$out_root"/prepared "$out_root"/paml_m0 "$out_root"/logs "$out_root"/status

echo "Building eligible gene list: $genes_file"
python3 - "$alignment_dir" "$min_sequences" "$genes_file" "$gene_counts_file" <<'PY'
import sys
from pathlib import Path

alignment_dir = Path(sys.argv[1])
min_sequences = int(sys.argv[2])
genes_file = Path(sys.argv[3])
gene_counts_file = Path(sys.argv[4])

eligible = []
counts = []
total = 0
for path in sorted(alignment_dir.glob("*.codon.dedup.fa")):
    total += 1
    nseq = sum(1 for line in path.open() if line.startswith(">"))
    gene = path.name.removesuffix(".codon.dedup.fa")
    counts.append((gene, nseq))
    if nseq >= min_sequences:
        eligible.append((nseq, gene))

genes_file.parent.mkdir(parents=True, exist_ok=True)
gene_counts_file.write_text(
    "gene\tdeduplicated_sequence_count\n"
    + "\n".join(f"{gene}\t{nseq}" for gene, nseq in sorted(counts, key=lambda x: (x[1], x[0])))
    + "\n"
)

# Run small alignments first. This gets most accessory genes completed quickly
# while leaving the largest core-gene alignments for later.
eligible_genes = [gene for nseq, gene in sorted(eligible, key=lambda x: (x[0], x[1]))]
genes_file.write_text("\n".join(eligible_genes) + ("\n" if eligible_genes else ""))
print(f"total_alignments\t{total}")
print(f"eligible_min_sequences_{min_sequences}\t{len(eligible)}")
PY

run_one_gene() {
  local gene="$1"
  local aln="$alignment_dir/$gene.codon.dedup.fa"
  local gene_out="$out_root/paml_m0/$gene"
  local summary="$gene_out/summary.tsv"
  local status_file="$out_root/status/$gene.status"
  local prep_log="$out_root/logs/$gene.prepare.log"
  local paml_log="$out_root/logs/$gene.paml_m0.log"

  if [[ -s "$summary" ]] && grep -q $'\tok\t' "$summary"; then
    printf "%s\t%s\n" "$gene" "skipped_existing_ok" > "$status_file"
    return 0
  fi

  if [[ ! -s "$aln" ]]; then
    printf "%s\t%s\n" "$gene" "missing_alignment" > "$status_file"
    return 0
  fi

  if ! Rscript "$work_dir/scripts/prepare_gene_inputs.R" \
    --alignment "$aln" \
    --metadata "$metadata" \
    --species-tree "$species_tree" \
    --out-fasta "$out_root/prepared/$gene.codon.fa" \
    --out-phylip "$out_root/prepared/$gene.phy" \
    --out-tree "$out_root/prepared/$gene.tree" \
    --manifest "$out_root/prepared/$gene.manifest.tsv" \
    --min-sequences "$min_sequences" > "$prep_log" 2>&1; then
    printf "%s\t%s\n" "$gene" "prepare_failed" > "$status_file"
    return 0
  fi

  if ! CODEML="$codeml_bin" Rscript "$work_dir/scripts/run_paml_m0.R" \
    --gene "$gene" \
    --phylip "$out_root/prepared/$gene.phy" \
    --tree "$out_root/prepared/$gene.tree" \
    --outdir "$gene_out" > "$paml_log" 2>&1; then
    printf "%s\t%s\n" "$gene" "paml_failed" > "$status_file"
    return 0
  fi

  printf "%s\t%s\n" "$gene" "ok" > "$status_file"
}

export -f run_one_gene
export alignment_dir metadata species_tree min_sequences out_root work_dir codeml_bin

echo "Running PAML M0 with JOBS=$jobs"
xargs -r -n 1 -P "$jobs" bash -c 'run_one_gene "$1"' _ < "$genes_file"

Rscript "$work_dir/scripts/collect_summaries.R" \
  "$out_root/paml_m0" \
  "$out_root/paml_m0_summary.tsv"

{
  echo -e "status\tn"
  awk -F '\t' '{count[$2]++} END {for (status in count) print status "\t" count[status]}' "$out_root"/status/*.status | sort
} > "$out_root/paml_m0_run_status.tsv"

echo "Wrote $out_root/paml_m0_summary.tsv"
echo "Wrote $out_root/paml_m0_run_status.tsv"
