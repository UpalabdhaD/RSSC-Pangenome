#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
work_dir="$repo_root/MS_REVISION/DNDS_ALT_MODELS"

genes_file="${1:-$work_dir/config/genes_pilot.txt}"
out_root="${OUT_ROOT:-$work_dir/results_shell}"
alignment_dir="${ALIGNMENT_DIR:-$repo_root/Results/DNDS/PAL2NAL_DEDUP}"
metadata="${METADATA:-$repo_root/Metadata/Full_RSSC_142_Metadata.tsv}"
species_tree="${SPECIES_TREE:-$repo_root/Results/UFbootstrap_iqtree/Midpoint_rooted/ufbootstrap_rssc_coresnpsites_midpointrooted.treefile}"
min_sequences="${MIN_SEQUENCES:-4}"

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

mkdir -p "$out_root"/prepared "$out_root"/paml_m0 "$out_root"/logs

while IFS= read -r gene || [[ -n "$gene" ]]; do
  [[ -z "$gene" || "$gene" =~ ^# ]] && continue

  aln="$alignment_dir/$gene.codon.dedup.fa"
  if [[ ! -s "$aln" ]]; then
    echo "WARN: missing alignment for $gene: $aln" >&2
    continue
  fi

  echo "Preparing $gene"
  Rscript "$work_dir/scripts/prepare_gene_inputs.R" \
    --alignment "$aln" \
    --metadata "$metadata" \
    --species-tree "$species_tree" \
    --out-fasta "$out_root/prepared/$gene.codon.fa" \
    --out-phylip "$out_root/prepared/$gene.phy" \
    --out-tree "$out_root/prepared/$gene.tree" \
    --manifest "$out_root/prepared/$gene.manifest.tsv" \
    --min-sequences "$min_sequences"

  echo "Running PAML M0 for $gene"
  CODEML="$codeml_bin" Rscript "$work_dir/scripts/run_paml_m0.R" \
    --gene "$gene" \
    --phylip "$out_root/prepared/$gene.phy" \
    --tree "$out_root/prepared/$gene.tree" \
    --outdir "$out_root/paml_m0/$gene" \
    > "$out_root/logs/$gene.paml_m0.log" 2>&1
done < "$genes_file"

Rscript "$work_dir/scripts/collect_summaries.R" \
  "$out_root/paml_m0" \
  "$out_root/paml_m0_summary.tsv"

echo "Wrote $out_root/paml_m0_summary.tsv"
