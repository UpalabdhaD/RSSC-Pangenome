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

hyphy_bin="${HYPHY_BIN:-$repo_root/.snakemake/conda/aab4c17d994da1523aea6cac8550a122_/bin/hyphy}"
hyphy_libpath="${HYPHY_LIBPATH:-$repo_root/.snakemake/conda/aab4c17d994da1523aea6cac8550a122_/share/hyphy}"
if [[ ! -x "$hyphy_bin" ]]; then
  hyphy_bin="$(command -v hyphy || true)"
fi
if [[ -z "$hyphy_bin" || ! -x "$hyphy_bin" ]]; then
  echo "ERROR: hyphy not found. Set HYPHY_BIN=/path/to/hyphy or CONDA_ENV=<env_with_hyphy>." >&2
  exit 1
fi
if [[ ! -d "$hyphy_libpath" ]]; then
  echo "ERROR: HyPhy library path not found. Set HYPHY_LIBPATH=/path/to/share/hyphy." >&2
  exit 1
fi

mkdir -p "$out_root"/prepared "$out_root"/hyphy_fel "$out_root"/logs

while IFS= read -r gene || [[ -n "$gene" ]]; do
  [[ -z "$gene" || "$gene" =~ ^# ]] && continue

  aln="$alignment_dir/$gene.codon.dedup.fa"
  if [[ ! -s "$aln" ]]; then
    echo "WARN: missing alignment for $gene: $aln" >&2
    continue
  fi

  if [[ ! -s "$out_root/prepared/$gene.codon.fa" || ! -s "$out_root/prepared/$gene.tree" ]]; then
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
  fi

  gene_out="$out_root/hyphy_fel/$gene"
  mkdir -p "$gene_out"
  echo "Running HyPhy FEL for $gene"
  "$hyphy_bin" "LIBPATH=$hyphy_libpath" FEL \
    --alignment "$out_root/prepared/$gene.codon.fa" \
    --tree "$out_root/prepared/$gene.tree" \
    --branches All \
    --srv Yes \
    --pvalue 0.1 \
    --ci No \
    --resample 0 \
    --output "$gene_out/fel.json" \
    > "$out_root/logs/$gene.hyphy_fel.log" 2>&1

  python "$work_dir/scripts/summarize_hyphy_fel.py" \
    "$gene" \
    "$gene_out/fel.json" \
    "$gene_out/summary.tsv"
done < "$genes_file"

python "$work_dir/scripts/collect_hyphy_fel_summaries.py" \
  "$out_root/hyphy_fel" \
  "$out_root/hyphy_fel_summary.tsv"

echo "Wrote $out_root/hyphy_fel_summary.tsv"
