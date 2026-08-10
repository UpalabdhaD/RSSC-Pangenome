#!/usr/bin/env bash
set -euo pipefail

# Run HyPhy BUSTED for the existing "Absent in large lineages" codon alignments.
# This is a data-analysis workflow only; it does not alter genome or alignment inputs.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
work_dir="$repo_root/MS_REVISION/DNDS_HYPHY"
prep_dir="$repo_root/MS_REVISION/DNDS_ALT_MODELS"

category_source="${CATEGORY_SOURCE:-$repo_root/MS_REVISION/DNDS_SENSITIVITY/results/kaks_sensitivity_joined.tsv}"
alignment_dir="${ALIGNMENT_DIR:-$repo_root/Results/DNDS/PAL2NAL_DEDUP}"
metadata="${METADATA:-$repo_root/Metadata/Full_RSSC_142_Metadata.tsv}"
species_tree="${SPECIES_TREE:-$repo_root/Results/UFbootstrap_iqtree/Midpoint_rooted/ufbootstrap_rssc_coresnpsites_midpointrooted.treefile}"
out_root="${OUT_ROOT:-$work_dir/results_absent_major_lineages_busted}"
hyphy_prefix="${HYPHY_PREFIX:-$work_dir/envs/hyphy_2.5.50}"
hyphy_bin="${HYPHY_BIN:-$hyphy_prefix/bin/hyphy}"
hyphy_libpath="${HYPHY_LIBPATH:-$hyphy_prefix/share/hyphy}"
min_sequences="${MIN_SEQUENCES:-4}"
timeout_seconds="${TIMEOUT_SECONDS:-3600}"
grid_size="${GRID_SIZE:-10}"

mkdir -p "$out_root"/{prepared,busted,logs,status}

if [[ ! -x "$hyphy_bin" ]]; then
  echo "ERROR: HyPhy executable not found: $hyphy_bin" >&2
  exit 1
fi
if [[ ! -d "$hyphy_libpath" ]]; then
  echo "ERROR: HyPhy library path not found: $hyphy_libpath" >&2
  exit 1
fi

python3 "$work_dir/scripts/build_category_manifest.py" \
  "$category_source" \
  "$alignment_dir" \
  "$out_root/category_manifest.tsv" \
  "$out_root/genes_min${min_sequences}.txt" \
  "$min_sequences"

run_one_gene() {
  local gene="$1"
  local aln="$alignment_dir/$gene.codon.dedup.fa"
  local prepared_fasta="$out_root/prepared/$gene.codon.fa"
  local prepared_phy="$out_root/prepared/$gene.phy"
  local prepared_tree="$out_root/prepared/$gene.tree"
  local manifest="$out_root/prepared/$gene.manifest.tsv"
  local gene_out="$out_root/busted/$gene"
  local json="$gene_out/busted.json"
  local summary="$gene_out/summary.tsv"
  local status_file="$out_root/status/$gene.status"

  mkdir -p "$gene_out"
  if [[ -s "$summary" ]] && grep -q $'\tok\t' "$summary"; then
    printf "%s\tskipped_existing_ok\n" "$gene" > "$status_file"
    return 0
  fi

  if ! Rscript "$prep_dir/scripts/prepare_gene_inputs.R" \
    --alignment "$aln" \
    --metadata "$metadata" \
    --species-tree "$species_tree" \
    --out-fasta "$prepared_fasta" \
    --out-phylip "$prepared_phy" \
    --out-tree "$prepared_tree" \
    --manifest "$manifest" \
    --min-sequences "$min_sequences" > "$out_root/logs/$gene.prepare.log" 2>&1; then
    printf "%s\tprepare_failed\n" "$gene" > "$status_file"
    return 0
  fi

  if timeout "$timeout_seconds" "$hyphy_bin" "LIBPATH=$hyphy_libpath" BUSTED \
    --alignment "$prepared_fasta" \
    --tree "$prepared_tree" \
    --branches All \
    --srv Yes \
    --rates 2 \
    --syn-rates 2 \
    --grid-size "$grid_size" \
    --starting-points 1 \
    --output "$json" > "$out_root/logs/$gene.busted.log" 2>&1; then
    if [[ -s "$json" ]]; then
      python3 "$work_dir/scripts/summarize_hyphy_busted.py" "$gene" "$json" "$out_root/logs/$gene.busted.log" "$summary"
      printf "%s\tok\n" "$gene" > "$status_file"
    else
      printf "%s\tjson_empty\n" "$gene" > "$status_file"
    fi
  else
    printf "%s\tbusted_failed_or_timeout\n" "$gene" > "$status_file"
  fi
}

while IFS= read -r gene || [[ -n "$gene" ]]; do
  [[ -z "$gene" ]] && continue
  echo "Running HyPhy BUSTED for $gene"
  run_one_gene "$gene"
done < "$out_root/genes_min${min_sequences}.txt"

python3 "$work_dir/scripts/collect_busted_results.py" \
  "$out_root/category_manifest.tsv" \
  "$out_root/status" \
  "$out_root/busted" \
  "$out_root/hyphy_busted_summary.tsv"

echo "Wrote $out_root/hyphy_busted_summary.tsv"
