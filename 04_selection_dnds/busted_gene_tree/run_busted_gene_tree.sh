#!/usr/bin/env bash
set -euo pipefail

# Rerun HyPhy BUSTED on FastTree gene trees for the 607 families that the original BUSTED run
# fitted on the backup ladder tree (families with paralogs). We ran it on a second lab server that
# mounts the analysis project at the same path. The alignment, HyPhy binary and options are the same as run_one_gene()
# in MS_REVISION/DNDS_HYPHY/scripts/run_busted_all_shell.sh, with three differences:
#   --tree     trees/<gene>.tree (FastTree gene tree) in place of the prepared backup tree
#   timeout    6 hours per family in place of 1 hour, so large families can finish
#   CPU=2      two threads per HyPhy process, to share the server's cores across parallel jobs
# Inputs are read from $PROJECT_ROOT; outputs go to results/ next to this script.
#
# Usage: PROJECT_ROOT=/path/to/analysis/project JOBS=48 nohup bash run_busted_gene_tree.sh > run.log 2>&1 &

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
hyphy_dir="${PROJECT_ROOT:?Set PROJECT_ROOT to the analysis project folder}/MS_REVISION/DNDS_HYPHY"
prepared="$hyphy_dir/results_all_busted/prepared"
hyphy_prefix="$hyphy_dir/envs/hyphy_2.5.50"
hyphy_bin="$hyphy_prefix/bin/hyphy"
hyphy_libpath="$hyphy_prefix/share/hyphy"
summarize="$hyphy_dir/scripts/summarize_hyphy_busted.py"
out_root="$here/results"
jobs="${JOBS:-48}"
threads="${THREADS:-2}"
timeout_seconds="${TIMEOUT_SECONDS:-21600}"
grid_size=50

mkdir -p "$out_root"/busted "$out_root"/logs "$out_root"/status

run_one_gene() {
  local gene="$1"
  local gene_out="$out_root/busted/$gene"
  local json="$gene_out/busted.json"
  local summary="$gene_out/summary.tsv"
  local status_file="$out_root/status/$gene.status"

  mkdir -p "$gene_out"
  if [[ -s "$summary" ]] && grep -q $'\tok\t' "$summary"; then
    printf '%s\tskipped_existing_ok\n' "$gene" > "$status_file"
    return 0
  fi

  local t0=$SECONDS
  if timeout "$timeout_seconds" "$hyphy_bin" "LIBPATH=$hyphy_libpath" "CPU=$threads" BUSTED \
      --alignment "$prepared/$gene.codon.fa" \
      --tree "$here/trees/$gene.tree" \
      --branches All \
      --srv Yes \
      --rates 2 \
      --syn-rates 2 \
      --grid-size "$grid_size" \
      --starting-points 1 \
      --output "$json" > "$out_root/logs/$gene.busted.log" 2>&1; then
    if [[ -s "$json" ]]; then
      python3 "$summarize" "$gene" "$json" "$out_root/logs/$gene.busted.log" "$summary"
      printf '%s\tok\t%s\n' "$gene" "$((SECONDS - t0))" > "$status_file"
    else
      printf '%s\tjson_empty\n' "$gene" > "$status_file"
    fi
  else
    printf '%s\tbusted_failed_or_timeout\t%s\n' "$gene" "$((SECONDS - t0))" > "$status_file"
  fi
}

export -f run_one_gene
export here prepared hyphy_bin hyphy_libpath summarize out_root threads timeout_seconds grid_size

echo "Started $(date -Is) on $(hostname): $(wc -l < "$here/genes_by_size.txt") families, JOBS=$jobs, CPU=$threads"
xargs -r -n 1 -P "$jobs" bash -c 'run_one_gene "$1"' _ < "$here/genes_by_size.txt"

{
  echo -e "status\tn"
  awk -F '\t' '{count[$2]++} END {for (s in count) print s "\t" count[s]}' "$out_root"/status/*.status | sort
} > "$out_root/busted_run_status.tsv"
echo "Finished $(date -Is)"
cat "$out_root/busted_run_status.tsv"
