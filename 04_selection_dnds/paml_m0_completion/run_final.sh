#!/usr/bin/env bash
set -euo pipefail

# Final PAML codeml M0 run for 697 gene families (queue.tsv: sequences, family, tree, origin).
#
#   species_tree  81 unfinished families with one sequence per strain. Same steps as
#                 run_one_gene() in run_paml_m0_all_shell.sh: the core-SNP species tree pruned
#                 to the family.
#   gene_tree     616 families with paralogs (two or more sequences from one strain): 356
#                 unfinished families and the 260 families that the original run fitted on the
#                 backup ladder tree. The species tree cannot hold paralogs, so we infer a gene
#                 tree from the prepared codon alignment with FastTree 2.2.0 (-nt -gtr), remove
#                 the support values, and give that tree to codeml.
#
# prepare_gene_inputs.R and run_paml_m0.R are the originals in ../dnds_alt_models/scripts/, so
# the codon alignment and the codeml settings are the same for both tree types. Inputs are read
# from the analysis project folder ($PROJECT_ROOT); nothing is written there.
#
# Usage: PROJECT_ROOT=/path/to/analysis/project JOBS=24 bash run_final.sh
# Optional: CODEML_BIN=/path/to/codeml (PAML 4.10.0), FASTTREE_BIN=/path/to/FastTree (2.2.0)

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_root="${PROJECT_ROOT:?Set PROJECT_ROOT to the analysis project folder}"

queue="${QUEUE:-$here/queue.tsv}"
out_root="${OUT_ROOT:-$here/results_final}"
alignment_dir="$analysis_root/Results/DNDS/PAL2NAL_DEDUP"
metadata="$analysis_root/Metadata/Full_RSSC_142_Metadata.tsv"
species_tree="$analysis_root/Results/UFbootstrap_iqtree/Midpoint_rooted/ufbootstrap_rssc_coresnpsites_midpointrooted.treefile"
codeml_bin="${CODEML_BIN:-$analysis_root/.snakemake/conda/283dac3d4c48743ab09cd5eb6bb0ab9e_/bin/codeml}"
fasttree_bin="${FASTTREE_BIN:-$(command -v FastTree || true)}"
script_dir="$here/../dnds_alt_models/scripts"
min_sequences=2
jobs="${JOBS:-24}"

for bin in "$codeml_bin" "$fasttree_bin"; do
  [[ -x "$bin" ]] || { echo "ERROR: not executable: $bin" >&2; exit 1; }
done

mkdir -p "$out_root"/prepared "$out_root"/gene_trees "$out_root"/paml_m0 "$out_root"/logs "$out_root"/status

run_one_gene() {
  local gene="$1" tree_type="$2"
  local aln="$alignment_dir/$gene.codon.dedup.fa"
  local gene_out="$out_root/paml_m0/$gene"
  local status_file="$out_root/status/$gene.status"
  local tree_file="$out_root/prepared/$gene.tree"

  if [[ -s "$gene_out/summary.tsv" ]] && grep -q $'\tok\t' "$gene_out/summary.tsv"; then
    printf "%s\t%s\t\t%s\n" "$gene" "skipped_existing_ok" "$tree_type" > "$status_file"
    return 0
  fi

  if [[ ! -s "$aln" ]]; then
    printf "%s\t%s\t\t%s\n" "$gene" "missing_alignment" "$tree_type" > "$status_file"
    return 0
  fi

  if ! Rscript "$script_dir/prepare_gene_inputs.R" \
    --alignment "$aln" \
    --metadata "$metadata" \
    --species-tree "$species_tree" \
    --out-fasta "$out_root/prepared/$gene.codon.fa" \
    --out-phylip "$out_root/prepared/$gene.phy" \
    --out-tree "$out_root/prepared/$gene.tree" \
    --manifest "$out_root/prepared/$gene.manifest.tsv" \
    --min-sequences "$min_sequences" > "$out_root/logs/$gene.prepare.log" 2>&1; then
    printf "%s\t%s\t\t%s\n" "$gene" "prepare_failed" "$tree_type" > "$status_file"
    return 0
  fi

  if [[ "$tree_type" == "gene_tree" ]]; then
    tree_file="$out_root/gene_trees/$gene.tree"
    if ! { "$fasttree_bin" -nt -gtr -quiet < "$out_root/prepared/$gene.codon.fa" \
             > "$out_root/gene_trees/$gene.fasttree.nwk" 2> "$out_root/logs/$gene.fasttree.log" &&
           Rscript -e 'suppressPackageStartupMessages(library(ape)); a <- commandArgs(TRUE); t <- read.tree(a[1]); t$node.label <- NULL; write.tree(t, a[2])' \
             "$out_root/gene_trees/$gene.fasttree.nwk" "$tree_file" >> "$out_root/logs/$gene.fasttree.log" 2>&1; }; then
      printf "%s\t%s\t\t%s\n" "$gene" "gene_tree_failed" "$tree_type" > "$status_file"
      return 0
    fi
  fi

  local t0=$SECONDS
  if ! CODEML="$codeml_bin" Rscript "$script_dir/run_paml_m0.R" \
    --gene "$gene" \
    --phylip "$out_root/prepared/$gene.phy" \
    --tree "$tree_file" \
    --outdir "$gene_out" > "$out_root/logs/$gene.paml_m0.log" 2>&1; then
    printf "%s\t%s\t\t%s\n" "$gene" "paml_failed" "$tree_type" > "$status_file"
    return 0
  fi

  printf "%s\t%s\t%s\t%s\n" "$gene" "ok" "$((SECONDS - t0))" "$tree_type" > "$status_file"
}

export -f run_one_gene
export alignment_dir metadata species_tree min_sequences out_root script_dir codeml_bin fasttree_bin

echo "Started $(date -Is): $(wc -l < "$queue") families, JOBS=$jobs"
cut -f2,3 "$queue" | xargs -r -n 2 -P "$jobs" bash -c 'run_one_gene "$1" "$2"' _

Rscript "$script_dir/collect_summaries.R" \
  "$out_root/paml_m0" \
  "$out_root/paml_m0_summary.tsv"

{
  echo -e "status\ttree\tn"
  awk -F '\t' '{count[$2 "\t" $4]++} END {for (k in count) print k "\t" count[k]}' "$out_root"/status/*.status | sort
} > "$out_root/paml_m0_run_status.tsv"

echo "Finished $(date -Is)"
cat "$out_root/paml_m0_run_status.tsv"
