#!/usr/bin/env bash
set -euo pipefail

# Test whether the backup tree changes the M0 omega. For 30 families that the original run
# fitted on the backup ladder tree (families with paralogs), infer a gene tree from the same
# prepared codon alignment with FastTree (GTR, nucleotides), remove support values, and rerun
# codeml M0 with the unchanged run_paml_m0.R. Inputs are read from the analysis folder.
#
# Usage: JOBS=6 bash run_tree_test.sh

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prepared="${PROJECT_ROOT:?Set PROJECT_ROOT to the analysis project folder}/MS_REVISION/DNDS_ALT_MODELS/results_all_paml_m0/prepared"
codeml_bin="${CODEML_BIN:-$PROJECT_ROOT/.snakemake/conda/283dac3d4c48743ab09cd5eb6bb0ab9e_/bin/codeml}"
fasttree_bin="${FASTTREE_BIN:-$(command -v FastTree || true)}"
script_dir="$here/../../dnds_alt_models/scripts"
jobs="${JOBS:-6}"

mkdir -p "$here"/gene_trees "$here"/paml_m0 "$here"/logs

run_one() {
  local gene="$1"
  "$fasttree_bin" -nt -gtr -quiet < "$prepared/$gene.codon.fa" > "$here/gene_trees/$gene.fasttree.nwk" 2> "$here/logs/$gene.fasttree.log"
  Rscript -e 'suppressPackageStartupMessages(library(ape)); a <- commandArgs(TRUE); t <- read.tree(a[1]); t$node.label <- NULL; write.tree(t, a[2])' \
    "$here/gene_trees/$gene.fasttree.nwk" "$here/gene_trees/$gene.tree"
  CODEML="$codeml_bin" Rscript "$script_dir/run_paml_m0.R" \
    --gene "$gene" \
    --phylip "$prepared/$gene.phy" \
    --tree "$here/gene_trees/$gene.tree" \
    --outdir "$here/paml_m0/$gene" > "$here/logs/$gene.paml_m0.log" 2>&1 || echo "$gene failed" >&2
}

export -f run_one
export here prepared codeml_bin fasttree_bin script_dir

cut -f2 "$here/sample30.tsv" | xargs -r -n 1 -P "$jobs" bash -c 'run_one "$1"' _

{
  echo -e "gene\tsequences\tomega_backup_tree\tomega_gene_tree"
  while IFS=$'\t' read -r n gene omega; do
    new=$(awk -F'\t' 'NR==2{print $3}' "$here/paml_m0/$gene/summary.tsv" 2>/dev/null || true)
    echo -e "$gene\t$n\t$omega\t${new:-NA}"
  done < "$here/sample30.tsv"
} > "$here/omega_backup_vs_gene_tree.tsv"
echo "Wrote $here/omega_backup_vs_gene_tree.tsv"
