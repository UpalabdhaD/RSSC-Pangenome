#!/usr/bin/env bash
set -euo pipefail
# Build a FastTree 2.2.0 gene tree (-nt -gtr) for each family in genes.txt from the prepared
# BUSTED codon alignment (read only), and remove support values so HyPhy reads plain branch
# lengths. The settings match run_final.sh in ../paml_m0_completion, so each family gets the
# same tree in the M0 and BUSTED analyses.
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prepared="${PROJECT_ROOT:?Set PROJECT_ROOT to the analysis project folder}/MS_REVISION/DNDS_HYPHY/results_all_busted/prepared"
fasttree_bin="${FASTTREE_BIN:-$(command -v FastTree || true)}"
mkdir -p "$here"/trees "$here"/logs
one() {
  local g="$1"
  "$fasttree_bin" -nt -gtr -quiet < "$prepared/$g.codon.fa" > "$here/trees/$g.fasttree.nwk" 2> "$here/logs/$g.fasttree.log"
  Rscript -e 'suppressPackageStartupMessages(library(ape)); a <- commandArgs(TRUE); t <- read.tree(a[1]); t$node.label <- NULL; write.tree(t, a[2])' \
    "$here/trees/$g.fasttree.nwk" "$here/trees/$g.tree"
}
export -f one; export here prepared fasttree_bin
xargs -r -n 1 -P "${JOBS:-6}" bash -c 'one "$1"' _ < "$here/genes.txt"
echo "trees: $(ls "$here"/trees/*.tree | wc -l)"
