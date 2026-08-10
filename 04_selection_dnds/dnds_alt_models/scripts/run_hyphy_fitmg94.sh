#!/usr/bin/env bash
set -euo pipefail

gene="$1"
alignment="$2"
tree="$3"
outdir="$4"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$outdir"
json="$outdir/fitmg94.json"
summary="$outdir/summary.tsv"

if ! command -v hyphy >/dev/null 2>&1; then
  printf 'gene\tstatus\tomega\tmessage\n' > "$summary"
  printf '%s\tmissing_hyphy\tNA\thyphy executable not found\n' "$gene" >> "$summary"
  exit 127
fi

if hyphy fitmg94 --help >/dev/null 2>&1; then
  hyphy fitmg94 --alignment "$alignment" --tree "$tree" --output "$json"
elif hyphy FitMG94 --help >/dev/null 2>&1; then
  hyphy FitMG94 --alignment "$alignment" --tree "$tree" --output "$json"
else
  printf 'gene\tstatus\tomega\tmessage\n' > "$summary"
  printf '%s\tunsupported_hyphy_cli\tNA\tFitMG94 CLI not available in this HyPhy build\n' "$gene" >> "$summary"
  exit 2
fi

python "$script_dir/summarize_hyphy_fitmg94.py" "$gene" "$json" "$summary"
