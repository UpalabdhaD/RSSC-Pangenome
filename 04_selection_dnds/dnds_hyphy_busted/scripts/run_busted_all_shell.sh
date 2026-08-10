#!/usr/bin/env bash
set -euo pipefail

# Run HyPhy BUSTED across all eligible deduplicated codon alignments.
# Existing successful summaries are skipped, so the run can be resumed safely.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
work_dir="$repo_root/MS_REVISION/DNDS_HYPHY"

out_root="${OUT_ROOT:-$work_dir/results_all_busted}"
alignment_dir="${ALIGNMENT_DIR:-$repo_root/Results/DNDS/PAL2NAL_DEDUP}"
metadata="${METADATA:-$repo_root/Metadata/Full_RSSC_142_Metadata.tsv}"
species_tree="${SPECIES_TREE:-$repo_root/Results/UFbootstrap_iqtree/Midpoint_rooted/ufbootstrap_rssc_coresnpsites_midpointrooted.treefile}"
hyphy_prefix="${HYPHY_PREFIX:-$work_dir/envs/hyphy_2.5.50}"
hyphy_bin="${HYPHY_BIN:-$hyphy_prefix/bin/hyphy}"
hyphy_libpath="${HYPHY_LIBPATH:-$hyphy_prefix/share/hyphy}"
min_sequences="${MIN_SEQUENCES:-4}"
jobs="${JOBS:-8}"
timeout_seconds="${TIMEOUT_SECONDS:-3600}"
grid_size="${GRID_SIZE:-50}"
genes_file="$out_root/genes_all_min${min_sequences}.txt"
sequence_counts_file="$out_root/genes_all_min${min_sequences}.sequence_counts.tsv"
manifest_file="$out_root/all_families_min${min_sequences}.manifest.tsv"

if [[ ! -x "$hyphy_bin" ]]; then
    echo "ERROR: HyPhy executable not found: $hyphy_bin" >&2
    exit 1
fi
if [[ ! -d "$hyphy_libpath" ]]; then
    echo "ERROR: HyPhy library path not found: $hyphy_libpath" >&2
    exit 1
fi

mkdir -p "$out_root"/{prepared,busted,logs,status}

python3 - "$alignment_dir" "$min_sequences" "$genes_file" "$sequence_counts_file" "$manifest_file" <<'PY'
import sys
from pathlib import Path

alignment_dir = Path(sys.argv[1])
min_sequences = int(sys.argv[2])
genes_file = Path(sys.argv[3])
counts_file = Path(sys.argv[4])
manifest_file = Path(sys.argv[5])

counts = []
eligible = []
for path in sorted(alignment_dir.glob("*.codon.dedup.fa")):
    gene = path.name.removesuffix(".codon.dedup.fa")
    n_sequences = sum(1 for line in path.open() if line.startswith(">"))
    counts.append((gene, n_sequences))
    if n_sequences >= min_sequences:
        eligible.append((n_sequences, gene))

counts_file.write_text(
    "gene\tdeduplicated_sequence_count\n"
    + "\n".join(f"{gene}\t{n_sequences}" for gene, n_sequences in counts)
    + "\n",
    encoding="utf-8",
)
manifest_file.write_text(
    "gene\tn_sequences\talignment\tstatus\n"
    + "\n".join(
        f"{gene}\t{n_sequences}\t{alignment_dir / (gene + '.codon.dedup.fa')}\t"
        f"{'eligible' if n_sequences >= min_sequences else 'too_few_sequences'}"
        for gene, n_sequences in counts
    )
    + "\n",
    encoding="utf-8",
)
# Finish smaller families first so progress is visible early in the run.
genes_file.write_text(
    "\n".join(gene for _, gene in sorted(eligible)) + "\n",
    encoding="utf-8",
)
print(f"total_alignments\t{len(counts)}")
print(f"eligible_min_sequences_{min_sequences}\t{len(eligible)}")
PY

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
        printf '%s\tskipped_existing_ok\n' "$gene" > "$status_file"
        return 0
    fi

    if ! Rscript "$work_dir/../DNDS_ALT_MODELS/scripts/prepare_gene_inputs.R" \
        --alignment "$aln" \
        --metadata "$metadata" \
        --species-tree "$species_tree" \
        --out-fasta "$prepared_fasta" \
        --out-phylip "$prepared_phy" \
        --out-tree "$prepared_tree" \
        --manifest "$manifest" \
        --min-sequences "$min_sequences" > "$out_root/logs/$gene.prepare.log" 2>&1; then
        printf '%s\tprepare_failed\n' "$gene" > "$status_file"
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
            printf '%s\tok\n' "$gene" > "$status_file"
        else
            printf '%s\tjson_empty\n' "$gene" > "$status_file"
        fi
    else
        printf '%s\tbusted_failed_or_timeout\n' "$gene" > "$status_file"
    fi
}

export -f run_one_gene
export alignment_dir metadata species_tree hyphy_bin hyphy_libpath min_sequences timeout_seconds grid_size out_root work_dir

echo "Running HyPhy BUSTED for all eligible families with JOBS=$jobs"
xargs -r -n 1 -P "$jobs" bash -c 'run_one_gene "$1"' _ < "$genes_file"

python3 "$work_dir/scripts/collect_busted_results.py" \
    "$manifest_file" \
    "$out_root/status" \
    "$out_root/busted" \
    "$out_root/hyphy_busted_summary.tsv"

{
    echo -e "status\tn"
    awk -F '\t' '{count[$2]++} END {for (status in count) print status "\t" count[status]}' "$out_root"/status/*.status | sort
} > "$out_root/busted_run_status.tsv"

echo "Wrote $out_root/hyphy_busted_summary.tsv"
echo "Wrote $out_root/busted_run_status.tsv"
