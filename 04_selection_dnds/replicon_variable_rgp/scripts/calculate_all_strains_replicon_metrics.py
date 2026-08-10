#!/usr/bin/env python3
"""All-strain chromosome vs non-chromosomal RGP/accessory metrics.

The script uses existing PanRGP/PPanGGOLiN outputs and RSSC FASTA files.
Outputs stay under MS_REVISION/REPLICON_VARIABLE_RGP/results.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
OUTDIR = REPO / "MS_REVISION" / "REPLICON_VARIABLE_RGP" / "results"
PANRGP_GENE_TABLE = REPO / "Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv"
RGP_TABLE = REPO / "Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv"
FASTA_DIR = REPO / "Data/PANGENOME_DATA/RSSC"

PANRGP_HEADER = [
    "genome",
    "gene",
    "contig",
    "start",
    "stop",
    "strand",
    "family",
    "nb_copy_in_genome",
    "partition",
    "persistent_neighbors",
    "shell_neighbors",
    "cloud_neighbors",
    "RGP",
    "Spot",
]


def contig_rank(contig: str) -> int:
    try:
        return int(contig.rsplit("_", 1)[1])
    except Exception:
        return 10**9


def replicon_class(contig: str) -> str:
    return "Chromosome" if contig == "contig_1" else "Non-chromosomal"


def read_fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    current = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
    if current:
        lengths.append(current)
    return lengths


def load_replicon_lengths(genome_contigs: dict[str, set[str]]) -> dict[tuple[str, str], int]:
    lengths: dict[tuple[str, str], int] = {}
    for genome, contigs in genome_contigs.items():
        fasta = FASTA_DIR / f"{genome}.fna"
        if not fasta.exists():
            continue
        seq_lengths = read_fasta_lengths(fasta)
        for contig in sorted(contigs, key=contig_rank):
            idx = contig_rank(contig) - 1
            if 0 <= idx < len(seq_lengths):
                lengths[(genome, contig)] = seq_lengths[idx]
    return lengths


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    def log_comb(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def log_hypergeom(x: int, row1: int, col1: int, n: int) -> float:
        return log_comb(col1, x) + log_comb(n - col1, row1 - x) - log_comb(n, row1)

    row1 = a + b
    col1 = a + c
    n = a + b + c + d
    observed_log = log_hypergeom(a, row1, col1, n)
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    included = [
        lp
        for x in range(lo, hi + 1)
        for lp in [log_hypergeom(x, row1, col1, n)]
        if lp <= observed_log + 1e-12
    ]
    max_log = max(included)
    log_p = max_log + math.log(sum(math.exp(lp - max_log) for lp in included))
    p_value = math.exp(log_p) if log_p > math.log(1e-300) else 0.0
    log10_p = log_p / math.log(10)
    odds = (a * d) / (b * c) if b and c else math.inf
    return odds, min(1.0, p_value), log10_p


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    gene_rows: list[dict] = []
    genome_contigs: dict[str, set[str]] = defaultdict(set)
    with PANRGP_GENE_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=PANRGP_HEADER)
        for row in reader:
            row["start"] = int(row["start"])
            row["stop"] = int(row["stop"])
            row["replicon_class"] = replicon_class(row["contig"])
            gene_rows.append(row)
            genome_contigs[row["genome"]].add(row["contig"])

    lengths = load_replicon_lengths(genome_contigs)

    retained_rgps: set[str] = set()
    rgp_by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    with RGP_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("wholeContig", "FALSE").upper() == "TRUE":
                continue
            retained_rgps.add(row["rgp_id"])
            row["length"] = int(float(row["length"]))
            row["genes"] = int(float(row["genes"]))
            row["replicon_class"] = replicon_class(row["contig"])
            rgp_by_key[(row["genome"], row["contig"])].append(row)

    per_replicon: list[dict] = []
    rows_by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in gene_rows:
        rows_by_key[(row["genome"], row["contig"])].append(row)

    for (genome, contig), rows in sorted(rows_by_key.items()):
        length_bp = lengths.get((genome, contig))
        if not length_bp:
            # Should not happen for this dataset; keep explicit if it does.
            length_bp = max(r["stop"] for r in rows)
        length_mb = length_bp / 1_000_000
        part = Counter(r["partition"] for r in rows)
        persistent = part["persistent"]
        shell = part["shell"]
        cloud = part["cloud"]
        variable = shell + cloud
        rgps = rgp_by_key.get((genome, contig), [])
        rgp_length = sum(r["length"] for r in rgps)
        rgp_genes = sum(1 for r in rows if r["RGP"] in retained_rgps)
        total = len(rows)
        per_replicon.append(
            {
                "genome": genome,
                "contig": contig,
                "replicon_class": replicon_class(contig),
                "replicon_length_bp": length_bp,
                "total_genes": total,
                "persistent_genes": persistent,
                "shell_genes": shell,
                "cloud_genes": cloud,
                "variable_genes": variable,
                "variable_gene_percent": 100 * variable / total,
                "accessory_genes_per_mb": variable / length_mb,
                "rgp_count": len(rgps),
                "rgp_length_bp": rgp_length,
                "rgp_coverage_percent": 100 * rgp_length / length_bp,
                "rgp_genes": rgp_genes,
                "rgp_genes_per_mb": rgp_genes / length_mb,
            }
        )

    # Aggregate by replicon class.
    agg: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in per_replicon:
        cls = row["replicon_class"]
        for key in [
            "replicon_length_bp",
            "total_genes",
            "persistent_genes",
            "shell_genes",
            "cloud_genes",
            "variable_genes",
            "rgp_count",
            "rgp_length_bp",
            "rgp_genes",
        ]:
            agg[cls][key] += row[key]

    summary_rows: list[dict] = []
    for cls in ["Chromosome", "Non-chromosomal"]:
        a = agg[cls]
        length_mb = a["replicon_length_bp"] / 1_000_000
        total = a["total_genes"]
        variable = a["variable_genes"]
        summary_rows.append(
            {
                "replicon_class": cls,
                "genome_count": len({r["genome"] for r in per_replicon if r["replicon_class"] == cls}),
                "replicon_count": sum(1 for r in per_replicon if r["replicon_class"] == cls),
                "total_length_bp": int(a["replicon_length_bp"]),
                "total_genes": int(total),
                "persistent_genes": int(a["persistent_genes"]),
                "shell_genes": int(a["shell_genes"]),
                "cloud_genes": int(a["cloud_genes"]),
                "variable_genes": int(variable),
                "variable_gene_percent": 100 * variable / total,
                "accessory_genes_per_mb": variable / length_mb,
                "rgp_count": int(a["rgp_count"]),
                "rgp_length_bp": int(a["rgp_length_bp"]),
                "rgp_coverage_percent": 100 * a["rgp_length_bp"] / a["replicon_length_bp"],
                "rgp_genes": int(a["rgp_genes"]),
                "rgp_genes_per_mb": a["rgp_genes"] / length_mb,
                "rgp_sequence_per_mb_bp": a["rgp_length_bp"] / length_mb,
            }
        )

    chrom = next(r for r in summary_rows if r["replicon_class"] == "Chromosome")
    plasmid = next(r for r in summary_rows if r["replicon_class"] == "Non-chromosomal")
    a = int(chrom["variable_genes"])
    b = int(chrom["persistent_genes"])
    c = int(plasmid["variable_genes"])
    d = int(plasmid["persistent_genes"])
    odds_chr, p_value, log10_p = fisher_exact_two_sided(a, b, c, d)
    odds_nonchrom = 1 / odds_chr
    p_text = f"{p_value:.8g}" if p_value > 0 else f"<1e-300; log10_p={log10_p:.3f}"

    with (OUTDIR / "all_strains_per_replicon_metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(per_replicon[0].keys()))
        writer.writeheader()
        writer.writerows(per_replicon)

    with (OUTDIR / "all_strains_replicon_class_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    with (OUTDIR / "all_strains_fisher_variable_vs_persistent.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "comparison",
                "chromosome_variable",
                "chromosome_persistent",
                "nonchromosomal_variable",
                "nonchromosomal_persistent",
                "odds_ratio_chromosome_vs_nonchromosomal",
                "odds_ratio_nonchromosomal_vs_chromosome",
                "fisher_exact_two_sided_p",
                "fisher_exact_two_sided_log10_p",
            ]
        )
        writer.writerow(
            [
                "shell_plus_cloud_vs_persistent_all_strains",
                a,
                b,
                c,
                d,
                f"{odds_chr:.8g}",
                f"{odds_nonchrom:.8g}",
                p_text,
                f"{log10_p:.8g}",
            ]
        )

    with (OUTDIR / "all_strains_replicon_class_summary.md").open("w") as handle:
        handle.write("# All-Strain Chromosome vs Non-chromosomal Summary\n\n")
        handle.write("| Metric | Chromosome | Non-chromosomal |\n")
        handle.write("|---|---:|---:|\n")
        metrics = [
            ("Genome count", "genome_count", "{:,.0f}"),
            ("Replicon count", "replicon_count", "{:,.0f}"),
            ("Total length", "total_length_bp", "{:,.0f} bp"),
            ("Total genes", "total_genes", "{:,.0f}"),
            ("Persistent genes", "persistent_genes", "{:,.0f}"),
            ("Shell genes", "shell_genes", "{:,.0f}"),
            ("Cloud genes", "cloud_genes", "{:,.0f}"),
            ("Shell + cloud genes", "variable_genes", "{:,.0f}"),
            ("Variable gene percentage", "variable_gene_percent", "{:.2f}%"),
            ("Accessory genes per Mb", "accessory_genes_per_mb", "{:.1f}"),
            ("RGP count", "rgp_count", "{:,.0f}"),
            ("RGP length", "rgp_length_bp", "{:,.0f} bp"),
            ("RGP coverage", "rgp_coverage_percent", "{:.2f}%"),
            ("RGP sequence per Mb", "rgp_sequence_per_mb_bp", "{:,.0f} bp/Mb"),
            ("RGP genes", "rgp_genes", "{:,.0f}"),
            ("RGP genes per Mb", "rgp_genes_per_mb", "{:.1f}"),
        ]
        for label, key, fmt in metrics:
            handle.write(f"| {label} | {fmt.format(chrom[key])} | {fmt.format(plasmid[key])} |\n")
        handle.write("\n## Fisher exact test\n\n")
        handle.write("| Replicon class | Variable genes | Persistent genes |\n")
        handle.write("|---|---:|---:|\n")
        handle.write(f"| Chromosome | {a:,} | {b:,} |\n")
        handle.write(f"| Non-chromosomal | {c:,} | {d:,} |\n\n")
        handle.write(f"- Odds ratio, non-chromosomal vs chromosome: {odds_nonchrom:.3f}\n")
        handle.write(f"- Fisher exact two-sided p-value: {p_text}\n")
        handle.write(f"- Fisher exact two-sided log10(p): {log10_p:.3f}\n\n")
        handle.write("Note: non-chromosomal corresponds to PanRGP `contig_2+` entries labelled `Plasmid` in the manuscript RGP table; in most complete RSSC genomes this is the megaplasmid-side replicon class.\n")

    readme = f"""# All-Strain Replicon Variable/RGP Analysis

This all-strain calculation extends the F1C1-only table to all {len(genome_contigs)} genomes represented in `Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv`.

## Source files

- Gene partitions and RGP membership: `{PANRGP_GENE_TABLE.relative_to(REPO)}`
- Retained RGP IDs, lengths, and coordinates: `{RGP_TABLE.relative_to(REPO)}`
- Replicon lengths: sequence lengths from `{FASTA_DIR.relative_to(REPO)}/*.fna`, mapped by FASTA order to PanRGP `contig_1`, `contig_2`, etc.

## Replicon classes

- `Chromosome`: PanRGP `contig_1`.
- `Non-chromosomal`: PanRGP `contig_2`, `contig_3`, or `contig_4`; this corresponds to the manuscript RGP table label `Plasmid`.

## Outputs

- `all_strains_per_replicon_metrics.tsv`: one row per genome/replicon.
- `all_strains_replicon_class_summary.tsv`: aggregate chromosome vs non-chromosomal metrics.
- `all_strains_replicon_class_summary.md`: manuscript-readable table.
- `all_strains_fisher_variable_vs_persistent.tsv`: Fisher exact test for shell+cloud vs persistent genes.

## Main test

The aggregate Fisher test compares:

| Replicon class | Variable genes | Persistent genes |
|---|---:|---:|
| Chromosome | {a:,} | {b:,} |
| Non-chromosomal | {c:,} | {d:,} |

- Odds ratio, non-chromosomal vs chromosome: {odds_nonchrom:.3f}
- Fisher exact two-sided p-value: {p_text}
- Fisher exact two-sided log10(p): {log10_p:.3f}
"""
    (OUTDIR / "README_all_strains.md").write_text(readme)

    print(f"Wrote {OUTDIR / 'all_strains_replicon_class_summary.md'}")
    print(f"Wrote {OUTDIR / 'all_strains_per_replicon_metrics.tsv'}")
    print(f"Wrote {OUTDIR / 'all_strains_fisher_variable_vs_persistent.tsv'}")


if __name__ == "__main__":
    main()
