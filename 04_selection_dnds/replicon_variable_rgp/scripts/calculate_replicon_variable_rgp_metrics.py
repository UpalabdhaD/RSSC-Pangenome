#!/usr/bin/env python3
"""Calculate chromosome vs megaplasmid variable-gene/RGP metrics for F1C1.

Inputs are existing repository outputs only. Results are written under
MS_REVISION/REPLICON_VARIABLE_RGP/results.
"""

from __future__ import annotations

import csv
import math
import re
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
OUTDIR = REPO / "MS_REVISION" / "REPLICON_VARIABLE_RGP" / "results"

PANRGP_GENE_TABLE = REPO / "Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv"
RGP_TABLE = REPO / "Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv"
F1C1_GFF = REPO / "Data/genbank_assemblies/gffs/F1C1.gff"
T3E_TABLE = REPO / "Results/Annotation/ABRICATE/F1C1__RSSC_T3EV3.tab"
RGP_DEFENSE_TABLE = REPO / "Manuscript_Figures/Figure_DefenseSystem/Data/SPOT_RGP_DEFENSESYS.tsv"

REPLICON_LABEL = {
    "contig_1": "Chromosome",
    "chromosome": "Chromosome",
    "contig_2": "Megaplasmid",
    "plasmid": "Megaplasmid",
}

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

MOBILE_KEYWORDS = [
    "transposase",
    "integrase",
    "recombinase",
    "resolvase",
    "insertion sequence",
    "is family",
    "mobile",
    "phage",
    "prophage",
    "conjug",
    "relaxase",
    "transposon",
]


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = urllib.parse.unquote(value)
    return attrs


def read_replicon_lengths_and_mobile() -> tuple[dict[str, int], dict[str, set[str]], list[dict]]:
    lengths: dict[str, int] = {}
    mobile_genes: dict[str, set[str]] = defaultdict(set)
    mobile_rows: list[dict] = []

    for line in F1C1_GFF.open():
        line = line.rstrip("\n")
        if line.startswith("##sequence-region"):
            _, seqid, _start, stop = line.split()
            label = REPLICON_LABEL.get(seqid)
            if label:
                lengths[label] = int(stop)
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        seqid, start, stop = parts[0], int(parts[3]), int(parts[4])
        label = REPLICON_LABEL.get(seqid)
        if not label:
            continue
        attrs = parse_attrs(parts[8])
        gene = attrs.get("locus_tag") or attrs.get("ID") or ""
        product = attrs.get("product", "")
        low = product.lower()
        if any(keyword in low for keyword in MOBILE_KEYWORDS):
            mobile_genes[label].add(gene)
            mobile_rows.append(
                {
                    "replicon": label,
                    "gene": gene,
                    "start": start,
                    "stop": stop,
                    "product": product,
                }
            )
    return lengths, mobile_genes, mobile_rows


def read_panrgp_f1c1() -> list[dict]:
    rows: list[dict] = []
    with PANRGP_GENE_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=PANRGP_HEADER)
        for row in reader:
            if row["genome"] != "F1C1":
                continue
            row["replicon"] = REPLICON_LABEL.get(row["contig"], row["contig"])
            row["start"] = int(row["start"])
            row["stop"] = int(row["stop"])
            row["in_rgp"] = bool(row["RGP"])
            rows.append(row)
    return rows


def read_f1c1_rgps() -> list[dict]:
    rows: list[dict] = []
    with RGP_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["genome"] != "F1C1":
                continue
            # The manuscript RGP code notes that whole-contig RGP calls are not
            # treated as spots. F1C1 currently has none, but keep the guard.
            if row.get("wholeContig", "FALSE") == "TRUE":
                continue
            row["replicon"] = REPLICON_LABEL.get(row["contig"], row.get("CONTIGS", row["contig"]))
            row["length"] = int(float(row["length"]))
            row["genes"] = int(float(row["genes"]))
            row["start"] = int(float(row["start"]))
            row["stop"] = int(float(row["stop"]))
            rows.append(row)
    return rows


def read_defense_genes() -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with RGP_DEFENSE_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["GENOME"] == "F1C1" and row["GENETYPE"] == "DEFENSE":
                label = REPLICON_LABEL.get(row["CONTIG"], row["CONTIGS"])
                out[label].add(row["GENE"])
    return out


def read_t3e_hits() -> tuple[dict[str, int], list[dict]]:
    counts: Counter[str] = Counter()
    rows: list[dict] = []
    with T3E_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seq = row["SEQUENCE"]
            label = REPLICON_LABEL.get(seq, seq)
            counts[label] += 1
            rows.append(
                {
                    "replicon": label,
                    "sequence": seq,
                    "start": row["START"],
                    "end": row["END"],
                    "gene": row["GENE"],
                    "identity": row["%IDENTITY"],
                    "coverage": row["%COVERAGE"],
                }
            )
    return dict(counts), rows


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Return odds ratio and two-sided Fisher exact p-value.

    Table orientation:
        [[a variable chromosome, b persistent chromosome],
         [c variable megaplasmid, d persistent megaplasmid]]
    """

    def log_comb(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def log_hypergeom(x: int, row1: int, col1: int, n: int) -> float:
        return log_comb(col1, x) + log_comb(n - col1, row1 - x) - log_comb(n, row1)

    row1 = a + b
    row2 = c + d
    col1 = a + c
    n = row1 + row2
    observed_log = log_hypergeom(a, row1, col1, n)
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    included_logs = [
        lp
        for x in range(lo, hi + 1)
        for lp in [log_hypergeom(x, row1, col1, n)]
        if lp <= observed_log + 1e-12
    ]
    max_log = max(included_logs)
    p_value = math.exp(max_log) * sum(math.exp(lp - max_log) for lp in included_logs)
    odds_ratio = (a * d) / (b * c) if b and c else math.inf
    return odds_ratio, min(1.0, p_value)


def fmt_int(x: int) -> str:
    return f"{x:,}"


def fmt_float(x: float, digits: int = 2) -> str:
    return f"{x:.{digits}f}"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    lengths, mobile_genes, mobile_rows = read_replicon_lengths_and_mobile()
    pan_rows = read_panrgp_f1c1()
    rgp_rows = read_f1c1_rgps()
    defense_genes = read_defense_genes()
    t3e_counts, t3e_rows = read_t3e_hits()
    included_rgp_ids = {row["rgp_id"] for row in rgp_rows}

    replicons = ["Chromosome", "Megaplasmid"]
    pan_by_rep = {rep: [r for r in pan_rows if r["replicon"] == rep] for rep in replicons}
    rgp_by_rep = {rep: [r for r in rgp_rows if r["replicon"] == rep] for rep in replicons}

    metrics: dict[str, dict[str, float | int | str]] = defaultdict(dict)
    for rep in replicons:
        rows = pan_by_rep[rep]
        partitions = Counter(r["partition"] for r in rows)
        total = len(rows)
        persistent = partitions["persistent"]
        shell = partitions["shell"]
        cloud = partitions["cloud"]
        variable = shell + cloud
        rgps = rgp_by_rep[rep]
        rgp_len = sum(r["length"] for r in rgps)
        rgp_gene_count = sum(1 for r in rows if r["RGP"] in included_rgp_ids)
        length_bp = lengths[rep]
        length_mb = length_bp / 1_000_000

        metrics["Replicon length"][rep] = f"{length_mb:.2f} Mb"
        metrics["Total genes"][rep] = total
        metrics["Persistent genes"][rep] = persistent
        metrics["Shell genes"][rep] = shell
        metrics["Cloud genes"][rep] = cloud
        metrics["Shell + cloud genes"][rep] = variable
        metrics["Percentage of genes classified as variable"][rep] = 100 * variable / total
        metrics["Number of F1C1 RGPs"][rep] = len(rgps)
        metrics["Total RGP length"][rep] = rgp_len
        metrics["Percentage of replicon covered by RGPs"][rep] = 100 * rgp_len / length_bp
        metrics["RGP sequence per Mb"][rep] = rgp_len / length_mb
        metrics["Genes within RGPs"][rep] = rgp_gene_count
        metrics["RGP genes per Mb"][rep] = rgp_gene_count / length_mb
        metrics["Accessory genes per Mb"][rep] = variable / length_mb
        metrics["Mobile-element-associated genes"][rep] = len(mobile_genes.get(rep, set()))
        metrics["Defense-system genes"][rep] = len(defense_genes.get(rep, set()))
        metrics["Type III effector genes"][rep] = t3e_counts.get(rep, 0)

    # Fisher exact test: variable vs persistent by replicon.
    a = int(metrics["Shell + cloud genes"]["Chromosome"])
    b = int(metrics["Persistent genes"]["Chromosome"])
    c = int(metrics["Shell + cloud genes"]["Megaplasmid"])
    d = int(metrics["Persistent genes"]["Megaplasmid"])
    odds_ratio_chr_over_mp, p_value = fisher_exact_two_sided(a, b, c, d)
    odds_ratio_mp_over_chr = 1 / odds_ratio_chr_over_mp

    metric_order = [
        "Replicon length",
        "Total genes",
        "Persistent genes",
        "Shell genes",
        "Cloud genes",
        "Shell + cloud genes",
        "Percentage of genes classified as variable",
        "Number of F1C1 RGPs",
        "Total RGP length",
        "Percentage of replicon covered by RGPs",
        "RGP sequence per Mb",
        "Genes within RGPs",
        "RGP genes per Mb",
        "Accessory genes per Mb",
        "Mobile-element-associated genes",
        "Defense-system genes",
        "Type III effector genes",
    ]

    with (OUTDIR / "chromosome_vs_megaplasmid_metrics.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Metric", "Chromosome", "Megaplasmid"])
        for metric in metric_order:
            chrom = metrics[metric]["Chromosome"]
            mega = metrics[metric]["Megaplasmid"]
            if isinstance(chrom, float):
                chrom = fmt_float(chrom)
            if isinstance(mega, float):
                mega = fmt_float(mega)
            writer.writerow([metric, chrom, mega])

    with (OUTDIR / "fisher_variable_vs_persistent.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "comparison",
                "chromosome_variable",
                "chromosome_persistent",
                "megaplasmid_variable",
                "megaplasmid_persistent",
                "odds_ratio_chromosome_vs_megaplasmid",
                "odds_ratio_megaplasmid_vs_chromosome",
                "fisher_exact_two_sided_p",
            ]
        )
        writer.writerow(
            [
                "shell_plus_cloud_vs_persistent",
                a,
                b,
                c,
                d,
                f"{odds_ratio_chr_over_mp:.6g}",
                f"{odds_ratio_mp_over_chr:.6g}",
                f"{p_value:.6g}",
            ]
        )

    with (OUTDIR / "mobile_keyword_genes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["replicon", "gene", "start", "stop", "product"])
        writer.writeheader()
        writer.writerows(mobile_rows)

    with (OUTDIR / "t3e_hits_by_replicon.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["replicon", "sequence", "start", "end", "gene", "identity", "coverage"],
        )
        writer.writeheader()
        writer.writerows(t3e_rows)

    readme = f"""# F1C1 Chromosome vs Megaplasmid Variable Gene and RGP Metrics

## Source files

- Gene partition and RGP membership: `{PANRGP_GENE_TABLE.relative_to(REPO)}`
- F1C1 RGP coordinates/lengths: `{RGP_TABLE.relative_to(REPO)}`
- Replicon lengths and product keywords: `{F1C1_GFF.relative_to(REPO)}`
- Defense genes: `{RGP_DEFENSE_TABLE.relative_to(REPO)}`
- Type III effector hits: `{T3E_TABLE.relative_to(REPO)}`

## Definitions

- Persistent genes: PanRGP `partition == persistent`.
- Variable/accessory genes: PanRGP `partition == shell` or `partition == cloud`.
- RGP genes per Mb: number of F1C1 PanRGP gene rows assigned to RGP IDs retained in `spot_merged_with_rgps.tsv`, divided by replicon length in Mb.
- Accessory genes per Mb: shell + cloud genes divided by replicon length in Mb.
- Mobile-element-associated genes: CDS in the F1C1 GFF whose product contains one of these keyword stems: {', '.join(MOBILE_KEYWORDS)}.
- Type III effector genes: ABRicate RSSC_T3Ev3 hits counted by contig; this is hit-based.

## Fisher exact test

Table orientation:

| Replicon | Variable genes | Persistent genes |
|---|---:|---:|
| Chromosome | {a} | {b} |
| Megaplasmid | {c} | {d} |

- Odds ratio, chromosome vs megaplasmid: {odds_ratio_chr_over_mp:.4g}
- Odds ratio, megaplasmid vs chromosome: {odds_ratio_mp_over_chr:.4g}
- Fisher exact two-sided p-value: {p_value:.4g}

Interpretation: the megaplasmid is enriched for variable genes relative to the chromosome when using shell+cloud versus persistent gene counts.
"""
    (OUTDIR / "README.md").write_text(readme)

    # Markdown report table for direct manuscript/revision use.
    with (OUTDIR / "chromosome_vs_megaplasmid_metrics.md").open("w") as handle:
        handle.write("# Chromosome vs Megaplasmid Metrics\n\n")
        handle.write("| Metric | Chromosome | Megaplasmid |\n")
        handle.write("|---|---:|---:|\n")
        for metric in metric_order:
            chrom = metrics[metric]["Chromosome"]
            mega = metrics[metric]["Megaplasmid"]
            if metric in {
                "Percentage of genes classified as variable",
                "Percentage of replicon covered by RGPs",
            }:
                chrom = f"{float(chrom):.2f}%"
                mega = f"{float(mega):.2f}%"
            elif metric == "Total RGP length":
                chrom = f"{int(chrom):,} bp"
                mega = f"{int(mega):,} bp"
            elif metric in {"RGP sequence per Mb"}:
                chrom = f"{float(chrom):,.0f} bp/Mb"
                mega = f"{float(mega):,.0f} bp/Mb"
            elif metric in {"RGP genes per Mb", "Accessory genes per Mb"}:
                chrom = f"{float(chrom):.1f}"
                mega = f"{float(mega):.1f}"
            elif isinstance(chrom, int):
                chrom = fmt_int(chrom)
                mega = fmt_int(int(mega))
            handle.write(f"| {metric} | {chrom} | {mega} |\n")

        handle.write("\n## Fisher exact test\n\n")
        handle.write("| Replicon | Variable genes | Persistent genes |\n")
        handle.write("|---|---:|---:|\n")
        handle.write(f"| Chromosome | {a:,} | {b:,} |\n")
        handle.write(f"| Megaplasmid | {c:,} | {d:,} |\n\n")
        handle.write(f"- Odds ratio, megaplasmid vs chromosome: {odds_ratio_mp_over_chr:.3f}\n")
        handle.write(f"- Fisher exact two-sided p-value: {p_value:.3e}\n")

    print(f"Wrote {OUTDIR / 'chromosome_vs_megaplasmid_metrics.tsv'}")
    print(f"Wrote {OUTDIR / 'chromosome_vs_megaplasmid_metrics.md'}")
    print(f"Wrote {OUTDIR / 'fisher_variable_vs_persistent.tsv'}")


if __name__ == "__main__":
    main()
