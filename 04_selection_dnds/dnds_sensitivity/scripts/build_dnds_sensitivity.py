#!/usr/bin/env python3
"""Build Ka/Ks sensitivity tables from existing pairwise KaKs outputs.

This script is intentionally read-only with respect to the original project
outputs. It parses Results/DNDS/PAL2NAL_DEDUP_KAKS/*.kaks.tsv, joins the
twilight class metadata, adds codon-alignment QC metrics, then reports how
the fraction of genes with Ka/Ks > 1 changes under conservative filters.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


FILTERS = {
    "baseline_numeric": lambda r: r["valid_numeric"],
    "ks_ge_0.01": lambda r: r["valid_numeric"] and r["ks"] >= 0.01,
    "ks_ge_0.05": lambda r: r["valid_numeric"] and r["ks"] >= 0.05,
    "alignment_len_ge300": lambda r: r["valid_numeric"]
    and r["alignment_nt_length"] >= 300,
    "equal_alignment_lengths": lambda r: r["valid_numeric"] and r["unequal_lengths"] == 0,
    "gap_le_0.20": lambda r: r["valid_numeric"] and r["gap_fraction"] <= 0.20,
    "ambiguous_le_0.05": lambda r: r["valid_numeric"]
    and r["ambiguous_fraction"] <= 0.05,
    "no_internal_stops": lambda r: r["valid_numeric"] and r["stop_codon_count"] == 0,
    "ks_le_2": lambda r: r["valid_numeric"] and r["ks"] <= 2.0,
    "qc_alignment": lambda r: (
        r["valid_numeric"]
        and r["alignment_nt_length"] >= 300
        and r["unequal_lengths"] == 0
        and r["gap_fraction"] <= 0.20
        and r["ambiguous_fraction"] <= 0.05
        and r["stop_codon_count"] == 0
    ),
    "min4_sequences": lambda r: r["valid_numeric"] and r["n_sequences"] >= 4,
    "strict_ks0.05_qc": lambda r: (
        r["valid_numeric"]
        and r["ks"] >= 0.05
        and r["alignment_nt_length"] >= 300
        and r["unequal_lengths"] == 0
        and r["gap_fraction"] <= 0.20
        and r["ambiguous_fraction"] <= 0.05
        and r["stop_codon_count"] == 0
    ),
    "strict_ks0.05_2_qc": lambda r: (
        r["valid_numeric"]
        and 0.05 <= r["ks"] <= 2.0
        and r["alignment_nt_length"] >= 300
        and r["unequal_lengths"] == 0
        and r["gap_fraction"] <= 0.20
        and r["ambiguous_fraction"] <= 0.05
        and r["stop_codon_count"] == 0
    ),
    "min4_strict_ks0.05_qc": lambda r: (
        r["valid_numeric"]
        and r["ks"] >= 0.05
        and r["n_sequences"] >= 4
        and r["alignment_nt_length"] >= 300
        and r["unequal_lengths"] == 0
        and r["gap_fraction"] <= 0.20
        and r["ambiguous_fraction"] <= 0.05
        and r["stop_codon_count"] == 0
    ),
}


def parse_float(value: str) -> float:
    value = (value or "").strip()
    if value in {"", "NA", "NaN", "nan", "Inf", "inf", "-Inf", "-inf"}:
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def is_finite(value: float) -> bool:
    return not math.isnan(value) and not math.isinf(value)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records


def alignment_qc(path: Path) -> dict[str, float | int | str]:
    if not path.exists():
        return {
            "alignment_file": str(path),
            "n_sequences": 0,
            "alignment_nt_length": 0,
            "unequal_lengths": 1,
            "gap_fraction": math.nan,
            "ambiguous_fraction": math.nan,
            "stop_codon_count": 0,
        }

    records = read_fasta(path)
    lengths = [len(seq) for _, seq in records]
    total = sum(lengths)
    gaps = sum(seq.count("-") for _, seq in records)
    ambiguous = sum(sum(1 for base in seq if base not in "ACGT-") for _, seq in records)
    stops = 0
    for _, seq in records:
        ungapped = seq.replace("-", "")
        for i in range(0, len(ungapped) - 2, 3):
            if ungapped[i : i + 3] in {"TAA", "TAG", "TGA"}:
                stops += 1
    return {
        "alignment_file": str(path),
        "n_sequences": len(records),
        "alignment_nt_length": max(lengths) if lengths else 0,
        "unequal_lengths": int(len(set(lengths)) > 1),
        "gap_fraction": gaps / total if total else math.nan,
        "ambiguous_fraction": ambiguous / total if total else math.nan,
        "stop_codon_count": stops,
    }


def load_classification(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["gene_name"]
            out[gene] = {
                "general_class": row.get("general_class", ""),
                "specific_class": row.get("specific_class", ""),
            }
    return out


def parse_kaks_outputs(kaks_dir: Path, align_dir: Path, classes: dict[str, dict[str, str]]) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(kaks_dir.glob("*.kaks.tsv")):
        raw_gene = path.name.removesuffix(".kaks.tsv")
        gene = raw_gene.replace("___", "~~~")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                # Most files have one pairwise/consensus line after deduplication.
                ka = parse_float(row.get("Ka", ""))
                ks = parse_float(row.get("Ks", ""))
                kaks = parse_float(row.get("Ka/Ks", ""))
                pvalue = parse_float(row.get("P-Value(Fisher)", ""))
                length = parse_float(row.get("Length", ""))
                s_sites = parse_float(row.get("S-Sites", ""))
                n_sites = parse_float(row.get("N-Sites", ""))
                s_substitutions = parse_float(row.get("S-Substitutions", ""))
                n_substitutions = parse_float(row.get("N-Substitutions", ""))
                pnps_proxy = math.nan
                if (
                    is_finite(s_substitutions)
                    and is_finite(n_substitutions)
                    and is_finite(s_sites)
                    and is_finite(n_sites)
                    and s_substitutions > 0
                    and s_sites > 0
                    and n_sites > 0
                ):
                    pnps_proxy = (n_substitutions / n_sites) / (s_substitutions / s_sites)
                qc = alignment_qc(align_dir / f"{raw_gene}.codon.dedup.fa")
                class_row = classes.get(gene, {})
                general = class_row.get("general_class", "Unclassified") or "Unclassified"
                specific = class_row.get("specific_class", "Unclassified") or "Unclassified"
                rows.append(
                    {
                        "gene_raw": raw_gene,
                        "gene": gene,
                        "sequence_pair": row.get("Sequence", ""),
                        "method": row.get("Method", ""),
                        "ka": ka,
                        "ks": ks,
                        "kaks": kaks,
                        "pvalue_fisher": pvalue,
                        "kaks_length": length,
                        "s_sites": s_sites,
                        "n_sites": n_sites,
                        "s_substitutions": s_substitutions,
                        "n_substitutions": n_substitutions,
                        "pnps_proxy": pnps_proxy,
                        "general_class": general,
                        "specific_class": specific,
                        "core_accessory": "Core" if "Core" in general else "Accessory",
                        "valid_numeric": is_finite(kaks) and is_finite(ka) and is_finite(ks),
                        "positive_kaks_gt1": is_finite(kaks) and kaks > 1.0,
                        "positive_pnps_proxy_gt1": is_finite(pnps_proxy) and pnps_proxy > 1.0,
                        **qc,
                    }
                )
    return rows


def write_table(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize(rows: list[dict], class_field: str) -> list[dict]:
    out: list[dict] = []
    for filter_name, filter_fn in FILTERS.items():
        groups: dict[str, list[dict]] = defaultdict(list)
        for row in rows:
            if filter_fn(row):
                groups[row[class_field]].append(row)
        for class_value, group_rows in sorted(groups.items()):
            total = len(group_rows)
            positive = sum(1 for row in group_rows if row["positive_kaks_gt1"])
            ks_low = sum(1 for row in group_rows if row["ks"] < 0.05)
            bad_qc = sum(
                1
                for row in group_rows
                if row["gap_fraction"] > 0.20
                or row["ambiguous_fraction"] > 0.05
                or row["stop_codon_count"] > 0
                or row["alignment_nt_length"] < 300
                or row["unequal_lengths"] != 0
            )
            proxy_total = sum(1 for row in group_rows if is_finite(row["pnps_proxy"]))
            proxy_positive = sum(1 for row in group_rows if row["positive_pnps_proxy_gt1"])
            out.append(
                {
                    "filter": filter_name,
                    "class_field": class_field,
                    "class": class_value,
                    "n": total,
                    "positive_n": positive,
                    "positive_fraction": positive / total if total else math.nan,
                    "positive_percent": 100 * positive / total if total else math.nan,
                    "low_ks_lt0.05_n": ks_low,
                    "qc_flag_n": bad_qc,
                    "pnps_proxy_n": proxy_total,
                    "pnps_proxy_gt1_n": proxy_positive,
                    "pnps_proxy_gt1_percent": 100 * proxy_positive / proxy_total
                    if proxy_total
                    else math.nan,
                }
            )
    return out


def write_markdown(path: Path, rows: list[dict], specific_summary: list[dict]) -> None:
    target = "Absent in large lineages"
    target_rows = [
        row for row in specific_summary if row["class"] == target and row["filter"] in FILTERS
    ]
    by_filter = {row["filter"]: row for row in target_rows}
    total_numeric = sum(1 for row in rows if row["valid_numeric"])
    total_positive = sum(1 for row in rows if row["valid_numeric"] and row["positive_kaks_gt1"])
    low_ks_positive = sum(
        1 for row in rows if row["valid_numeric"] and row["positive_kaks_gt1"] and row["ks"] < 0.05
    )

    def fmt(filter_name: str) -> str:
        row = by_filter.get(filter_name)
        if not row:
            return "not available"
        return f"{row['positive_percent']:.1f}% ({row['positive_n']}/{row['n']})"

    def retained(filter_name: str) -> str:
        row = by_filter.get(filter_name)
        if not row:
            return "not available"
        return f"{row['n']} retained; {row['positive_percent']:.1f}% positive"

    text = f"""# DNDS Sensitivity Summary

This analysis starts from the existing KaKs_Calculator NG outputs in `Results/DNDS/PAL2NAL_DEDUP_KAKS/` and joins pangenome class labels from `Manuscript_Figures/Figure_2/Data/twilight_classification_updated.tab`.

## Main sensitivity checks

- Baseline numeric Ka/Ks rows: {total_numeric}
- Baseline positive rows (`Ka/Ks > 1`): {total_positive}
- Positive rows with low synonymous divergence (`Ks < 0.05`): {low_ks_positive}
- `Absent in large lineages` baseline positive fraction: {fmt('baseline_numeric')}
- `Absent in large lineages` after `Ks >= 0.05`: {fmt('ks_ge_0.05')}
- `Absent in large lineages` after strict `Ks >= 0.05` plus alignment QC: {fmt('strict_ks0.05_qc')}
- `Absent in large lineages` after strict QC and excluding saturated `Ks > 2`: {fmt('strict_ks0.05_2_qc')}
- `Absent in large lineages` after additionally requiring at least 4 deduplicated sequences: {fmt('min4_strict_ks0.05_qc')}

## Individual artifact checks for `Absent in large lineages`

- Low synonymous divergence filter (`Ks >= 0.05`): {retained('ks_ge_0.05')}
- Short-alignment filter (`alignment length >= 300 nt`): {retained('alignment_len_ge300')}
- Unequal-length filter: {retained('equal_alignment_lengths')}
- High-gap filter (`gap fraction <= 20%`): {retained('gap_le_0.20')}
- High-ambiguity filter (`ambiguous fraction <= 5%`): {retained('ambiguous_le_0.05')}
- Internal-stop filter: {retained('no_internal_stops')}
- Saturated-divergence filter (`Ks <= 2`): {retained('ks_le_2')}

## Filters used

- `ks_ge_0.05`: removes genes where omega is most unstable due to near-zero synonymous divergence.
- `alignment_len_ge300`, `equal_alignment_lengths`, `gap_le_0.20`, `ambiguous_le_0.05`, and `no_internal_stops`: report each alignment artifact check separately.
- `qc_alignment`: keeps alignments with nucleotide alignment length >= 300, equal sequence lengths, gap fraction <= 20%, ambiguous base fraction <= 5%, and no internal stop codons.
- `strict_ks0.05_qc`: combines the low-Ks and alignment-quality filters.
- `strict_ks0.05_2_qc`: additionally removes high-Ks cases that may be saturated.
- `min4_strict_ks0.05_qc`: adds a prevalence/sample-size check and should be interpreted separately for rare categories.

Use `results/sensitivity_by_specific_class.tsv` and `plots/positive_fraction_by_specific_class.png` for manuscript-facing sensitivity reporting.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kaks-dir", required=True, type=Path)
    parser.add_argument("--align-dir", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    classes = load_classification(args.classification)
    rows = parse_kaks_outputs(args.kaks_dir, args.align_dir, classes)

    row_fields = [
        "gene_raw",
        "gene",
        "sequence_pair",
        "method",
        "ka",
        "ks",
        "kaks",
        "pvalue_fisher",
        "kaks_length",
        "s_sites",
        "n_sites",
        "s_substitutions",
        "n_substitutions",
        "pnps_proxy",
        "general_class",
        "specific_class",
        "core_accessory",
        "valid_numeric",
        "positive_kaks_gt1",
        "positive_pnps_proxy_gt1",
        "alignment_file",
        "n_sequences",
        "alignment_nt_length",
        "unequal_lengths",
        "gap_fraction",
        "ambiguous_fraction",
        "stop_codon_count",
    ]
    write_table(args.outdir / "kaks_sensitivity_joined.tsv", rows, row_fields)
    write_table(args.outdir / "sensitivity_by_specific_class.tsv", summarize(rows, "specific_class"), [
        "filter",
        "class_field",
        "class",
        "n",
        "positive_n",
        "positive_fraction",
        "positive_percent",
        "low_ks_lt0.05_n",
        "qc_flag_n",
        "pnps_proxy_n",
        "pnps_proxy_gt1_n",
        "pnps_proxy_gt1_percent",
    ])
    write_table(args.outdir / "sensitivity_by_general_class.tsv", summarize(rows, "general_class"), [
        "filter",
        "class_field",
        "class",
        "n",
        "positive_n",
        "positive_fraction",
        "positive_percent",
        "low_ks_lt0.05_n",
        "qc_flag_n",
        "pnps_proxy_n",
        "pnps_proxy_gt1_n",
        "pnps_proxy_gt1_percent",
    ])
    write_table(args.outdir / "sensitivity_by_core_accessory.tsv", summarize(rows, "core_accessory"), [
        "filter",
        "class_field",
        "class",
        "n",
        "positive_n",
        "positive_fraction",
        "positive_percent",
        "low_ks_lt0.05_n",
        "qc_flag_n",
        "pnps_proxy_n",
        "pnps_proxy_gt1_n",
        "pnps_proxy_gt1_percent",
    ])
    write_markdown(args.outdir / "DNDS_sensitivity_summary.md", rows, summarize(rows, "specific_class"))


if __name__ == "__main__":
    main()
