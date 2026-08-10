# DNDS Sensitivity Analysis

This directory contains a shell-only sensitivity check for the manuscript claim
that accessory/integration-spot genes show elevated Ka/Ks.

The analysis reuses existing outputs:

- `Results/DNDS/PAL2NAL_DEDUP_KAKS/*.kaks.tsv`
- `Results/DNDS/PAL2NAL_DEDUP/*.codon.dedup.fa`
- `Manuscript_Figures/Figure_2/Data/twilight_classification_updated.tab`

## Run

Use the current environment if Python 3 and R/ggplot2 are already available:

```bash
bash MS_REVISION/DNDS_SENSITIVITY/scripts/run_dnds_sensitivity.sh
```

Or create and use the minimal conda environment:

```bash
mamba env create -f MS_REVISION/DNDS_SENSITIVITY/envs/dnds_sensitivity.yml
bash MS_REVISION/DNDS_SENSITIVITY/scripts/run_dnds_sensitivity.sh dnds_sensitivity
```

## Outputs

- `results/kaks_sensitivity_joined.tsv`: full joined Ka/Ks, metadata, and alignment QC table.
- `results/sensitivity_by_specific_class.tsv`: positive Ka/Ks fraction by pangenome class and filter.
- `results/sensitivity_by_general_class.tsv`: positive Ka/Ks fraction by general class and filter.
- `results/DNDS_sensitivity_summary.md`: short reviewer-facing numerical summary.
- `plots/absent_large_lineages_sensitivity.png`: focused check for the high-risk class.
- `plots/absent_large_lineages_individual_artifact_checks.png`: focused check showing each artifact filter separately.
- `plots/positive_fraction_by_specific_class.png`: all specific-class sensitivity comparison.
- `plots/strict_kaks_distribution_by_specific_class.png`: distribution after strict filtering.

## Interpretation

The key manuscript-safe comparison is the baseline versus `strict_ks0.05_qc`
or `strict_ks0.05_2_qc` result. These remove low-Ks genes, short/poor
alignments, and optionally saturated Ks values. The `min4_*` filters are
reported separately because requiring at least four deduplicated sequences can
exclude genuinely rare accessory genes and should not be treated as a pure
alignment-artifact filter.

The main focused plot groups the alignment checks into `Alignment QC` for
readability. The individual-filter plot separates all checks: `Ks >= 0.05`,
alignment length >=300 nt, equal alignment lengths, gap fraction <=20%,
ambiguous-base fraction <=5%, no internal stops, and `Ks <= 2`.
