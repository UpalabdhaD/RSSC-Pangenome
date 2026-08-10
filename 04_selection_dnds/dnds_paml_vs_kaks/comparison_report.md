# PAML M0 versus prior KaKs_Calculator NG comparison

## Technical summary

The PAML run is not fully complete: 9,187 of 9,624 eligible gene families have valid completed summaries (95.5%). Four families failed preparation and 433 were not reached. The comparison below therefore uses only the 6,271 families present in both the valid PAML set and the prior Ka/Ks table.

The two methods do not produce interchangeable estimates. PAML M0 gives a lower median omega (0.0955) than the prior Ka/Ks values (0.2437), and identifies far fewer families as positively selected: 121/6,271 (1.93%) versus 1,374/6,271 (21.9%). The rank relationship is moderate (Spearman rho = 0.481), but raw Pearson correlation is close to zero because of extreme values and boundary estimates. Selection-class agreement is 78.9%, while Cohen's kappa is only 0.084, indicating weak agreement beyond the heavily imbalanced purifying class.

## Scope and metric definitions

- PAML input: `MS_REVISION/DNDS_ALT_MODELS/results_all_paml_m0/`, codeml one-ratio M0 (`model = 0`, `NSsites = 0`).
- Prior method: `Results/DNDS/RSSC_kaks2_df.csv`, generated from the first data row of each `PAL2NAL_DEDUP_KAKS/*.kaks.tsv` file using KaKs_Calculator NG.
- Comparison unit: one gene family, joined by the exact family identifier.
- Positive selection threshold: PAML `omega > 1`; prior method `Ka/Ks > 1`.
- Matched analysis set: 6,271 gene families with finite, non-negative values from both methods.

## Completion and data-quality checks

| Check | Result |
|---|---:|
| Eligible PAML families | 9,624 |
| Completed PAML statuses | 9,187 |
| Valid PAML summaries | 9,187 |
| Preparation failures | 4 |
| Not run / no status | 433 |
| Prior Ka/Ks valid unique families | 8,876 |
| Matched families | 6,271 |
| PAML omega lower-boundary values (<=0.0001) | 277 |
| PAML omega upper-boundary values (>=999) | 7 |

The PAML summary validation found no invalid summaries among completed statuses. Boundary values should be treated as censored/optimization-bound results rather than precise point estimates, especially when interpreting the largest discrepancies.

## Method comparison

| Measure | Result |
|---|---:|
| PAML median omega | 0.0955 |
| Prior Ka/Ks median | 0.2437 |
| PAML positive fraction | 1.93% |
| Prior Ka/Ks positive fraction | 21.9% |
| Spearman correlation | 0.481 |
| Pearson correlation on raw values | 0.059 |
| Pearson correlation on log10 values, positive pairs only | 0.314 |
| Median absolute raw difference | 0.1015 |
| Threshold agreement | 78.9% |
| Cohen's kappa | 0.084 |

The confusion table is:

| | Prior Ka/Ks <= 1 | Prior Ka/Ks > 1 |
|---|---:|---:|
| PAML omega <= 1 | 4,863 | 1,287 |
| PAML omega > 1 | 34 | 87 |

Thus, most disagreements are prior-method positive / PAML purifying calls (1,287 families), not the reverse (34 families). This supports a conservative interpretation of PAML M0 for positive-selection screening in this dataset.

## Why the estimates differ

PAML M0 is a phylogeny-aware one-ratio estimate fitted across the family using a tree and codon model. The prior workflow used KaKs_Calculator NG on an AXT record assembled by concatenating all deduplicated family sequences into one record. This difference in model, tree usage, sequence representation, and handling of among-lineage variation is sufficient to produce systematic shifts; the two columns should not be treated as replicates of the same estimator.

## Recommended interpretation

Use the PAML M0 results for the primary family-level comparison only after the remaining 433 families are completed or explicitly excluded from the target gene universe. For manuscript claims about positive selection, report the method and threshold explicitly and avoid combining the two methods into one positive-selection list. Treat PAML values at 0.0001 or 999 as boundary estimates and review their `mlc`/alignment diagnostics before making gene-specific claims.

## Reproducibility outputs

- [comparison metrics](comparison_metrics.tsv)
- [PAML completion audit](paml_completion_audit.tsv)
- [matched family-level data](paml_vs_kaks_matched.tsv)
- [selection confusion table](selection_confusion.tsv)
- [largest discrepancies](largest_method_discrepancies.tsv)
- [comparison figure](plots/paml_vs_kaks_comparison.png)
- [reproducible analysis script](scripts/compare_paml_kaks.py)
- [analysis notebook](notebooks/paml_vs_kaks_analysis.ipynb)
