# “Absent in large lineages”: KaKs_Calculator NG versus PAML M0

## Scope

This is the category-specific comparison requested for the gene families labelled **“Absent in large lineages”** in the existing Ka/Ks sensitivity data. “NG” refers to the previously used KaKs_Calculator NG estimate, while “M0” refers to the PAML codeml one-ratio estimate (`omega`). The two quantities are related but are not methodologically interchangeable: NG is the prior family-level Ka/Ks calculation, whereas PAML M0 estimates one omega across the family using the supplied phylogenetic model.

## Coverage and data-quality audit

The category contains 151 rows in the sensitivity input. Of these, 106 have valid numeric NG Ka/Ks values, which is the denominator underlying the previously reported 90/106 (84.9%) result. Most category members cannot be assessed by the current PAML run because they are singletons.

| Population | Families | Share |
|---|---:|---:|
| Category rows in sensitivity input | 151 | 100.0% |
| Valid NG Ka/Ks values | 106 | 70.2% of category |
| Singleton families (`n = 1`) | 144 | 95.4% of category |
| Families with at least two sequences | 7 | 4.6% of category |
| PAML-eligible families | 7 | 100.0% of multi-sequence families |
| Completed and usable PAML M0 summaries | 6 | 85.7% of PAML-eligible families |
| Matched valid NG + usable M0 values | 5 | 4.7% of valid NG families |

One PAML-eligible family (`pilV`) has not yet been run. One completed PAML family (`fabF`) has an invalid NG value and is therefore excluded from the direct NG–M0 comparison. The 144 singleton families are not represented in the matched comparison because the current PAML input was restricted to families with at least two sequences.

## Direct comparison of matched families

| Measure | Previous NG | PAML M0 | Matched families |
|---|---:|---:|---:|
| Median ratio | 0.602499 | 0.054480 | 5 |
| Ratio > 1 | 1/5 (20.0%) | 0/5 (0.0%) | 5 |
| Ratio <= 1 | 4/5 (80.0%) | 5/5 (100.0%) | 5 |
| Median absolute raw difference | 0.540139 | — | 5 |
| Median log10 difference (M0 − NG) | -1.617125 | — | 5 |

The family-level values are:

| Gene family | Sequences | NG Ka/Ks | PAML M0 omega | NG > 1 | M0 > 1 |
|---|---:|---:|---:|---:|---:|
| `ecfA2` | 27 | 2.256160 | 0.054480 | yes | no |
| `fabD` | 38 | 0.521248 | 0.054790 | no | no |
| `merP` | 2 | 0.441198 | 0.002650 | no | no |
| `merT` | 2 | 0.602530 | 0.002830 | no | no |
| `relE` | 26 | 0.602499 | 0.062360 | no | no |

The two methods agree on the >1 classification for 4/5 families (80.0%). Cohen’s kappa is 0.0 because the matched sample contains no M0-positive family and is too small and imbalanced for the agreement statistic to provide useful evidence. The raw Pearson correlation is 0.381 (`p = 0.526`), and the Spearman correlation is 0.200 (`p = 0.747`); neither is informative at this sample size.

## Interpretation

For the small matched subset, the prior NG calculation identifies one positive family (`ecfA2`), whereas PAML M0 identifies none. PAML values are lower for all five matched families. This is consistent with the broader observation that the previous NG estimates can produce more Ka/Ks > 1 calls than PAML M0, but the category-specific comparison is **not sufficiently complete or representative to validate or refute the 84.9% category-level estimate**. In particular, 5 matched families represent only 4.7% of the 106 valid NG families, and 95.4% of the category consists of singletons that are not PAML-assessable in the current run.

Accordingly, the defensible conclusion is that the high NG fraction in this category should not be described as PAML-confirmed evidence of an arms-race process. It remains a result consistent with diversifying selection, but recent acquisition, low synonymous divergence, alignment/model sensitivity, and the extreme singleton composition remain important alternative explanations. A category-wide PAML comparison would require a model appropriate for singleton or otherwise unreplicated families; the current M0 output cannot supply that validation.

## Reproducible outputs

- [Category coverage audit](category_coverage.tsv)
- [Matched NG–M0 values](paml_m0_vs_ng_matched.tsv)
- [Selection-threshold confusion table](selection_confusion.tsv)
- [Comparison metrics](comparison_metrics.tsv)
- [Analysis script](../scripts/compare_absent_major_lineages.py)

The underlying inputs are [`kaks_sensitivity_joined.tsv`](../../DNDS_SENSITIVITY/results/kaks_sensitivity_joined.tsv) and [`paml_m0_aggregated.tsv`](../paml_m0_aggregated.tsv).
