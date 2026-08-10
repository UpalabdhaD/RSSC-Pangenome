# DNDS Sensitivity Summary

This analysis starts from the existing KaKs_Calculator NG outputs in `Results/DNDS/PAL2NAL_DEDUP_KAKS/` and joins pangenome class labels from `Manuscript_Figures/Figure_2/Data/twilight_classification_updated.tab`.

## Main sensitivity checks

- Baseline numeric Ka/Ks rows: 8599
- Baseline positive rows (`Ka/Ks > 1`): 3492
- Positive rows with low synonymous divergence (`Ks < 0.05`): 82
- `Absent in large lineages` baseline positive fraction: 84.9% (90/106)
- `Absent in large lineages` after `Ks >= 0.05`: 84.9% (90/106)
- `Absent in large lineages` after strict `Ks >= 0.05` plus alignment QC: 85.5% (71/83)
- `Absent in large lineages` after strict QC and excluding saturated `Ks > 2`: 94.3% (66/70)
- `Absent in large lineages` after additionally requiring at least 4 deduplicated sequences: 33.3% (1/3)

## Individual artifact checks for `Absent in large lineages`

- Low synonymous divergence filter (`Ks >= 0.05`): 106 retained; 84.9% positive
- Short-alignment filter (`alignment length >= 300 nt`): 83 retained; 85.5% positive
- Unequal-length filter: 106 retained; 84.9% positive
- High-gap filter (`gap fraction <= 20%`): 106 retained; 84.9% positive
- High-ambiguity filter (`ambiguous fraction <= 5%`): 106 retained; 84.9% positive
- Internal-stop filter: 106 retained; 84.9% positive
- Saturated-divergence filter (`Ks <= 2`): 89 retained; 93.3% positive

## Filters used

- `ks_ge_0.05`: removes genes where omega is most unstable due to near-zero synonymous divergence.
- `alignment_len_ge300`, `equal_alignment_lengths`, `gap_le_0.20`, `ambiguous_le_0.05`, and `no_internal_stops`: report each alignment artifact check separately.
- `qc_alignment`: keeps alignments with nucleotide alignment length >= 300, equal sequence lengths, gap fraction <= 20%, ambiguous base fraction <= 5%, and no internal stop codons.
- `strict_ks0.05_qc`: combines the low-Ks and alignment-quality filters.
- `strict_ks0.05_2_qc`: additionally removes high-Ks cases that may be saturated.
- `min4_strict_ks0.05_qc`: adds a prevalence/sample-size check and should be interpreted separately for rare categories.

Use `results/sensitivity_by_specific_class.tsv` and `plots/positive_fraction_by_specific_class.png` for manuscript-facing sensitivity reporting.
