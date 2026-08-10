# Reviewer Response: Ka/Ks Interpretation in Accessory Genes

We thank the reviewer for raising this important point. We agree that elevated
Ka/Ks values in recently acquired or highly variable accessory genes must be
interpreted cautiously, because low synonymous divergence and poor codon
alignments can inflate omega estimates. We have therefore added a sensitivity
analysis to test whether the observed enrichment of genes with Ka/Ks > 1 in the
`Absent in large lineages` category is robust to these potential artifacts.

Specifically, we re-analysed the existing KaKs_Calculator NG outputs after
applying filters designed to remove unstable or low-quality estimates. These
filters excluded genes with low synonymous divergence (`Ks < 0.05`), short
codon alignments (<300 nt), unequal alignment lengths, high gap content (>20%),
high ambiguous-base content (>5%), and alignments containing internal stop
codons. We also performed an additional saturated-divergence check by excluding
genes with `Ks > 2`. The full analysis and plots are provided in
`MS_REVISION/DNDS_SENSITIVITY/`. For readability, the main summary plot groups
the five codon-alignment checks into a single `Alignment QC` category; a
separate detailed plot and table report each artifact check individually.

The elevated signal in the `Absent in large lineages` category remained stable
under these artifact-control filters. In the original analysis, 90/106 genes
(84.9%) in this category had Ka/Ks > 1. After removing genes with `Ks < 0.05`,
the proportion was unchanged at 90/106 (84.9%). Considering each alignment
artifact check separately, the retained positive fractions remained high:
71/83 (85.5%) after excluding short alignments, 90/106 (84.9%) after requiring
equal alignment lengths, 90/106 (84.9%) after applying the gap-content filter,
90/106 (84.9%) after applying the ambiguous-base filter, and 90/106 (84.9%)
after excluding alignments with internal stop codons. After applying the
combined `Ks >= 0.05` and codon-alignment QC filters, 71/83 genes (85.5%)
still had Ka/Ks > 1. After additionally excluding potentially saturated
comparisons (`Ks > 2`), 66/70 genes (94.3%) had Ka/Ks > 1. These results
indicate that the signal is not driven primarily by low synonymous divergence
or obvious alignment-quality artifacts.

We also evaluated a stricter prevalence-based filter requiring at least four
deduplicated sequences per gene. This filter reduced the `Absent in large
lineages` set to only three genes, so we do not treat it as a primary artifact
filter for rare accessory genes; instead, we report it separately as a
sample-size caveat. We have revised the text to state that the Ka/Ks pattern is
consistent with diversifying selection in these variable compartments, while
retaining caution that omega estimates for rare and recently acquired genes
should not be interpreted as definitive evidence for positive selection at
individual loci without further model-based confirmation.

In response to the reviewer's suggestion, we have also prepared an alternative
model-based framework using HyPhy FEL and PAML M0 for targeted validation. These
models address a different question from the pairwise NG estimates: FEL tests
site-wise synonymous and nonsynonymous rates, whereas PAML M0 estimates a
single gene-wide omega under an explicit codon-substitution model. We consider
these appropriate complementary checks, but we have kept the revised manuscript
claim focused on the robust class-level sensitivity result rather than
overstating individual-gene inference.
