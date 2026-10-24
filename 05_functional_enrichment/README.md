# 05_functional_enrichment

Corresponds to STAR Methods > Method details, "Functional enrichment analysis."

## What's here

- `subworkflows/annotate_RSSC_eggnog.smk`, `envs/emapper.yaml` — pangenome-wide
  eggNOG-mapper annotation, run per-genome across all 142 RSSC assemblies (output
  path `EGGNOG_RSSC_142/{genome}`), providing the COG category calls this stage's
  enrichment analysis is built on.

## Where the enrichment analyses live

- GO enrichment for Figure 4: `../04_selection_dnds/figure4_paml/Figure4_PAML_M0.Rmd`
  (clusterProfiler `enricher`, Benjamini–Hochberg-adjusted P < 0.05). The original notebook,
  `../04_selection_dnds/scripts/RSSC_dnds.Rmd`, ran the same GO enrichment for the first
  submission.
- COG over-representation across the pangenome classes (Figure 3F and Figure S5): not in this
  repository. The notebook is `Manuscript_Figures/Figure_2/COG_enrichment.Rmd` in the
  analysis project.
