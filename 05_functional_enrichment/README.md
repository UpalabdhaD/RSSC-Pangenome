# 05_functional_enrichment

Corresponds to Methods §5.12, "Functional annotation and enrichment analysis."

## What's here

- `subworkflows/annotate_RSSC_eggnog.smk`, `envs/emapper.yaml` — pangenome-wide
  eggNOG-mapper annotation, run per-genome across all 142 RSSC assemblies (output
  path `EGGNOG_RSSC_142/{genome}`), providing the COG category calls this stage's
  enrichment analysis is built on.

## Where the enrichment analysis itself lives

The over-representation analysis (ClusterProfiler's `enrichR`, comparing core,
intermediate, and rare gene classes) is implemented inside
`../04_selection_dnds/scripts/RSSC_dnds.Rmd`, not as a separate script — the
manuscript's Figure 4 presents Ka/Ks and functional-category enrichment together, and
the underlying notebook does too. This directory documents the annotation half of
§5.12; the analysis half is cross-referenced from `04_selection_dnds`.
