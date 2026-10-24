# 03_genomic_plasticity_and_defense

Corresponds to STAR Methods > Method details, "Genomic plasticity regions, integration spots, and
specialized gene systems" and "Co-occurrence analysis of accessory gene families."

## What's here

- `subworkflows/RGP_annotation.smk`, `subworkflows/ppanggolin.smk` — PPanGGOLiN's
  PanRGP module (v1.2.105), detecting regions of genomic plasticity (RGPs) across the
  142-genome pangenome and collapsing them into the 651 conserved integration spots.
- `scripts/process_RGP_SPOTS.sh`, `scripts/process_rgps.sh`, `scripts/rgp_getfasta.sh`,
  `scripts/merge_ppanggolin_tables.sh` — RGP/spot post-processing and FASTA extraction.
- `subworkflows/abricate_virulence_finder.smk`, `envs/abricate.yaml` — Abricate v1.0.1
  against CARD, VFDB, the custom RalstoT3E (T3SS) database, and (via
  `annotate_secretion_systems.smk`) the SecReT6 (T6SS) database.
- `subworkflows/annotate_secretion_systems.smk` — T3SS/T6SS-specific secretion-system
  annotation.
- `subworkflows/antismash.smk`, `envs/antismash.yaml`, `scripts/antismash.sh`,
  `scripts/combine_antismash.py` — antiSMASH v7.0.0 secondary-metabolite BGC detection
  (`--taxon bacteria --minimal`).
- `subworkflows/*` + `envs/padloc.yaml`, `scripts/merge_PADLOC.sh`,
  `scripts/padloc_stats.py`, `padloc.patch`, `rssc_defense.qmd` — PADLOC antiviral
  defense-system detection and per-spot defense-island summarization (Wadjet,
  CRISPR-Cas, etc).
- `scripts/plot_spot_gene_family_distribution.R` — the per-spot gene-family-count
  distribution (Figure 5-adjacent).
- `scripts/RGP.Rmd` — RGP/spot exploratory notebook.
- `subworkflows/cooccurance.smk`, `envs/coinfinder.tested.yaml` — Coinfinder v1.2.1 in
  association mode: gene families that co-occur across the phylogeny more often than
  expected, and the connected components of the resulting network.
- `figure_s11/redraw_supp_fig11.R` (added for the revision) — Figure S11: spots on each
  replicon class, genes per RGP and RGP length, for widespread spots (an RGP in more than
  100 of the 142 genomes) and other spots, with two-sided Wilcoxon rank-sum tests (ggpubr
  `stat_compare_means`).
- `figure_s17/redraw_supp_fig17.R` (added for the revision) — Figure S17: distinct gene
  families per spot against distinct defense-associated gene families per spot, with
  Spearman's rank correlation and a least-squares line.

