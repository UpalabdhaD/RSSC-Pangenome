# 03_genomic_plasticity_and_defense

Corresponds to Methods §5.10, "Genomic plasticity regions, integration spots, and
specialized gene systems."

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

## Flag before wider release

`subworkflows/cooccurance.smk`, `envs/coinfinder.yaml` — this runs the gene
co-occurrence / connected-components analysis behind Results §3.4 and Figure 4D
("Identification of coincident accessory genes revealed a non-random modular
structure... Network reconstruction identified connected components"). The
manuscript's Methods text never names Coinfinder or describes this step, so a reader
can't currently reproduce Figure 4D from the paper alone — worth adding a Methods
sentence if this repo ships alongside the paper.
