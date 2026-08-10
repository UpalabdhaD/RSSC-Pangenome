# 01_pangenome

Corresponds to Methods §5.8, "Pangenome construction and population-informed
partitioning."

## What's here

- `scripts/download_and_parse_ncbi_genome.sh` — retrieves the 158 "complete"-status
  Ralstonia genome assemblies from NCBI GenBank via the `datasets` CLI (142 were
  retained after filtering to RSSC species).
- `scripts/raname_gbk_assemblies.sh` — renames the downloaded assemblies by strain ID
  ahead of annotation and pangenome construction.
- `subworkflows/panaroo.smk` — Panaroo pangenome construction (95% protein similarity
  threshold, `--clean-mode sensitive`, MAFFT v7.505 for core-gene alignment).
- `subworkflows/postprocess_panaroo.smk` — pangenome accumulation curves and PCA via
  the Pagoo R package; alpha parameter via micropan.
- `scripts/RSSC-pagoo.Rmd` — the Pagoo accumulation-curve/PCA notebook.
- `subworkflows/partitionPangenome.smk` and `scripts/twilight.R` — the
  population-structure-informed partition using `twilight` (core/rare thresholds
  0.95/0.15).
- `envs/fastbaps.yaml` — FastBAPS hierarchical Bayesian clustering, used to define the
  population groups that twilight partitions against.
- `scripts/twilight_stat_anslysis.Rmd`, `scripts/twilight_fastbaps_cluster_sankey.Rmd`,
  `scripts/backstrace_twilight.Rmd` — downstream twilight/FastBAPS analysis and
  reporting notebooks.

## Not included

Roary, PEPPAN, panX, PPanGGOLiN-as-pangenome-caller (PPanGGOLiN is used later for RGP
detection, §5.10, not pangenome construction itself), PopPUNK, and the pangenome
graph tools (PGGB, Minigraph, PanGraph) — none of these are named as the pangenome
method in the manuscript, which specifies Panaroo.
