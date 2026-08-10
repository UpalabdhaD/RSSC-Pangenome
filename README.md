# RSSC-completeGenome-analysis

Companion code for the comparative-genomics pipeline described in the F1C1 genome
manuscript (*Integration spots organize genome plasticity and shape accessory-genome
architecture in Ralstonia solanacearum*), Methods §5.8–§5.12. Organized by pipeline
stage, starting from pangenome construction through to functional enrichment.

## Scope

This directory covers the population-genomic analysis (142 RSSC genomes), not the
F1C1 assembly itself. Excluded on purpose:

- **Methods §5.1–§5.7** (culture, DNA extraction, sequencing, hybrid assembly,
  assembly QC, single-genome annotation, species-level ANI/Mash placement) — these
  precede pangenome construction and concern the F1C1 genome alone, not the 142-genome
  comparative analysis this directory documents.
- **Tools tried during method development but not named in the manuscript**: PEPPAN,
  Roary, PopPUNK, PGGB, Minigraph, PanGraph, panX, FastTree, Panstripe. These are
  active or present elsewhere in the lab's working `workflow/` directory but aren't
  part of the published pipeline.
- **Orphaned draft subworkflows** (`pangenome2.smk`, `pangenome_3.smk`, etc.) with no
  live reference anywhere in the project's Snakefile.
- **The RNASeq/differential-expression scripts** (`DEG.Rmd`, `Downstream_analysis.Rmd`)
  — a separate assay, unrelated to this manuscript.

## Stages

| Directory | Methods section | Covers |
|---|---|---|
| `01_pangenome/` | §5.8 | NCBI genome retrieval, Panaroo pangenome construction, Pagoo/micropan accumulation curves, twilight lineage-aware partitioning, FastBAPS clustering |
| `02_phylogeny/` | §5.9 | FastANI-based population clustering (networkx/Cytoscape), core-gene SNP-sites + IQ-TREE/ModelFinder phylogeny, Grapple accessory-distance clustering |
| `03_genomic_plasticity_and_defense/` | §5.10 | PPanGGOLiN PanRGP integration spots, Abricate (CARD/VFDB/RalstoT3E T3SS/SecReT6 T6SS), antiSMASH BGCs, PADLOC antiviral defense systems |
| `04_selection_dnds/` | §5.11 | KaKs_Calculator2 (NG method) Ka/Ks estimation pipeline: seqtk, MAFFT L-INS-i, PAL2NAL, seqkit, custom AXT conversion. Also holds peer-review-driven extensions not yet in the manuscript's Methods text: alternative dN/dS models (PAML M0, HyPhy FEL/FITMG94), a Ka/Ks sensitivity analysis, a PAML-vs-KaKs_Calculator2 comparison, HyPhy BUSTED episodic-selection tests, and replicon-variable-RGP metrics (chromosome vs. megaplasmid) — see the README in that directory |
| `05_functional_enrichment/` | §5.12 | Pangenome-wide eggNOG-mapper annotation feeding the COG over-representation analysis (the analysis itself is carried out in `04_selection_dnds/scripts/RSSC_dnds.Rmd`, alongside Ka/Ks, since Figure 4 combines both) |
