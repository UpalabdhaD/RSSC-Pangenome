# RSSC-completeGenome-analysis

Companion code for the comparative-genomics analysis in the F1C1 genome manuscript
(*Integration spots organize genome plasticity and shape accessory-genome architecture in
Ralstonia solanacearum*, iScience). The folders follow the pipeline from pangenome construction
to functional enrichment. The table under "Stages" links each folder to its subsection under
STAR Methods > Method details.

## Scope

This repository covers the population-genomic analysis of 142 RSSC genomes, not the F1C1
assembly itself. Excluded on purpose:

- **The F1C1 genome steps** (culture, DNA extraction, sequencing, hybrid assembly, assembly
  quality control, single-genome annotation and species-level ANI/Mash placement). These come
  before pangenome construction and concern the F1C1 genome alone.
- **Tools tried during method development but not named in the manuscript**: PEPPAN, Roary,
  PopPUNK, PGGB, Minigraph, PanGraph, panX and Panstripe. They remain in the lab's working
  `workflow/` directory but are not part of the published pipeline.
- **Orphaned draft subworkflows** (`pangenome2.smk`, `pangenome_3.smk`, etc.) with no live
  reference in the project's Snakefile.
- **The RNA-seq and differential-expression scripts** (`DEG.Rmd`, `Downstream_analysis.Rmd`),
  which belong to a separate assay.

## Stages

| Directory | Method details subsection | Covers |
|---|---|---|
| `01_pangenome/` | Pangenome construction and population-informed partitioning | NCBI genome retrieval, Panaroo pangenome construction, Pagoo/micropan accumulation curves, Twilight lineage-aware partitioning, FastBAPS clustering, and the Table S6 accession list |
| `02_phylogeny/` | Phylogenetic reconstruction and population structure analysis | FastANI-based population clustering (networkx/Cytoscape), core-gene SNP-sites + IQ-TREE/ModelFinder phylogeny, Grapple accessory-distance clustering |
| `03_genomic_plasticity_and_defense/` | Genomic plasticity regions, integration spots, and specialized gene systems; Co-occurrence analysis of accessory gene families | PPanGGOLiN panRGP integration spots, ABRicate (CARD/VFDB/RalstoT3E T3SS/SecReT6 T6SS), antiSMASH BGCs, PADLOC antiviral defense systems, Coinfinder co-occurrence, and Figures S11 and S17 |
| `04_selection_dnds/` | Codon-model estimation of ω and episodic diversifying selection | Codon alignments (seqtk, MAFFT L-INS-i, PAL2NAL, seqkit), PAML codeml M0 ω per gene family, HyPhy BUSTED, FastTree gene trees for families with paralogs, Figure 4, Tables S7 and S8, and the chromosome-versus-megaplasmid comparison. Also keeps the KaKs_Calculator2 pipeline of the first submission |
| `05_functional_enrichment/` | Functional enrichment analysis | Pangenome-wide eggNOG-mapper annotation for the COG and GO enrichment analyses |

## Changes for the iScience revision (October 2026)

The revised manuscript replaces pairwise Ka/Ks from KaKs_Calculator2 with PAML codeml M0 ω, and
it reports HyPhy BUSTED tests for episodic diversifying selection. This version adds the code
behind the revised numbers, figures and tables:

| Folder | Produces |
|---|---|
| `04_selection_dnds/paml_m0_completion/` | M0 ω for all 9,624 gene families with at least two distinct sequences (`paml_m0_complete.tsv`), with FastTree gene trees for families with paralogs |
| `04_selection_dnds/busted_gene_tree/` | BUSTED rerun on FastTree gene trees for the 607 families with paralogs (`busted_complete.tsv`) |
| `04_selection_dnds/figure4_paml/` | Figure 4 |
| `04_selection_dnds/dnds_summary/` | The Results numbers on coding-sequence evolution, Table S7 and Table S8 |
| `03_genomic_plasticity_and_defense/figure_s11/` | Figure S11 and its Wilcoxon rank-sum tests |
| `03_genomic_plasticity_and_defense/figure_s17/` | Figure S17 and its Spearman correlation |
| `01_pangenome/genome_accessions/` | Table S6 (GenBank accessions of the 142 genomes) |

Tool versions for the new steps: PAML 4.10.0 (codeml), HyPhy 2.5.48 and FastTree 2.2.0.

The COG over-representation analysis behind Figure 3F and Figure S5 is not in this repository
(see `05_functional_enrichment/README.md`).

## Running the scripts

The scripts read their inputs from the analysis project folder, which holds the genomes, the
Panaroo and PPanGGOLiN outputs, the codon alignments and the trees. The scripts added for the
revision take its location from the `PROJECT_ROOT` environment variable, for example:

```bash
PROJECT_ROOT=/path/to/analysis/project Rscript --vanilla redraw_supp_fig17.R
```

Set `CODEML_BIN` and `FASTTREE_BIN` when codeml or FastTree is not on your `PATH`. Each of these
scripts writes its outputs next to itself and writes nothing into the analysis project folder.
