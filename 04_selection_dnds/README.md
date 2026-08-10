# 04_selection_dnds

Corresponds to Methods §5.11, "Ka/Ks estimation."

## What's here

- `subworkflows/dnds.smk` — the Ka/Ks pipeline: per-Panaroo-gene-family protein
  extraction (`seqtk`), MAFFT L-INS-i alignment, PAL2NAL codon back-translation,
  `seqkit`-based deduplication, custom multi-FASTA→pairwise-AXT conversion, then
  KaKs_Calculator2 with the bacterial genetic code (`-c 11`) and the Nei–Gojobori
  method (`-m NG`).
- `scripts/parseFastaIntoAXT.pl` — the custom AXT-block parser named explicitly in
  §5.11 ("a custom parser").
- `envs/kaks2.yaml`, `envs/pal2nal.yaml`, `envs/seqtk.yaml`,
  `envs/protein_alignment.yaml` — per-tool conda environments for the steps above.
- `scripts/RSSC_dnds.Rmd` — parses and summarizes the Ka/Ks output tables in R; this
  notebook also carries out the COG over-representation analysis described in §5.12
  (Figure 4 presents Ka/Ks and functional-category enrichment together), so it's the
  canonical location for both — see `../05_functional_enrichment/README.md`.
- `scripts/dnds.Rmd`, `scripts/Filter_valid_kaks.ipynb` — supporting exploration and
  filtering of gene families with unstable/low-divergence Ka/Ks estimates (see
  Methods §5.11's caution on this).

## Peer-review-driven extensions

None of the tools or analyses below are currently described in the manuscript's
Methods text — this is a snapshot of the reviewer rebuttal work, folded into this
directory since it all extends the same Ka/Ks pipeline above.

- `dnds_alt_models/` — alternative selection-model comparison against the primary
  KaKs_Calculator2/NG result: PAML M0 (`run_paml_m0*`) and HyPhy FEL/FITMG94
  (`run_hyphy_fel*`, `run_hyphy_fitmg94.sh`), with summarization and collection
  scripts for each.
- `dnds_sensitivity/` — sensitivity analysis testing how stable the Ka/Ks positive-
  selection calls are (see `reviewer_response_dnds_sensitivity.md` and
  `reviewer_response_comment_5.md` for the reviewer comment this responds to).
- `dnds_paml_vs_kaks/` — direct comparison of PAML M0 vs. KaKs_Calculator2 estimates,
  including a per-family confusion analysis and the "absent from major lineages" gene
  subset comparison flagged as needing cautious interpretation in Methods §5.11.
- `dnds_hyphy_busted/` — HyPhy BUSTED episodic diversifying-selection tests, run both
  on the "absent from major lineages" subset and genome-wide (`all_busted/`). The
  installed HyPhy conda environment used to run these (`DNDS_HYPHY/envs/hyphy_2.5.50/`,
  ~580 MB) is intentionally excluded — rebuild from `envs/hyphy.yaml` in
  `dnds_alt_models/` instead of copying the installed environment.
- `replicon_variable_rgp/` — chromosome-vs-megaplasmid metrics for variable
  (RGP-associated) vs. persistent genome fractions, including a Fisher's exact test
  comparing the two replicons.

Result tables (`.tsv`/`.csv`/`.json`), plots (`.png`/`.pdf`), and the `.zip` archive
under `REPLICON_VARIABLE_RGP/` were left out of the directories above to keep this a
code/write-up snapshot rather than a data dump — regenerate them by running the
scripts above.
