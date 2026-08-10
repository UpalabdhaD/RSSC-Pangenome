# Alternative DNDS Models

This folder adds a reproducible pilot workflow for estimating gene-family-level dN/dS with models beyond KaKs_Calculator.

Inputs are the existing codon alignments:

```text
Results/DNDS/PAL2NAL_DEDUP/*.codon.dedup.fa
```

The PAML path runs `codeml` one-ratio M0 (`model = 0`, `NSsites = 0`) and reports one omega estimate per gene family. The helper script maps gene IDs back to strains using `Metadata/Full_RSSC_142_Metadata.tsv`, prunes the core-SNP species tree, and renames tree tips to match the codon alignment. If a family contains duplicated or unmapped strains, it falls back to a deterministic balanced sequence tree and records that in the manifest.

Shell-only pilot run:

```bash
bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_paml_m0_shell.sh
bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_hyphy_fel_shell.sh
```

To use a named conda environment instead of the default discovered binaries:

```bash
CONDA_ENV=my_paml_env bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_paml_m0_shell.sh
CONDA_ENV=my_hyphy_env bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_hyphy_fel_shell.sh
```

To run a custom gene set, put one gene family per line in a text file and pass it as the first argument:

```bash
bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_paml_m0_shell.sh my_genes.txt
bash MS_REVISION/DNDS_ALT_MODELS/scripts/run_hyphy_fel_shell.sh my_genes.txt
```

The shell wrappers default to existing materialized conda binaries under `.snakemake/conda`. They do not create new environments.
