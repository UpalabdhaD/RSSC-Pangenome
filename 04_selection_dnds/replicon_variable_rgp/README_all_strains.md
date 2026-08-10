# All-Strain Replicon Variable/RGP Analysis

This all-strain calculation extends the F1C1-only table to all 142 genomes represented in `Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv`.

## Source files

- Gene partitions and RGP membership: `Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv`
- Retained RGP IDs, lengths, and coordinates: `Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv`
- Replicon lengths: sequence lengths from `Data/PANGENOME_DATA/RSSC/*.fna`, mapped by FASTA order to PanRGP `contig_1`, `contig_2`, etc.

## Replicon classes

- `Chromosome`: PanRGP `contig_1`.
- `Non-chromosomal`: PanRGP `contig_2`, `contig_3`, or `contig_4`; this corresponds to the manuscript RGP table label `Plasmid`.

## Outputs

- `all_strains_per_replicon_metrics.tsv`: one row per genome/replicon.
- `all_strains_replicon_class_summary.tsv`: aggregate chromosome vs non-chromosomal metrics.
- `all_strains_replicon_class_summary.md`: manuscript-readable table.
- `all_strains_fisher_variable_vs_persistent.tsv`: Fisher exact test for shell+cloud vs persistent genes.

## Main test

The aggregate Fisher test compares:

| Replicon class | Variable genes | Persistent genes |
|---|---:|---:|
| Chromosome | 128,533 | 355,865 |
| Non-chromosomal | 100,795 | 126,710 |

- Odds ratio, non-chromosomal vs chromosome: 2.202
- Fisher exact two-sided p-value: <1e-300; log10_p=-4749.917
- Fisher exact two-sided log10(p): -4749.917
