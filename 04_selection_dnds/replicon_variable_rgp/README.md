# F1C1 Chromosome vs Megaplasmid Variable Gene and RGP Metrics

## Source files

- Gene partition and RGP membership: `Results/Pangenome/PPANGGOLIN/PANRGP/tables_merged.tsv`
- F1C1 RGP coordinates/lengths: `Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv`
- Replicon lengths and product keywords: `Data/genbank_assemblies/gffs/F1C1.gff`
- Defense genes: `Manuscript_Figures/Figure_DefenseSystem/Data/SPOT_RGP_DEFENSESYS.tsv`
- Type III effector hits: `Results/Annotation/ABRICATE/F1C1__RSSC_T3EV3.tab`

## Definitions

- Persistent genes: PanRGP `partition == persistent`.
- Variable/accessory genes: PanRGP `partition == shell` or `partition == cloud`.
- RGP genes per Mb: number of F1C1 PanRGP gene rows assigned to RGP IDs retained in `spot_merged_with_rgps.tsv`, divided by replicon length in Mb.
- Accessory genes per Mb: shell + cloud genes divided by replicon length in Mb.
- Mobile-element-associated genes: CDS in the F1C1 GFF whose product contains one of these keyword stems: transposase, integrase, recombinase, resolvase, insertion sequence, is family, mobile, phage, prophage, conjug, relaxase, transposon.
- Type III effector genes: ABRicate RSSC_T3Ev3 hits counted by contig; this is hit-based.

## Fisher exact test

Table orientation:

| Replicon | Variable genes | Persistent genes |
|---|---:|---:|
| Chromosome | 825 | 2527 |
| Megaplasmid | 720 | 892 |

- Odds ratio, chromosome vs megaplasmid: 0.4045
- Odds ratio, megaplasmid vs chromosome: 2.472
- Fisher exact two-sided p-value: 4.719e-45

Interpretation: the megaplasmid is enriched for variable genes relative to the chromosome when using shell+cloud versus persistent gene counts.
