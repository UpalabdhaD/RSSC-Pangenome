# 02_phylogeny

Corresponds to STAR Methods > Method details, "Phylogenetic reconstruction and population structure
analysis."

## What's here

- `subworkflows/ANI.smk`, `scripts/ANI_clustermap.py`, `scripts/genarate_nw.sh` —
  pairwise FastANI across all 142 genomes, filtered at a 95% ANI cutoff, clustered by
  connected components (networkx) and visualized in Cytoscape / UPGMA (pheatmap).
  `genarate_nw.sh` is the wrapper that chains ANI distances through to the clustered
  network output; it hardcodes an old pre-migration absolute path to
  `ANI_clustermap.py` — fix before running outside the original machine.
- `scripts/create_anicluster.Rmd`, `scripts/jaccard_distance_isolates.Rmd` — notebooks
  building on the ANI/accessory-distance clustering above.
- `envs/snpsites.yaml` — SNP-sites, used to filter the 3,348-gene Panaroo core-gene
  alignment (3,263,944 bp) down to polymorphic sites (591,669 bp) ahead of tree
  building.
- `scripts/iq_modelfinder.sh` — IQ-TREE v2.3.6 maximum-likelihood phylogeny with 1,000
  ultrafast bootstraps, model selected by ModelFinder (GTR+F+I+R2).
- `scripts/midpoint_root.R`, `envs/ape_phytools.yaml` — midpoint-roots the IQ-TREE
  output for display (ape/phytools) ahead of the iTOL visualization shown in Figure 2.
- `scripts/grapple/pw_similarity.py` — pairwise accessory-genome distance calculation
  from the Grapple suite, feeding the igraph/ggraph clustering described in Methods.

## Not included

FastTree — an env spec exists in the working `workflow/envs/`, but the manuscript's
published ML phylogeny uses IQ-TREE, not FastTree, so it's left out here.
