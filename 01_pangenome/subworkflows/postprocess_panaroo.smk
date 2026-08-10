"""
    Post process panaroo pangnome outputs for evolutionary analysis
"""

 rule filter_aln:
     "Filter the monomorphic sites and keep the snp sites from panaroo core gene alignments, generates a multifasta of alignments"
     input:
         os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', '{pangenome}', 'core_gene_alignment.aln'),
     output:
         os.path.join(config['RESULTS'], 'Pangenome', 'SNP-sites', '{pangenome}_core_gene_aln_snps.aln'),
     conda:
         "../envs/snpsites.yaml",
     shell:
         """
             snp-sites {input} > {output}
         """
 rule fasttree:
     input:
         os.path.join(config['RESULTS'], 'Pangenome', 'SNP-sites', '{pangenome}_core_gene_aln_snps.aln'),
     output:
         os.path.join(config['RESULTS'], 'Pangenome', 'SNP-tree_fasttree', '{pangenome}_core_snpsites_tree.tree'),
     conda:
         "../envs/fasttree.yaml",
     log:
         "logs/snp_fasttree/{pangenome}_snpsites_fasttree.log",
     shell:
         """
             fasttree -nt -gtr -gamma {input} > {output} \
             2> {log}
         """

 rule midpointroot:
     input:
         os.path.join(config['RESULTS'], 'Pangenome', 'SNP-tree_fasttree', '{pangenome}_core_snpsites_tree.tree'),
     output:
         os.path.join(config['RESULTS'], 'Pangenome', 'SNP-tree_fasttree', '{pangenome}_coresnpsites_midpointrooted.tree'),
     params:
         script = 'scripts/midpoint_root.R'
     conda:
         "../envs/ape_phytools.yaml",
     shell:
         "Rscript {params.script} {input} {output}"
 rule fastbaps_from_snpaln:
     "The SNP sites-filtered core-gene alignment was used as input to FastBaps, alongside the rooted phylogenetic tree to provide a guide for the hierarchical partitioning."
     input:
         aln = os.path.join(config['RESULTS'], 'Pangenome', 'SNP-sites', '{pangenome}_core_gene_aln_snps.aln'),
         rooted_tree = os.path.join(config['RESULTS'], 'Pangenome', 'SNP-tree_fasttree', "{pangenome}_coresnpsites_midpointrooted.tree"), 
     output:
         os.path.join(config['RESULTS'], 'Pangenome', 'FASTBAPS_CLUSTERS', "{pangenome}.fastbapsClusters.csv"),
     conda:
         "../envs/fastbaps.yaml",
     threads:
         20,
     shell:
         """
             fastbapspath=`Rscript -e 'cat(system.file("run_fastbaps", package = "fastbaps"))'` && \
             Rscript $fastbapspath \
             -i  {input.aln} \
             -o {output} \
             -l 4 \
             --phylogeny  {input.rooted_tree} \
             -t {threads}
         """

 rule midpoint_iqtree:
     input:
         os.path.join(config['RESULTS'], 'UFbootstrap_iqtree', '{treefile}.treefile'),
     output:
         os.path.join(config['RESULTS'], 'UFbootstrap_iqtree', 'Midpoint_rooted', '{treefile}_coresnpsites_midpointrooted.treefile'),
     params:
         script = 'scripts/midpoint_root.R'
     conda:
         "../envs/ape_phytools.yaml",
     shell:
         "Rscript {params.script} {input} {output}"

rule fastbaps_with_iqtree:
    input:
        aln = os.path.join(config['RESULTS'], 'Pangenome', 'SNP-sites', 'RSSC-sensitive_core_gene_aln_snps.aln'),
        rooted_tree = os.path.join(config['RESULTS'],'UFbootstrap_iqtree', 'Midpoint_rooted', 'ufbootstrap_rssc_coresnpsites_midpointrooted.treefile'), 
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'FASTBAPS_CLUSTERS_IQTREEML_l4', "RSSC_iqtree_coreSNP.fastbapsClusters_L4_{prior}.csv"),
    conda:
        "../envs/fastbaps.yaml",
    params: "{prior}",
    threads:
        20,
    shell:
        """
            fastbapspath=`Rscript -e 'cat(system.file("run_fastbaps", package = "fastbaps"))'` && \
            Rscript $fastbapspath \
            -i  {input.aln} \
            -o {output} -p {params} \
            -l 4 \
            --phylogeny  {input.rooted_tree} \
            -t {threads}
        """

rule twilight_partitions:
    input:
        fastbaps = os.path.join(config['RESULTS'], 'Pangenome', 'FASTBAPS_CLUSTERS_IQTREEML_l4', "RSSC_iqtree_coreSNP.fastbapsClusters_L4.csv"),
        gpa = os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', 'RSSC-sensitive', 'gene_presence_absence.Rtab'),
    output:
        directory(os.path.join(config['RESULTS'], 'TWILIGHT_ANALYSIS_IQTREE_FASTBAPS')),
    params:
        bapsl4 = os.path.join(config['RESULTS'], 'Pangenome', 'FASTBAPS_CLUSTERS_IQTREEML_l4', "baps_l4.tsv"),
        script = "scripts/twilight/twilight.R",
    shell:
        """
            cat {input.fastbaps} | cut -d, -f 1,5 | awk 'NR>1' | sed 's/,/\t/' > {params.bapsl4} && \
            Rscript {params.script} \
                -p {input.gpa} \
                -g {params.bapsl4} \
                -o {output} \
                -s 2
        """

rule annotate_pangenome:
    "translate the CDS into protein and search with diamond"
    input:
        seq = os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', 'RSSC-sensitive', 'pan_genome_reference.fa'),
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', 'RSSC-sensitive', 'PANGENOME_EGGNOG', 'process.done'),
    params:
        db = '<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/Eggnog_database',
        mode = 'diamond',
        itype = 'CDS',
        # genepred = 'prodigal',
        outdir = os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', 'RSSC-sensitive', 'PANGENOME_EGGNOG'),
    threads:
        20,
    conda:
        "../envs/emapper.yaml",
    shell:
        """
            emapper.py \
                --data_dir {params.db} \
                -m {params.mode} \
                --itype {params.itype} \
                --translate \
                -i {input.seq} \
                -o "RSSCsensitivePANAROO" \
                --output_dir {params.outdir} \
                --excel \
                --cpu {threads} && \
                touch {output}
        """
