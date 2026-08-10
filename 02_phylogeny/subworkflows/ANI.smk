rule ani:
    input:
        os.path.join(config['input'], 'Rso_{db}_complete'),
    output:
        os.path.join(config['RESULTS'], 'ANI', 'Rso_{db}_complete', 'ANIclustermap_dendrogram.nwk'),
    params:
        outdir = os.path.join(config['RESULTS'], 'ANI', 'Rso_{db}_complete'),
        mode = "fastani",
        figwidth = 20,
        figheight = 20,
        # color = "white,orange,red",
        cmap_range = "80,85,90,92.5,95,97.5,100",
    threads:
        10,
    conda:
        "../envs/aniclustermap.yaml",
    shell:
        """
            ANIclustermap \
            -m {params.mode} \
            -i {input} \
            -o {params.outdir} \
            --fig_width {params.figwidth} \
            --fig_height {params.figheight} \
            --cmap_range {params.cmap_range}
        """