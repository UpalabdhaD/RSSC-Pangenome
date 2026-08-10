"Granular partitioning of the pangenome with population aware fashion"

rule twilight:
    input:
        gpa = os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', '{groups}', "gene_presence_absence.Rtab"),
        baps = os.path.join(config['RESULTS'], 'Pangenome', 'FASTBAPS_CLUSTERS', "{groups}.fastbapsClusters.csv"),
    output:
        bapsfile = os.path.join(config['RESULTS'], 'Pangenome', 'twilight_partitions', "{groups}", "{groups}.fastbapsClusters.toTwilight.tab"),
        flag = os.path.join(config['RESULTS'], 'Pangenome', 'twilight_partitions', "{groups}", "process.done"),
    params:
        min_size = 1,
        core = float(0.95),
        rare = float(0.15),
        bapsL2 = int(3), #here using the second level
        script = "scripts/twilight/twilight.R",
        outdir = os.path.join(config['RESULTS'], 'Pangenome', 'twilight_partitions', "{groups}"),
    shell:
        """
            cat {input.baps} | awk 'NR>1' | cut -d "," -f 1,{params.bapsL2} | sed 's/,/\t/' > {output.bapsfile} && \
            Rscript {params.script} \
                -p {input.gpa} \
                -g {output.bapsfile} \
                -o {params.outdir} \
                --min_size {params.min_size} \
                --core_threshold {params.core} \
                --rare_threshold {params.rare} && \
                touch {output.flag}
        """