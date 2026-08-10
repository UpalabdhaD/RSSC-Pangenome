"pangenome co-occurance"

rule coinfinder:
    input:
        pa = os.path.join(config["RESULTS"], "Pangenome", "CoinfinderInput", "{groups}.roaryFormatted.csv"),
        phylo = os.path.join(config["RESULTS"], "Pangenome", "SNP-tree_fasttree", "{groups}_coresnpsites_midpointrooted_nonZero_e-6.tree"),
    output:
        flag = os.path.join(config["RESULTS"], "Cooccurance", "{groups}", "coinfinder.done"),
    params:
        maxmode = "--associate",
        signif = int(0.05),
        correction = "--bonferroni",
        prefix = os.path.join(config["RESULTS"], "Cooccurance", "{groups}"),
    threads:
        6,
    conda:
        "../envs/coinfinder.yaml",
    log:
        "logs/Coinfinder/{groups}.log",
    shell:
        """
            coinfinder \
                --input {input.pa} \
                --phylogeny {input.phylo} \
                --inputroary \
                --associate \
                --level {params.signif} \
                {params.correction} \
                --num_cores {threads} \
                --output {params.prefix} && \
                touch {output.flag}
        """
