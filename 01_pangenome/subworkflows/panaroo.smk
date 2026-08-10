"""
helps to identify possible issues with samples prior to running the full pangenome algorithm
https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
"""


rule panaroo_qc:
    input:
        os.path.join(config['RESULTS'], 'Pangenome', 'GFFS', "{collection}"),
    output:
        os.path.join(config['RESULTS'], 'EVALUATION', 'PANAROO', '{collection}', "mash_dist.txt"),
    params:
        db=config["PANAROO_QC_DB"],
        out=os.path.join(config['RESULTS'], 'EVALUATION', 'PANAROO', '{collection}'),
    threads:
        8,
    conda:
        "../envs/panaroo.yaml",
    shell:
        """
            panaroo-qc -i {input}/*.gff \
                -o {params.out} \
                -t {threads} \
                --graph_type all \
                --ref_db {params.db} \
                && touch {output}
        """

rule panaroo:
    input:
        os.path.join(config['RESULTS'], 'Pangenome', 'GFFS', "{collection}"),
    output:
        directory(os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', '{collection}-{mode}')),
    params:
        cleanmode = "{mode}",
        iden_thr = config['Panaroo']['iden_thr'],
        alignment = config['Panaroo']['alignment'],
        aligner = config['Panaroo']['aligner'],
    threads:
        20,
    conda:
        "../envs/panaroo.yaml",
    log:
        "logs/Pangenome/Panaroo/{collection}-{mode}.log",
    shell:
        """
            panaroo \
                -t {threads} \
                -i {input}/*.gff \
                -o {output} \
                --clean-mode {params.cleanmode} \
                --threshold {params.iden_thr} \
                --alignment {params.alignment} \
                --aligner {params.aligner} \
                &> {log}
        """
