
rule eggnogmapper:
    input:
        seq = os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'PANRGP_SEQ_FROM_BED', "{genome}.seq"),
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', "EMAPPER", "{genome}", "process.done"),
    params:
        db = '<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/Eggnog_database',
        mode = 'diamond',
        itype = 'genome',
        genepred = 'prodigal',
        outdir = os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', "EMAPPER", "{genome}"),
    threads:
        8,
    conda:
        "../envs/emapper.yaml",
    shell:
        """
            emapper.py \
                --data_dir {params.db} \
                -m {params.mode} \
                --itype {params.itype} \
                --genepred {params.genepred} \
                -i {input.seq} \
                -o {wildcards.genome} \
                --output_dir {params.outdir} \
                --excel \
                --cpu {threads} && \
                touch {output}
        """