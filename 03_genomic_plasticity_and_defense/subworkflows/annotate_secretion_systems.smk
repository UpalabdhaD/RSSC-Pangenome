
rule abricate:
    input:
        os.path.join(config['RESULTS'], 'Annotation', 'Ralstonia-159-complete-BAKTA', '{genome}_BAKTA', '{genome}.gbff'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'ABRICATE', "{genome}__{db}.tab"),
    params:
        db = "{db}",
        dbdir = "<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/PROKARYOTIC_GENOME_ANNOTATION/ABRICATE_DB/"
    threads:
        10,
    conda:
        "../envs/abricate.yaml",
    shell:
        """
            abricate \
                --db {params.db} \
                --datadir {params.dbdir} \
                {input} > {output}
        """

rule concat_tab_per_db:
    input:
        os.path.join(config['RESULTS'], 'Annotation', 'ABRICATE'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'ABRICATE_CONCAT_perDB', "{db}.tab"),
    params:
        db = "{db}",
    shell:
        """
            cat {input}/*__{params.db}* | grep -v "#" > {output}  
        """

rule summarize_abricate_perDB:
    input:
        os.path.join(config['RESULTS'], 'Annotation', 'ABRICATE'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'ABRICATE_SUMMARY_perDB', "{db}.tab"),
    params:
        db = '{db}',
    conda:
        "../envs/abricate.yaml",
    shell:
        """
            abricate \
                --summary \
                --nopath \
                {input}/*__{params.db}* \
                > {output}
        """
