

rule annotate_rssc_eggnog:
    "use protein faa"
    input:
        faa = os.path.join(config['RESULTS'], 'Annotation', 'GenBank-complete-BAKTA', '{genome}_BAKTA', '{genome}.faa'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'EGGNOG_RSSC_142', "{genome}", 'process.done'),
    params:
        db = '<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/Eggnog_database',
        mode = 'diamond',
        itype = 'proteins',
        outdir = os.path.join(config['RESULTS'], 'Annotation', 'EGGNOG_RSSC_142', "{genome}"),
    threads:
        16,
    conda:
        "../envs/emapper.yaml",
    shell:
        """
            emapper.py \
                --data_dir {params.db} \
                -m {params.mode} \
                --itype {params.itype} \
                -i {input.faa} \
                -o "{wildcards.genome}" \
                --output_dir {params.outdir} \
                --excel \
                --cpu {threads} && \
                touch {output}
        """
