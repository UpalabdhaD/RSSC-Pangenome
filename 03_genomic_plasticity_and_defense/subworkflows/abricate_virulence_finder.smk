rule removefastafromgff:
    input:
        os.path.join(config['RESULTS'], 'Annotation', 'Ralstonia-159-complete-BAKTA', '{genome}_BAKTA', '{genome}.gff3'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'DEFENSE_SYSTEMS', 'PADLOC', "GFFS_WO_FASTA", '{genome}.gff'),
    shell:
        """
            sed '/^##FASTA/Q' {input} >  {output}
        """

rule padloc:
    input:
        gff = os.path.join(config['RESULTS'],'Annotation', 'DEFENSE_SYSTEMS', 'PADLOC', "GFFS_WO_FASTA", '{genome}.gff'),
        faa = os.path.join(config['RESULTS'], 'Annotation', 'Ralstonia-159-complete-BAKTA', '{genome}_BAKTA', '{genome}.faa'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'DEFENSE_SYSTEMS', 'PADLOC', "{genome}", "{genome}.done"),
    params:
        outdir = os.path.join(config['RESULTS'], 'Annotation', 'DEFENSE_SYSTEMS', 'PADLOC', "{genome}"),
        padlocdb = "<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/PADLOC",
    # conda:
    #     "../envs/padloc.yaml",
    threads:
        4,
    shell:
        """
            padloc \
                --faa {input.faa} \
                --gff {input.gff} \
                --outdir {params.outdir} \
                --data {params.padlocdb} \
                --cpu {threads} \
                --force && \
                touch {output}
        """

rule defensefinder:
    input:
        faa = os.path.join(config['RESULTS'], 'Annotation', 'Ralstonia-159-complete-BAKTA', '{genome}_BAKTA', '{genome}.faa'),
    output:
        os.path.join(config['RESULTS'], 'Annotation', 'DEFENSE_SYSTEMS', "{genome}", "{genome}_defense_finder_systems.tsv"),
    threads:
        4,
    params:
        outdir = os.path.join(config['RESULTS'], 'Annotation', 'DEFENSE_SYSTEMS', "{genome}"),
        db = "<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/DEFENSE_FINDER",
    conda:
        "../envs/defensefinder.yaml",
    shell:
        """
            defense-finder \
                run \
                -o {params.outdir} \
                -w {threads} \
                --models-dir {params.db} \
                {input}
        """

# rule genomad:
#     input:
#         os.path.join(config['RESULTS'], 'Annotation', 'GenBank-complete-BAKTA', '{genome}_BAKTA', '{genome}.fna'),
#     output:
#         directory(os.path.join(config['RESULTS'], 'Annotation', 'GENOMAD', "{genome}")),
#     params:
#         db = "<LOCAL_DB_ROOT>/BIOINFO_LOCAL_DBS/PROKARYOTIC_GENOME_ANNOTATION/GENOMAD",
#         # outdir = os.path.join(config['RESULTS'], 'Annotation', 'GENOMAD', "{genome}"),
#     threads:
#         10,
#     conda:
#         "../envs/genomad.yaml",
#     shell:
#         """
#             genomad \
#                 end-to-end \
#                 --threads {threads} \
#                 --cleanup \
#                 --splits 8 \
#                 {input} \
#                 {output} \
#                 {params.db}
#         """
    