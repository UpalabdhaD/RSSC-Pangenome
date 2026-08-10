"Genarate Ppanggolin Pangenome"

rule panaroo_to_pang:
    input:
        os.path.join(config['RESULTS'], 'Pangenome', 'PANAROO', 'RSSC-sensitive', 'gene_presence_absence_roary.csv'),
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'panaroo_cluster.csv'),
    params:
        script = os.path.join(config['SCRIPTS'], 'roary_to_ppanggolin.py'),
    shell:
        """
            python {params.script} -r {input} -p {output}
        """

rule remove_non_ASCII:
    input:
        os.path.join(config['RESULTS'], 'Pangenome', 'GFFS', 'RSSC', '{sample}.gff'),
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'INPUT_GFFS', '{sample}.gff'),
    shell:
        """
            sed 's/’//g'  {input} > {output}
        """


rule integrate:
    input:
        os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'genome.txt'),
    output:
        directory(os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'ANNOTATE')),
    threads:
        30,
    conda:
        "../envs/ppanggolin.yaml",
    shell:
        """
            ppanggolin annotate \
                --anno {input} \
                --output {output} \
                -c {threads} \
                --basename "RSSC" 
        """

rule cluster:
    input:
        h5 = os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'ANNOTATE', 'RSSC.h5'),
        cluster = os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'panaroo_cluster.csv'),
    output:
        directory(os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'CLUSTER'))
    threads:
        40,
    shell:
        """
            ppanggolin cluster \
                -p {input.h5} \
                --clusters {input.cluster} \
                --infer_singleton \
                -c {threads}
        """


rule p_workflow:
    input:
        anno=os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'genome.txt'),
        clusters=os.path.join(os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'panaroo_cluster.csv')),
    output:
        os.path.join(config['RESULTS'], 'Pangenome', 'PPANGGOLIN', 'workflow'),
    
    params:
        identity=95,
        
    conda:
        "../envs/ppanggolin.yaml",
    threads:
        40,
    shell:

