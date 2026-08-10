"""
https://carpentries-incubator.github.io/genome-mining/02-antismash/index.html
"""

rule antismash:
    input:
        os.path.join(config['RESULTS'], 'Annotation', 'Ralstonia-159-complete-BAKTA', '{genome}_BAKTA', '{genome}.gbff'),
    output:
        directory(os.path.join(config['RESULTS'], 'Annotation', 'AntiSMASH', '{genome}')),
    params:
        db=config['ANTISMASH_DB'],
    threads:
        4,
    conda:
        '../envs/antismash.yaml',
    log:
        'logs/antismash/{genome}.log',
    shell:
        """
            antismash \
                --minimal \
                {input} \
                -t bacteria \
                -c {threads} \
                --databases {params.db} \
                --output-dir {output} \
                --output-basename {wildcards.genome}
        """