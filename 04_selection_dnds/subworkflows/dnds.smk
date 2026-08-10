
# rule align:
#     input:
#         os.path.join(config['RESULTS'], "Pangenome", "PANAROO", "RSSC-sensitive", "rssc_pan_genefam_proteins", "{genefamily}.list.fa"),
#     output:
#         os.path.join(config['RESULTS'], "DNDS", "GENE_FAM_ALN", "{genefamily}.prot.aln"),
#     threads:
#         1,
#     conda:
#         "../envs/protein_alignment.yaml",
#     shell:
#         """
#             mafft-linsi {input} > {output}
#         """

# rule pal2nal:
#     input:
#         aa_aln = os.path.join(config['RESULTS'], "DNDS", "GENE_FAM_ALN", "{genefamily}.prot.aln"),
#         cds = os.path.join(config['RESULTS'], "DNDS", "CDS2", "{genefamily}.nt.fa"),
#     output:
#         os.path.join(config['RESULTS'], "DNDS", "PAL2NAL", "{genefamily}.codon.fa"),
#     params:
#         remove_stop = '-nogap',
#         output_fmt = 'fasta',
#         codontable = 11,
#     conda:
#         "../envs/pal2nal.yaml",
#     threads:
#         1,
#     shell:
#         """
#             pal2nal.pl \
#                 {input.aa_aln} {input.cds} \
#                 {params.remove_stop} \
#                 -codontable {params.codontable} \
#                 -output {params.output_fmt} \
#                 > {output}
#         """


rule remove_duplicates:
    input:
        aln=os.path.join(config['RESULTS'], "DNDS", "PAL2NAL", "{genefamily}.codon.fa"),
    output:
        aln_unique=os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP", "{genefamily}.codon.dedup.fa"),
    conda:
        "../envs/hyphy.yaml",
    shell:
        """
            seqkit rmdup -s {input.aln} -o {output.aln_unique}
        """

rule faToaxt:
    input: 
        os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP", "{genefamily}.codon.dedup.fa"),
    output:
        os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP_AXT", "{genefamily}.codon.dedup.axt"),
    params:
        prefix=os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP_AXT", "{genefamily}.codon.dedup"),
    shell:
        "perl scripts/parseFastaIntoAXT.pl {input} {params.prefix}"

rule kaks2:
    input:
        os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP_AXT", "{genefamily}.codon.dedup.axt"),
    output:
        os.path.join(config['RESULTS'], "DNDS", "PAL2NAL_DEDUP_KAKS", "{genefamily}.kaks.tsv"),
    params:
        codon="11",
        method="NG"
    conda:
        "../envs/kaks2.yaml",
    shell:
        """
            KaKs_Calculator -i {input} -o {output} -c {params.codon} -m {params.method}
        """