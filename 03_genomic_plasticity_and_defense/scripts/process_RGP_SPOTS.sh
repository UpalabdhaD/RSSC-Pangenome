#!/bin/bash


SPOT_DB="<PROJECT_ROOT>/Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv"
SPOTID="$1"
OUTDIR=$(realpath "$2")
GFFDIR=$(realpath "$3")


mkdir -p "$OUTDIR"

cat "$SPOT_DB" |  awk 'NR>1' | awk '$11!="TRUE" && $12!="TRUE"'| grep "$SPOTID" | cut -f 1,2,3,7,8 | \
while read RGPID GENOME CONTIG START END;
    do  
        # get rgp from that spot
        INTERVAL=$(echo -e "${CONTIG}\t${START}\t${END}")
        
        # grab the gff file from BAKTA annotation
        GFFFILE="${GFFDIR}/${GENOME}_BAKTA/${GENOME}.gff3"

        # remove the sequence from the gff3 end, to use in bedtools
        GFFCLEAN=$(cat "$GFFFILE" | sed '/^##FASTA$/,$d')

        echo "$INTERVAL" | \
        bedtools intersect \
            -f 1 \
            -a <(echo "$GFFCLEAN") -b - \
            > "${OUTDIR}/${RGPID}.gff"
    done