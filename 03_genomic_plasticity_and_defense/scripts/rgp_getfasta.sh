#!/bin/bash

BEDDIR=$(realpath $1)
GENOMEDIR=$(realpath $2)
OUTDIR=$(realpath $3)

mkdir -p $OUTDIR

for bedfilepath in `ls "$BEDDIR"`;
    do
        genome=$(basename "$bedfilepath" | sed 's/.bed//')    
        echo "Hold on! getting seq from $genome using $bedfilepath"
        
        fnafile="$GENOMEDIR/${genome}_BAKTA/${genome}.fna"

        bedtools getfasta -nameOnly -fi "$fnafile" -bed "$BEDDIR/$bedfilepath" | fold -w 60 \
        > "$OUTDIR/$genome.seq"
    done


