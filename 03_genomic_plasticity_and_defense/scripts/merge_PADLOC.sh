#!/bin/bash

# create a combined df for PADLOC predictions accross the RSSC genome

PADLOC_ROOT=$(realpath $1)
OUTDIR=$(realpath $2)


mkdir -p "$OUTDIR"


for filepath in `ls $PADLOC_ROOT/*/*.csv`
do
    STRAIN_NAME=$(basename "$filepath" | sed 's/_padloc.csv//')
    echo "Hold On! Processing $STRAIN_NAME"

    cat "$filepath" | awk 'NR>1' | awk -F, -v STRAIN_NAME="$STRAIN_NAME" 'OFS="," {print STRAIN_NAME, $0}' \
    > "$OUTDIR/$STRAIN_NAME.tmp"
done

HEADER="strain,system.number,seqid,system,target.name,hmm.accession,hmm.name,protein.name,full.seq.E.value,domain.iE.value,target.coverage,hmm.coverage,start,end,strand,target.description,relative.position,contig.end,all.domains,best.hits"
echo "$HEADER" > "$OUTDIR/PADLOC_AGGREGATED.csv"

# Merge all strains
cat "$OUTDIR"/*.tmp >> "$OUTDIR/PADLOC_AGGREGATED.csv"

echo "Merged file has been written to $OUTDIR/PADLOC_AGGREGATED.csv"
