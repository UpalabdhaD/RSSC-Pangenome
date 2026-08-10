#!/bin/bash

INPUT=$(realpath $1)
OUTPUT_DIR=$(realpath $2)

mkdir -p "$OUTPUT_DIR"

cat "$INPUT" |awk 'NR>1' | cut -f 1,2,3,7,8,11,12 | awk '$11!="TRUE" && $12!="TRUE"'| cut -f 1-5 | \
while read RGPID GENOME CONTIG START END;
do
    echo "Hold on! Processing $RGPID" 
    # mkdir -p "${OUTPUT_DIR}/${GENOME}"
    echo -e "${CONTIG}\t${START}\t${END}\t${RGPID}" >> "$OUTPUT_DIR/${GENOME}.bed"
done