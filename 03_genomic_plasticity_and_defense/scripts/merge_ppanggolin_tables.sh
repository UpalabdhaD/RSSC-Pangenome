#!/bin/bash

TABLE_DIR=$(realpath $1)
OUTPUTFILEPATH=$(realpath $2)

for tsvfilepath in `ls $TABLE_DIR/*.tsv`;
do
    strain = 
    cat tsvfilepath 
done
