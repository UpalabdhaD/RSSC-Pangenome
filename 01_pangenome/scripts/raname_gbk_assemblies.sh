#!/bin/bash
# strain to path 
for i in genbank/bacteria/Ralstonia/*/*/*.fna.gz; do strain=$(echo "$i" | cut -d '/' -f 5); echo $strain,$i >> Ralstonia_strain_assembly_rename.csv;done

# then copy by strain name

