#!/bin/bash

INDIR=$(realpath "$1")
OUTDIR=$(realpath "$2")
DB=$(realpath "$3")
THREADS="$4"


mkdir -p $OUTDIR
for gbkpath in `ls $INDIR/*.gbff`;
	do
		prefix=$(basename "${gbkpath}" | sed 's/.gbff//')
		echo "Antismash running for ${prefix}"

		antismash \
			"$gbkpath" \
			-c "$THREADS" \
			--databases "$DB" \
			--output-dir "${OUTDIR}/${prefix}"
	done
