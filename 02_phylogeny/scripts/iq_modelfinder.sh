#!/bin/bash

seed=1028


ufbootstrap="1000"
aln=$(realpath $1)
model="MFP"
outprefix=$(realpath $2)
maxthreads=40
maxmem="80%"
threads="AUTO"

mkdir -p "$outprefix"

time \
	iqtree \
	-s "$aln" \
	-m "$model" \
	--prefix "$outprefix/ufbootstrap_rssc" \
	-T "$threads" \
	--threads-max "$maxthreads" \
	--mem "$maxmem" \
	--seed "$seed" \
	-B "$ufbootstrap"
