#!/usr/bin/env Rscript
library(phytools)
library(ape)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the required arguments are provided
if (length(args) < 2) {
  stop("Two arguments must be supplied: input file and output file.", call. = FALSE)
}

input_file <- args[1]
output_file <- args[2]

tree_rssc <- read.tree(input_file)
# Perform midpoint rooting
midpoint <- midpoint.root(tree_rssc)
write.tree(midpoint, file = output_file)

cat("Successfully rooted the tree and saved to", output_file, "\n")