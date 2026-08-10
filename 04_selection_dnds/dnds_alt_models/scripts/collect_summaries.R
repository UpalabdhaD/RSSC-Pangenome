#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
indir <- if (length(args) >= 1) args[[1]] else "MS_REVISION/DNDS_ALT_MODELS/results/paml_m0"
outfile <- if (length(args) >= 2) args[[2]] else "MS_REVISION/DNDS_ALT_MODELS/results/paml_m0_summary.tsv"

files <- list.files(indir, pattern = "^summary\\.tsv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No summary.tsv files found under ", indir)

dfs <- lapply(files, function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
})
out <- do.call(rbind, dfs)
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
write.table(out, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
