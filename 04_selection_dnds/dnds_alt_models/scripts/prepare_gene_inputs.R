#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  args[[idx + 1]]
}

alignment <- get_arg("--alignment")
metadata_file <- get_arg("--metadata")
species_tree_file <- get_arg("--species-tree")
out_fasta <- get_arg("--out-fasta")
out_phylip <- get_arg("--out-phylip")
out_tree <- get_arg("--out-tree")
manifest_file <- get_arg("--manifest")
min_sequences <- as.integer(get_arg("--min-sequences", "4"))

if (any(vapply(list(alignment, metadata_file, species_tree_file, out_fasta, out_phylip, out_tree, manifest_file), is.null, logical(1)))) {
  stop("Missing required argument")
}

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- grep("^>", lines)
  if (!length(headers)) stop("No FASTA records in ", path)
  ends <- c(headers[-1] - 1, length(lines))
  names <- sub("^>", "", lines[headers])
  seqs <- mapply(function(start, end) {
    paste(gsub("\\s+", "", lines[(start + 1):end]), collapse = "")
  }, headers, ends, USE.NAMES = FALSE)
  names(seqs) <- names
  toupper(seqs)
}

write_fasta <- function(seqs, path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    writeLines(gsub("(.{1,80})", "\\1\n", seqs[[nm]], perl = TRUE), con)
  }
}

write_phylip <- function(seqs, path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(sprintf("%d %d", length(seqs), nchar(seqs[[1]])), con)
  for (nm in names(seqs)) {
    writeLines(sprintf("%s  %s", nm, seqs[[nm]]), con)
  }
}

map_to_strain <- function(ids, strains) {
  strains <- strains[order(nchar(strains), decreasing = TRUE)]
  vapply(ids, function(id) {
    hit <- strains[id == strains | startsWith(id, paste0(strains, "_"))]
    if (length(hit)) hit[[1]] else NA_character_
  }, character(1))
}

dir.create(dirname(out_fasta), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_phylip), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tree), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)

seqs <- read_fasta(alignment)
seqs <- seqs[nchar(seqs) > 0]
seq_lengths <- unique(nchar(seqs))
if (length(seq_lengths) != 1) stop("Alignment has unequal sequence lengths: ", alignment)
if (seq_lengths %% 3 != 0) stop("Alignment length is not divisible by 3: ", alignment)
if (length(seqs) < min_sequences) stop("Too few sequences after filtering: ", alignment)

metadata <- read.delim(metadata_file, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
strain_map <- map_to_strain(names(seqs), metadata$STRAIN)
tree_source <- "balanced_sequence_tree"

if (!anyNA(strain_map) && !anyDuplicated(strain_map)) {
  species_tree <- read.tree(species_tree_file)
  if (all(strain_map %in% species_tree$tip.label)) {
    drop_tips <- setdiff(species_tree$tip.label, strain_map)
    gene_tree <- drop.tip(species_tree, drop_tips)
    gene_tree$tip.label <- names(strain_map)[match(gene_tree$tip.label, strain_map)]
    tree_source <- "pruned_species_tree"
  } else {
    gene_tree <- stree(length(seqs), type = "left")
    gene_tree$tip.label <- names(seqs)
    gene_tree$edge.length <- rep(0.1, nrow(gene_tree$edge))
  }
} else {
  gene_tree <- stree(length(seqs), type = "left")
  gene_tree$tip.label <- names(seqs)
  gene_tree$edge.length <- rep(0.1, nrow(gene_tree$edge))
}

write_fasta(seqs, out_fasta)
write_phylip(seqs, out_phylip)
write.tree(gene_tree, file = out_tree)

manifest <- data.frame(
  sequence_id = names(seqs),
  strain = strain_map,
  sequence_length = nchar(seqs),
  tree_source = tree_source,
  stringsAsFactors = FALSE
)
write.table(manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
