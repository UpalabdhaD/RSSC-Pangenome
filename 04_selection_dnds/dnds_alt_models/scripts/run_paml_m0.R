#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  args[[idx + 1]]
}

gene <- get_arg("--gene")
phylip <- normalizePath(get_arg("--phylip"), mustWork = TRUE)
tree <- normalizePath(get_arg("--tree"), mustWork = TRUE)
outdir <- get_arg("--outdir")
codeml <- Sys.getenv("CODEML", unset = "codeml")

if (any(vapply(list(gene, phylip, tree, outdir), is.null, logical(1)))) {
  stop("Missing required argument")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ctl <- file.path(outdir, "codeml.ctl")
mlc <- file.path(outdir, "mlc")
summary <- file.path(outdir, "summary.tsv")

ctl_lines <- c(
  sprintf("seqfile = %s", phylip),
  sprintf("treefile = %s", tree),
  sprintf("outfile = %s", mlc),
  "noisy = 0",
  "verbose = 0",
  "runmode = 0",
  "seqtype = 1",
  "CodonFreq = 2",
  "clock = 0",
  "model = 0",
  "NSsites = 0",
  "icode = 0",
  "fix_kappa = 0",
  "kappa = 2",
  "fix_omega = 0",
  "omega = 0.2",
  "getSE = 0",
  "RateAncestor = 0",
  "Small_Diff = 1e-6",
  "cleandata = 1",
  "method = 0"
)
writeLines(ctl_lines, ctl)

old_wd <- getwd()
setwd(outdir)
on.exit(setwd(old_wd), add = TRUE)

status <- system2(codeml, "codeml.ctl", stdout = TRUE, stderr = TRUE)
writeLines(status, file.path(outdir, "codeml.stdout"))

if (!file.exists(mlc)) {
  write.table(
    data.frame(gene = gene, status = "failed", omega = NA, kappa = NA, lnL = NA, np = NA),
    summary, sep = "\t", row.names = FALSE, quote = FALSE
  )
  quit(status = 1)
}

txt <- readLines(mlc, warn = FALSE)

extract_first <- function(pattern) {
  hit <- regmatches(txt, regexec(pattern, txt, perl = TRUE))
  hit <- hit[vapply(hit, length, integer(1)) > 1]
  if (!length(hit)) return(NA_character_)
  hit[[1]][[2]]
}

omega <- extract_first("omega \\(dN/dS\\) =\\s*([0-9.Ee+-]+)")
kappa <- extract_first("kappa \\(ts/tv\\) =\\s*([0-9.Ee+-]+)")
lnl_match <- regmatches(txt, regexec("lnL\\(ntime:\\s*\\d+\\s+np:\\s*(\\d+)\\):\\s*([-0-9.Ee+]+)", txt, perl = TRUE))
lnl_match <- lnl_match[vapply(lnl_match, length, integer(1)) > 2]
np <- if (length(lnl_match)) lnl_match[[1]][[2]] else NA_character_
lnL <- if (length(lnl_match)) lnl_match[[1]][[3]] else NA_character_

write.table(
  data.frame(
    gene = gene,
    status = "ok",
    omega = as.numeric(omega),
    kappa = as.numeric(kappa),
    lnL = as.numeric(lnL),
    np = as.integer(np)
  ),
  summary, sep = "\t", row.names = FALSE, quote = FALSE
)
