#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_dnds_sensitivity.R <results_dir> <plots_dir>", call. = FALSE)
}

results_dir <- args[[1]]
plots_dir <- args[[2]]
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

specific <- read_tsv_base(file.path(results_dir, "sensitivity_by_specific_class.tsv"))
general <- read_tsv_base(file.path(results_dir, "sensitivity_by_general_class.tsv"))
joined <- read_tsv_base(file.path(results_dir, "kaks_sensitivity_joined.tsv"))

main_filter_levels <- c(
  "baseline_numeric",
  "ks_ge_0.05",
  "qc_alignment",
  "strict_ks0.05_qc",
  "strict_ks0.05_2_qc"
)
detailed_filter_levels <- c(
  "baseline_numeric",
  "ks_ge_0.05",
  "alignment_len_ge300",
  "equal_alignment_lengths",
  "gap_le_0.20",
  "ambiguous_le_0.05",
  "no_internal_stops",
  "ks_le_2",
  "strict_ks0.05_qc",
  "strict_ks0.05_2_qc"
)
filter_levels <- c(
  "baseline_numeric",
  "ks_ge_0.01",
  "ks_ge_0.05",
  "alignment_len_ge300",
  "equal_alignment_lengths",
  "gap_le_0.20",
  "ambiguous_le_0.05",
  "no_internal_stops",
  "ks_le_2",
  "qc_alignment",
  "min4_sequences",
  "strict_ks0.05_qc",
  "strict_ks0.05_2_qc",
  "min4_strict_ks0.05_qc"
)
filter_labels <- c(
  baseline_numeric = "Baseline",
  ks_ge_0.01 = "Ks >= 0.01",
  ks_ge_0.05 = "Ks >= 0.05",
  alignment_len_ge300 = "Length >= 300 nt",
  equal_alignment_lengths = "Equal lengths",
  gap_le_0.20 = "Gap <= 20%",
  ambiguous_le_0.05 = "Ambiguous <= 5%",
  no_internal_stops = "No internal stops",
  ks_le_2 = "Ks <= 2",
  qc_alignment = "Alignment QC",
  min4_sequences = ">=4 dedup seqs",
  strict_ks0.05_qc = "Ks >= 0.05 + QC",
  strict_ks0.05_2_qc = "0.05 <= Ks <= 2 + QC",
  min4_strict_ks0.05_qc = ">=4 seqs + Ks >= 0.05 + QC"
)

specific$filter_raw <- specific$filter
general$filter_raw <- general$filter
specific$filter <- factor(specific$filter, levels = filter_levels, labels = filter_labels)
general$filter <- factor(general$filter, levels = filter_levels, labels = filter_labels)

plot_specific <- specific[specific$filter_raw %in% main_filter_levels, ]
plot_general <- general[general$filter_raw %in% main_filter_levels, ]

baseline <- plot_specific[plot_specific$filter == "Baseline", c("class", "positive_percent")]
baseline <- baseline[order(baseline$positive_percent, decreasing = TRUE), ]
plot_specific$class <- factor(plot_specific$class, levels = baseline$class)

p1 <- ggplot(plot_specific, aes(x = class, y = positive_percent, fill = filter)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Genes with Ka/Ks > 1 (%)",
    fill = "Sensitivity filter",
    title = "Ka/Ks > 1 sensitivity by pangenome class"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(plots_dir, "positive_fraction_by_specific_class.png"),
  plot = p1,
  width = 10,
  height = 6.5,
  dpi = 300
)

target <- plot_specific[plot_specific$class == "Absent in large lineages", ]
p2 <- ggplot(target, aes(x = filter, y = positive_percent, group = 1)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = paste0(round(positive_percent, 1), "%\n", positive_n, "/", n)),
            vjust = -0.55, size = 3) +
  ylim(0, max(100, max(target$positive_percent, na.rm = TRUE) + 8)) +
  labs(
    x = NULL,
    y = "Genes with Ka/Ks > 1 (%)",
    title = "Sensitivity for genes absent in large lineages"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  filename = file.path(plots_dir, "absent_large_lineages_sensitivity.png"),
  plot = p2,
  width = 8,
  height = 4.8,
  dpi = 300
)

detailed_target <- specific[
  specific$filter_raw %in% detailed_filter_levels &
    specific$class == "Absent in large lineages",
]
detailed_target$filter <- factor(
  detailed_target$filter,
  levels = filter_labels[detailed_filter_levels]
)

p2b <- ggplot(detailed_target, aes(x = filter, y = positive_percent, group = 1)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = paste0(round(positive_percent, 1), "%\n", positive_n, "/", n)),
            vjust = -0.55, size = 3) +
  ylim(0, max(100, max(detailed_target$positive_percent, na.rm = TRUE) + 8)) +
  labs(
    x = NULL,
    y = "Genes with Ka/Ks > 1 (%)",
    title = "Detailed artifact checks for genes absent in large lineages"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  filename = file.path(plots_dir, "absent_large_lineages_individual_artifact_checks.png"),
  plot = p2b,
  width = 10,
  height = 5.2,
  dpi = 300
)

plot_general$class <- factor(plot_general$class, levels = unique(plot_general$class[plot_general$filter == "Baseline"]))
p3 <- ggplot(plot_general, aes(x = class, y = positive_percent, fill = filter)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76) +
  labs(
    x = NULL,
    y = "Genes with Ka/Ks > 1 (%)",
    fill = "Sensitivity filter",
    title = "Ka/Ks > 1 sensitivity by general class"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(plots_dir, "positive_fraction_by_general_class.png"),
  plot = p3,
  width = 8,
  height = 4.8,
  dpi = 300
)

joined$kaks <- as.numeric(joined$kaks)
joined$ks <- as.numeric(joined$ks)
strict <- joined[
  is.finite(joined$kaks) &
    joined$ks >= 0.05 &
    joined$ks <= 2 &
    joined$alignment_nt_length >= 300 &
    joined$unequal_lengths == 0 &
    joined$gap_fraction <= 0.20 &
    joined$ambiguous_fraction <= 0.05 &
    joined$stop_codon_count == 0,
]
strict$kaks_capped <- pmin(strict$kaks, 5)
strict$specific_class <- factor(strict$specific_class, levels = levels(plot_specific$class))

p4 <- ggplot(strict, aes(x = specific_class, y = kaks_capped)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Ka/Ks, capped at 5",
    title = "Ka/Ks distribution after strict Ks and alignment QC"
  ) +
  theme_bw(base_size = 11)

ggsave(
  filename = file.path(plots_dir, "strict_kaks_distribution_by_specific_class.png"),
  plot = p4,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

message("Wrote plots to ", plots_dir)
