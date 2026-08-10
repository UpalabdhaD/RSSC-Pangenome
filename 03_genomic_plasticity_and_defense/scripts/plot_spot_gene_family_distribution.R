#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: plot_spot_gene_family_distribution.R INPUT_TSV OUTPUT_DIR", call. = FALSE)
}

input_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spots <- read_tsv(input_path, show_col_types = FALSE) %>%
    transmute(spot, gene_family_count = nb_families) %>%
    filter(is.finite(gene_family_count), gene_family_count >= 0)

if (nrow(spots) == 0) {
    stop("No non-negative gene-family counts were found.", call. = FALSE)
}

summary <- spots %>%
    summarise(
        spots = n(),
        minimum = min(gene_family_count),
        first_quartile = quantile(gene_family_count, 0.25),
        median = median(gene_family_count),
        mean = mean(gene_family_count),
        third_quartile = quantile(gene_family_count, 0.75),
        maximum = max(gene_family_count)
    )
write_tsv(summary, file.path(output_dir, "spot_gene_family_count_summary.tsv"))

median_count <- summary$median[[1]]
histogram <- ggplot(spots, aes(x = gene_family_count)) +
    geom_histogram(binwidth = 25, boundary = 0, colour = "white", fill = "#277da1") +
    geom_vline(xintercept = median_count, colour = "#d1495b", linewidth = 0.7, linetype = "dashed") +
    annotate(
        "text",
        x = median_count,
        y = Inf,
        label = sprintf("Median = %.0f", median_count),
        vjust = 1.5,
        hjust = -0.05,
        colour = "#a22c3a",
        size = 3.4
    ) +
    labs(
        x = "Gene families per spot",
        y = "Number of spots"
    ) +
    theme_classic(base_size = 12)

ecdf_plot <- ggplot(spots, aes(x = gene_family_count)) +
    stat_ecdf(geom = "step", linewidth = 0.8, colour = "#277da1") +
    geom_vline(xintercept = median_count, colour = "#d1495b", linewidth = 0.7, linetype = "dashed") +
    scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(
        x = "Gene families per spot",
        y = "Cumulative proportion of spots"
    ) +
    theme_classic(base_size = 12)

ggsave(file.path(output_dir, "spot_gene_family_distribution.pdf"), histogram, width = 7.2, height = 4.8)
ggsave(file.path(output_dir, "spot_gene_family_distribution.png"), histogram, width = 7.2, height = 4.8, dpi = 300)
ggsave(file.path(output_dir, "spot_gene_family_ecdf.pdf"), ecdf_plot, width = 7.2, height = 4.8)
ggsave(file.path(output_dir, "spot_gene_family_ecdf.png"), ecdf_plot, width = 7.2, height = 4.8, dpi = 300)
