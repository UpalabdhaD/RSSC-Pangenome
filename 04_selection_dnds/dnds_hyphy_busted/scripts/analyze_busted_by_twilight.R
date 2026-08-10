#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(rstatix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop(
        "Usage: analyze_busted_by_twilight.R BUSTED_SUMMARY TWILIGHT_CLASSIFICATION OUTPUT_DIR",
        call. = FALSE
    )
}

busted_path <- args[[1]]
twilight_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

busted <- read_tsv(busted_path, show_col_types = FALSE, na = c("", "NA", "NaN"))
twilight <- read_tsv(twilight_path, show_col_types = FALSE)

analysis_data <- busted %>%
    filter(
        input_status == "eligible",
        busted_status == "ok",
        is.finite(global_omega_from_log),
        global_omega_from_log > 0,
        is.finite(busted_pvalue),
        busted_pvalue >= 0
    ) %>%
    inner_join(
        twilight %>% select(gene_name, general_class, specific_class),
        by = c("gene" = "gene_name")
    ) %>%
    mutate(
        general_class = factor(
            general_class,
            levels = c("Core", "Intermediate", "Rare", "Varied", "Absent in large lineages")
        ),
        busted_log10_pvalue = -log10(pmax(busted_pvalue, .Machine$double.xmin)),
        busted_log10_pvalue_capped = pmin(busted_log10_pvalue, 10)
    ) %>%
    filter(!is.na(general_class))

if (n_distinct(analysis_data$general_class) < 2) {
    stop("Fewer than two Twilight categories have valid HyPhy BUSTED results.", call. = FALSE)
}

write_tsv(analysis_data, file.path(out_dir, "busted_twilight_analysis_data.tsv"))

metrics <- tibble::tribble(
    ~metric, ~value_column,
    "BUSTED global omega", "global_omega_from_log",
    "BUSTED p-value", "busted_log10_pvalue"
)

kruskal_rows <- list()
pairwise_rows <- list()
group_summary_rows <- list()

for (i in seq_len(nrow(metrics))) {
    metric <- metrics$metric[[i]]
    value_column <- metrics$value_column[[i]]
    metric_data <- analysis_data %>%
        transmute(general_class, value = .data[[value_column]]) %>%
        filter(is.finite(value))

    kw <- kruskal_test(metric_data, value ~ general_class) %>%
        mutate(metric = metric, n = nrow(metric_data), groups = n_distinct(metric_data$general_class))
    kruskal_rows[[i]] <- kw

    pairwise_rows[[i]] <- pairwise_wilcox_test(
        metric_data,
        value ~ general_class,
        p.adjust.method = "BH",
        exact = FALSE
    ) %>%
        mutate(metric = metric, n = nrow(metric_data))

    group_summary_rows[[i]] <- metric_data %>%
        group_by(general_class) %>%
        summarise(
            n = n(),
            median = median(value),
            iqr = IQR(value),
            .groups = "drop"
        ) %>%
        mutate(metric = metric, .before = 1)
}

kruskal_results <- bind_rows(kruskal_rows) %>%
    mutate(p_fdr = p.adjust(p, method = "BH")) %>%
    select(metric, n, groups, statistic, df, p, p_fdr, method)
pairwise_results <- bind_rows(pairwise_rows) %>%
    select(metric, n, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)
group_summary <- bind_rows(group_summary_rows) %>%
    arrange(metric, general_class)

write_tsv(kruskal_results, file.path(out_dir, "busted_twilight_kruskal_wallis.tsv"))
write_tsv(pairwise_results, file.path(out_dir, "busted_twilight_pairwise_wilcoxon_bh.tsv"))
write_tsv(group_summary, file.path(out_dir, "busted_twilight_group_summary.tsv"))

plot_metric <- function(metric_name, value_column, y_label, log_scale = FALSE) {
    plot_data <- analysis_data %>%
        transmute(general_class, value = .data[[value_column]])
    kw_p <- kruskal_results %>% filter(metric == metric_name) %>% pull(p_fdr)
    subtitle <- sprintf("Kruskal-Wallis BH-adjusted p = %.3g", kw_p)

    plot <- ggplot(plot_data, aes(x = general_class, y = value, fill = general_class)) +
        geom_violin(trim = FALSE, alpha = 0.72, colour = "grey25", linewidth = 0.25) +
        geom_boxplot(width = 0.14, outlier.shape = NA, colour = "grey15", linewidth = 0.28) +
        stat_summary(fun = median, geom = "point", size = 1.6, colour = "black") +
        scale_fill_brewer(palette = "Set2", guide = "none") +
        scale_x_discrete(labels = c("Absent in large lineages" = "Absent in\nlarge lineages")) +
        labs(x = "Twilight group", y = y_label, subtitle = subtitle) +
        theme_classic(base_size = 11) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    if (log_scale) {
        plot <- plot + scale_y_log10()
    }
    plot
}

omega_plot <- plot_metric(
    "BUSTED global omega",
    "global_omega_from_log",
    "BUSTED global omega (log10 scale)",
    log_scale = TRUE
)
pvalue_plot <- plot_metric(
    "BUSTED p-value",
    "busted_log10_pvalue_capped",
    "-log10(BUSTED p-value), capped at 10"
)

ggsave(file.path(out_dir, "busted_global_omega_by_twilight.pdf"), omega_plot, width = 8.6, height = 5.4)
ggsave(file.path(out_dir, "busted_global_omega_by_twilight.png"), omega_plot, width = 8.6, height = 5.4, dpi = 300)
ggsave(file.path(out_dir, "busted_pvalue_by_twilight.pdf"), pvalue_plot, width = 8.6, height = 5.4)
ggsave(file.path(out_dir, "busted_pvalue_by_twilight.png"), pvalue_plot, width = 8.6, height = 5.4, dpi = 300)
