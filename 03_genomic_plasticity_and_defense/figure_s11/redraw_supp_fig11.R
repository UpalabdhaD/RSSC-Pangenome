# Redraw Supplementary Figure 11 with occupancy-based group labels.
#
# This script repeats the data steps of the original figure code
# ($PROJECT_ROOT/Manuscript_Figures/Figure_RGP/RGP.Rmd)
# on the same input files, which it only reads. It then replaces the labels
# "Hotspot" and "Spot" with labels that state the grouping rule the original
# code applied in practice: whether more than 100 of the 142 genomes carry an
# RGP at the spot. The script stops if the new grouping differs from the
# original one in any row.
#
# Run from this folder: PROJECT_ROOT=/path/to/analysis/project Rscript --vanilla redraw_supp_fig11.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
})

data_dir <- paste0(Sys.getenv("PROJECT_ROOT", "."), "/Manuscript_Figures/Figure_RGP/Data/")   # set PROJECT_ROOT to the analysis project folder
out_png <- "Supp_Figure_11_spot_occupancy.png"
out_pdf <- "Supp_Figure_11_spot_occupancy.pdf"

N_GENOMES <- 142
OCC_CUT <- 100

# ---- Data steps, as in RGP.Rmd ----------------------------------------------
rgps <- read_tsv(paste0(data_dir, "regions_of_genomic_plasticity.tsv"), show_col_types = FALSE) %>%
  rename("rgp_id" = "region")
spots <- read_tsv(paste0(data_dir, "spots.tsv"), show_col_types = FALSE)
spot_summary <- read_tsv(paste0(data_dir, "summarize_spots.tsv"), show_col_types = FALSE) %>%
  rename("spot_id" = "spot")

spot_characterization <- rgps %>%
  left_join(spots, by = "rgp_id") %>%
  left_join(spot_summary, by = "spot_id") %>%
  mutate(SPOT_CHAR = ifelse(nb_unique_family_sets > 100, "Hotspot", "Spot")) %>%
  mutate(CONTIGS = case_when(contig == "contig_1" ~ "Chromosome",
                             TRUE ~ "Non-chromosomal")) %>%
  filter(SPOT_CHAR != "NA")

# ---- New grouping: number of genomes with an RGP at the spot ------------------
label_wide  <- "Widespread spots\n(>100 of 142 genomes)"
label_other <- "Other spots\n(≤100 genomes)"

genomes_per_spot <- spot_characterization %>%
  group_by(spot_id) %>%
  summarise(genomes = n_distinct(genome), .groups = "drop")

spot_characterization <- spot_characterization %>%
  left_join(genomes_per_spot, by = "spot_id") %>%
  mutate(GROUP = factor(ifelse(genomes > OCC_CUT, label_wide, label_other),
                        levels = c(label_wide, label_other)))

# Check 1: the new grouping matches the original one in every row.
same_group <- (spot_characterization$SPOT_CHAR == "Hotspot") ==
              (spot_characterization$GROUP == label_wide)
stopifnot(all(same_group))

# Check 2: the expected totals from the original figure.
stopifnot(nrow(spot_characterization) == 10110,
          n_distinct(spot_characterization$spot_id) == 651,
          n_distinct(spot_characterization$spot_id[spot_characterization$GROUP == label_wide]) == 44)

cat("Grouping check passed: new groups match the original Hotspot/Spot groups in all",
    nrow(spot_characterization), "RGP rows.\n\n")

# ---- Numbers for the legend and the statistics section -----------------------
cat("Spots per group:\n")
print(spot_characterization %>% group_by(GROUP) %>% summarise(spots = n_distinct(spot_id), rgps = n(), .groups = "drop"))
cat("\nSpots per group and replicon class (top-left bars):\n")
print(spot_characterization %>% group_by(CONTIGS, GROUP) %>% summarise(spots = n_distinct(spot_id), .groups = "drop"))
cat("\nRGPs per group and replicon class (right panels):\n")
print(spot_characterization %>% count(GROUP, CONTIGS))
cat("\nGenes per RGP, widespread vs other spots (bottom-left panel):\n")
print(compare_means(genes ~ GROUP, data = spot_characterization))
cat("\nRGP length, chromosomal vs non-chromosomal, within each group (right panels):\n")
print(compare_means(length ~ CONTIGS, data = spot_characterization, group.by = "GROUP"))
cat("\nMedian and interquartile range of RGP length:\n")
print(spot_characterization %>% group_by(GROUP, CONTIGS) %>%
        summarise(median = median(length), IQR = IQR(length), .groups = "drop"))

# ---- Plots, as in RGP.Rmd except for labels ----------------------------------
spot_chars <- spot_characterization %>%
  group_by(CONTIGS, GROUP) %>%
  summarise(n = n_distinct(spot_id), .groups = "drop") %>%
  ggplot(aes(x = GROUP, y = n, fill = CONTIGS)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_d3(name = "Replicon") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("") +
  ylab("Number of spots") +
  coord_flip()

spot_genecount <- spot_characterization %>%
  ggplot(aes(x = GROUP, y = genes, fill = GROUP)) +
  geom_jitter(alpha = 1/20, size = 1) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point", size = 2, color = "white") +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  scale_fill_aaas() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Number of genes in RGP") +
  stat_compare_means(aes(label = after_stat(p.signif)), label.x = 1.5) +
  coord_flip()

RGP_len <- spot_characterization %>%
  ggplot(aes(x = CONTIGS, y = length, fill = CONTIGS)) +
  geom_jitter(alpha = 1/20, size = 0.2) +
  geom_violin(trim = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "white") +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  facet_wrap(~ GROUP) +
  scale_fill_aaas() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("RGP length (bp)") +
  stat_compare_means(aes(label = after_stat(p.signif)), label.x = 1.5)

fig <- (spot_chars / spot_genecount | RGP_len)

# Fixed seed so the jittered points are reproducible. Same size as the original (10 x 7 in, 300 dpi).
set.seed(2026)
ggsave(out_png, fig, width = 10, height = 7, dpi = 300)
set.seed(2026)
ggsave(out_pdf, fig, width = 10, height = 7, device = cairo_pdf)
cat("\nSaved", out_png, "and", out_pdf, "\n")
