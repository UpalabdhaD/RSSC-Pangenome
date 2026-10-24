# Redraw Figure S17: gene-family richness of each integration spot against the number of
# defense-associated gene families in the same spot.
#
# The original plot ($PROJECT_ROOT/Manuscript_Figures/Figure_RGP/
# Spot_with_defense.Rmd, chunk after "No of gene family vs no of defense genes") put the number of
# distinct RGPs per spot on the x axis (n_distinct(rgp_id)) and the number of distinct defense-system
# types on the y axis (n_distinct(SYSTEM)), while its axis labels, the legend and the Results text
# describe gene families. Its colours came from SPOT_CHAR, which marks spots with more than 100 RGPs
# (44 spots), not the main-text hotspot rule of more than 100 gene families (51 spots).
#
# This script repeats the original join on the same input files, which it only reads, and plots what
# the text describes:
#   x  distinct gene families in the spot (panRGP nb_families, the value behind the hotspot rule)
#   y  distinct gene families among the PADLOC-annotated genes in the spot's RGPs
#   colour  hotspot = more than 100 distinct gene families
# As in the original, the plot keeps spots with at least one defense gene, and it draws a
# least-squares line with its 95% confidence band. It adds Spearman's rank correlation.
#
# Run from this folder: PROJECT_ROOT=/path/to/analysis/project Rscript --vanilla redraw_supp_fig17.R > redraw_log.txt

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
})

prj_root <- paste0(Sys.getenv("PROJECT_ROOT", "."), "/")   # set PROJECT_ROOT to the analysis project folder
out_png <- "Supp_Figure_17_defense_gene_families.png"
out_pdf <- "Supp_Figure_17_defense_gene_families.pdf"
FONT <- "Liberation Sans"   # Arial letter widths
HOT_CUT <- 100              # hotspot: more than 100 distinct gene families (main text)

# ---- Data steps, as in Spot_with_defense.Rmd -----------------------------------------------------
spot_rgp <- read_tsv(paste0(prj_root, "Manuscript_Figures/Figure_RGP/Data/spot_merged_with_rgps.tsv"),
                     show_col_types = FALSE)
rgp_genes <- read_tsv(paste0(prj_root, "Results/Pangenome/PPANGGOLIN/PANRGP/SPOT_RGP_GENEID_TABLE.tsv"),
                      show_col_types = FALSE)
padloc <- read_tsv(paste0(prj_root, "Manuscript_Figures/Figure_DefenseSystem/Data/PADLOC_DFinder_merged.tsv"),
                   show_col_types = FALSE)

spot_rgp_defense <- spot_rgp %>%
  left_join(rgp_genes, by = c("rgp_id" = "RGP")) %>%
  left_join(padloc, by = c("gene" = "target.name")) %>%
  mutate(GENETYPE = ifelse(PREDICTION %in% c("PADLOC", "DF_PADLOC", "DF"), "DEFENSE", "OTHER"))

# Check: the totals reported in the manuscript.
spots <- spot_rgp %>% distinct(spot_id, nb_rgp, nb_families)
stopifnot(n_distinct(spot_rgp$rgp_id) == 10110,
          nrow(spots) == 651,
          sum(spots$nb_families > HOT_CUT) == 51,
          sum(spot_rgp_defense$GENETYPE == "DEFENSE") == 3486,
          nrow(padloc) == 4785)
cat("Input checks passed: 651 spots, 10,110 RGPs, 51 spots with more than 100 gene families,",
    "3,486 of 4,785 PADLOC defense-gene annotations inside spot RGPs.\n\n")

# ---- Per-spot values -------------------------------------------------------------------------------
families_from_genes <- spot_rgp_defense %>%
  group_by(spot_id) %>%
  summarise(families_counted = n_distinct(family), .groups = "drop")

defense_per_spot <- spot_rgp_defense %>%
  filter(GENETYPE == "DEFENSE") %>%
  group_by(spot_id) %>%
  summarise(defense_families = n_distinct(family),
            defense_systems = n_distinct(system),
            defense_rows = n(), .groups = "drop")

plot_df <- spots %>%
  left_join(families_from_genes, by = "spot_id") %>%
  inner_join(defense_per_spot, by = "spot_id") %>%   # spots with at least one defense gene
  mutate(GROUP = factor(ifelse(nb_families > HOT_CUT, "Hotspot", "Other spot"),
                        levels = c("Hotspot", "Other spot")))

cat("panRGP nb_families equals the distinct families counted from the gene table in",
    sum(spots$nb_families == families_from_genes$families_counted[match(spots$spot_id, families_from_genes$spot_id)]),
    "of 651 spots\n")
cat("Spots with at least one defense gene:", nrow(plot_df), "\n")
print(plot_df %>% count(GROUP, name = "spots"))
cat("\nSpots with the most defense-associated gene families:\n")
print(plot_df %>% arrange(desc(defense_families)) %>%
        select(spot_id, nb_rgp, nb_families, defense_families, defense_systems, defense_rows) %>% head(12))

# ---- Statistics for the legend --------------------------------------------------------------------
rho <- cor.test(plot_df$nb_families, plot_df$defense_families, method = "spearman", exact = FALSE)
fit <- lm(defense_families ~ nb_families, data = plot_df)
cat(sprintf("\nSpearman rho = %.3f, P = %.2e, n = %d spots\n", rho$estimate, rho$p.value, nrow(plot_df)))
cat(sprintf("Least-squares line: defense families = %.3f + %.4f x gene families, R^2 = %.3f\n",
            coef(fit)[1], coef(fit)[2], summary(fit)$r.squared))
rho_sys <- cor.test(plot_df$nb_families, plot_df$defense_systems, method = "spearman", exact = FALSE)
cat(sprintf("For comparison, gene families against defense-system types: Spearman rho = %.3f, P = %.2e\n",
            rho_sys$estimate, rho_sys$p.value))

# ---- Plot -------------------------------------------------------------------------------------------
LABEL_TOP <- 8   # name the spots with the most defense-associated gene families
plot_df <- plot_df %>%
  mutate(label = ifelse(rank(-defense_families, ties.method = "first") <= LABEL_TOP,
                        sub("spot_", "spot ", spot_id), ""))

fig <- ggplot(plot_df, aes(x = nb_families, y = defense_families)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
              colour = "black", fill = "grey70", linewidth = 0.4) +
  geom_rug(aes(colour = GROUP), linewidth = 0.25, alpha = 0.7) +
  geom_point(aes(colour = GROUP), size = 1.3, alpha = 0.85) +
  geom_text_repel(aes(label = label), family = FONT, size = 6 / .pt, min.segment.length = 0,
                  segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf, seed = 2026) +
  scale_colour_d3(name = NULL,
                  labels = c(sprintf("Hotspot (>%d gene families; n = %d)", HOT_CUT, sum(plot_df$GROUP == "Hotspot")),
                             sprintf("Other spot (≤%d gene families; n = %d)", HOT_CUT, sum(plot_df$GROUP == "Other spot")))) +
  labs(x = "Distinct gene families in spot",
       y = "Distinct defense-associated gene families in spot") +
  theme_bw(base_size = 7, base_family = FONT) +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        panel.grid.minor = element_blank())

ggsave(out_png, fig, width = 114, height = 100, units = "mm", dpi = 600, device = ragg::agg_png)
ggsave(out_pdf, fig, width = 114, height = 100, units = "mm", device = cairo_pdf)
cat("\nSaved", out_png, "and", out_pdf, "\n")
