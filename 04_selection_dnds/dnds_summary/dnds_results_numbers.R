# Numbers for the Results section on coding-sequence evolution, from the final PAML M0 table and
# the final BUSTED table. Writes:
#   numbers.txt                  every number the Results draft cites
#   selection_candidates.tsv     families with M0 omega > 1 or BUSTED q <= 0.05, with annotation
#   busted_by_family_size.tsv    BUSTED detection by number of distinct sequences
#   coverage_by_class.tsv        PAML estimate coverage per lineage-aware class
#   busted_by_cog.tsv            BUSTED q <= 0.05 per COG category, Fisher's exact test with BH
#   busted_by_class.tsv          BUSTED q <= 0.05 per lineage-aware class
#   busted_function_groups.tsv   tested families in the function groups the Results text names
# Inputs are read only. Run from this folder: PROJECT_ROOT=/path/to/analysis/project Rscript dnds_results_numbers.R

suppressPackageStartupMessages(library(tidyverse))

prj_root <- paste0(Sys.getenv("PROJECT_ROOT", "."), "/")   # set PROJECT_ROOT to the analysis project folder
paml <- read_tsv("../paml_m0_completion/paml_m0_complete.tsv", show_col_types = FALSE)
busted <- read_tsv("../busted_gene_tree/busted_complete.tsv", show_col_types = FALSE)
twilight <- read_tsv(paste0(prj_root, "Manuscript_Figures/Figure_2/Data/twilight_classification_updated.tab"),
                     show_col_types = FALSE)
eggnog <- read_tsv(paste0(prj_root, "Results/Pangenome/PANAROO/RSSC-sensitive/PANGENOME_EGGNOG/RSSCsensitivePANAROO.emapper.annotations"),
                   comment = "##", show_col_types = FALSE) %>%
  rename(query = `#query`) %>%
  select(query, COG_category, Preferred_name, Description)

# Pangenome tables write "~~~" where the alignment file names have "___"
fix <- function(x) str_replace_all(x, "___", "~~~")
classes <- twilight %>% select(gene = gene_name, general_class, specific_class)

out <- character()
say <- function(...) { line <- paste0(...); out <<- c(out, line); cat(line, "\n") }

# PAML M0 coverage --------------------------------------------------------------------------
n_families <- nrow(twilight)
eligible <- nrow(paml)
estimated <- paml %>% filter(run_status %in% c("ok", "skipped_existing_ok"), !is.na(omega))
boundary <- estimated %>% filter(omega >= 999)
used <- estimated %>% filter(omega < 999) %>% mutate(gene = fix(gene)) %>%
  left_join(classes, by = "gene")
stopifnot(!anyNA(used$specific_class))

say("Pangenome families: ", n_families)
say("Eligible for M0 (>= 2 distinct sequences): ", eligible, "; single distinct sequence: ", n_families - eligible)
say("M0 estimates: ", nrow(estimated), " of ", eligible, " (", sprintf("%.1f", 100 * nrow(estimated) / eligible), "%)")
say("Not estimated: ", paste(names(table(paml$run_status[!paml$run_status %in% c("ok", "skipped_existing_ok")])),
                             table(paml$run_status[!paml$run_status %in% c("ok", "skipped_existing_ok")]), collapse = "; "))
say("Boundary omega = 999 excluded: ", nrow(boundary), "; families used: ", nrow(used))
say("Tree source among used: ", paste(names(table(used$tree_source)), table(used$tree_source), collapse = "; "))

core_acc <- used %>% mutate(fraction = ifelse(grepl("Core", general_class), "core", "accessory")) %>%
  group_by(fraction) %>%
  summarise(n = n(), below1 = sum(omega < 1), above1 = sum(omega > 1), median = median(omega), .groups = "drop")
for (i in seq_len(nrow(core_acc))) with(core_acc[i, ], say(
  fraction, ": n = ", n, ", omega < 1 in ", below1, " (", sprintf("%.1f", 100 * below1 / n), "%), omega > 1 in ",
  above1, ", median ", sprintf("%.3f", median)))
say("All used: omega <= 1 in ", sum(used$omega <= 1), " of ", nrow(used), " (", sprintf("%.1f", 100 * mean(used$omega <= 1)), "%)")

by_class <- used %>% group_by(specific_class) %>%
  summarise(n = n(), median = median(omega), above1 = sum(omega > 1),
            pct_above1 = 100 * mean(omega > 1), .groups = "drop") %>% arrange(median)
kw <- kruskal.test(omega ~ specific_class, data = used)
say("Kruskal-Wallis across ", n_distinct(used$specific_class), " classes: chi2 = ", sprintf("%.1f", kw$statistic),
    ", df = ", kw$parameter, ", P = ", format.pval(kw$p.value))
core_ref <- used$omega[used$specific_class == "Collection core"]
wil <- by_class %>% filter(specific_class != "Collection core") %>%
  mutate(p = map_dbl(specific_class, ~ wilcox.test(used$omega[used$specific_class == .x], core_ref)$p.value),
         p_bh = p.adjust(p, "BH"))
by_class <- by_class %>% left_join(wil %>% select(specific_class, p_bh), by = "specific_class")
for (i in seq_len(nrow(by_class))) with(by_class[i, ], say(
  "  ", specific_class, ": n = ", n, ", median ", sprintf("%.3f", median), ", omega > 1 in ", above1,
  " (", sprintf("%.1f", pct_above1), "%)", ifelse(is.na(p_bh), "", paste0(", BH P vs Collection core = ", format.pval(p_bh, digits = 2)))))
big <- wil %>% filter(n > 100)
say("Classes with n > 100: ", nrow(big), "; all higher median than Collection core: ",
    all(big$median > median(core_ref)), "; largest BH P = ", format.pval(max(big$p_bh), digits = 2))

coverage <- twilight %>% select(gene = gene_name, specific_class) %>%
  left_join(paml %>% mutate(gene = fix(gene)) %>% select(gene, omega, run_status), by = "gene") %>%
  mutate(status = case_when(is.na(run_status) ~ "single distinct sequence",
                            !run_status %in% c("ok", "skipped_existing_ok") ~ "failed",
                            omega >= 999 ~ "boundary (omega = 999)",
                            TRUE ~ "estimated")) %>%
  count(specific_class, status) %>% pivot_wider(names_from = status, values_from = n, values_fill = 0)
write_tsv(coverage, "coverage_by_class.tsv")
single_rare <- twilight %>% filter(general_class == "Rare") %>% pull(gene_name) %>%
  setdiff(fix(paml$gene)) %>% length()
say("Single-distinct-sequence families in the Rare classes: ", single_rare)
absent <- coverage %>% filter(specific_class == "Absent in large lineages")
say("Absent in large lineages: ", paste(names(absent)[-1], unlist(absent[1, -1]), collapse = "; "))

# BUSTED ------------------------------------------------------------------------------------
b_ok <- busted %>% filter(!is.na(busted_pvalue))
say("BUSTED eligible (>= 4 distinct sequences): ", nrow(busted), "; with a result: ", nrow(b_ok))
say("BUSTED p <= 0.05: ", sum(b_ok$busted_pvalue <= 0.05), "; q <= 0.05: ", sum(b_ok$busted_qvalue <= 0.05))
say("BUSTED global omega < 1: ", sum(b_ok$global_omega < 1, na.rm = TRUE), " of ", nrow(b_ok),
    " (", sprintf("%.1f", 100 * mean(b_ok$global_omega < 1, na.rm = TRUE)), "%)")

size_bins <- b_ok %>% mutate(bin = cut(n_sequences, c(3, 9, 19, 49, 99, Inf),
                                       labels = c("4-9", "10-19", "20-49", "50-99", ">=100"))) %>%
  group_by(bin) %>% summarise(families = n(), positive_q005 = sum(busted_qvalue <= 0.05),
                              pct = 100 * positive_q005 / families, .groups = "drop")
write_tsv(size_bins, "busted_by_family_size.tsv")
for (i in seq_len(nrow(size_bins))) with(size_bins[i, ], say(
  "  BUSTED ", bin, " sequences: ", positive_q005, " of ", families, " (", sprintf("%.1f", pct), "%) at q <= 0.05"))

# PAML against BUSTED -----------------------------------------------------------------------
both <- used %>% select(gene, omega, specific_class) %>%
  inner_join(b_ok %>% mutate(gene = fix(gene)), by = "gene")
rho <- cor(both$omega, both$global_omega, method = "spearman", use = "complete.obs")
say("Families with M0 and BUSTED: ", nrow(both), "; Spearman rho (M0 omega, BUSTED global omega) = ", sprintf("%.3f", rho))
say("  omega > 1 and q <= 0.05: ", sum(both$omega > 1 & both$busted_qvalue <= 0.05),
    "; omega > 1 only: ", sum(both$omega > 1 & both$busted_qvalue > 0.05),
    "; q <= 0.05 only: ", sum(both$omega <= 1 & both$busted_qvalue <= 0.05),
    "; neither: ", sum(both$omega <= 1 & both$busted_qvalue > 0.05))

# Selection candidates ----------------------------------------------------------------------
candidates <- used %>% select(gene, omega, tree_source, sequence_count, general_class, specific_class) %>%
  full_join(b_ok %>% mutate(gene = fix(gene)) %>%
              select(gene, busted_pvalue, busted_qvalue, global_omega, busted_tree = tree_source), by = "gene") %>%
  filter(omega > 1 | busted_qvalue <= 0.05) %>%
  left_join(classes %>% rename(gc = general_class, sc = specific_class), by = "gene") %>%
  mutate(general_class = coalesce(general_class, gc), specific_class = coalesce(specific_class, sc)) %>%
  select(-gc, -sc) %>%
  left_join(eggnog %>% mutate(query = fix(query)), by = c("gene" = "query")) %>%
  mutate(evidence = case_when(omega > 1 & busted_qvalue <= 0.05 ~ "M0 omega > 1 and BUSTED",
                              omega > 1 ~ "M0 omega > 1",
                              TRUE ~ "BUSTED")) %>%
  arrange(factor(evidence, c("M0 omega > 1 and BUSTED", "M0 omega > 1", "BUSTED")), desc(omega))
write_tsv(candidates, "selection_candidates.tsv")
say("Selection candidates: ", nrow(candidates), " (", paste(names(table(candidates$evidence)),
    table(candidates$evidence), collapse = "; "), ")")

# What the selection candidates are ---------------------------------------------------------
# M0 omega > 1 against the number of distinct sequences in a family
seq_bins <- used %>%
  mutate(seqs = factor(ifelse(sequence_count >= 5, "5 or more", sequence_count), c("2", "3", "4", "5 or more"))) %>%
  group_by(seqs) %>% summarise(families = n(), above1 = sum(omega > 1), .groups = "drop")
for (i in seq_len(nrow(seq_bins))) with(seq_bins[i, ], say(
  "  M0 omega > 1 with ", seqs, " distinct sequences: ", above1, " of ", families,
  " (", sprintf("%.1f", 100 * above1 / families), "%)"))
m0_pos <- candidates %>% filter(omega > 1)
no_function <- is.na(m0_pos$COG_category) | m0_pos$COG_category %in% c("-", "S")
say("M0 omega > 1: ", nrow(m0_pos), "; in the three Lineage-specific classes: ",
    sum(grepl("Lineage specific", m0_pos$specific_class)), "; with two distinct sequences: ",
    sum(m0_pos$sequence_count == 2), "; no eggNOG hit, no COG category or COG S only: ", sum(no_function))
for (i in which(candidates$evidence == "M0 omega > 1 and BUSTED")) with(candidates[i, ], say(
  "  Both criteria: ", gene, ", omega ", sprintf("%.2f", omega), ", q ", signif(busted_qvalue, 2), ", ",
  specific_class, ", ", sequence_count, " distinct sequences, ", Preferred_name, ", COG ", COG_category, ", ", Description))

# BUSTED-positive families by COG category: Fisher's exact test against all tested families, BH
b_ann <- b_ok %>% mutate(gene = fix(gene), positive = busted_qvalue <= 0.05) %>%
  left_join(eggnog %>% mutate(gene = fix(query)) %>% select(-query), by = "gene")
cog <- b_ann %>% mutate(cog = coalesce(COG_category, "?")) %>%
  separate_longer_position(cog, 1) %>%
  mutate(cog = ifelse(cog == "?", "no eggNOG hit", cog)) %>%
  distinct(gene, cog, positive)
n_pos <- sum(b_ann$positive)
cog_test <- cog %>% group_by(cog) %>% summarise(families = n(), positive = sum(positive), .groups = "drop") %>%
  mutate(pct = 100 * positive / families,
         p = map2_dbl(positive, families, ~ fisher.test(matrix(c(.x, .y - .x, n_pos - .x,
                                                                   nrow(b_ann) - .y - (n_pos - .x)), 2))$p.value),
         p_bh = p.adjust(p, "BH")) %>% arrange(p)
write_tsv(cog_test, "busted_by_cog.tsv")
busted_class <- b_ann %>% left_join(classes, by = "gene") %>% group_by(specific_class) %>%
  summarise(families = n(), positive = sum(positive), pct = 100 * positive / families, .groups = "drop") %>%
  arrange(pct)
write_tsv(busted_class, "busted_by_class.tsv")
for (i in seq_len(nrow(busted_class))) with(busted_class[i, ], say(
  "  BUSTED q <= 0.05 in ", specific_class, ": ", positive, " of ", families, " (", sprintf("%.1f", pct), "%)"))
say("BUSTED q <= 0.05 by COG category: smallest BH P = ", format.pval(min(cog_test$p_bh), digits = 2),
    " (", cog_test$cog[1], ", ", cog_test$positive[1], " of ", cog_test$families[1], ")")

# Function groups named in the Results text. We chose these groups after reading the candidate
# list, so the shares describe the list and are not an enrichment test.
effector_name <- "^(rip|pop[BCP]|hop|avr|yopJ)"
groups <- list(
  "type III effector" = grepl(effector_name, b_ann$gene, ignore.case = TRUE) |
    grepl("type III effector", b_ann$Description, ignore.case = TRUE),
  "Rhs, VgrG or toxin" = grepl("Rhs|toxin|hemagglutinin", b_ann$Description, ignore.case = TRUE),
  "c-di-GMP (GGDEF or EAL)" = grepl("diguanylate|GGDEF|EAL", b_ann$Description))
examples <- imap_dfr(groups, ~ b_ann[.x, ] %>% mutate(group = .y)) %>%
  select(group, gene, n_sequences, busted_qvalue, positive, Preferred_name, Description)
write_tsv(examples, "busted_function_groups.tsv")
for (g in names(groups)) say("  BUSTED q <= 0.05 among ", g, " families: ", sum(groups[[g]] & b_ann$positive),
                            " of ", sum(groups[[g]]), "; positive: ",
                            paste(b_ann$gene[groups[[g]] & b_ann$positive], collapse = ", "))
say("  BUSTED q <= 0.05 among all tested families: ", n_pos, " of ", nrow(b_ann),
    " (", sprintf("%.1f", 100 * n_pos / nrow(b_ann)), "%)")

writeLines(out, "numbers.txt")
