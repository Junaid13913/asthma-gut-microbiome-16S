#############################################################
# 04_differential_abundance.R
# Adjusted differential abundance, patients vs controls (MaAsLin2)
#
#   Each visit compared with controls in a separate model:
#     CLR abundance ~ group + age + sex      (primary, Model 1)
#     CLR abundance ~ group + age + sex + BMI (sensitivity, Model 2)
#   Pre-computed CLR values as input (normalization = NONE, transform = NONE)
#   3 levels x 2 thresholds x 2 comparisons x 2 models = 24 runs
#   FDR: Benjamini-Hochberg over taxa for the group term, within each run
#
#   Figure 1D   genus/family, 0.1%/10%, Model 1, taxa with P < 0.05
#   Figure E4B-D phylum/family/genus, both thresholds
#   Table S2    DA columns (genus, 0.1%/10%, Model 1)
#############################################################

source("R/00_setup.R")
suppressPackageStartupMessages({ library(Maaslin2); library(ggtext); library(grid) })

maaslin_dir <- file.path(out_dir, "maaslin2", "group_DA")
metadata_clean <- readRDS(f_metadata)
rd <- function(f) readRDS(file.path(obj_dir, f))
genus_count_01  <- rd("genus_count_01pct.rds");  genus_count_05  <- rd("genus_count_05pct.rds")
family_count_01 <- rd("family_count_01pct.rds"); family_count_05 <- rd("family_count_05pct.rds")
phylum_count_01 <- rd("phylum_count_01pct.rds"); phylum_count_05 <- rd("phylum_count_05pct.rds")
genus_clr_01    <- rd("genus_clr_01pct.rds");    genus_clr_05    <- rd("genus_clr_05pct.rds")
family_clr_01   <- rd("family_clr_01pct.rds");   family_clr_05   <- rd("family_clr_05pct.rds")
phylum_clr_01   <- rd("phylum_clr_01pct.rds");   phylum_clr_05   <- rd("phylum_clr_05pct.rds")
genus_count_all  <- readRDS(f_unfiltered[["Genus"]])
family_count_all <- readRDS(f_unfiltered[["Family"]])
phylum_count_all <- readRDS(f_unfiltered[["Phylum"]])

#############################################################
# 1. MaAsLin2 runs
#############################################################
meta_all <- metadata_clean %>%
  mutate(sample_id = as.character(sample_id),
         age_num = suppressWarnings(as.numeric(age)),
         sex_num = as.numeric(factor(sex)) - 1,
         bmi_num = suppressWarnings(as.numeric(bmi))) %>%
  as.data.frame()
rownames(meta_all) <- meta_all$sample_id

count_lookup <- list(Genus_01pct = genus_count_01, Genus_05pct = genus_count_05,
                     Family_01pct = family_count_01, Family_05pct = family_count_05,
                     Phylum_01pct = phylum_count_01, Phylum_05pct = phylum_count_05)
clr_objects <- list(
  list(clr = genus_clr_01,  count_key = "Genus_01pct",  level = "Genus",  thr = "01pct"),
  list(clr = genus_clr_05,  count_key = "Genus_05pct",  level = "Genus",  thr = "05pct"),
  list(clr = family_clr_01, count_key = "Family_01pct", level = "Family", thr = "01pct"),
  list(clr = family_clr_05, count_key = "Family_05pct", level = "Family", thr = "05pct"),
  list(clr = phylum_clr_01, count_key = "Phylum_01pct", level = "Phylum", thr = "01pct"),
  list(clr = phylum_clr_05, count_key = "Phylum_05pct", level = "Phylum", thr = "05pct"))
comparisons <- c("Control vs Asthma Time1", "Control vs Asthma Time2")
models <- list(list(covars = c("age_num", "sex_num"),            label = "Model1_age_sex"),
               list(covars = c("age_num", "sex_num", "bmi_num"), label = "Model2_age_sex_BMI"))

run_maaslin_group <- function(clr_matrix, count_matrix, comparison, covariates, out_dir,
                              tax_level, threshold_label, model_label) {
  groups <- if (comparison == "Control vs Asthma Time1") c("Control", "Asthma Time1") else c("Control", "Asthma Time2")
  meta_sub <- meta_all %>%
    filter(analysis_group %in% groups) %>%
    mutate(group = ifelse(analysis_group == "Control", "Control", "Asthma")) %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  clr_sub <- clr_matrix[meta_sub$sample_id, , drop = FALSE]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- tryCatch(
    Maaslin2(input_data = clr_sub, input_metadata = meta_sub, output = out_dir,
             fixed_effects = c("group", covariates), normalization = "NONE", transform = "NONE",
             analysis_method = "LM", min_prevalence = 0, min_abundance = 0, standardize = FALSE,
             max_significance = 0.25, reference = c("group,Control"),
             plot_heatmap = FALSE, plot_scatter = FALSE),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  all_terms <- read.delim(file.path(out_dir, "all_results.tsv"), stringsAsFactors = FALSE)
  
  ra_sub <- sweep(count_matrix[meta_sub$sample_id, , drop = FALSE], 1,
                  rowSums(count_matrix[meta_sub$sample_id, , drop = FALSE]), "/") * 100
  colnames(ra_sub) <- make.names(colnames(ra_sub), unique = TRUE)
  mean_ra_control <- colMeans(ra_sub[meta_sub$group == "Control", , drop = FALSE])
  mean_ra_asthma  <- colMeans(ra_sub[meta_sub$group == "Asthma",  , drop = FALSE])
  
  all_terms %>%
    filter(metadata == "group") %>%
    rename(qval_maaslin_all_terms = qval) %>%
    mutate(n_tests_maaslin_all_terms = nrow(all_terms),
           qval = p.adjust(pval, method = "BH"),        # BH over taxa, group term only
           n_tests_fdr = sum(!is.na(pval)),
           mean_RA_control = round(mean_ra_control[feature], 4),
           mean_RA_asthma  = round(mean_ra_asthma[feature], 4),
           comparison = comparison, tax_level = tax_level, threshold = threshold_label,
           model_name = model_label, n_used = nrow(meta_sub))
}

all_results <- list(); k <- 1
for (obj in clr_objects) for (comp in comparisons) for (m in models) {
  key <- paste(obj$level, obj$thr, gsub(" ", "_", comp), m$label)
  cat(sprintf("[%d/24] %s\n", k, key)); k <- k + 1
  res <- run_maaslin_group(obj$clr, count_lookup[[obj$count_key]], comp, m$covars,
                           file.path(maaslin_dir, key), obj$level, obj$thr, m$label)
  if (!is.null(res)) all_results[[key]] <- res
}
stopifnot(length(all_results) == 24)

#############################################################
# 2. Combined results
#############################################################
name_lookup_all <- bind_rows(
  lapply(list(genus_count_all, family_count_all, phylum_count_all), function(m)
    data.frame(feature = make.names(colnames(m), unique = TRUE),
               original_name = colnames(m), stringsAsFactors = FALSE))) %>%
  distinct(feature, .keep_all = TRUE)

all_df <- bind_rows(all_results) %>%
  left_join(name_lookup_all, by = "feature") %>%
  mutate(
    feature_clean = ifelse(is.na(original_name), feature, original_name),
    tax_level  = factor(tax_level, levels = c("Phylum", "Family", "Genus")),
    threshold  = factor(threshold, levels = c("01pct", "05pct")),
    comparison = factor(comparison, levels = comparisons),
    model_name = factor(model_name, levels = c("Model1_age_sex", "Model2_age_sex_BMI")),
    significance = case_when(qval < 0.05 ~ "q < 0.05", qval < 0.10 ~ "q < 0.10",
                             pval < 0.05 ~ "p < 0.05", TRUE ~ "NS"),
    significance = factor(significance, levels = c("q < 0.05", "q < 0.10", "p < 0.05", "NS")),
    conf_low  = coef - 1.96 * stderr,
    conf_high = coef + 1.96 * stderr,
    axis_label = case_when(
      tax_level == "Genus" & !grepl("group|Unassigned|ND3007|R-7", feature_clean, ignore.case = TRUE) ~
        paste0("*", feature_clean, "*"),
      TRUE ~ feature_clean),
    ra_annotation = paste0(round(mean_RA_control, 2), "/", round(mean_RA_asthma, 2), "%"))
stopifnot(sum(is.na(all_df$original_name)) == 0)
write_csv(all_df, file.path(tables_dir, "group_DA_MaAsLin2_all_results.csv"))

#############################################################
# 3. Figure 1D (0.1%/10%, age + sex)
#############################################################
status_levels <- c("q < 0.05", "q < 0.10", "p < 0.05", "NS")
make_forest_plot <- function(data, facet_cols_var, title, subtitle) {
  status_colors <- c("q < 0.05" = "#D1495B", "q < 0.10" = "#7296BD", "p < 0.05" = "#C49351", "NS" = "#B8BEC6")
  data <- data %>% mutate(significance = factor(significance, levels = status_levels))
  ggplot(data, aes(x = coef, y = axis_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#929AA5", linewidth = 0.45) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high, color = significance), height = 0.18, linewidth = 0.65) +
    geom_point(aes(shape = significance, fill = significance), size = 3.2, stroke = 0.55, color = "#253247") +
    facet_grid(rows = vars(tax_level), cols = vars(.data[[facet_cols_var]]), scales = "free_y", space = "free_y") +
    scale_shape_manual(name = "Statistical status",
                       values = c("q < 0.05" = 21, "q < 0.10" = 23, "p < 0.05" = 22, "NS" = 21),
                       limits = status_levels, drop = FALSE) +
    scale_fill_manual(name = "Statistical status",
                      values = c("q < 0.05" = "#D1495B", "q < 0.10" = "#7296BD", "p < 0.05" = "#C49351", "NS" = "white"),
                      limits = status_levels, drop = FALSE) +
    scale_color_manual(values = status_colors, limits = status_levels, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_discrete(expand = expansion(add = 0.65)) +
    labs(title = title, subtitle = subtitle, x = "MaAsLin2 coefficient (CLR scale)", y = NULL) +
    guides(color = "none", shape = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = "#253247", hjust = 0, margin = margin(b = 7)),
      plot.subtitle = element_text(size = 9, color = "#626C79", hjust = 0, margin = margin(b = 14)),
      axis.title.x = element_text(size = 10, margin = margin(t = 10)),
      axis.text.x = element_text(color = "#394454", size = 9),
      axis.text.y = ggtext::element_markdown(color = "#253247", size = 10, margin = margin(r = 7)),
      axis.ticks.y = element_blank(),
      strip.background = element_rect(fill = "#F1F4F8", color = "black", linewidth = 0.9),
      strip.text = element_text(face = "bold", size = 11, color = "#253247", margin = margin(8, 8, 8, 8)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.grid.major.x = element_line(color = "#EDF0F4", linewidth = 0.3),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 9),
      legend.key = element_rect(fill = "white", color = NA), legend.key.width = grid::unit(1.1, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(14, 16, 10, 14))
}

main_data <- all_df %>% filter(threshold == "01pct", model_name == "Model1_age_sex")
sig_features_main <- main_data %>% filter(pval < 0.05) %>% pull(feature) %>% unique()
order_key <- main_data %>% filter(feature %in% sig_features_main) %>%
  group_by(tax_level, axis_label) %>% summarise(mean_coef = mean(coef, na.rm = TRUE), .groups = "drop") %>%
  arrange(tax_level, mean_coef)
plot_main <- main_data %>% filter(feature %in% sig_features_main) %>%
  mutate(axis_label = factor(axis_label, levels = unique(order_key$axis_label)))

p_fig1d_main <- make_forest_plot(
  plot_main, "comparison", "Figure 1D · Adjusted microbial group comparisons",
  paste("Adjusted for age and sex", "Precomputed CLR", "Positive coefficients indicate enrichment in asthma", sep = "  |  "))
print(p_fig1d_main)
ggsave(file.path(fig_dir, "Figure1D_group_DA.png"),
       p_fig1d_main, width = 8.08, height = 6.07, dpi = 300, bg = "white")

fig1d_taxa <- plot_main %>% distinct(tax_level, feature_clean) %>% arrange(tax_level, feature_clean)
cat("\nFigure 1D taxa (p < 0.05 in >=1 comparison, primary model):\n")
print(table(fig1d_taxa$tax_level))

#############################################################
# 4. Figure E4 B-D
#############################################################
make_level_panel <- function(level_name, panel_letter) {
  data <- all_df %>% filter(tax_level == level_name, model_name == "Model1_age_sex")
  sig_feat <- data %>% filter(pval < 0.05) %>% pull(feature) %>% unique()
  if (length(sig_feat) == 0) return(NULL)
  order_key <- data %>% filter(feature %in% sig_feat, threshold == "01pct") %>%
    group_by(axis_label) %>% summarise(mean_coef = mean(coef, na.rm = TRUE), .groups = "drop") %>% arrange(mean_coef)
  label_order <- if (nrow(order_key) > 0) order_key$axis_label else unique(data$axis_label[data$feature %in% sig_feat])
  plot_data <- data %>% filter(feature %in% sig_feat) %>%
    mutate(axis_label = factor(axis_label, levels = unique(label_order)),
           threshold_label = ifelse(threshold == "01pct", "0.1%/10%", "0.5%/10%"))
  ggplot(plot_data, aes(x = coef, y = axis_label, color = significance, shape = significance, fill = significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25, linewidth = 0.5) +
    geom_point(size = 2.6, stroke = 0.5, color = "black") +
    facet_grid(rows = vars(comparison), cols = vars(threshold_label), scales = "free_y", space = "free_y") +
    scale_shape_manual(values = c("q < 0.05" = 21, "q < 0.10" = 23, "p < 0.05" = 22, "NS" = 24), drop = FALSE) +
    scale_fill_manual(values = c("q < 0.05" = "#D1495B", "q < 0.10" = "#00798C", "p < 0.05" = "#EDAE49", "NS" = "#FFFFFF"), drop = FALSE) +
    scale_color_manual(values = c("q < 0.05" = "#D1495B", "q < 0.10" = "#00798C", "p < 0.05" = "#EDAE49", "NS" = "gray40"), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.2))) +
    theme_bw(base_size = 9) +
    labs(title = panel_letter, subtitle = paste0(level_name, "-level · age+sex model"),
         x = "MaAsLin2 coefficient (CLR scale)", y = NULL, shape = "Status", fill = "Status", color = "Status") +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8, color = "gray40"),
          axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_markdown(color = "black", size = 7.5),
          strip.background = element_rect(fill = "gray95", color = "gray80"),
          strip.text = element_text(face = "bold", size = 7.5),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.text = element_text(size = 7),
          legend.title = element_text(size = 7.5, face = "bold"),
          panel.spacing = unit(0.4, "lines"))
}
for (lv in list(c("Phylum", "B"), c("Family", "C"), c("Genus", "D"))) {
  p <- make_level_panel(lv[1], lv[2])
  if (!is.null(p)) {
    print(p)
    ggsave(file.path(fig_dir, paste0("FigureE4", lv[2], "_", lv[1], "_group_DA.png")),
           p, width = 230, height = 130, units = "mm", dpi = 300, bg = "white")
  }
}

#############################################################
# 5. Table S2: differential-abundance columns
#############################################################
da_wide <- all_df %>%
  filter(tax_level == "Genus", threshold == "01pct", model_name == "Model1_age_sex") %>%
  mutate(visit = ifelse(comparison == "Control vs Asthma Time1", "T1", "T2")) %>%
  select(original_name = feature_clean, visit, coef, pval, qval) %>%
  pivot_wider(names_from = visit, values_from = c(coef, pval, qval)) %>%
  mutate(DA_direction_T1 = ifelse(coef_T1 > 0, "higher", "lower"),
         DA_direction_T2 = ifelse(coef_T2 > 0, "higher", "lower"),
         DA_nominal_T1 = pval_T1 < 0.05, DA_nominal_T2 = pval_T2 < 0.05,
         DA_FDR_T1 = qval_T1 < 0.05,    DA_FDR_T2 = qval_T2 < 0.05)
stopifnot(nrow(da_wide) == 71, !anyDuplicated(da_wide$original_name))

write_csv(da_wide, file.path(tables_dir, "TableS2_DA_columns.csv"))
cat("\nGenera with P < 0.05 at either visit:", sum(da_wide$DA_nominal_T1 | da_wide$DA_nominal_T2),
    "| FDR q < 0.05:", paste(da_wide$original_name[da_wide$DA_FDR_T1 | da_wide$DA_FDR_T2], collapse = ", "), "\n")

#############################################################
# Reproducibility check against the original results
#############################################################
f_orig <- file.path(repo_dir, "..", "02_result_tables", "v4_20261117",
                    "Fig1D_MaAsLin2_Group_AllResults_v4_20261117.csv")
if (file.exists(f_orig)) {
  o <- read.csv(f_orig)
  m <- inner_join(all_df %>% mutate(across(c(tax_level, threshold, comparison, model_name), as.character)),
                  o, by = c("feature", "tax_level", "threshold", "comparison", "model_name"))
  cat(sprintf("\nComparison with original: %d/%d rows matched | max |diff| coef = %.2e, P = %.2e, q = %.2e\n",
              nrow(m), nrow(all_df), max(abs(m$coef.x - m$coef.y)),
              max(abs(m$pval.x - m$pval.y)), max(abs(m$qval.x - m$qval.y))))
}
cat("Manuscript: Odoribacter (Time 1) is the only taxon with q < 0.05 (q = 0.008)\n")