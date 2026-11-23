#############################################################
# 05_clinical_associations.R
# Associations between genus abundance and clinical outcomes in asthma
#
# Outcomes: ACT, FEV1 % predicted, FEV1/FVC, -log(FeNO + 1),
#           -log(blood eosinophils + 1)  (higher = better for all)
#
#   Figure 4A  pooled mixed-effects model (MaAsLin2, both visits, n = 88)
#              CLR ~ clinical + time + age + sex, random effect = subject
#   Figure 4B  baseline-adjusted change model (lm, 44 patients)
#              dCLR ~ dclinical + baseline CLR + baseline clinical + age + sex
#   Figure E7A/B  same models at the 0.5%/10% threshold
#   Figure E7C    change model without baseline adjustment
#   Figure E8A/B  treatment subgroups: excluding OCS/biologics users (n = 34);
#                 stable GINA step 2-4 (n = 20)
#   Tables S6 (pooled) and S7 (change) are built from these results
#
# FDR: Benjamini-Hochberg over genera within each clinical outcome
# Symbols: * P<0.05, ** P<0.01, *** P<0.001 (unadjusted); dagger = q<0.05
#############################################################

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Maaslin2); library(ggtext); library(compositions); library(tibble)
})

maaslin_root <- file.path(out_dir, "maaslin2")
fig4_dir <- fig_dir; supp_dir <- fig_dir
SIG_LEGEND <- "* p < 0.05   ** p < 0.01   *** p < 0.001   † FDR q < 0.05"

#############################################################
# 0. Data
#############################################################
metadata_clean <- readRDS(f_metadata)
crf_raw        <- readRDS(f_crf)

crf_delta_clean <- crf_raw %>%
  mutate(
    neg_log_FeNO_t1        = -log(FeNO_t1 + 1),
    neg_log_FeNO_t2        = -log(FeNO_t2 + 1),
    neg_log_BloodEos_t1    = -log(Blood_eos_t1 + 1),
    neg_log_BloodEos_t2    = -log(Blood_eos_t2 + 1),
    delta_neg_log_FeNO     = neg_log_FeNO_t2 - neg_log_FeNO_t1,
    delta_neg_log_BloodEos = neg_log_BloodEos_t2 - neg_log_BloodEos_t1
  )

genus_count_01  <- readRDS(file.path(obj_dir, "genus_count_01pct.rds"))
genus_count_05  <- readRDS(file.path(obj_dir, "genus_count_05pct.rds"))
genus_clr_01    <- readRDS(file.path(obj_dir, "genus_clr_01pct.rds"))
genus_clr_05    <- readRDS(file.path(obj_dir, "genus_clr_05pct.rds"))
genus_count_all <- readRDS(file.path(input_dir, "genus_count_all_UNFILTERED_20260922.rds"))

family_count_01  <- readRDS(file.path(obj_dir, "family_count_01pct.rds"))
family_count_05  <- readRDS(file.path(obj_dir, "family_count_05pct.rds"))
family_clr_01    <- readRDS(file.path(obj_dir, "family_clr_01pct.rds"))
family_clr_05    <- readRDS(file.path(obj_dir, "family_clr_05pct.rds"))
family_count_all <- readRDS(file.path(input_dir, "family_count_all_UNFILTERED_20260922.rds"))

phylum_count_01  <- readRDS(file.path(obj_dir, "phylum_count_01pct.rds"))
phylum_count_05  <- readRDS(file.path(obj_dir, "phylum_count_05pct.rds"))
phylum_clr_01    <- readRDS(file.path(obj_dir, "phylum_clr_01pct.rds"))
phylum_clr_05    <- readRDS(file.path(obj_dir, "phylum_clr_05pct.rds"))
phylum_count_all <- readRDS(file.path(input_dir, "phylum_count_all_UNFILTERED_20260922.rds"))

## MaAsLin2 renames features with make.names(); this maps them back
build_name_lookup <- function(count_matrix) {
  data.frame(feature = make.names(colnames(count_matrix), unique = TRUE),
             original_name = colnames(count_matrix), stringsAsFactors = FALSE)
}
name_lookup_all <- bind_rows(build_name_lookup(genus_count_all),
                             build_name_lookup(family_count_all),
                             build_name_lookup(phylum_count_all)) %>%
  distinct(feature, .keep_all = TRUE)
genus_name_map <- setNames(colnames(genus_clr_01), make.names(colnames(genus_clr_01)))

cat("Data loaded. metadata:", nrow(metadata_clean), "| CRF:", nrow(crf_delta_clean), "\n")

#############################################################
# SECTION 1: FIGURE 4A — POOLED MIXED-EFFECTS MaAsLin2
#############################################################
maaslin_4a_dir <- file.path(maaslin_root, "Fig4A")
dir.create(maaslin_4a_dir, recursive = TRUE, showWarnings = FALSE)

meta_pooled <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  mutate(
    sample_id = as.character(sample_id),
    age_num = suppressWarnings(as.numeric(age)),
    sex_num = as.numeric(factor(sex)) - 1,
    bmi_num = suppressWarnings(as.numeric(bmi)),
    time = factor(ifelse(analysis_group == "Asthma Time1", "T1", "T2"), levels = c("T1", "T2"))
  ) %>%
  left_join(
    crf_delta_clean %>%
      select(subject_name_cohort, ACT_t1, ACT_t2, FEV1_pct_t1, FEV1_pct_t2,
             FEV1_FVC_t1, FEV1_FVC_t2, neg_log_FeNO_t1, neg_log_FeNO_t2,
             neg_log_BloodEos_t1, neg_log_BloodEos_t2),
    by = "subject_name_cohort"
  ) %>%
  mutate(
    ACT = ifelse(time == "T1", ACT_t1, ACT_t2),
    FEV1_pct = ifelse(time == "T1", FEV1_pct_t1, FEV1_pct_t2),
    FEV1_FVC = ifelse(time == "T1", FEV1_FVC_t1, FEV1_FVC_t2),
    neg_log_FeNO = ifelse(time == "T1", neg_log_FeNO_t1, neg_log_FeNO_t2),
    neg_log_BloodEos = ifelse(time == "T1", neg_log_BloodEos_t1, neg_log_BloodEos_t2),
    subject_id = subject_name_cohort
  ) %>%
  as.data.frame()
rownames(meta_pooled) <- meta_pooled$sample_id
cat("Pooled metadata:", nrow(meta_pooled), "rows (expect 88)\n")
stopifnot(nrow(meta_pooled) == 88)

count_lookup <- list(
  Genus_01pct = genus_count_01, Genus_05pct = genus_count_05,
  Family_01pct = family_count_01, Family_05pct = family_count_05,
  Phylum_01pct = phylum_count_01, Phylum_05pct = phylum_count_05
)
clr_lookup <- list(
  Genus_01pct = genus_clr_01, Genus_05pct = genus_clr_05,
  Family_01pct = family_clr_01, Family_05pct = family_clr_05,
  Phylum_01pct = phylum_clr_01, Phylum_05pct = phylum_clr_05
)

clinical_vars <- c("ACT", "FEV1_pct", "FEV1_FVC", "neg_log_FeNO", "neg_log_BloodEos")
covariate_models <- list(primary     = c("age_num", "sex_num"),
                         sensitivity = c("age_num", "sex_num", "bmi_num"))

run_maaslin_4a <- function(level, threshold, clin_var, model_name, covariates) {
  key <- paste0(level, "_", threshold)
  clr_data <- clr_lookup[[key]]
  count_data <- count_lookup[[key]]
  
  meta_sub <- meta_pooled %>%
    filter(sample_id %in% rownames(clr_data), !is.na(.data[[clin_var]])) %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  
  clr_sub <- clr_data[meta_sub$sample_id, , drop = FALSE]
  out_dir <- file.path(maaslin_4a_dir, paste(key, clin_var, model_name, sep = "_"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fit <- tryCatch(
    Maaslin2(
      input_data = clr_sub, input_metadata = meta_sub, output = out_dir,
      fixed_effects = c(clin_var, "time", covariates),
      random_effects = c("subject_id"),
      normalization = "NONE", transform = "NONE",
      analysis_method = "LM", min_prevalence = 0, min_abundance = 0,
      reference = c("time,T1"),
      plot_heatmap = FALSE, plot_scatter = FALSE
    ),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  res_path <- file.path(out_dir, "all_results.tsv")
  if (!file.exists(res_path)) return(NULL)
  all_terms <- read.delim(res_path, stringsAsFactors = FALSE)
  
  ra_sub <- sweep(count_data[meta_sub$sample_id, , drop = FALSE], 1,
                  rowSums(count_data[meta_sub$sample_id, , drop = FALSE]), "/") * 100
  colnames(ra_sub) <- make.names(colnames(ra_sub), unique = TRUE)
  mean_ra_all <- colMeans(ra_sub)
  
  all_terms %>%
    filter(metadata == clin_var) %>%
    rename(qval_maaslin_all_terms = qval) %>%
    mutate(
      n_tests_maaslin_all_terms = nrow(all_terms),
      mean_RA_pct = round(mean_ra_all[feature], 4),
      level = level, threshold = threshold, clinical_var = clin_var,
      model = model_name, n_used = nrow(meta_sub)
    )
}

all_results_4a <- list()
run_counter <- 1
total_runs <- length(count_lookup) * length(clinical_vars) * length(covariate_models)
for (key in names(count_lookup)) {
  level_thresh <- strsplit(key, "_")[[1]]
  for (cv in clinical_vars) {
    for (model_name in names(covariate_models)) {
      cat(sprintf("[%d/%d] %s | %s | %s\n", run_counter, total_runs, key, cv, model_name))
      res <- run_maaslin_4a(level_thresh[1], level_thresh[2], cv, model_name,
                            covariate_models[[model_name]])
      if (!is.null(res)) all_results_4a[[length(all_results_4a) + 1]] <- res
      run_counter <- run_counter + 1
    }
  }
}

## q within each clinical outcome (BH over genera), not over all terms
all_df_4a <- bind_rows(all_results_4a) %>%
  left_join(name_lookup_all, by = "feature") %>%
  mutate(feature_clean = ifelse(is.na(original_name), feature, original_name)) %>%
  group_by(level, threshold, clinical_var, model) %>%
  mutate(qval = p.adjust(pval, method = "BH"), n_tests_fdr = sum(!is.na(pval))) %>%
  ungroup()

write.csv(all_df_4a, file.path(tables_dir, "Figure4A_PooledMixedEffects_AllResults.csv"),
          row.names = FALSE)

cat("\n=== Megamonas ~ ACT in 4A (primary, genus 0.1%) ===\n")
all_df_4a %>%
  filter(level == "Genus", threshold == "01pct", model == "primary",
         feature_clean == "Megamonas", clinical_var == "ACT") %>%
  select(feature_clean, clinical_var, coef, stderr, pval, qval, qval_maaslin_all_terms) %>%
  print()
cat("\n>>> Figure 4A statistical runs complete.\n")

#############################################################
# SECTION 2: FIGURE 4B — BASELINE-ADJUSTED LONGITUDINAL CHANGE
#############################################################
paired_subjects <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  group_by(subject_name_cohort) %>%
  filter(n_distinct(analysis_group) == 2) %>%
  pull(subject_name_cohort) %>% unique()
stopifnot(length(paired_subjects) == 44)

meta_t1 <- metadata_clean %>%
  filter(analysis_group == "Asthma Time1", subject_name_cohort %in% paired_subjects) %>%
  mutate(sample_id = as.character(sample_id),
         age_num = suppressWarnings(as.numeric(age)),
         sex_num = as.numeric(factor(sex)) - 1,
         bmi_num = suppressWarnings(as.numeric(bmi))) %>%
  arrange(subject_name_cohort)
meta_t2 <- metadata_clean %>%
  filter(analysis_group == "Asthma Time2", subject_name_cohort %in% paired_subjects) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  arrange(subject_name_cohort)
stopifnot(identical(meta_t1$subject_name_cohort, meta_t2$subject_name_cohort))

clinical_pairs <- list(
  ACT              = list(delta = "delta_ACT",              baseline = "ACT_t1"),
  FEV1_pct         = list(delta = "delta_FEV1_pct",         baseline = "FEV1_pct_t1"),
  FEV1_FVC         = list(delta = "delta_FEV1_FVC",         baseline = "FEV1_FVC_t1"),
  neg_log_FeNO     = list(delta = "delta_neg_log_FeNO",     baseline = "neg_log_FeNO_t1"),
  neg_log_BloodEos = list(delta = "delta_neg_log_BloodEos", baseline = "neg_log_BloodEos_t1")
)

run_delta_models_one_combo <- function(level, threshold, clin_var, model_name, covariates) {
  key <- paste0(level, "_", threshold)
  clr_data <- clr_lookup[[key]]
  count_data <- count_lookup[[key]]
  
  stopifnot(all(c(meta_t1$sample_id, meta_t2$sample_id) %in% rownames(clr_data)))
  t1_clr <- clr_data[meta_t1$sample_id, , drop = FALSE]
  t2_clr <- clr_data[meta_t2$sample_id, , drop = FALSE]
  delta_clr <- t2_clr - t1_clr
  rownames(delta_clr) <- meta_t1$sample_id
  stopifnot(sum(is.na(delta_clr)) == 0)
  
  meta_use <- meta_t1 %>%
    select(subject_name_cohort, sample_id, age_num, sex_num, bmi_num) %>%
    left_join(crf_delta_clean %>%
                select(subject_name_cohort,
                       delta_col    = all_of(clinical_pairs[[clin_var]]$delta),
                       baseline_col = all_of(clinical_pairs[[clin_var]]$baseline)),
              by = "subject_name_cohort") %>%
    filter(!is.na(delta_col), !is.na(baseline_col),
           if_all(all_of(covariates), ~ !is.na(.)))
  n_used <- nrow(meta_use)
  
  genera <- colnames(clr_data)
  results_primary <- vector("list", length(genera))
  results_sens    <- vector("list", length(genera))
  get_term <- function(fit, term) {
    if (is.null(fit) || !(term %in% rownames(summary(fit)$coefficients)))
      return(data.frame(coef = NA, se = NA, pval = NA))
    co <- summary(fit)$coefficients[term, ]
    data.frame(coef = co["Estimate"], se = co["Std. Error"], pval = co["Pr(>|t|)"])
  }
  for (i in seq_along(genera)) {
    g <- genera[i]
    df <- meta_use %>%
      mutate(delta_taxon    = delta_clr[match(sample_id, rownames(delta_clr)), g],
             baseline_taxon = t1_clr[match(sample_id, rownames(t1_clr)), g])
    fit_p <- tryCatch(lm(as.formula(paste("delta_taxon ~ delta_col + baseline_taxon + baseline_col +",
                                          paste(covariates, collapse = " + "))), data = df),
                      error = function(e) NULL)
    fit_s <- tryCatch(lm(as.formula(paste("delta_taxon ~ delta_col +",
                                          paste(covariates, collapse = " + "))), data = df),
                      error = function(e) NULL)
    results_primary[[i]] <- cbind(data.frame(feature = g), get_term(fit_p, "delta_col"))
    results_sens[[i]]    <- cbind(data.frame(feature = g), get_term(fit_s, "delta_col"))
  }
  
  ra_sub <- sweep(count_data, 1, rowSums(count_data), "/") * 100
  colnames(ra_sub) <- make.names(colnames(ra_sub), unique = TRUE)
  mean_ra_all <- colMeans(ra_sub)
  finalize <- function(res_list, model_type) {
    bind_rows(res_list) %>%
      mutate(mean_RA_pct = round(mean_ra_all[make.names(feature)], 4),
             level = level, threshold = threshold, clinical_var = clin_var,
             model = model_name, model_type = model_type, n_used = n_used)
  }
  list(primary     = finalize(results_primary, "baseline_adjusted"),
       sensitivity = finalize(results_sens,    "unadjusted_sensitivity"))
}

all_results_4b_baseline <- list(); all_results_4b_sens <- list()
run_counter <- 1
for (key in names(count_lookup)) {
  level_thresh <- strsplit(key, "_")[[1]]
  for (cv in names(clinical_pairs)) {
    for (model_name in names(covariate_models)) {
      cat(sprintf("[%d/%d] %s | %s | %s\n", run_counter, total_runs, key, cv, model_name))
      out <- run_delta_models_one_combo(level_thresh[1], level_thresh[2], cv, model_name,
                                        covariate_models[[model_name]])
      all_results_4b_baseline[[length(all_results_4b_baseline) + 1]] <- out$primary
      all_results_4b_sens[[length(all_results_4b_sens) + 1]]         <- out$sensitivity
      run_counter <- run_counter + 1
    }
  }
}

clean_names <- function(df) {
  df %>% left_join(name_lookup_all, by = "feature") %>%
    mutate(feature_clean = ifelse(is.na(original_name), feature, original_name))
}
add_fdr_by_combo <- function(df) {
  df %>% group_by(level, threshold, clinical_var, model) %>%
    mutate(qval = p.adjust(pval, method = "BH"), n_tests_fdr = sum(!is.na(pval))) %>%
    ungroup()
}
all_df_4b_baseline <- bind_rows(all_results_4b_baseline) %>% clean_names() %>% add_fdr_by_combo()
all_df_4b_sens     <- bind_rows(all_results_4b_sens)     %>% clean_names() %>% add_fdr_by_combo()

write.csv(all_df_4b_baseline,
          file.path(tables_dir, "Figure4B_BaselineAdjusted_AllResults.csv"), row.names = FALSE)
write.csv(all_df_4b_sens,
          file.path(tables_dir, "Figure4B_Sensitivity_NoBaseline_AllResults.csv"), row.names = FALSE)

cat("\n=== Megamonas ~ ACT, change model (0.1% and 0.5%) ===\n")
all_df_4b_baseline %>%
  filter(level == "Genus", model == "primary", feature_clean == "Megamonas", clinical_var == "ACT") %>%
  select(threshold, coef, se, pval, qval, n_used, n_tests_fdr) %>% print()
cat("\n>>> Figure 4B statistical runs complete.\n")

#############################################################
# SECTION 3: UNIFIED GENUS LIST  (rule: p<0.05 in 4A or 4B + 7 named taxa)
#############################################################
named_genera_of_interest <- c(
  "Lachnospiraceae ND3007 group", "[Eubacterium] eligens group", "Segatella",
  "Megamonas", "Anaerostipes", "Anaerobutyricum", "Blautia"
)
sig_4a <- all_df_4a %>%
  filter(level == "Genus", threshold == "01pct", model == "primary", pval < 0.05) %>%
  pull(feature_clean) %>% unique()
sig_4b <- all_df_4b_baseline %>%
  filter(level == "Genus", threshold == "01pct", model == "primary", pval < 0.05) %>%
  pull(feature_clean) %>% unique()

unified_genus_list <- union(union(sig_4a, sig_4b), named_genera_of_interest)
unified_genus_list <- unified_genus_list[!grepl("^Unassigned", unified_genus_list, ignore.case = TRUE)]
cat("\nUnified genus list:", length(unified_genus_list), "genera\n")
stopifnot(length(setdiff(named_genera_of_interest, colnames(genus_count_01))) == 0)

write.csv(data.frame(Genus = unified_genus_list),
          file.path(tables_dir, "Figure4_UnifiedGenusList.csv"), row.names = FALSE)

fig4a_display          <- all_df_4a          %>% filter(level == "Genus", feature_clean %in% unified_genus_list)
fig4b_baseline_display <- all_df_4b_baseline %>% filter(level == "Genus", feature_clean %in% unified_genus_list)
fig4b_sens_display     <- all_df_4b_sens     %>% filter(level == "Genus", feature_clean %in% unified_genus_list)
write.csv(fig4a_display,          file.path(tables_dir, "Figure4A_UnifiedList_Display.csv"), row.names = FALSE)
write.csv(fig4b_baseline_display, file.path(tables_dir, "Figure4B_UnifiedList_Display.csv"), row.names = FALSE)

#############################################################
# SECTION 4: FIGURE 4 MAIN VISUALIZATION (4A + 4B)
#############################################################
clin_var_labels <- c(
  ACT = "ACT", FEV1_pct = "FEV1%", FEV1_FVC = "FEV1/FVC",
  neg_log_FeNO = "−log(FeNO)", neg_log_BloodEos = "−log(Blood Eos)"
)
ordered_genus_list <- c(
  intersect(named_genera_of_interest, unified_genus_list),
  sort(setdiff(unified_genus_list, named_genera_of_interest))
)
ra_01 <- sweep(genus_count_01, 1, rowSums(genus_count_01), "/") * 100
mean_ra_lookup <- colMeans(ra_01)
make_axis_label <- function(genus) {
  ra_val <- if (genus %in% names(mean_ra_lookup)) round(mean_ra_lookup[genus], 2) else NA
  marker <- if (genus %in% named_genera_of_interest) "● " else ""
  paste0(marker, genus, " (", ra_val, "%)")
}
axis_labels <- setNames(sapply(ordered_genus_list, make_axis_label), ordered_genus_list)

COL_NEG <- "#5C3D6E"; COL_MID <- "#F7F5F0"; COL_POS <- "#A68A56"; NA_TILE <- "#EAE6E1"

build_heatmap_data <- function(df, panel_label) {
  df %>%
    filter(threshold == "01pct", model == "primary", feature_clean %in% ordered_genus_list) %>%
    mutate(
      feature_clean = factor(feature_clean, levels = rev(ordered_genus_list),
                             labels = rev(axis_labels[ordered_genus_list])),
      clinical_label = factor(clin_var_labels[clinical_var], levels = clin_var_labels),
      sig_label = make_sig(pval, qval),
      panel = panel_label
    )
}
heat_4a <- build_heatmap_data(fig4a_display,          "4A: Pooled mixed-effects")
heat_4b <- build_heatmap_data(fig4b_baseline_display, "4B: Baseline-adjusted Δ")

max_abs_4a <- max(abs(heat_4a$coef), na.rm = TRUE)
max_abs_4b <- max(abs(heat_4b$coef), na.rm = TRUE)

make_heatmap_panel <- function(data, title_text, subtitle_text = NULL,
                               max_abs_coef, show_y_axis = TRUE) {
  present_genera <- data %>%
    group_by(feature_clean) %>%
    summarise(any_data = any(!is.na(coef)), .groups = "drop") %>%
    filter(any_data) %>% pull(feature_clean)
  data <- data %>% filter(feature_clean %in% present_genera)
  data$feature_clean <- droplevels(data$feature_clean)
  
  ggplot(data, aes(x = clinical_label, y = feature_clean, fill = coef)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label,
                  color = !is.na(coef) & abs(coef) > 0.55 * max_abs_coef),
              size = 3.8, fontface = "bold", vjust = 0.6) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#222222"), guide = "none") +
    scale_fill_gradient2(
      low = COL_NEG, mid = COL_MID, high = COL_POS, midpoint = 0,
      limits = c(-max_abs_coef, max_abs_coef), na.value = NA_TILE,
      name = "Regression\ncoefficient",
      guide = guide_colorbar(barwidth = unit(3.8, "cm"), barheight = unit(0.45, "cm"),
                             title.position = "top", title.hjust = 0.5,
                             ticks.colour = "grey50", frame.colour = "grey50",
                             frame.linewidth = 0.5)
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme_bw(base_size = 9.5) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white"),
      panel.border     = element_rect(color = "#333333", fill = NA, linewidth = 0.9),
      panel.grid       = element_blank(),
      plot.title       = element_text(face = "bold", size = 10.5, hjust = 0.5,
                                      color = "#111111", margin = margin(b = 3)),
      plot.subtitle    = element_text(size = 7.5, hjust = 0.5, color = "grey50",
                                      margin = margin(b = 4)),
      plot.margin      = margin(8, 12, 6, 12),
      axis.text.x      = element_text(angle = 30, hjust = 1, color = "#111111",
                                      size = 8.5, face = "bold"),
      axis.text.y      = if (show_y_axis) element_text(color = "#111111", size = 8.2,
                                                       face = "plain", margin = margin(r = 3))
      else element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.margin    = margin(t = 6, b = 6),
      legend.title     = element_text(size = 8.5, face = "bold", color = "#111111"),
      legend.text      = element_text(size = 8, color = "#333333")
    )
}

p_4a <- make_heatmap_panel(heat_4a, title_text = "4A: Pooled Mixed-Effects Model",
                           subtitle_text = "CLR ~ clinical + time + age + sex  |  random = subject",
                           max_abs_coef = max_abs_4a, show_y_axis = TRUE)
p_4b <- make_heatmap_panel(heat_4b, title_text = "4B: Baseline-Adjusted Longitudinal Change",
                           subtitle_text = "ΔCLR ~ Δclinical + baseline CLR + baseline clinical + age + sex",
                           max_abs_coef = max_abs_4b, show_y_axis = TRUE)

p_fig4_main <- (p_4a | p_4b) + plot_layout(guides = "keep", widths = c(1, 1)) +
  plot_annotation(
    title    = "Figure 4 — Gut Microbiome Associations with Asthma Clinical Markers",
    subtitle = paste0(
      "Genus level · 0.1% prevalence / 10% abundance thresholds · primary model (age + sex)\n",
      "Unified genus list: 4A/4B nominal hits (p < 0.05) ∪ 7 pre-specified taxa (●)\n",
      "BRONZE/GOLD = taxon higher with better clinical status (favorable)  |  ",
      "PURPLE = taxon higher with worse clinical status (unfavorable)\n",
      SIG_LEGEND, "   |   Colour scales are panel-specific (independent symmetric limits)"),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                     color = "#1A1A2E", margin = margin(b = 4)),
      plot.subtitle   = element_text(size = 8, color = "grey40", hjust = 0.5,
                                     lineheight = 1.4, margin = margin(b = 8)),
      plot.margin     = margin(14, 14, 10, 14)))
print(p_fig4_main)

ggsave(file.path(fig4_dir, "Figure4A_Pooled.png"),
       p_4a, width = 4.84, height = 6.04, dpi = 300, bg = "white")
ggsave(file.path(fig4_dir, "Figure4B_Change.png"),
       p_4b, width = 4.84, height = 6.04, dpi = 300, bg = "white")
cat("✓ Figure 4A and 4B saved\n")

#############################################################
# SECTION 5: FIGURE E7A/E7B (0.5% threshold)
#############################################################
ra_05 <- sweep(genus_count_05, 1, rowSums(genus_count_05), "/") * 100
mean_ra_lookup_05 <- colMeans(ra_05)
make_axis_label_05 <- function(genus) {
  ra_val <- if (genus %in% names(mean_ra_lookup_05)) round(mean_ra_lookup_05[genus], 2) else NA
  marker <- if (genus %in% named_genera_of_interest) "● " else ""
  label_text <- if (is.na(ra_val)) paste0(genus, " (not in 0.5% table)") else paste0(genus, " (", ra_val, "%)")
  paste0(marker, label_text)
}
axis_labels_05 <- setNames(sapply(ordered_genus_list, make_axis_label_05), ordered_genus_list)
n_dropped_at_05 <- sum(!ordered_genus_list %in% colnames(genus_count_05))

build_heatmap_data_05 <- function(df) {
  df %>%
    filter(threshold == "05pct", model == "primary", feature_clean %in% ordered_genus_list) %>%
    mutate(
      feature_clean = factor(feature_clean, levels = rev(ordered_genus_list),
                             labels = rev(axis_labels_05[ordered_genus_list])),
      clinical_label = factor(clin_var_labels[clinical_var], levels = clin_var_labels),
      sig_label = make_sig(pval, qval)
    )
}
heat_4a_05 <- build_heatmap_data_05(fig4a_display)
heat_4b_05 <- build_heatmap_data_05(fig4b_baseline_display)
max_abs_4a_05 <- max(abs(heat_4a_05$coef), na.rm = TRUE)
max_abs_4b_05 <- max(abs(heat_4b_05$coef), na.rm = TRUE)

p_4a_05 <- make_heatmap_panel(heat_4a_05, title_text = "4A: Pooled Mixed-Effects (0.5%/10%)",
                              max_abs_coef = max_abs_4a_05, show_y_axis = TRUE)
p_4b_05 <- make_heatmap_panel(heat_4b_05, title_text = "4B: Baseline-Adjusted Δ (0.5%/10%)",
                              max_abs_coef = max_abs_4b_05, show_y_axis = TRUE)
print((p_4a_05 | p_4b_05) + plot_layout(guides = "keep", widths = c(1, 1)))

ggsave(file.path(supp_dir, "FigureE7A_Pooled_05pct.png"),
       p_4a_05, width = 4.37, height = 4.69, dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureE7B_Change_05pct.png"),
       p_4b_05, width = 4.37, height = 4.69, dpi = 300, bg = "white")
cat("✓ Figure E7A and E7B saved\n")

#############################################################
# SECTION 7: FIGURE E7C — change model with vs without baseline adjustment
#############################################################
comparison_table <- all_df_4b_baseline %>%
  filter(level == "Genus", threshold == "01pct", model == "primary", feature_clean %in% ordered_genus_list) %>%
  select(feature_clean, clinical_var, coef_adjusted = coef, se_adjusted = se,
         pval_adjusted = pval, qval_adjusted = qval) %>%
  inner_join(all_df_4b_sens %>%
               filter(level == "Genus", threshold == "01pct", model == "primary",
                      feature_clean %in% ordered_genus_list) %>%
               select(feature_clean, clinical_var, coef_unadjusted = coef, se_unadjusted = se,
                      pval_unadjusted = pval, qval_unadjusted = qval),
             by = c("feature_clean", "clinical_var")) %>%
  mutate(
    sig_adjusted = pval_adjusted < 0.05, sig_unadjusted = pval_unadjusted < 0.05,
    agreement = case_when(
      sig_adjusted & sig_unadjusted  ~ "Significant in both",
      sig_adjusted & !sig_unadjusted ~ "Only significant WITH baseline adjustment",
      !sig_adjusted & sig_unadjusted ~ "Only significant WITHOUT baseline adjustment",
      TRUE ~ "NS in both"),
    coef_direction_flip = sign(coef_adjusted) != sign(coef_unadjusted))
write.csv(comparison_table, file.path(tables_dir, "Figure4B_SensitivityComparison.csv"),
          row.names = FALSE)
print(table(comparison_table$agreement))

build_sens_heatmap_data <- function(coef_col, pval_col, qval_col) {
  comparison_table %>%
    transmute(
      feature_clean = factor(feature_clean, levels = rev(ordered_genus_list),
                             labels = rev(axis_labels[ordered_genus_list])),
      clinical_label = factor(clin_var_labels[clinical_var], levels = clin_var_labels),
      coef = .data[[coef_col]], pval = .data[[pval_col]], qval = .data[[qval_col]],
      sig_label = make_sig(pval, qval))
}
heat_adj   <- build_sens_heatmap_data("coef_adjusted",   "pval_adjusted",   "qval_adjusted")
heat_unadj <- build_sens_heatmap_data("coef_unadjusted", "pval_unadjusted", "qval_unadjusted")
max_abs_adj   <- max(abs(heat_adj$coef),   na.rm = TRUE)
max_abs_unadj <- max(abs(heat_unadj$coef), na.rm = TRUE)

p_adj <- make_heatmap_panel(heat_adj, title_text = "Baseline-Adjusted (primary model)",
                            subtitle_text = "ΔCLR ~ Δclinical + baseline CLR + baseline clinical + age + sex",
                            max_abs_coef = max_abs_adj, show_y_axis = TRUE)
p_unadj <- make_heatmap_panel(heat_unadj, title_text = "Unadjusted (sensitivity model)",
                              subtitle_text = "ΔCLR ~ Δclinical + age + sex",
                              max_abs_coef = max_abs_unadj, show_y_axis = TRUE)
print((p_adj | p_unadj) + plot_layout(guides = "keep", widths = c(1, 1)))
ggsave(file.path(supp_dir, "FigureE7C_Change_NoBaseline.png"),
       p_unadj, width = 4.17, height = 5.14, dpi = 300, bg = "white")
cat("✓ Figure E7C saved\n")

#############################################################
# SECTION 8: FIGURE E8 — TREATMENT SENSITIVITY
#############################################################
covariates <- c("age_num", "sex_num")

## Part A: samples of patients who received OCS or biologics at either visit
## (CRF variables regular_OCS and Biologics); 10 of the 44 paired patients are
## excluded, 34 remain.
ocs_bio_samples_hyphens <- c(
  "A-6-0",   "A-6-24",
  "PS-7-0",  "PS-7-12",
  "A-13-12", "A-13-24",
  "PS-4-0",  "PS-5-0",
  "A-21-12", "PS-8-12",
  "A-23-0",  "A-23-12",
  "A-27-0",  "A-27-24",
  "A-40-0",  "A-40-12",
  "PS-6-0",  "PS-6-12",
  "MA-10-0", "PS-9-1"
)
ocs_bio_samples <- gsub("-", "_", ocs_bio_samples_hyphens)
stopifnot("an OCS/biologics sample ID is not in metadata" =
            all(ocs_bio_samples %in% as.character(metadata_clean$sample_id)))

## Part B: stable treatment step 2-4 (step unchanged T1 -> T2)
crf_stable <- crf_delta_clean %>%
  filter(!is.na(treatment_step_t1), !is.na(treatment_step_t2),
         treatment_step_t1 == treatment_step_t2, treatment_step_t1 %in% 2:4)
stable_t1_ids <- gsub("-", "_", crf_stable$sample_id_t1)
stable_t2_ids <- gsub("-", "_", crf_stable$sample_id_t2)
stopifnot(nrow(crf_stable) == 20,
          all(c(stable_t1_ids, stable_t2_ids) %in% as.character(metadata_clean$sample_id)))

run_sensitivity_4A <- function(meta_filtered, label, out_subdir) {
  out_dir <- file.path(maaslin_root, out_subdir, "Fig4A")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta_pooled_sub <- meta_filtered %>%
    filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
    mutate(sample_id = as.character(sample_id),
           age_num = suppressWarnings(as.numeric(age)),
           sex_num = as.numeric(factor(sex)) - 1,
           time = factor(ifelse(analysis_group == "Asthma Time1", "T1", "T2"), levels = c("T1", "T2"))) %>%
    left_join(crf_delta_clean %>%
                select(subject_name_cohort, ACT_t1, ACT_t2, FEV1_pct_t1, FEV1_pct_t2,
                       FEV1_FVC_t1, FEV1_FVC_t2, neg_log_FeNO_t1, neg_log_FeNO_t2,
                       neg_log_BloodEos_t1, neg_log_BloodEos_t2),
              by = "subject_name_cohort") %>%
    mutate(ACT = ifelse(time == "T1", ACT_t1, ACT_t2),
           FEV1_pct = ifelse(time == "T1", FEV1_pct_t1, FEV1_pct_t2),
           FEV1_FVC = ifelse(time == "T1", FEV1_FVC_t1, FEV1_FVC_t2),
           neg_log_FeNO = ifelse(time == "T1", neg_log_FeNO_t1, neg_log_FeNO_t2),
           neg_log_BloodEos = ifelse(time == "T1", neg_log_BloodEos_t1, neg_log_BloodEos_t2),
           subject_id = subject_name_cohort) %>%
    as.data.frame()
  rownames(meta_pooled_sub) <- meta_pooled_sub$sample_id
  
  all_res <- list()
  for (cv in clinical_vars) {
    meta_sub <- meta_pooled_sub %>%
      filter(sample_id %in% rownames(genus_clr_01), !is.na(.data[[cv]]), !is.na(age_num), !is.na(sex_num))
    clr_sub <- genus_clr_01[meta_sub$sample_id, , drop = FALSE]
    fit <- tryCatch(
      Maaslin2(input_data = clr_sub, input_metadata = meta_sub, output = file.path(out_dir, cv),
               fixed_effects = c(cv, "time", covariates), random_effects = c("subject_id"),
               normalization = "NONE", transform = "NONE", analysis_method = "LM",
               min_prevalence = 0, min_abundance = 0, reference = c("time,T1"),
               plot_heatmap = FALSE, plot_scatter = FALSE),
      error = function(e) { cat("ERROR in 4A", cv, ":", conditionMessage(e), "\n"); NULL })
    if (!is.null(fit)) {
      all_terms <- read.delim(file.path(out_dir, cv, "all_results.tsv"), stringsAsFactors = FALSE)
      all_res[[cv]] <- all_terms %>% filter(metadata == cv) %>%
        rename(qval_maaslin_all_terms = qval) %>%
        mutate(qval = p.adjust(pval, method = "BH"),
               n_tests_fdr = sum(!is.na(pval)),
               clinical_var = cv, sensitivity = label, model_type = "Fig4A",
               n_used = nrow(meta_sub),
               n_subjects = n_distinct(meta_sub$subject_id))
    }
  }
  bind_rows(all_res)
}

run_sensitivity_4B <- function(meta_filtered, crf_filtered, label) {
  paired <- meta_filtered %>%
    filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
    group_by(subject_name_cohort) %>% filter(n_distinct(analysis_group) == 2) %>%
    pull(subject_name_cohort) %>% unique()
  m_t1 <- meta_filtered %>%
    filter(analysis_group == "Asthma Time1", subject_name_cohort %in% paired) %>%
    mutate(sample_id = as.character(sample_id),
           age_num = suppressWarnings(as.numeric(age)), sex_num = as.numeric(factor(sex)) - 1) %>%
    arrange(subject_name_cohort)
  m_t2 <- meta_filtered %>%
    filter(analysis_group == "Asthma Time2", subject_name_cohort %in% paired) %>%
    mutate(sample_id = as.character(sample_id)) %>% arrange(subject_name_cohort)
  stopifnot(identical(m_t1$subject_name_cohort, m_t2$subject_name_cohort))
  cat(sprintf("\n[%s] Fig4B: N=%d paired subjects\n", label, nrow(m_t1)))
  
  t1_clr <- genus_clr_01[m_t1$sample_id, , drop = FALSE]
  t2_clr <- genus_clr_01[m_t2$sample_id, , drop = FALSE]
  delta_clr <- t2_clr - t1_clr
  rownames(delta_clr) <- m_t1$sample_id
  genera <- colnames(genus_clr_01)
  
  all_res <- list()
  for (cv in clinical_vars) {
    meta_use <- m_t1 %>%
      select(subject_name_cohort, sample_id, age_num, sex_num) %>%
      left_join(crf_filtered %>%
                  select(subject_name_cohort,
                         delta_col = all_of(clinical_pairs[[cv]]$delta),
                         baseline_col = all_of(clinical_pairs[[cv]]$baseline)),
                by = "subject_name_cohort") %>%
      filter(!is.na(delta_col), !is.na(baseline_col), !is.na(age_num), !is.na(sex_num))
    n_used <- nrow(meta_use)
    res <- lapply(genera, function(g) {
      df <- meta_use %>% mutate(delta_taxon = delta_clr[match(sample_id, rownames(delta_clr)), g],
                                baseline_taxon = t1_clr[match(sample_id, rownames(t1_clr)), g])
      fit <- tryCatch(lm(as.formula(paste("delta_taxon ~ delta_col + baseline_taxon + baseline_col +",
                                          paste(covariates, collapse = " + "))), data = df),
                      error = function(e) NULL)
      if (!is.null(fit) && "delta_col" %in% rownames(summary(fit)$coefficients)) {
        co <- summary(fit)$coefficients["delta_col", ]
        data.frame(feature = g, coef = co[["Estimate"]], se = co[["Std. Error"]], pval = co[["Pr(>|t|)"]])
      } else data.frame(feature = g, coef = NA, se = NA, pval = NA)
    })
    all_res[[cv]] <- bind_rows(res) %>%
      mutate(n_used = n_used, clinical_var = cv, sensitivity = label, model_type = "Fig4B",
             qval = p.adjust(pval, method = "BH"), n_tests_fdr = sum(!is.na(pval)))
  }
  bind_rows(all_res)
}

meta_A <- metadata_clean %>% filter(!sample_id %in% ocs_bio_samples)
crf_A  <- crf_delta_clean %>%
  filter(!(gsub("-", "_", sample_id_t1) %in% ocs_bio_samples) &
           !(gsub("-", "_", sample_id_t2) %in% ocs_bio_samples))
meta_B <- metadata_clean %>% filter(sample_id %in% c(stable_t1_ids, stable_t2_ids))
crf_B  <- crf_stable

res_A_4a <- run_sensitivity_4A(meta_A, "excl_OCS_Bio",  "sensitivity_excl_ocs_bio")
res_A_4b <- run_sensitivity_4B(meta_A, crf_A, "excl_OCS_Bio")
res_B_4a <- run_sensitivity_4A(meta_B, "stable_step24", "sensitivity_stable_step24")
res_B_4b <- run_sensitivity_4B(meta_B, crf_B, "stable_step24")

## map MaAsLin2 names back with make.names() (gsub(".", " ") lost 6 genera)
clean_feature <- function(df) {
  df %>% mutate(feature_clean = ifelse(feature %in% names(genus_name_map),
                                       unname(genus_name_map[feature]), feature))
}
all_4a <- bind_rows(res_A_4a, res_B_4a) %>% clean_feature()
all_4b <- bind_rows(res_A_4b, res_B_4b) %>% clean_feature()
stopifnot(all(all_4a$feature_clean %in% colnames(genus_clr_01)),
          all(all_4b$feature_clean %in% colnames(genus_clr_01)))
write.csv(all_4a, file.path(tables_dir, "SupplementaryTable_Sensitivity2_Fig4A.csv"), row.names = FALSE)
write.csv(all_4b, file.path(tables_dir, "SupplementaryTable_Sensitivity2_Fig4B.csv"), row.names = FALSE)

## n taken from the data (34 / 20 paired patients)
nA <- max(res_A_4b$n_used); nB <- max(res_B_4b$n_used)
stopifnot(nA == 34, nB == 20)
LAB_MAIN <- "Main Analysis\n(N=44)"
LAB_A    <- sprintf("Excl OCS/Biologics\n(N=%d)", nA)
LAB_B    <- sprintf("Stable Step 2–4\n(N=%d)", nB)

main_4a <- all_df_4a %>% filter(level == "Genus", threshold == "01pct", model == "primary",
                                feature_clean %in% ordered_genus_list)
main_4b <- all_df_4b_baseline %>% filter(level == "Genus", threshold == "01pct", model == "primary",
                                         feature_clean %in% ordered_genus_list)
reshape_sens <- function(df, sens_label)
  df %>% filter(sensitivity == sens_label, feature_clean %in% ordered_genus_list)
sensA_4a <- reshape_sens(all_4a, "excl_OCS_Bio");  sensA_4b <- reshape_sens(all_4b, "excl_OCS_Bio")
sensB_4a <- reshape_sens(all_4a, "stable_step24"); sensB_4b <- reshape_sens(all_4b, "stable_step24")

build_sens_heatmap <- function(df, panel_label) {
  df %>%
    filter(feature_clean %in% ordered_genus_list) %>%
    mutate(feature_clean = factor(feature_clean, levels = rev(ordered_genus_list),
                                  labels = rev(axis_labels[ordered_genus_list])),
           clinical_label = factor(clin_var_labels[clinical_var], levels = clin_var_labels),
           sig_label = make_sig(pval, qval),
           panel = panel_label)
}

e8_panels <- function(d_main, d_A, d_B) {
  h_main <- build_sens_heatmap(d_main, "Main"); h_A <- build_sens_heatmap(d_A, "A")
  h_B    <- build_sens_heatmap(d_B, "B")
  m <- max(abs(c(h_main$coef, h_A$coef, h_B$coef)), na.rm = TRUE)
  (make_heatmap_panel(h_main, title_text = LAB_MAIN, max_abs_coef = m, show_y_axis = TRUE) |
      make_heatmap_panel(h_A,    title_text = LAB_A,    max_abs_coef = m, show_y_axis = FALSE) |
      make_heatmap_panel(h_B,    title_text = LAB_B,    max_abs_coef = m, show_y_axis = FALSE)) +
    plot_layout(guides = "collect", widths = c(1.3, 1, 1)) &
    theme(legend.position = "bottom")
}
e8_annot <- function(title, line2) plot_annotation(
  title = title,
  subtitle = paste0(line2,
                    "BRONZE/GOLD = taxon higher with better clinical status (favorable)  |  ",
                    "PURPLE = taxon higher with worse clinical status (unfavorable)\n",
                    "Shared color scale across all three panels\n", SIG_LEGEND),
  theme = theme(plot.background = element_rect(fill = "white", color = NA),
                plot.title = element_text(face = "bold", size = 12, hjust = 0.5,
                                          color = "#1A1A2E", margin = margin(b = 4)),
                plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5,
                                             lineheight = 1.4, margin = margin(b = 8)),
                plot.margin = margin(12, 12, 10, 12)))

p_supp_4a <- e8_panels(main_4a, sensA_4a, sensB_4a) +
  e8_annot("Supplementary Figure — Treatment Sensitivity: Figure 4A (Pooled Mixed-Effects)",
           "Genus level · 0.1%/10% threshold · primary model (age + sex) · random effect = subject\n")
ggsave(file.path(supp_dir, "FigureE8A_TreatmentSensitivity_Pooled.png"),
       p_supp_4a, width = 12, height = 7.5, dpi = 300, bg = "white")

p_supp_4b <- e8_panels(main_4b, sensA_4b, sensB_4b) +
  e8_annot("Supplementary Figure — Treatment Sensitivity: Figure 4B (Baseline-Adjusted Δ)",
           paste0("Genus level · 0.1%/10% threshold · primary model (age + sex)\n",
                  "ΔCLR ~ Δclinical + baseline CLR + baseline clinical + age + sex\n"))
ggsave(file.path(supp_dir, "FigureE8B_TreatmentSensitivity_Change.png"),
       p_supp_4b, width = 12, height = 7.5, dpi = 300, bg = "white")
cat("✓ Figure E8A and E8B saved\n")

## key-taxa table (same layout as before)
extract_key <- function(df, model_type, sens_label) {
  df %>% filter(feature_clean %in% named_genera_of_interest) %>%
    transmute(Genus = feature_clean, Clinical_var = clinical_var, Model = model_type,
              Sensitivity = sens_label, Coef = round(coef, 4), P_value = signif(pval, 3),
              Q_value = signif(qval, 3), Sig = make_sig(pval, qval))
}
supp_table <- bind_rows(
  extract_key(main_4a,  "Fig4A (Pooled MixedFX)", "Main (N=44)"),
  extract_key(sensA_4a, "Fig4A (Pooled MixedFX)", sprintf("Excl OCS/Bio (N=%d)", nA)),
  extract_key(sensB_4a, "Fig4A (Pooled MixedFX)", sprintf("Stable Step2-4 (N=%d)", nB)),
  extract_key(main_4b,  "Fig4B (Delta model)",    "Main (N=44)"),
  extract_key(sensA_4b, "Fig4B (Delta model)",    sprintf("Excl OCS/Bio (N=%d)", nA)),
  extract_key(sensB_4b, "Fig4B (Delta model)",    sprintf("Stable Step2-4 (N=%d)", nB))
) %>% arrange(Genus, Model, Clinical_var, Sensitivity)
write.csv(supp_table, file.path(tables_dir, "SupplementaryTable_Sensitivity2_Treatment_KeyTaxa.csv"),
          row.names = FALSE)

#############################################################
# Reproducibility check against the original results
#############################################################
orig_dir <- file.path(repo_dir, "..", "02_result_tables", "v4_20261117")
compare_results <- function(new, file, keys, label) {
  f <- file.path(orig_dir, file)
  if (!file.exists(f)) { cat(sprintf("  %-34s original not found\n", label)); return(invisible()) }
  o <- read.csv(f, stringsAsFactors = FALSE)
  new <- new %>% mutate(across(all_of(keys), as.character))
  o   <- o   %>% mutate(across(all_of(keys), as.character))
  m <- inner_join(new, o, by = keys)
  cat(sprintf("  %-34s %d/%d rows | max |diff| coef %.1e, P %.1e, q %.1e\n", label,
              nrow(m), nrow(new), max(abs(m$coef.x - m$coef.y), na.rm = TRUE),
              max(abs(m$pval.x - m$pval.y), na.rm = TRUE), max(abs(m$qval.x - m$qval.y), na.rm = TRUE)))
}
cat("\nComparison with original results:\n")
k <- c("level", "threshold", "clinical_var", "model", "feature")
compare_results(all_df_4a,          "Figure4A_PooledMixedEffects_AllResults_v4_20261117.csv",    k, "Pooled model (4A, E7A)")
compare_results(all_df_4b_baseline, "Figure4B_BaselineAdjusted_AllResults_v4_20261117.csv",      k, "Change model (4B, E7B)")
compare_results(all_df_4b_sens,     "Figure4B_Sensitivity_NoBaseline_AllResults_v4_20261117.csv", k, "No baseline adjustment (E7C)")
k8 <- c("sensitivity", "clinical_var", "feature")
compare_results(all_4a, "SupplementaryTable_Sensitivity2_Fig4A_v4_20261117.csv", k8, "Treatment subgroups, pooled (E8A)")
compare_results(all_4b, "SupplementaryTable_Sensitivity2_Fig4B_v4_20261117.csv", k8, "Treatment subgroups, change (E8B)")