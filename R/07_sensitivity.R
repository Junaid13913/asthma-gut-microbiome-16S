#############################################################
# 07_sensitivity.R
# Sensitivity analyses for the patient-control comparison (Time 1 / Time 2)
#
#   A. Excluded-9: 9 Time 1 samples reported in a previous study removed
#      (35 patients vs 41 controls). LEfSe step 1 (Kruskal-Wallis on relative
#      abundance) for Lachnospiraceae, and MaAsLin2 with the Figure 1D model
#      (family + genus)                                   -> Table S8
#   B. Sample map: sample, de-identified subject code, visit, sequencing run
#                                                         -> Table S9
#   C. Run-adjusted differential abundance: Figure 1D model in Runs 2-3 only,
#      with and without run as a covariate                -> Table S10
#
# MaAsLin2 specification as Figure 1D: group factor (reference = Control),
# age + sex, pre-computed CLR, standardize = FALSE, BH over taxa (group term)
#############################################################

source("R/00_setup.R")
suppressPackageStartupMessages(library(Maaslin2))

metadata_clean   <- readRDS(f_metadata)
family_count_all <- readRDS(f_unfiltered[["Family"]])
family_count_01  <- readRDS(file.path(obj_dir, "family_count_01pct.rds"))
family_clr_01    <- readRDS(file.path(obj_dir, "family_clr_01pct.rds"))
genus_count_01   <- readRDS(file.path(obj_dir, "genus_count_01pct.rds"))
genus_clr_01     <- readRDS(file.path(obj_dir, "genus_clr_01pct.rds"))

## Figure 1D results from 04_differential_abundance.R (used for validation)
f1d_all <- read.csv(file.path(tables_dir, "group_DA_MaAsLin2_all_results.csv"), stringsAsFactors = FALSE)
fmtp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                           ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p)))

#############################################################
# A. Excluded-9
#############################################################
maaslin_dir <- file.path(out_dir, "maaslin2", "excluded9"); dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)
exclude_samples <- c("PS_1_0", "A_6_0", "A_7_0", "PS_4_0", "A_19_0",
                     "A_22_0", "A_23_0", "A_27_0", "A_32_0")
stopifnot(all(exclude_samples %in% as.character(metadata_clean$sample_id)))

prep <- function(m) m %>%
  mutate(sample_id = as.character(sample_id),
         group = factor(ifelse(analysis_group == "Control", "Control", "Asthma"), levels = c("Control", "Asthma")),
         age_num = suppressWarnings(as.numeric(age)),
         sex_num = as.numeric(factor(sex)) - 1) %>% as.data.frame()
sets <- list(
  Full_T1  = prep(metadata_clean %>% filter(analysis_group %in% c("Control", "Asthma Time1"))),
  Excl9_T1 = prep(metadata_clean %>% filter(analysis_group %in% c("Control", "Asthma Time1"),
                                            !sample_id %in% exclude_samples)))
for (s in names(sets)) rownames(sets[[s]]) <- sets[[s]]$sample_id
print(sapply(sets, function(m) table(m$group)))
stopifnot(sum(sets$Excl9_T1$group == "Asthma") == 35, sum(sets$Excl9_T1$group == "Control") == 41)

## ---------------------------------------------------------------------------
## 1. LEfSe step 1 in R (Kruskal-Wallis on relative abundance), Lachnospiraceae
## ---------------------------------------------------------------------------
lefse_kw <- function(meta, counts, taxon, label) {
  ra <- counts[meta$sample_id, , drop = FALSE]
  ra <- sweep(ra, 1, rowSums(ra), "/") * 100
  x <- ra[, taxon]
  data.frame(set = label, taxon = taxon,
             n_control = sum(meta$group == "Control"), n_asthma = sum(meta$group == "Asthma"),
             mean_RA_control = round(mean(x[meta$group == "Control"]), 2),
             mean_RA_asthma  = round(mean(x[meta$group == "Asthma"]), 2),
             P_KW = signif(kruskal.test(x ~ meta$group)$p.value, 3))
}
lefse_res <- bind_rows(
  lefse_kw(sets$Full_T1,  family_count_all, "Lachnospiraceae", "Full T1 (unfiltered table)"),
  lefse_kw(sets$Full_T1,  family_count_01,  "Lachnospiraceae", "Full T1 (0.1% table)"),
  lefse_kw(sets$Excl9_T1, family_count_all, "Lachnospiraceae", "Excl-9 T1 (unfiltered table)"),
  lefse_kw(sets$Excl9_T1, family_count_01,  "Lachnospiraceae", "Excl-9 T1 (0.1% table)"))
cat("\n=== LEfSe step 1 (KW). VALIDATION: Full T1 should match Figure 1C (24% vs 29%, P = 0.002) ===\n")
print(lefse_res, row.names = FALSE)

## LEfSe input (Galaxy) for the excluded-9 set, family level
ra_ex9 <- family_count_all[sets$Excl9_T1$sample_id, , drop = FALSE]
ra_ex9 <- sweep(ra_ex9, 1, rowSums(ra_ex9), "/")
lefse_in <- rbind(class = as.character(sets$Excl9_T1$group),
                  subject = sets$Excl9_T1$sample_id,
                  t(round(ra_ex9, 8)))
write.table(lefse_in, file.path(tables_dir, "LEfSe_input_Family_Excl9_T1.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)

run_ma <- function(meta, clr, key) {
  out <- file.path(maaslin_dir, key)
  fit <- Maaslin2(input_data = clr[meta$sample_id, , drop = FALSE], input_metadata = meta, output = out,
                  fixed_effects = c("group", "age_num", "sex_num"), normalization = "NONE",
                  transform = "NONE", analysis_method = "LM", min_prevalence = 0, min_abundance = 0,
                  standardize = FALSE, max_significance = 0.25, reference = c("group,Control"),
                  plot_heatmap = FALSE, plot_scatter = FALSE)
  fit$results %>% filter(metadata == "group") %>%
    rename(qval_maaslin_all_terms = qval) %>%
    mutate(qval = p.adjust(pval, "BH"), n_tests_fdr = n(), set = key)
}
fam_map <- setNames(colnames(family_clr_01), make.names(colnames(family_clr_01)))
gen_map <- setNames(colnames(genus_clr_01),  make.names(colnames(genus_clr_01)))
ma <- bind_rows(
  run_ma(sets$Full_T1,  family_clr_01, "Family_Full_T1")  %>% mutate(level = "Family", taxon = fam_map[feature]),
  run_ma(sets$Excl9_T1, family_clr_01, "Family_Excl9_T1") %>% mutate(level = "Family", taxon = fam_map[feature]),
  run_ma(sets$Full_T1,  genus_clr_01,  "Genus_Full_T1")   %>% mutate(level = "Genus",  taxon = gen_map[feature]),
  run_ma(sets$Excl9_T1, genus_clr_01,  "Genus_Excl9_T1")  %>% mutate(level = "Genus",  taxon = gen_map[feature]))
stopifnot(!anyNA(ma$taxon))
write_csv(ma, file.path(tables_dir, "excluded9_MaAsLin2_family_genus.csv"))

## validation: full-cohort runs must reproduce Figure 1D
f1d <- f1d_all %>%
  filter(model_name == "Model1_age_sex", threshold == "01pct", comparison == "Control vs Asthma Time1",
         tax_level %in% c("Family", "Genus"))
chk <- ma %>% filter(grepl("Full", set)) %>%
  inner_join(f1d, by = c("feature", "level" = "tax_level"), suffix = c("", ".f1d"))
cat(sprintf("\nVALIDATION vs Figure 1D: matched %d rows | max|d coef| = %.2e | max|d p| = %.2e | max|d q| = %.2e\n",
            nrow(chk), max(abs(chk$coef - chk$coef.f1d)), max(abs(chk$pval - chk$pval.f1d)),
            max(abs(chk$qval - chk$qval.f1d))))


key_taxa <- c("Lachnospiraceae", "Lachnospiraceae ND3007 group", "Odoribacter", "[Eubacterium] eligens group", "Anaerobutyricum")
cat("\n=== MaAsLin2 (age + sex), Time 1 vs controls: full cohort vs excluded-9 ===\n")
ma %>% filter(taxon %in% key_taxa) %>%
  transmute(level, taxon, set, beta = round(coef, 3), SE = round(stderr, 3), P = signif(pval, 3),
            q = signif(qval, 3), n_tests = n_tests_fdr, N) %>%
  arrange(level, taxon, set) %>% print(row.names = FALSE)

cat("\nExcluded-9: all q < 0.05 results:\n")
ma %>% filter(grepl("Excl9", set), qval < 0.05) %>% select(level, taxon, coef, pval, qval) %>% print(row.names = FALSE)


S8 <- ma %>% filter(grepl("Excl9", set)) %>%
  arrange(factor(level, levels = c("Family", "Genus")), pval) %>%
  transmute(Level = level, Taxon = taxon, b = sprintf("%.3f", coef), SE = sprintf("%.3f", stderr),
            P = fmtp(pval), q = fmtp(qval), Significance = make_sig(pval, qval))
names(S8)[3] <- intToUtf8(946)
stopifnot(sum(S8$Level == "Genus") == 71)
write_excel_csv(S8, file.path(tables_dir, "TableS8_excluded9.csv"), na = "")

#############################################################
# B. Sample map (Table S9)
#############################################################
maaslin_dir <- file.path(out_dir, "maaslin2", "run_adjusted"); dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)
RUN_LAB <- c(HN00165901 = "Run 1", HN00223676 = "Run 2", last = "Run 3")
RUN_ID  <- c(HN00165901 = "HN00165901", HN00223676 = "HN00223676", last = "HN00256115")
stopifnot(all(metadata_clean$sequencing_time %in% names(RUN_LAB)))

meta <- metadata_clean %>%
  mutate(sample_id = as.character(sample_id),
         run = unname(RUN_LAB[sequencing_time]), run_id = unname(RUN_ID[sequencing_time]),
         visit = case_when(analysis_group == "Asthma Time1" ~ "Time 1",
                           analysis_group == "Asthma Time2" ~ "Time 2",
                           TRUE ~ "Single visit"),
         grp = ifelse(analysis_group == "Control", "Control", "Asthma"))

paired_subjects <- meta %>% filter(grp == "Asthma") %>% group_by(subject_name_cohort) %>%
  filter(n_distinct(analysis_group) == 2) %>% pull(subject_name_cohort) %>% unique()
stopifnot(length(paired_subjects) == 44)

## de-identified subject codes (P01-P44 patients, C01-C41 controls), assigned in sorted order
subj_codes <- bind_rows(
  tibble(subject_name_cohort = sort(unique(meta$subject_name_cohort[meta$grp == "Asthma"])), prefix = "P"),
  tibble(subject_name_cohort = sort(unique(meta$subject_name_cohort[meta$grp == "Control"])), prefix = "C")) %>%
  group_by(prefix) %>% mutate(subject_code = sprintf("%s%02d", prefix, row_number())) %>% ungroup()
stopifnot(sum(subj_codes$prefix == "P") == 44, sum(subj_codes$prefix == "C") == 41)

S9 <- meta %>% left_join(subj_codes %>% select(subject_name_cohort, subject_code), by = "subject_name_cohort") %>%
  arrange(grp, subject_code, visit) %>%
  transmute(`Sample ID (SRA sample name)` = sample_id, `Subject code` = subject_code,
            Group = grp, Visit = visit, `Sequencing run` = run, `Run ID` = run_id,
            `Previously reported sample` = ifelse(run == "Run 1", "Yes", "No"))
stopifnot(nrow(S9) == 129, !anyNA(S9$`Subject code`),
          all(table(S9$`Subject code`[S9$Group == "Asthma"]) == 2))
write_excel_csv(S9, file.path(tables_dir, "TableS9_sample_map.csv"), na = "")

cat("\n=== Samples per run x group/visit ===\n")
print(table(S9$`Sequencing run`, paste(S9$Group, S9$Visit)))

## ---------------------------------------------------------------------------
## 2. PAIRS ACROSS RUNS
## ---------------------------------------------------------------------------
pairs <- S9 %>% filter(Group == "Asthma") %>%
  select(`Subject code`, Visit, `Sequencing run`) %>%
  pivot_wider(names_from = Visit, values_from = `Sequencing run`)
cat("\n=== 44 pairs: run at Time 1 (rows) x run at Time 2 (columns) ===\n")
print(table(T1 = pairs$`Time 1`, T2 = pairs$`Time 2`))
n_diff <- sum(pairs$`Time 1` != pairs$`Time 2`)
cat(sprintf("\nPairs with T1 and T2 in DIFFERENT runs: %d of 44 | same run: %d\n", n_diff, 44 - n_diff))

#############################################################
# C. Run-adjusted differential abundance (Table S10)
#############################################################
run_ma <- function(m, clr, covs, key) {
  m <- m %>% mutate(group = factor(grp, levels = c("Control", "Asthma")),
                    run = factor(run), age_num = suppressWarnings(as.numeric(age)),
                    sex_num = as.numeric(factor(sex)) - 1) %>% as.data.frame()
  rownames(m) <- m$sample_id
  fit <- Maaslin2(input_data = clr[m$sample_id, , drop = FALSE], input_metadata = m,
                  output = file.path(maaslin_dir, key), fixed_effects = c("group", covs),
                  normalization = "NONE", transform = "NONE", analysis_method = "LM",
                  min_prevalence = 0, min_abundance = 0, standardize = FALSE, max_significance = 0.25,
                  reference = c("group,Control", if ("run" %in% covs) "run,Run 2"),
                  plot_heatmap = FALSE, plot_scatter = FALSE)
  fit$results %>% filter(metadata == "group") %>%
    mutate(qval = p.adjust(pval, "BH"), n_tests = n(), N_samples = nrow(m), analysis = key)
}

fam_map <- setNames(colnames(family_clr_01), make.names(colnames(family_clr_01)))
gen_map <- setNames(colnames(genus_clr_01),  make.names(colnames(genus_clr_01)))
res <- list()
for (cmp in c("Asthma Time1", "Asthma Time2")) {
  full <- meta %>% filter(analysis_group %in% c("Control", cmp))
  r23  <- full %>% filter(run != "Run 1")
  stopifnot(length(unique(r23$run)) == 2,
            all(table(r23$run, r23$grp) > 0))          # both groups present in both runs
  tag <- ifelse(cmp == "Asthma Time1", "T1", "T2")
  for (lv in c("Family", "Genus")) {
    clr <- if (lv == "Family") family_clr_01 else genus_clr_01
    mp  <- if (lv == "Family") fam_map else gen_map
    for (a in list(list(d = full, c = c("age_num", "sex_num"),        k = "All runs"),
                   list(d = r23,  c = c("age_num", "sex_num"),        k = "Runs 2-3"),
                   list(d = r23,  c = c("age_num", "sex_num", "run"), k = "Runs 2-3, run-adjusted"))) {
      key <- paste(lv, tag, a$k)
      res[[key]] <- run_ma(a$d, clr, a$c, gsub("[^A-Za-z0-9]+", "_", key)) %>%
        mutate(level = lv, comparison = tag, model = a$k, taxon = unname(mp[feature]))
    }
  }
}
da <- bind_rows(res); stopifnot(!anyNA(da$taxon))
write_csv(da, file.path(tables_dir, "run_adjusted_DA.csv"))

## validation: "All runs" must reproduce Figure 1D
f1d <- f1d_all %>%
  filter(model_name == "Model1_age_sex", threshold == "01pct", tax_level %in% c("Family", "Genus")) %>%
  mutate(comparison = ifelse(grepl("Time1", comparison), "T1", "T2"))
v <- da %>% filter(model == "All runs") %>%
  inner_join(f1d, by = c("feature", "level" = "tax_level", "comparison"), suffix = c("", ".f1d"))
cat(sprintf("\nVALIDATION vs Figure 1D: matched %d | max|d coef| = %.2e | max|d p| = %.2e\n",
            nrow(v), max(abs(v$coef - v$coef.f1d)), max(abs(v$pval - v$pval.f1d))))

## Figure 1D taxa: main vs Runs 2-3 run-adjusted
fig1d_taxa <- da %>% filter(model == "All runs", pval < 0.05) %>% distinct(level, taxon)
wide <- da %>% semi_join(fig1d_taxa, by = c("level", "taxon")) %>%
  mutate(v = sprintf("%.2f (P=%s)", coef, signif(pval, 2))) %>%
  select(level, taxon, comparison, model, v) %>%
  pivot_wider(names_from = model, values_from = v) %>% arrange(level, taxon, comparison)
cat("\n=== Figure 1D taxa: beta (P) — all runs vs Runs 2-3 vs Runs 2-3 run-adjusted ===\n")
print(as.data.frame(wide), row.names = FALSE)

cat("\n=== Direction kept / P < 0.05 kept in run-adjusted Runs 2-3 (Figure 1D cells with P < 0.05) ===\n")
cmp_tab <- da %>% filter(model %in% c("All runs", "Runs 2-3, run-adjusted")) %>%
  select(level, taxon, comparison, model, coef, pval) %>%
  pivot_wider(names_from = model, values_from = c(coef, pval)) %>%
  filter(`pval_All runs` < 0.05)
cat(sprintf("Cells: %d | same direction: %d | still P < 0.05: %d\n", nrow(cmp_tab),
            sum(sign(cmp_tab$`coef_All runs`) == sign(cmp_tab$`coef_Runs 2-3, run-adjusted`)),
            sum(cmp_tab$`pval_Runs 2-3, run-adjusted` < 0.05)))

cat("\n=== q < 0.05 in any analysis ===\n")
da %>% filter(qval < 0.05) %>% select(level, comparison, model, taxon, coef, pval, qval, N_samples) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== Lachnospiraceae family, all models ===\n")
da %>% filter(taxon == "Lachnospiraceae") %>% select(comparison, model, coef, pval, qval, N_samples) %>%
  as.data.frame() %>% print(row.names = FALSE)

S10 <- da %>% filter(model %in% c("All runs", "Runs 2-3, run-adjusted")) %>%
  mutate(model = ifelse(model == "All runs", "All runs", "Runs 2-3 run-adj"),
         b = sprintf("%.3f", coef), P = fmtp(pval), q = fmtp(qval)) %>%
  select(Level = level, Taxon = taxon, Comparison = comparison, model, b, P, q) %>%
  pivot_wider(names_from = model, values_from = c(b, P, q), names_glue = "{model}: {.value}") %>%
  select(Level, Taxon, Comparison, starts_with("All runs"), starts_with("Runs 2-3")) %>%
  arrange(Level, Comparison, Taxon)
names(S10) <- sub(": b$", paste0(": ", intToUtf8(946)), names(S10))
write_excel_csv(S10, file.path(tables_dir, "TableS10_run_adjusted_DA.csv"), na = "")

#############################################################
# Reproducibility check against the final Supplementary Tables S8-S10
#############################################################
ms_tab <- file.path(repo_dir, "..", "Manuscript", "Supplementary_Tables")
same_table <- function(new, file, label) {
  f <- file.path(ms_tab, file)
  if (!file.exists(f)) { cat(sprintf("  %-4s final table not found\n", label)); return(invisible()) }
  old <- read.csv(f, check.names = FALSE, colClasses = "character", encoding = "UTF-8")
  new <- as.data.frame(lapply(new, as.character), check.names = FALSE)
  key <- intersect(names(new), names(old))
  m <- inner_join(new, old, by = key)
  cat(sprintf("  %-4s %d/%d rows identical in all %d columns\n", label, nrow(m), nrow(new), length(key)))
}
cat("\nComparison with final Supplementary Tables:\n")
same_table(S8,  "SupplementaryTable_S8_Excl9_GroupDifferences.csv", "S8")
same_table(S9,  "SupplementaryTable_S9_SampleMap.csv",              "S9")
same_table(S10, "SupplementaryTable_S10_RunAdjustedDA.csv",         "S10")