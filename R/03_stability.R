#############################################################
# 03_stability.R
# Within-patient stability of the gut microbiome
#   Figure 3A  ICC per genus (one-way, agreement, single; CLR values)
#   Figure 3B  mean CLR change per genus, Time 2 - Time 1 (one-sample t-test)
#   Figure 2E  community types: PAM on Aitchison distance (Euclidean on CLR),
#              k chosen by maximum mean silhouette (k = 2-10), all 129 samples
#   Figure 2F  community-type transitions, Time 1 -> Time 2
#############################################################

source("R/00_setup.R")
library(irr)
library(cluster)

metadata_clean <- readRDS(f_metadata) %>% mutate(sample_id = as.character(sample_id))
clr   <- list("0.1%/10%" = readRDS(file.path(obj_dir, "genus_clr_01pct.rds")),
              "0.5%/10%" = readRDS(file.path(obj_dir, "genus_clr_05pct.rds")))
count <- list("0.1%/10%" = readRDS(file.path(obj_dir, "genus_count_01pct.rds")),
              "0.5%/10%" = readRDS(file.path(obj_dir, "genus_count_05pct.rds")))

meta_asthma <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  group_by(subject_name_cohort) %>% filter(n_distinct(analysis_group) == 2) %>% ungroup()
meta_t1 <- meta_asthma %>% filter(analysis_group == "Asthma Time1") %>% arrange(subject_name_cohort)
meta_t2 <- meta_asthma %>% filter(analysis_group == "Asthma Time2") %>% arrange(subject_name_cohort)
stopifnot(identical(meta_t1$subject_name_cohort, meta_t2$subject_name_cohort), nrow(meta_t1) == 44)

## Figure 3A: ICC per genus
compute_icc <- function(clr_matrix, count_matrix, threshold) {
  t1 <- clr_matrix[meta_t1$sample_id, , drop = FALSE]
  t2 <- clr_matrix[meta_t2$sample_id, , drop = FALSE]
  cp <- count_matrix[c(meta_t1$sample_id, meta_t2$sample_id), , drop = FALSE]
  mean_ra <- colMeans(sweep(cp, 1, rowSums(cp), "/") * 100)
  bind_rows(lapply(colnames(clr_matrix), function(g) {
    fit <- irr::icc(cbind(t1[, g], t2[, g]), model = "oneway", type = "agreement", unit = "single")
    data.frame(Genus = g, Threshold = threshold,
               ICC = round(fit$value, 4),
               CI_lower = round(max(fit$lbound, 0), 4),
               CI_upper = round(min(fit$ubound, 1), 4),
               Mean_RA_pct = round(mean_ra[g], 3))
  }))
}
icc_all <- bind_rows(lapply(names(clr), function(th) compute_icc(clr[[th]], count[[th]], th)))
write_csv(icc_all, file.path(tables_dir, "ICC_per_genus.csv"))

icc_01 <- icc_all %>% filter(Threshold == "0.1%/10%")
cat(sprintf("\nFigure 3A (0.1%%/10%%): n = %d | median ICC = %.3f | ICC >= 0.50: %.1f%%\n",
            nrow(icc_01), median(icc_01$ICC), 100 * mean(icc_01$ICC >= 0.50)))
print(icc_01 %>% arrange(desc(ICC)) %>% head(5) %>% select(Genus, ICC, CI_lower, CI_upper))

## Figure 3B: CLR change per genus (T2 - T1)
compute_delta <- function(clr_matrix, threshold) {
  delta <- clr_matrix[meta_t2$sample_id, , drop = FALSE] - clr_matrix[meta_t1$sample_id, , drop = FALSE]
  bind_rows(lapply(colnames(delta), function(g) {
    tt <- t.test(delta[, g])
    data.frame(Genus = g, Threshold = threshold, Mean_change = round(mean(delta[, g]), 4),
               CI_lower = round(tt$conf.int[1], 4), CI_upper = round(tt$conf.int[2], 4),
               P = round(tt$p.value, 4))
  }))
}
delta_all <- bind_rows(lapply(names(clr), function(th) compute_delta(clr[[th]], th)))
write_csv(delta_all, file.path(tables_dir, "CLR_change_per_genus.csv"))
cat("\nFigure 3B (0.1%/10%): genera with P < 0.05\n")
print(delta_all %>% filter(Threshold == "0.1%/10%", P < 0.05) %>% arrange(Mean_change), row.names = FALSE)

## Figure 2E-F: community types
run_clustering <- function(clr_matrix, threshold, k_range = 2:10) {
  d <- stats::dist(clr_matrix, method = "euclidean")
  set.seed(SEED)
  sil <- sapply(k_range, function(k) mean(silhouette(pam(d, k = k, diss = TRUE))[, 3]))
  k_opt <- k_range[which.max(round(sil, 4))]
  set.seed(SEED)
  fit <- pam(d, k = k_opt, diss = TRUE)
  list(silhouette = data.frame(Threshold = threshold, k = k_range, silhouette = round(sil, 4)),
       assignment = data.frame(sample_id = rownames(clr_matrix), Threshold = threshold,
                               cluster = paste0("CT", fit$clustering),
                               silhouette_width = round(silhouette(fit)[, 3], 4)) %>%
         left_join(metadata_clean %>% select(sample_id, analysis_group, subject_name_cohort),
                   by = "sample_id"))
}
ct <- lapply(names(clr), function(th) run_clustering(clr[[th]], th)); names(ct) <- names(clr)
assign_all <- bind_rows(lapply(ct, `[[`, "assignment"))
write_csv(assign_all, file.path(tables_dir, "community_type_assignments.csv"))
write_csv(bind_rows(lapply(ct, `[[`, "silhouette")), file.path(tables_dir, "community_type_silhouette.csv"))

for (th in names(ct)) {
  a  <- ct[[th]]$assignment
  ra <- sweep(count[[th]], 1, rowSums(count[[th]]), "/")
  for (cl in sort(unique(a$cluster))) {
    m <- colMeans(ra[a$sample_id[a$cluster == cl], , drop = FALSE])
    cat(sprintf("  %s %s: top genus %s\n", th, cl, names(which.max(m))))
  }
}

ct_summary <- assign_all %>% count(Threshold, analysis_group, cluster) %>%
  group_by(Threshold, analysis_group) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
write_csv(ct_summary, file.path(tables_dir, "community_type_summary.csv"))

trans <- assign_all %>% filter(subject_name_cohort %in% meta_t1$subject_name_cohort,
                               analysis_group != "Control") %>%
  select(Threshold, subject_name_cohort, analysis_group, cluster) %>%
  pivot_wider(names_from = analysis_group, values_from = cluster)
cat("\nFigure 2E (0.1%/10%): % per group\n")
print(ct_summary %>% filter(Threshold == "0.1%/10%"), row.names = FALSE)
cat("\nFigure 2F: patients keeping the same community type\n")
print(trans %>% group_by(Threshold) %>%
        summarise(same = sum(`Asthma Time1` == `Asthma Time2`), n = n(), .groups = "drop"))

## Reproducibility check against the original result tables
orig <- file.path(repo_dir, "..", "02_result_tables")
cat("\nComparison with original tables:\n")
f_icc <- file.path(orig, "Figure2C_ICC_both_thresholds_20260922.csv")
if (file.exists(f_icc)) {
  o <- read.csv(f_icc)
  m <- inner_join(icc_all, o, by = c("Genus", "Threshold" = "threshold"))
  cat(sprintf("  ICC: %d/%d genera matched, max |difference| = %.6f\n",
              nrow(m), nrow(icc_all), max(abs(m$ICC.x - m$ICC.y))))
}
f_ct <- file.path(orig, "Figure3B_ClusterAssignments_20260923.csv")
if (file.exists(f_ct)) {
  o <- read.csv(f_ct)
  m <- inner_join(assign_all, o, by = c("sample_id", "Threshold" = "threshold"))
  cat(sprintf("  Community types: %d/%d samples with identical cluster\n",
              sum(m$cluster.x == m$cluster.y), nrow(m)))
}
f_d <- file.path(orig, "SupplementaryFig3_PanelC_DeltaCLR_20261019.csv")
if (file.exists(f_d)) {
  o <- read.csv(f_d)
  m <- inner_join(delta_all, o, by = c("Genus", "Threshold" = "threshold"))
  cat(sprintf("  CLR change: %d/%d genera matched, max |difference| in P = %.6f\n",
              nrow(m), nrow(delta_all), max(abs(m$P - m$p_value))))
}

cat("\nValues in the manuscript:\n",
    " 3A  n = 71, median ICC = 0.514, ICC >= 0.50: 53.5%\n",
    " 3B  P < 0.05: Clostridium, Veillonella, Phascolarctobacterium (down); Unassigned, Blautia (up)\n",
    " 2E  CT1: Control 36.6%, Time 1 47.7%, Time 2 43.2%\n",
    " 2F  40 of 44 kept their community type\n")
