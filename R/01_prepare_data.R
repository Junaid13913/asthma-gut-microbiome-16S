#############################################################
# 01_prepare_data.R
# Load data, filter taxa at two thresholds, CLR-transform.
#   Filtering: mean relative abundance >= 0.1% (primary) or 0.5%
#              (sensitivity) AND prevalence >= 10%, per level
#   CLR      : pseudo-count 0.5, after filtering, per level
# Outputs: results/data_objects/<level>_count_<thr>.rds, <level>_clr_<thr>.rds
#############################################################

source("R/00_setup.R")
library(compositions)

## ---- 1. Load ----
metadata_clean <- readRDS(f_metadata)
seqtab_final   <- readRDS(f_asv)
crf_raw        <- readRDS(f_crf)
count_all      <- lapply(f_unfiltered, readRDS)

cat("\nSamples:", nrow(metadata_clean), "| ASVs:", ncol(seqtab_final), "\n")
print(table(metadata_clean$analysis_group))
stopifnot(nrow(metadata_clean) == 129, nrow(seqtab_final) == 129)

paired_subjects <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  group_by(subject_name_cohort) %>%
  filter(n_distinct(analysis_group) == 2) %>%
  pull(subject_name_cohort) %>% unique()
stopifnot(length(paired_subjects) == 44)

## ---- 2. Filter ----
filter_taxa <- function(count_matrix, min_ra, min_prev) {
  ra   <- sweep(count_matrix, 1, rowSums(count_matrix), "/")
  keep <- colMeans(ra) >= min_ra & colMeans(ra > 0) >= min_prev
  count_matrix[, keep, drop = FALSE]
}

## ---- 3. CLR ----
apply_clr <- function(count_matrix) as.matrix(compositions::clr(count_matrix + PSEUDOCOUNT))

cat("\nTaxa retained:\n")
for (lev in names(count_all)) {
  for (thr in names(THRESHOLDS)) {
    cnt <- filter_taxa(count_all[[lev]], THRESHOLDS[[thr]], MIN_PREV)
    clr <- apply_clr(cnt)
    stopifnot(max(abs(rowMeans(clr))) < 1e-8)
    nm <- tolower(lev)
    saveRDS(cnt, file.path(obj_dir, sprintf("%s_count_%s.rds", nm, thr)))
    saveRDS(clr, file.path(obj_dir, sprintf("%s_clr_%s.rds",   nm, thr)))
    cat(sprintf("  %-7s %s : %d taxa\n", lev, thr, ncol(cnt)))
  }
}
## expected: Genus 71 / 33, Family 34 / 22, Phylum 7 / 4

## ---- 4. Reproducibility check against the original objects ----
cat("\nComparison with original objects in 01_data_objects:\n")
for (lev in tolower(names(count_all))) {
  for (thr in names(THRESHOLDS)) {
    for (type in c("count", "clr")) {
      new <- readRDS(file.path(obj_dir, sprintf("%s_%s_%s.rds", lev, type, thr)))
      check_same(new, file.path(input_dir, sprintf("%s_%s_%s_20260922.rds", lev, type, thr)),
                 sprintf("%s %s %s", lev, type, thr))
    }
  }
}
