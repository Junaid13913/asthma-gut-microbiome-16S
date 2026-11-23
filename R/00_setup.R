#############################################################
# 00_setup.R
# Paths, packages and shared settings for all analysis scripts.
# Run every script from the github_release folder.
# Input data (not in the repository): ../01_data_objects/
# All outputs are written to results/ inside this folder.
#############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(vegan)
})

## ---- Paths ----
repo_dir   <- normalizePath(getwd())
input_dir  <- normalizePath(file.path(repo_dir, "..", "01_data_objects"), mustWork = TRUE)
out_dir    <- file.path(repo_dir, "results")
obj_dir    <- file.path(out_dir, "data_objects")
tables_dir <- file.path(out_dir, "tables")
fig_dir    <- file.path(out_dir, "figures")
for (d in c(obj_dir, tables_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- Input files ----
f_metadata <- file.path(input_dir, "metadata_verified_UPDATED_20260617.rds")
f_asv      <- file.path(input_dir, "ASV_table_verified_UPDATED_20260617.rds")
f_taxonomy <- file.path(input_dir, "taxonomy_verified_UPDATED_20260617.rds")
f_crf      <- file.path(input_dir, "CRF_clinical_data_clean_UPDATED_20260617.rds")
f_unfiltered <- c(
  Genus  = file.path(input_dir, "genus_count_all_UNFILTERED_20260922.rds"),
  Family = file.path(input_dir, "family_count_all_UNFILTERED_20260922.rds"),
  Phylum = file.path(input_dir, "phylum_count_all_UNFILTERED_20260922.rds")
)

## ---- Shared settings ----
THRESHOLDS   <- c("01pct" = 0.001, "05pct" = 0.005)
MIN_PREV     <- 0.10
PSEUDOCOUNT  <- 0.5
SEED         <- 42
GROUP_LEVELS <- c("Control", "Asthma Time1", "Asthma Time2")

## ---- Helpers ----
make_sig <- function(p, q = NA) {
  stars <- dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                            p < 0.05 ~ "*", TRUE ~ "")
  paste0(stars, ifelse(!is.na(q) & q < 0.05, "\u2020", ""))
}

check_same <- function(new, original_file, label, tol = 1e-8) {
  if (!file.exists(original_file)) {
    cat(sprintf("  %-28s original not found (skipped)\n", label)); return(invisible(NA))
  }
  old <- readRDS(original_file)
  same <- identical(dim(new), dim(old)) &&
    identical(rownames(new), rownames(old)) &&
    identical(colnames(new), colnames(old)) &&
    isTRUE(all.equal(unname(as.matrix(new)), unname(as.matrix(old)), tolerance = tol))
  cat(sprintf("  %-28s %s\n", label, ifelse(same, "IDENTICAL", "*** DIFFERENT ***")))
  invisible(same)
}

cat("Setup done.\n  input :", input_dir, "\n  output:", out_dir, "\n")
