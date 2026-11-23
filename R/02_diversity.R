#############################################################
# 02_diversity.R
# Alpha and beta diversity (ASV table after DADA2 QC, no taxon filtering)
#   Figure 1A  alpha diversity, controls vs asthma at each visit
#   Figure 1B  Bray-Curtis PCoA + PERMANOVA; distance to controls
#   Figure 2A  paired alpha diversity, Time 1 vs Time 2
#   Figure 2B  paired PCoA + PERMANOVA blocked by subject
#   Figure 2C  within-patient vs between-patient dissimilarity
#   Figure 2D  change in distance to controls, Time 1 -> Time 2
# Output: results/tables/diversity_statistics.csv (+ sample-level tables)
#############################################################

source("R/00_setup.R")
library(permute)

metadata_clean <- readRDS(f_metadata) %>% mutate(sample_id = as.character(sample_id))
seqtab_final   <- readRDS(f_asv)

## ---- Paired asthma samples (44 subjects, same order at both visits) ----
meta_asthma <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  group_by(subject_name_cohort) %>% filter(n_distinct(analysis_group) == 2) %>% ungroup()
meta_t1 <- meta_asthma %>% filter(analysis_group == "Asthma Time1") %>% arrange(subject_name_cohort)
meta_t2 <- meta_asthma %>% filter(analysis_group == "Asthma Time2") %>% arrange(subject_name_cohort)
stopifnot(identical(meta_t1$subject_name_cohort, meta_t2$subject_name_cohort), nrow(meta_t1) == 44)
control_samples <- metadata_clean$sample_id[metadata_clean$analysis_group == "Control"]
stopifnot(length(control_samples) == 41)

## ---- Shared objects ----
alpha_div <- data.frame(
  sample_id     = rownames(seqtab_final),
  observed_ASVs = rowSums(seqtab_final > 0),
  shannon       = vegan::diversity(seqtab_final, index = "shannon")
) %>% left_join(metadata_clean %>% select(sample_id, analysis_group, subject_name_cohort),
                by = "sample_id")

to_ra <- function(m) { r <- sweep(m, 1, rowSums(m), "/"); r[!is.finite(r)] <- 0; r }
bc_matrix <- as.matrix(vegdist(to_ra(seqtab_final), method = "bray"))

stats <- list()
add_stat <- function(figure, test, comparison, statistic, p) {
  stats[[length(stats) + 1]] <<- data.frame(Figure = figure, Test = test, Comparison = comparison,
                                            Statistic = statistic, P = p)
}

## Figure 1A: alpha diversity, controls vs asthma (Mann-Whitney)
for (metric in c("observed_ASVs", "shannon")) {
  for (grp in c("Asthma Time1", "Asthma Time2")) {
    d <- alpha_div %>% filter(analysis_group %in% c("Control", grp))
    t <- wilcox.test(as.formula(paste(metric, "~ analysis_group")), data = d)
    add_stat("1A", "Wilcoxon rank-sum", paste(metric, ": Control vs", grp), NA, t$p.value)
  }
}

## Figure 1B: PERMANOVA and distance to controls
for (grp in c("Asthma Time1", "Asthma Time2")) {
  m  <- metadata_clean %>% filter(analysis_group %in% c("Control", grp))
  dd <- as.dist(bc_matrix[m$sample_id, m$sample_id])
  set.seed(SEED)
  pm <- adonis2(dd ~ analysis_group, data = m, permutations = 999)
  add_stat("1B", "PERMANOVA (R2)", paste("Control vs", grp), pm$R2[1], pm$`Pr(>F)`[1])
}
pcoa_all <- cmdscale(as.dist(bc_matrix), k = 2, eig = TRUE)

bc_to_ctrl <- data.frame(sample_id = metadata_clean$sample_id,
                         analysis_group = metadata_clean$analysis_group,
                         bc_to_ctrl = apply(bc_matrix[metadata_clean$sample_id, control_samples], 1, median))
for (grp in c("Asthma Time1", "Asthma Time2")) {
  t <- wilcox.test(bc_to_ctrl$bc_to_ctrl[bc_to_ctrl$analysis_group == grp],
                   bc_to_ctrl$bc_to_ctrl[bc_to_ctrl$analysis_group == "Control"], exact = FALSE)
  add_stat("1B", "Wilcoxon rank-sum", paste("Distance to controls:", grp, "vs Control"), NA, t$p.value)
}

## Figure 2A: paired alpha diversity (Wilcoxon signed-rank)
a1 <- alpha_div[match(meta_t1$sample_id, alpha_div$sample_id), ]
a2 <- alpha_div[match(meta_t2$sample_id, alpha_div$sample_id), ]
for (metric in c("observed_ASVs", "shannon")) {
  t <- wilcox.test(a1[[metric]], a2[[metric]], paired = TRUE, exact = FALSE)
  add_stat("2A", "Wilcoxon signed-rank", paste(metric, ": Time 1 vs Time 2"), NA, t$p.value)
}

## Figure 2B: paired PCoA, PERMANOVA with permutations within subject
paired_ids <- c(meta_t1$sample_id, meta_t2$sample_id)
bc_paired  <- vegdist(to_ra(seqtab_final[paired_ids, , drop = FALSE]), method = "bray")
meta_paired <- metadata_clean[match(paired_ids, metadata_clean$sample_id), ]
perm_scheme <- how(within = Within(type = "free"),
                   blocks = meta_paired$subject_name_cohort, nperm = 999)
set.seed(SEED)
pm <- adonis2(bc_paired ~ analysis_group, data = meta_paired, permutations = perm_scheme)
add_stat("2B", "PERMANOVA blocked by subject (R2)", "Time 1 vs Time 2", pm$R2[1], pm$`Pr(>F)`[1])
pcoa_paired <- cmdscale(bc_paired, k = 2, eig = TRUE)

## Figure 2C: within-patient (T1-T2) vs between-patient (T1-T1)
within_dist  <- sapply(seq_len(44), function(i) bc_matrix[meta_t1$sample_id[i], meta_t2$sample_id[i]])
bc_t1        <- as.matrix(vegdist(to_ra(seqtab_final[meta_t1$sample_id, , drop = FALSE]), method = "bray"))
between_dist <- sapply(seq_len(44), function(i) median(bc_t1[i, -i]))
t <- wilcox.test(within_dist, between_dist, paired = TRUE, exact = FALSE)
add_stat("2C", "Wilcoxon signed-rank",
         sprintf("Within (median %.3f) vs between (median %.3f)", median(within_dist), median(between_dist)),
         NA, t$p.value)

## Figure 2D: distance to controls, Time 1 vs Time 2
t1_dist <- bc_to_ctrl$bc_to_ctrl[match(meta_t1$sample_id, bc_to_ctrl$sample_id)]
t2_dist <- bc_to_ctrl$bc_to_ctrl[match(meta_t2$sample_id, bc_to_ctrl$sample_id)]
t <- wilcox.test(t1_dist, t2_dist, paired = TRUE, exact = FALSE)
add_stat("2D", "Wilcoxon signed-rank",
         sprintf("Distance to controls, T1 vs T2 (%d moved closer)", sum(t2_dist < t1_dist)), NA, t$p.value)

## Save
stats_df <- bind_rows(stats) %>% mutate(Statistic = round(Statistic, 4), P = signif(P, 3))
write_csv(stats_df, file.path(tables_dir, "diversity_statistics.csv"))
write_csv(alpha_div %>% select(sample_id, analysis_group, observed_ASVs, shannon),
          file.path(tables_dir, "alpha_diversity_per_sample.csv"))
write_csv(data.frame(sample_id = rownames(pcoa_all$points), PCoA1 = pcoa_all$points[, 1],
                     PCoA2 = pcoa_all$points[, 2]),
          file.path(tables_dir, "pcoa_braycurtis_all_samples.csv"))

cat("\nPCoA variance explained (all samples):",
    round(100 * pcoa_all$eig[1:2] / sum(pcoa_all$eig[pcoa_all$eig > 0]), 1), "%\n")
cat("PCoA variance explained (paired):     ",
    round(100 * pcoa_paired$eig[1:2] / sum(pcoa_paired$eig[pcoa_paired$eig > 0]), 1), "%\n\n")
print(stats_df, row.names = FALSE)

cat("\nValues in the manuscript figures:\n",
    " 1A  observed P = 0.496 / 0.894; Shannon P = 0.454 / 0.892\n",
    " 1B  PCoA 9.7% / 5.5%; distance to controls P = 0.27 / 0.094\n",
    " 2A  observed P = 0.246; Shannon P = 0.156\n",
    " 2B  PCoA 10.2% / 6.8%; R2 = 0.006, P = 0.255\n",
    " 2C  medians 0.591 vs 0.835, P < 0.001\n",
    " 2D  P = 0.50\n")
