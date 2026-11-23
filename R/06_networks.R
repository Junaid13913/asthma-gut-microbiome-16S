#############################################################
# 06_networks.R
# Co-abundance and co-change networks in asthma (44 paired patients)
#
#   Co-abundance (Time 1, Time 2): SPIEC-EASI (Meinshausen-Buhlmann, StARS)
#     input = filtered genus COUNTS; SPIEC-EASI applies its own CLR (clr(x + 1)),
#     so the data are CLR-transformed once
#   Co-change (Time 2 - Time 1): GeneNet partial correlations on delta-CLR,
#     edges with tail-area FDR q < 0.2
#   Edges with |partial correlation| < 1e-6 are listed but not used
#   Modules: Louvain communities (>= 2 genera) in the positive-edge subgraph
#   Module score = mean CLR of member genera; lm(score ~ clinical + age + sex);
#     BH over all module tests within a threshold
#
#   Figure 5, Figure E9, Tables S3 (edges), S4 (modules), S5 (module-clinical)
#############################################################

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(SpiecEasi); library(igraph); library(ggraph); library(tidygraph); library(ggforce)
  library(corpcor); library(GeneNet); library(ggalluvial)
})

## ---- Settings ----
EDGE_TOL        <- 1e-6    # |pcor| below this = numerically zero = no edge
COCHANGE_Q      <- 0.2     # GeneNet edge threshold
MIN_MODULE_SIZE <- 2
COVARIATES      <- c("age_num", "sex_num")
COL_POS <- "#A68A56"; COL_NEG <- "#5C3D6E"; COL_MID <- "#F7F5F0"
fig5_dir <- fig_dir; supp_dir <- fig_dir; out_tab <- tables_dir

## ---- Data ----
metadata_clean  <- readRDS(f_metadata)
crf_delta_clean <- readRDS(f_crf)
genus_clr_01    <- readRDS(file.path(obj_dir, "genus_clr_01pct.rds"))
genus_clr_05    <- readRDS(file.path(obj_dir, "genus_clr_05pct.rds"))
genus_count_01  <- readRDS(file.path(obj_dir, "genus_count_01pct.rds"))
genus_count_05  <- readRDS(file.path(obj_dir, "genus_count_05pct.rds"))
stopifnot(identical(colnames(genus_count_01), colnames(genus_clr_01)),
          identical(colnames(genus_count_05), colnames(genus_clr_05)))
paired_subjects <- metadata_clean %>%
  filter(analysis_group %in% c("Asthma Time1", "Asthma Time2")) %>%
  group_by(subject_name_cohort) %>% filter(n_distinct(analysis_group) == 2) %>%
  pull(subject_name_cohort) %>% unique()
meta_t1 <- metadata_clean %>%
  filter(analysis_group == "Asthma Time1", subject_name_cohort %in% paired_subjects) %>%
  mutate(sample_id = as.character(sample_id), age_num = suppressWarnings(as.numeric(age)),
         sex_num = as.numeric(factor(sex)) - 1) %>% arrange(subject_name_cohort)
meta_t2 <- metadata_clean %>%
  filter(analysis_group == "Asthma Time2", subject_name_cohort %in% paired_subjects) %>%
  mutate(sample_id = as.character(sample_id)) %>% arrange(subject_name_cohort)
stopifnot(length(paired_subjects) == 44)

stopifnot(exists("meta_t1"), exists("meta_t2"), exists("genus_clr_01"),
          exists("genus_clr_05"), exists("crf_delta_clean"),
          identical(meta_t1$subject_name_cohort, meta_t2$subject_name_cohort))

## Build all clinical columns from the raw CRF values (works whether crf_delta_clean
## is the raw CRF file or the Figure 4 version)
need_raw <- c("ACT_t1", "ACT_t2", "FEV1_pct_t1", "FEV1_pct_t2", "FEV1_FVC_t1", "FEV1_FVC_t2",
              "FeNO_t1", "FeNO_t2", "Blood_eos_t1", "Blood_eos_t2")
miss <- setdiff(need_raw, names(crf_delta_clean))
if (length(miss)) stop("crf_delta_clean is missing raw columns: ", paste(miss, collapse = ", "))
crf_delta_clean <- crf_delta_clean %>%
  mutate(
    neg_log_FeNO_t1     = -log(FeNO_t1 + 1),      neg_log_FeNO_t2     = -log(FeNO_t2 + 1),
    neg_log_BloodEos_t1 = -log(Blood_eos_t1 + 1), neg_log_BloodEos_t2 = -log(Blood_eos_t2 + 1),
    delta_ACT              = ACT_t2 - ACT_t1,
    delta_FEV1_pct         = FEV1_pct_t2 - FEV1_pct_t1,
    delta_FEV1_FVC         = FEV1_FVC_t2 - FEV1_FVC_t1,
    delta_neg_log_FeNO     = neg_log_FeNO_t2 - neg_log_FeNO_t1,
    delta_neg_log_BloodEos = neg_log_BloodEos_t2 - neg_log_BloodEos_t1)

CLIN_VARS   <- c("ACT", "FEV1_pct", "FEV1_FVC", "neg_log_FeNO", "neg_log_BloodEos")
CLIN_LABELS <- c(ACT = "ACT", FEV1_pct = "FEV1%", FEV1_FVC = "FEV1/FVC",
                 neg_log_FeNO = "−log(FeNO)", neg_log_BloodEos = "−log(blood Eos)")


#############################################################
# 1. Networks, modules and module-clinical regressions
#############################################################
## CLR matrices (rows = subject, identical order at T1, T2, Δ)
sub_ids <- meta_t1$subject_name_cohort
mk <- function(clr, ids) { m <- clr[ids, , drop = FALSE]; rownames(m) <- sub_ids; m }
clr_in <- list(
  "0.1%/10%" = list(T1 = mk(genus_clr_01, meta_t1$sample_id), T2 = mk(genus_clr_01, meta_t2$sample_id)),
  "0.5%/10%" = list(T1 = mk(genus_clr_05, meta_t1$sample_id), T2 = mk(genus_clr_05, meta_t2$sample_id))
)
for (th in names(clr_in)) clr_in[[th]]$Delta <- clr_in[[th]]$T2 - clr_in[[th]]$T1

## filtered COUNT matrices for SPIEC-EASI (same subjects, same genera, same order)
count_in <- list(
  "0.1%/10%" = list(T1 = mk(genus_count_01, meta_t1$sample_id), T2 = mk(genus_count_01, meta_t2$sample_id)),
  "0.5%/10%" = list(T1 = mk(genus_count_05, meta_t1$sample_id), T2 = mk(genus_count_05, meta_t2$sample_id))
)
for (th in names(count_in)) for (tp in c("T1", "T2"))
  stopifnot(identical(dimnames(count_in[[th]][[tp]]), dimnames(clr_in[[th]][[tp]])),
            all(count_in[[th]][[tp]] >= 0))

## ── Graph + module helpers ───────────────────────────────────────────────────
make_graph <- function(edges_all, genera) {
  e <- edges_all %>% filter(status == "retained")
  g <- make_empty_graph(n = length(genera), directed = FALSE)
  V(g)$name <- genera
  if (nrow(e) > 0) {
    g <- add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))),
                   attr = list(pcor = e$pcor,
                               sign = ifelse(e$pcor > 0, "positive", "negative"),
                               weight = abs(e$pcor)))
  }
  V(g)$degree <- as.numeric(igraph::degree(g))
  g
}

assign_modules <- function(g, prefix) {
  V(g)$module_label <- "None"
  if (ecount(g) == 0) return(g)
  pos <- which(E(g)$sign == "positive")
  if (length(pos) == 0) return(g)
  gp <- subgraph_from_edges(g, pos, delete.vertices = FALSE)
  set.seed(42)
  mem   <- as.character(membership(cluster_louvain(gp, weights = E(gp)$weight)))
  sizes <- table(mem)
  keep  <- names(sizes)[sizes >= MIN_MODULE_SIZE]
  if (length(keep) == 0) return(g)
  first <- sapply(keep, function(k) min(V(gp)$name[mem == k]))
  ord   <- keep[order(-as.integer(sizes[keep]), first)]
  lab   <- setNames(paste0(prefix, "-", LETTERS[seq_along(ord)]), ord)
  V(g)$module_label <- ifelse(mem %in% keep, lab[mem], "None")
  g
}

summarise_net <- function(id, threshold, timepoint, g, edges_all) {
  list(id = id, threshold = threshold, timepoint = timepoint, graph = g,
       edges_all = edges_all,
       n_edges = ecount(g),
       n_pos = if (ecount(g)) sum(E(g)$sign == "positive") else 0L,
       n_neg = if (ecount(g)) sum(E(g)$sign == "negative") else 0L,
       modules = sort(setdiff(unique(V(g)$module_label), "None")))
}

## ── SPIEC-EASI (T1, T2) ──────────────────────────────────────────────────────
## ---- SPIEC-EASI (T1, T2): input = filtered counts ----
spiec_rds <- function(tp, th) file.path(obj_dir, sprintf("SpiecEasi_%s_%s_fit.rds", tp,
                                                         ifelse(th == "0.1%/10%", "01pct", "05pct")))
get_fit <- function(tp, th, counts) {
  set.seed(42)
  fit <- spiec.easi(as.matrix(counts), method = "mb", lambda.min.ratio = 1e-2, nlambda = 20,
                    pulsar.params = list(rep.num = 50, seed = 42))
  saveRDS(list(fit = fit, input = "filtered counts (SPIEC-EASI applies clr(x + 1))"), spiec_rds(tp, th))
  fit
}

build_spiec_net <- function(tp, th) {
  clr    <- clr_in[[th]][[tp]]
  genera <- colnames(clr)
  fit    <- get_fit(tp, th, count_in[[th]][[tp]])
  adj    <- as.matrix(getRefit(fit))
  beta   <- as.matrix(symBeta(getOptBeta(fit), mode = "maxabs"))
  idx    <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  edges_all <- data.frame(from = genera[idx[, 1]], to = genera[idx[, 2]],
                          pcor = beta[idx], pval = NA_real_, qval = NA_real_,
                          stringsAsFactors = FALSE)
  edges_all$status <- ifelse(abs(edges_all$pcor) < EDGE_TOL,
                             "removed: |pcor| < 1e-6 (numerically zero)", "retained")
  g <- assign_modules(make_graph(edges_all, genera), tp)
  summarise_net(paste(tp, th), th, tp, g, edges_all)
}

## ── Co-change (ΔCLR) ─────────────────────────────────────────────────────────
build_cochange_net <- function(th) {
  delta  <- clr_in[[th]]$Delta
  genera <- colnames(delta)
  et  <- GeneNet::network.test.edges(corpcor::pcor.shrink(delta, verbose = FALSE),
                                     plot = FALSE, verbose = FALSE)
  sig <- et[et$qval < COCHANGE_Q, ]
  edges_all <- data.frame(from = genera[sig$node1], to = genera[sig$node2],
                          pcor = sig$pcor, pval = sig$pval, qval = sig$qval,
                          stringsAsFactors = FALSE)
  edges_all$status <- ifelse(abs(edges_all$pcor) < EDGE_TOL,
                             "removed: |pcor| < 1e-6 (numerically zero)", "retained")
  g <- assign_modules(make_graph(edges_all, genera), "Δ")
  summarise_net(paste("Delta", th), th, "Delta", g, edges_all)
}

nets <- list()
for (th in c("0.1%/10%", "0.5%/10%")) {
  nets[[paste("T1", th)]]    <- build_spiec_net("T1", th)
  nets[[paste("T2", th)]]    <- build_spiec_net("T2", th)
  nets[[paste("Delta", th)]] <- build_cochange_net(th)
}

## ── Module scores + regressions ──────────────────────────────────────────────
base_meta <- meta_t1 %>%
  transmute(subject_name_cohort,
            age_num = suppressWarnings(as.numeric(age)),
            sex_num = as.numeric(factor(sex)) - 1)
clin_meta <- list(
  T1 = base_meta %>% left_join(crf_delta_clean %>% select(subject_name_cohort,
                                                          ACT = ACT_t1, FEV1_pct = FEV1_pct_t1, FEV1_FVC = FEV1_FVC_t1,
                                                          neg_log_FeNO = neg_log_FeNO_t1, neg_log_BloodEos = neg_log_BloodEos_t1),
                               by = "subject_name_cohort"),
  T2 = base_meta %>% left_join(crf_delta_clean %>% select(subject_name_cohort,
                                                          ACT = ACT_t2, FEV1_pct = FEV1_pct_t2, FEV1_FVC = FEV1_FVC_t2,
                                                          neg_log_FeNO = neg_log_FeNO_t2, neg_log_BloodEos = neg_log_BloodEos_t2),
                               by = "subject_name_cohort"),
  Delta = base_meta %>% left_join(crf_delta_clean %>% select(subject_name_cohort,
                                                             ACT = delta_ACT, FEV1_pct = delta_FEV1_pct, FEV1_FVC = delta_FEV1_FVC,
                                                             neg_log_FeNO = delta_neg_log_FeNO, neg_log_BloodEos = delta_neg_log_BloodEos),
                                  by = "subject_name_cohort")
)
BLOCK_NAME <- c(T1 = "Time 1", T2 = "Time 2", Delta = "Δ (co-change)")

module_regressions <- function(net) {
  tp  <- net$timepoint
  clr <- clr_in[[net$threshold]][[tp]]
  g   <- net$graph
  if (length(net$modules) == 0) {
    return(data.frame(threshold = net$threshold, block = BLOCK_NAME[[tp]],
                      module = "no module", members = NA_character_, clinical_var = NA,
                      coef = NA_real_, se = NA_real_, pval = NA_real_, n = NA_integer_,
                      note = sprintf("No module of >=%d genera (%d retained edges, %d positive)",
                                     MIN_MODULE_SIZE, net$n_edges, net$n_pos)))
  }
  out <- list()
  for (m in net$modules) {
    mem   <- V(g)$name[V(g)$module_label == m]
    score <- rowMeans(clr[, mem, drop = FALSE])
    for (cv in CLIN_VARS) {
      df <- clin_meta[[tp]] %>%
        mutate(mod_score = score[match(subject_name_cohort, names(score))],
               clin_val  = .data[[cv]]) %>%
        filter(!is.na(mod_score), !is.na(clin_val),
               if_all(all_of(COVARIATES), ~ !is.na(.)))
      fit <- lm(as.formula(paste("mod_score ~ clin_val +", paste(COVARIATES, collapse = " + "))), data = df)
      co  <- summary(fit)$coefficients["clin_val", ]
      out[[length(out) + 1]] <- data.frame(
        threshold = net$threshold, block = BLOCK_NAME[[tp]], module = m,
        members = paste(sort(mem), collapse = "; "), clinical_var = cv,
        coef = co[["Estimate"]], se = co[["Std. Error"]], pval = co[["Pr(>|t|)"]],
        n = nrow(df), note = "")
    }
  }
  bind_rows(out)
}

## BH across all real tests within a threshold (Time 1 + Time 2 + Δ blocks);
## placeholder "no module" rows keep qval = NA
module_clinical <- bind_rows(lapply(nets, module_regressions)) %>% group_by(threshold) %>%
  mutate(qval = { q <- rep(NA_real_, n()); ok <- !is.na(pval)
  q[ok] <- p.adjust(pval[ok], method = "BH"); q }) %>% ungroup()

fig5_obj <- list(settings = list(EDGE_TOL = EDGE_TOL, COCHANGE_Q = COCHANGE_Q,
                                 MIN_MODULE_SIZE = MIN_MODULE_SIZE, covariates = COVARIATES),
                 networks = nets, module_clinical = module_clinical)
saveRDS(fig5_obj, file.path(obj_dir, "network_analysis_object.rds"))

#############################################################
# 2. Tables S3, S4, S5
#############################################################
S3 <- bind_rows(lapply(fig5_obj$networks, function(n)
  n$edges_all %>% mutate(network = ifelse(n$timepoint == "Delta", "Co-change (ΔCLR)",
                                          paste0("Co-abundance ", n$timepoint)),
                         threshold = n$threshold,
                         sign = ifelse(pcor > 0, "positive", "negative"),
                         .before = 1)))
S4 <- bind_rows(lapply(fig5_obj$networks, function(n)
  data.frame(network = ifelse(n$timepoint == "Delta", "Co-change (ΔCLR)",
                              paste0("Co-abundance ", n$timepoint)),
             threshold = n$threshold, genus = V(n$graph)$name,
             degree = V(n$graph)$degree, module = V(n$graph)$module_label)))
S5 <- fig5_obj$module_clinical %>%
  mutate(clinical_var = unname(CLIN_LABELS[clinical_var])) %>%
  mutate(across(c(coef, se), ~ round(., 4)), pval = signif(pval, 3), qval = signif(qval, 3))

write.csv(S3, file.path(out_tab, "SupplementaryTable_S3_EdgeLists_AllNetworks.csv"), row.names = FALSE)
write.csv(S4, file.path(out_tab, "SupplementaryTable_S4_ModuleMembership.csv"),     row.names = FALSE)
write.csv(S5, file.path(out_tab, "SupplementaryTable_S5_ModuleClinicalAssociations.csv"), row.names = FALSE)
cat("✓ S3, S4, S5 written\n\n")

cat("=== NETWORK SUMMARY (after removing |pcor| < 1e-6) ===\n")
for (n in fig5_obj$networks) {
  cat(sprintf("%-16s edges=%2d (+%d/-%d)  removed=%d  modules: %s\n", n$id, n$n_edges, n$n_pos,
              n$n_neg, sum(n$edges_all$status != "retained"),
              if (length(n$modules)) paste(n$modules, collapse = ", ") else "none"))
  for (m in n$modules)
    cat("     ", m, ":", paste(sort(V(n$graph)$name[V(n$graph)$module_label == m]), collapse = ", "), "\n")
}
cat("\n=== NOMINAL MODULE–CLINICAL HITS (p < 0.05) ===\n")
print(as.data.frame(S5 %>% filter(!is.na(pval), pval < 0.05) %>%
                      select(threshold, block, module, members, clinical_var, coef, pval, qval)),
      row.names = FALSE)


#############################################################
# 3. Consistency checks (figure = S3 = S4 = S5)
#############################################################
for (nw in fig5_obj$networks) {
  net_name <- ifelse(nw$timepoint == "Delta", "Co-change (ΔCLR)",
                     paste0("Co-abundance ", nw$timepoint))
  s3_n <- S3 %>% filter(threshold == nw$threshold, network == net_name, status == "retained")
  s4_mods <- S4 %>% filter(threshold == nw$threshold, network == net_name, module != "None") %>%
    pull(module) %>% unique() %>% sort()
  s5_mods <- S5 %>% filter(threshold == nw$threshold, block == BLOCK_NAME[[nw$timepoint]],
                           module != "no module") %>% pull(module) %>% unique() %>% sort()
  stopifnot(nrow(s3_n) == nw$n_edges,
            identical(s4_mods, nw$modules),
            identical(s5_mods, nw$modules))
}
cat("✓ Consistency verified: edge counts (Fig 5 subtitles = S3) and module labels\n",
    "  (Fig 5 = S4 = S5) agree for all 6 networks.\n")


#############################################################
# 4. Figure 5 and Figure E9 panels
#############################################################
stars <- function(p) dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                                      p < 0.05 ~ "*", TRUE ~ "")

block_heatmap <- function(th, block, lims) {
  d <- fig5_obj$module_clinical %>% filter(threshold == th, block == !!block)
  if (nrow(d) == 0 || all(d$module == "no module")) {
    return(ggplot() +
             annotate("text", x = 0, y = 0, label = "No supported module", size = 3.2,
                      colour = "grey40", fontface = "italic") +
             labs(title = block) + theme_void(base_size = 9.5) +
             theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                   panel.border = element_rect(colour = "#333333", fill = NA)))
  }
  d <- d %>% mutate(clinical_label = factor(CLIN_LABELS[clinical_var], levels = CLIN_LABELS),
                    module = factor(module, levels = rev(sort(unique(module)))),
                    star = stars(pval))
  ggplot(d, aes(clinical_label, module, fill = coef)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = star), fontface = "bold", size = 3.6, vjust = 0.75) +
    scale_fill_gradient2(low = COL_NEG, mid = COL_MID, high = COL_POS, midpoint = 0,
                         limits = lims, name = "Regression\ncoefficient") +
    labs(title = block, x = NULL, y = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"))
}


suppressPackageStartupMessages({ library(ggraph); library(tidygraph); library(igraph)
  library(ggforce); library(ggplot2); library(patchwork); library(dplyr); library(ggalluvial) })
COL_POS <- "#B39B6F"; COL_NEG <- "#9C88AE"
nets    <- fig5_obj$networks
panel_dir <- file.path(fig_dir, "network_panels"); dir.create(panel_dir, showWarnings = FALSE)

pal_fun <- function(m) {
  real <- sort(setdiff(unique(m), "Unconnected"))
  cols <- c("#5C3D6E","#A68A56","#4A7FA5","#6B9E78","#C4734A","#8E6B9E","#5A8A72",
            "#B5956A","#7A5C3E","#3D6B8E","#9E6B5A","#6B7A3E","#5A6B9E","#9E5A6B","#3E7A6B")
  c(setNames(cols[seq_along(real)], real), Unconnected = "#D9D6D0")
}
TITLE <- c("T1 0.1%/10%" = "Asthma Time 1 co-abundance network · 0.1%/10%",
           "T2 0.1%/10%" = "Asthma Time 2 co-abundance network · 0.1%/10%",
           "T1 0.5%/10%" = "Asthma Time 1 co-abundance network · 0.5%/10%",
           "T2 0.5%/10%" = "Asthma Time 2 co-abundance network · 0.5%/10%",
           "Delta 0.1%/10%" = "Co-change (ΔCLR) network · 0.1%/10%",
           "Delta 0.5%/10%" = "Co-change (ΔCLR) network · 0.5%/10%")

## ---------------------------------------------------------------- NETWORKS ----
net_panel <- function(id, genus_names = FALSE, keep_isolated = FALSE) {
  n   <- nets[[id]]
  sub <- sprintf("%d edges (%d positive, %d negative) · %d modules",
                 n$n_edges, n$n_pos, n$n_neg, length(n$modules))
  base_theme <- theme_void(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.title    = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 3)),
          plot.subtitle = element_text(size = 8, colour = "grey30", hjust = 0.5, margin = margin(b = 4)),
          legend.position = "bottom", legend.box = "horizontal",
          legend.title = element_text(face = "bold", size = 8), legend.text = element_text(size = 7.5),
          plot.margin = margin(8, 8, 8, 8))
  if (n$n_edges == 0)
    return(ggplot() + annotate("text", x = 0, y = 0, size = 3.2, colour = "grey50", fontface = "italic",
                               label = paste0("No supported edges or modules\n(1 selected edge with |partial correlation| < 10⁻⁶ removed)")) +
             labs(title = TITLE[[id]], subtitle = sub) + base_theme)
  g <- n$graph
  V(g)$module_label <- ifelse(V(g)$module_label == "None", "Unconnected", V(g)$module_label)
  V(g)$degree <- igraph::degree(g)
  E(g)$weight <- abs(E(g)$pcor)
  E(g)$assoc  <- factor(ifelse(E(g)$sign == "positive", "Positive", "Negative"), levels = c("Negative", "Positive"))
  if (!keep_isolated) g <- delete_vertices(g, which(V(g)$degree == 0))
  set.seed(42)
  p <- ggraph(as_tbl_graph(g), layout = "fr") +
    geom_mark_hull(aes(x, y, group = module_label, fill = module_label,
                       filter = module_label != "Unconnected"),
                   alpha = 0.12, colour = NA, expand = unit(3, "mm"), radius = unit(3, "mm"),
                   con.type = "none", show.legend = FALSE) +
    geom_edge_link(aes(edge_colour = assoc, edge_linetype = assoc, width = weight), alpha = 0.9) +
    geom_node_point(aes(fill = module_label, size = pmax(degree, 0.5)), shape = 21, colour = "white", stroke = 0.6) +
    scale_edge_colour_manual(values = c(Negative = COL_NEG, Positive = COL_POS), name = "Association", drop = FALSE) +
    scale_edge_linetype_manual(values = c(Negative = "dashed", Positive = "solid"), name = "Association", drop = FALSE) +
    scale_edge_width(range = c(0.6, 2.4), guide = "none") +
    scale_size(range = c(4, 10), breaks = c(1, 3), name = "Connection degree") +
    scale_fill_manual(values = pal_fun(V(g)$module_label), guide = "none") +
    labs(title = TITLE[[id]], subtitle = sub) + base_theme
  if (genus_names)
    p <- p + geom_node_text(aes(label = ifelse(degree > 0, name, "")), repel = TRUE, size = 2.4,
                            fontface = "italic", colour = "#1A1A1A", segment.colour = "#1A1A1A",
                            segment.size = 0.25, min.segment.length = 0.25, box.padding = 0.25,
                            max.overlaps = Inf, seed = 42)
  p
}
save_net <- function(id, stem, keep_isolated = FALSE) {
  ggsave(file.path(panel_dir, paste0(stem, ".png")), net_panel(id, FALSE, keep_isolated),
         width = 120, height = 115, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(panel_dir, paste0(stem, "_withGenusNames.png")), net_panel(id, TRUE, keep_isolated),
         width = 120, height = 115, units = "mm", dpi = 300, bg = "white")
}
save_net("T1 0.1%/10%",    "Fig5A_T1_network_01pct")
save_net("T2 0.1%/10%",    "Fig5B_T2_network_01pct")
save_net("T1 0.5%/10%",    "E9_T1_network_05pct", keep_isolated = TRUE)
save_net("T2 0.5%/10%",    "E9_T2_network_05pct", keep_isolated = TRUE)
save_net("Delta 0.1%/10%", "E9_Delta_network_01pct")
save_net("Delta 0.5%/10%", "E9_Delta_network_05pct")

## ---------------------------------------------------------------- HEATMAPS ----
## one combined panel per threshold, all labels kept
heat_all <- function(th) {
  d <- fig5_obj$module_clinical %>% filter(threshold == th, !is.na(coef))
  lim <- max(abs(d$coef))
  (block_heatmap(th, "Time 1", c(-lim, lim)) | block_heatmap(th, "Time 2", c(-lim, lim)) |
      block_heatmap(th, "Δ (co-change)", c(-lim, lim))) +
    plot_layout(guides = "collect", widths = c(1, 1.15, 0.45))
}
ggsave(file.path(panel_dir, "Fig5C_Heatmap_01pct.png"), heat_all("0.1%/10%"),
       width = 240, height = 120, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(panel_dir, "E9_Heatmap_05pct.png"), heat_all("0.5%/10%"),
       width = 200, height = 90, units = "mm", dpi = 300, bg = "white")

## ---------------------------------------------------------------- E9A ALLUVIAL ----
memb <- function(id) { g <- nets[[id]]$graph
data.frame(genus = V(g)$name, module = V(g)$module_label, stringsAsFactors = FALSE) }
trans <- inner_join(memb("T1 0.1%/10%"), memb("T2 0.1%/10%"), by = "genus", suffix = c("_T1", "_T2")) %>%
  filter(module_T1 != "None", module_T2 != "None")
alluv <- dplyr::count(trans, module_T1, module_T2, name = "n")
cat("E9A: genera in a module at both visits:", nrow(trans), "\n")
lv  <- c(sort(unique(alluv$module_T1)), sort(unique(alluv$module_T2)))
pal <- setNames(rep(c("#5C3D6E","#A68A56","#4A7FA5","#6B9E78","#C4734A","#8E6B9E","#5A8A72",
                      "#B5956A","#7A5C3E","#3D6B8E","#9E6B5A","#6B7A3E"), length.out = length(lv)), lv)
e9a <- function(stratum_labels) {
  p <- ggplot(alluv, aes(axis1 = module_T1, axis2 = module_T2, y = n)) +
    geom_flow(aes(fill = module_T1), alpha = 0.45, colour = "white", linewidth = 0.3, curve_type = "sigmoid") +
    geom_stratum(aes(fill = after_stat(stratum)), width = 0.38, colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(limits = c("Asthma\nTime 1", "Asthma\nTime 2"), expand = c(0.15, 0.15)) +
    labs(x = NULL, y = "Number of genera",
         title = sprintf("Module membership, Time 1 → Time 2 (0.1%%/10%%; %d genera)", nrow(trans))) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
          axis.text = element_text(colour = "black"), plot.background = element_rect(fill = "white", colour = NA))
  if (stratum_labels) p <- p + geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6)
  p
}
ggsave(file.path(panel_dir, "E9A_ModuleRewiring_01pct.png"), e9a(FALSE),
       width = 90, height = 100, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(panel_dir, "E9A_ModuleRewiring_01pct_withModuleNames.png"), e9a(TRUE),
       width = 90, height = 100, units = "mm", dpi = 300, bg = "white")

cat("✓ All panels saved in:", panel_dir, "\n")
print(list.files(panel_dir))

#############################################################
# Reproducibility check against the final Supplementary Tables S3 and S5
#############################################################
ms_tab <- file.path(repo_dir, "..", "Manuscript", "Supplementary_Tables")
f3 <- file.path(ms_tab, "SupplementaryTable_S3_EdgeLists_AllNetworks.csv")
f5 <- file.path(ms_tab, "SupplementaryTable_S5_ModuleClinicalAssociations.csv")
cat("\nComparison with final Supplementary Tables:\n")
if (file.exists(f3)) {
  o <- read.csv(f3, check.names = FALSE)
  key_new <- with(S3 %>% filter(status == "retained"),
                  paste(threshold, pmin(from, to), pmax(from, to), sprintf("%.4f", pcor)))
  key_old <- with(o %>% filter(Status == "Retained"),
                  paste(Threshold, pmin(`Genus 1`, `Genus 2`), pmax(`Genus 1`, `Genus 2`),
                        sprintf("%.4f", `Partial correlation`)))
  cat(sprintf("  S3 retained edges: new %d, final table %d, identical (pair + pcor): %d\n",
              length(key_new), length(key_old), length(intersect(key_new, key_old))))
}
if (file.exists(f5)) {
  o <- read.csv(f5, check.names = FALSE) %>% filter(!is.na(Coefficient))
  n <- S5 %>% filter(!is.na(coef))
  norm <- function(x) tolower(gsub("\u2212", "-", x))   # "−" vs "-" in labels
  m <- inner_join(n %>% mutate(k = paste(threshold, module, norm(clinical_var))),
                  o %>% mutate(k = paste(Threshold, Module, norm(`Clinical variable`))), by = "k")
  cat(sprintf("  S5 module-clinical tests: %d/%d matched | max |diff| coef %.1e, P %.1e\n",
              nrow(m), nrow(n), max(abs(m$coef - m$Coefficient)), max(abs(m$pval - m$`P value`))))
}