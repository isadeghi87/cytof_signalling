#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(uwot)
  library(mgcv)
})

set.seed(1)

out_txt_sample <- "check_sample_table.txt"
out_txt_phos   <- "check_phos_medians.txt"

sink(out_txt_sample, split = FALSE)
cat("=== CHECK: sample_table.csv ===\n")

st <- fread("sample_table.csv")

if (!"file_base" %in% names(st)) {
  cat("INFO: 'file_base' missing; deriving from 'file'.\n")
  if ("file" %in% names(st)) {
    st[, file_base := basename(file)]
  } else {
    stop("sample_table.csv must contain 'file' or 'file_base'.")
  }
}

if (!"cytokine" %in% names(st)) {
  cat("INFO: 'cytokine' missing; deriving from 'stim'.\n")
  if (!"stim" %in% names(st)) {
    stop("sample_table.csv must contain 'cytokine' or 'stim'.")
  }
  stim_to_cytokine <- function(stim_vec) {
    v <- tolower(as.character(stim_vec))
    res <- ifelse(grepl("il21", v), "IL21",
           ifelse(grepl("il15", v), "IL15",
           ifelse(grepl("il2",  v), "IL2",
           ifelse(grepl("pac", v), "PAC",
           ifelse(grepl("cd3", v) & !grepl("il", v), "CD3",
           ifelse(grepl("unstim", v) | v == "" | is.na(v),
                  "UNSTIM", "OTHER"))))))
    res
  }
  st[, cytokine := stim_to_cytokine(stim)]
}

req_cols_st <- c("file_base","batch","is_anchor","time_min","cytokine")
missing_st <- setdiff(req_cols_st, names(st))
if (length(missing_st)) {
  cat("WARNING: sample_table missing columns:", paste(missing_st, collapse=", "), "\n")
}

cat("\nHead(sample_table):\n")
print(head(st))
cat("\nStr(sample_table):\n")
print(str(st))

cat("\nUnique cytokine values:\n")
print(sort(unique(st$cytokine)))

cat("\nUnique time_min values:\n")
print(sort(unique(st$time_min)))

cat("\nCheck is_anchor values:\n")
print(table(st$is_anchor, useNA = "ifany"))

cat("\nCytokine × time_min coverage:\n")
tab_cov <- as.data.table(table(cytokine = st$cytokine,
                               time_min = st$time_min,
                               useNA = "no"))
print(tab_cov)

sink()

## Coverage heatmap
cov_mat <- dcast(tab_cov, cytokine ~ time_min, value.var = "N", fill = 0)
cov_mat_rownames <- cov_mat$cytokine
cov_mat <- as.matrix(cov_mat[ , -1])
rownames(cov_mat) <- cov_mat_rownames

dir.create("diagnostic_plots", showWarnings = FALSE)

png("diagnostic_plots/cytokine_time_coverage.png", width = 1200, height = 800)
pheatmap(cov_mat, main = "Cytokine × time_min coverage")
dev.off()
pdf("diagnostic_plots/cytokine_time_coverage.pdf", width = 10, height = 7)
pheatmap(cov_mat, main = "Cytokine × time_min coverage")
dev.off()

###############################################################################
## PHOSPHO MEDIANS CHECK
###############################################################################

sink(out_txt_phos, split = FALSE)
cat("=== CHECK: phos_medians_by_sample_cluster.csv ===\n")

phos_path <- if (file.exists("diff_phos_results/phos_medians_by_sample_cluster.csv")) {
  "diff_phos_results/phos_medians_by_sample_cluster.csv"
} else {
  "phospho_medians/phos_medians_by_sample_cluster.csv"
}

pm <- fread(phos_path)

req_cols_pm <- c("file_base","cluster","time_min","cytokine")
missing_pm <- setdiff(req_cols_pm, names(pm))
if (length(missing_pm)) {
  cat("WARNING: phos_medians missing columns:", paste(missing_pm, collapse=", "), "\n")
}

cat("\nSource:", phos_path, "\n")
cat("\nHead(phos_medians):\n")
print(head(pm))
cat("\nStr(phos_medians):\n")
print(str(pm))

phos_markers <- setdiff(names(pm), req_cols_pm)
cat("\nPhospho markers detected:\n")
print(phos_markers)

cat("\nCluster coverage per cytokine × time_min:\n")
cov_clusters <- pm[, .(n_clusters = uniqueN(cluster)),
                   by = .(cytokine, time_min)]
print(cov_clusters)

cat("\nCheck all 20 clusters per cytokine × time_min (n_clusters == 20):\n")
print(cov_clusters[n_clusters < 20])

cat("\nDuplicate rows (file_base × cluster × cytokine × time_min):\n")
pm[, dup_key := paste(file_base, cluster, cytokine, time_min, sep = "_")]
dup_tab <- pm[, .N, by = dup_key][N > 1]
print(dup_tab)

cat("\nRange check for phospho markers (min, max):\n")
rng <- pm[, lapply(.SD, function(x) c(min = min(x, na.rm = TRUE),
                                      max = max(x, na.rm = TRUE))),
          .SDcols = phos_markers]
print(rng)

sink()

###############################################################################
## METACLUSTER MEDIANS CHECK
###############################################################################

mc <- fread("flowsom_results/metacluster_medians_asinh.csv")
mc_clusters <- unique(mc$cluster)

mc_msg <- paste0("Metacluster_medians: n clusters = ", length(mc_clusters),
                 " (expected ~20).")
cat(mc_msg, "\n", file = out_txt_phos, append = TRUE)

lineage_markers <- setdiff(names(mc), "cluster")
sd_by_marker <- mc[, lapply(.SD, sd), .SDcols = lineage_markers]
flat_markers <- names(sd_by_marker)[which(as.numeric(sd_by_marker[1, ]) < 1e-3)]

cat("\nLineage marker SD across metaclusters:\n",
    file = out_txt_phos, append = TRUE)
capture.output(print(sd_by_marker),
               file = out_txt_phos, append = TRUE)

cat("\nPotential flat / mis-annotated markers (SD < 1e-3):\n",
    file = out_txt_phos, append = TRUE)
capture.output(print(flat_markers),
               file = out_txt_phos, append = TRUE)

###############################################################################
## LIMMA RESULTS CHECK
###############################################################################

if (file.exists("diff_phos_results/diff_phospho_limma_results.csv")) {
  lim <- fread("diff_phos_results/diff_phospho_limma_results.csv")
  cat("\n=== CHECK: diff_phospho_limma_results.csv ===\n",
      file = out_txt_phos, append = TRUE)
  contr <- sort(unique(lim$contrast))
  cat("Contrasts present:\n", file = out_txt_phos, append = TRUE)
  capture.output(print(contr), file = out_txt_phos, append = TRUE)

  sig <- lim[adj.P.Val < 0.05]
  sig_counts <- sig[, .N, by = .(cluster, cytokine)]
  cat("\nSignificant markers per cluster × cytokine (FDR < 0.05):\n",
      file = out_txt_phos, append = TRUE)
  capture.output(print(sig_counts), file = out_txt_phos, append = TRUE)
}

###############################################################################
## PHENOGRAPH LABELS CHECK
###############################################################################

if (file.exists("phenograph_results/phenograph_labels.csv")) {
  pl <- fread("phenograph_results/phenograph_labels.csv")
  cat("\n=== CHECK: phenograph_labels.csv ===\n",
      file = out_txt_phos, append = TRUE)
  capture.output(print(head(pl)), file = out_txt_phos, append = TRUE)

  ## assume same order as dimred_metadata
  dim_meta <- readRDS("dimred_flowsom/dimred_metadata.rds")
  if (nrow(pl) == nrow(dim_meta)) {
    dim_dt <- as.data.table(dim_meta)
    dim_dt[, flowsom_cluster := as.integer(as.character(cluster))]
    pl[, phenograph_cluster := cluster]
    pl[, cluster := NULL]
    full_lab <- cbind(dim_dt, pl)

    ## contingency table FlowSOM vs PhenoGraph
    tab_fp <- full_lab[, .N, by = .(flowsom_cluster, phenograph_cluster)]
    tab_wide <- dcast(tab_fp, flowsom_cluster ~ phenograph_cluster,
                      value.var = "N", fill = 0)

    cat("\nContingency table FlowSOM vs PhenoGraph (head):\n",
        file = out_txt_phos, append = TRUE)
    capture.output(print(head(tab_wide)),
                   file = out_txt_phos, append = TRUE)

    ## save full contingency
    fwrite(tab_fp, "FlowSOM_vs_Phenograph_contingency.csv")

    ## heatmap
    mat_fp <- as.matrix(tab_wide[, -1])
    rownames(mat_fp) <- tab_wide$flowsom_cluster
    png("diagnostic_plots/flowsom_vs_phenograph_heatmap.png",
        width = 1000, height = 800)
    pheatmap(mat_fp, main = "FlowSOM vs PhenoGraph")
    dev.off()
    pdf("diagnostic_plots/flowsom_vs_phenograph_heatmap.pdf",
        width = 9, height = 7)
    pheatmap(mat_fp, main = "FlowSOM vs PhenoGraph")
    dev.off()

    ## crude ARI / NMI approximations (using base R)
    ## convert to vectors of labels
    fs_lab <- full_lab$flowsom_cluster
    ph_lab <- full_lab$phenograph_cluster

    ## Adjusted Rand Index (simple implementation)
    ari_simple <- function(x, y) {
      tab <- table(x, y)
      n <- sum(tab)
      sum_comb <- sum(choose(tab, 2))
      sum_row  <- sum(choose(rowSums(tab), 2))
      sum_col  <- sum(choose(colSums(tab), 2))
      expected <- sum_row * sum_col / choose(n, 2)
      max_val  <- (sum_row + sum_col) / 2
      (sum_comb - expected) / (max_val - expected)
    }

    ari_val <- ari_simple(fs_lab, ph_lab)
    cat("\nApproximate ARI FlowSOM vs PhenoGraph:", ari_val, "\n",
        file = "cluster_stability_report.txt")
  } else {
    cat("\nWARNING: phenograph_labels and dimred_metadata have different row counts; cannot compute contingency.\n",
        file = out_txt_phos, append = TRUE)
  }

  cl_col <- if ("phenograph_cluster" %in% names(pl)) "phenograph_cluster" else "cluster"
  dist_tab <- pl[, .N, by = cl_col]
  cat("\nPhenograph cluster distribution:\n",
      file = out_txt_phos, append = TRUE)
  capture.output(print(dist_tab), file = out_txt_phos, append = TRUE)
}

###############################################################################
## DIAGNOSTIC QA PLOTS FROM PHOSPHO MEDIANS
###############################################################################

## Cluster "counts" per condition (here: number of rows per comb)
cl_counts <- pm[, .N, by = .(cluster, cytokine, time_min)]

gg_cl <- ggplot(cl_counts,
                aes(x = time_min, y = N, fill = factor(cluster))) +
  geom_col(position = "stack") +
  facet_wrap(~ cytokine, scales = "free_y") +
  theme_bw() +
  labs(title = "Cluster representation per cytokine × time",
       x = "time_min", y = "count (rows)")

ggsave("diagnostic_plots/cluster_counts_per_condition.png", gg_cl,
       width = 10, height = 6)
ggsave("diagnostic_plots/cluster_counts_per_condition.pdf", gg_cl,
       width = 10, height = 6)

## Marker distribution plots
dir.create("diagnostic_plots/marker_distributions", showWarnings = FALSE)

pm_long <- melt(pm,
                id.vars = req_cols_pm,
                measure.vars = phos_markers,
                variable.name = "marker",
                value.name = "value")

for (mk in phos_markers) {
  sub <- pm_long[marker == mk]

  p_hist <- ggplot(sub, aes(x = value)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "black") +
    theme_bw() + labs(title = paste("Histogram:", mk), x = "asinh value")

  p_box <- ggplot(sub, aes(x = cytokine, y = value)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() + labs(title = paste("Boxplot by cytokine:", mk),
                      x = "cytokine", y = "asinh value")

  p_time <- ggplot(sub[, .(mean = mean(value, na.rm = TRUE),
                           sd = sd(value, na.rm = TRUE),
                           n = .N),
                      by = .(cytokine, time_min)],
                   aes(x = time_min, y = mean, colour = cytokine, group = cytokine)) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    labs(title = paste("Time trend:", mk),
         x = "time_min", y = "mean asinh value")

  fname_base <- file.path("diagnostic_plots/marker_distributions", mk)
  ggsave(paste0(fname_base, "_hist.png"), p_hist, width = 6, height = 4)
  ggsave(paste0(fname_base, "_hist.pdf"), p_hist, width = 6, height = 4)
  ggsave(paste0(fname_base, "_box.png"), p_box, width = 6, height = 4)
  ggsave(paste0(fname_base, "_box.pdf"), p_box, width = 6, height = 4)
  ggsave(paste0(fname_base, "_time.png"), p_time, width = 6, height = 4)
  ggsave(paste0(fname_base, "_time.pdf"), p_time, width = 6, height = 4)
}

###############################################################################
## PCA & UMAP ON PHOSPHO MEDIANS
###############################################################################

pm_mat <- as.matrix(pm[, ..phos_markers])
pm_scaled <- scale(pm_mat)

pca <- prcomp(pm_scaled, center = FALSE, scale. = FALSE)
pca_scores <- as.data.table(pca$x[, 1:2])
setnames(pca_scores, c("PC1","PC2"))
pca_scores <- cbind(pm[, .(file_base, cluster, time_min, cytokine)], pca_scores)

fwrite(pca_scores, "PCA_phospho_coordinates.csv")

gp_pca <- ggplot(pca_scores, aes(PC1, PC2, colour = cytokine)) +
  geom_point(alpha = 0.7, size = 0.8) +
  theme_bw() +
  labs(title = "PCA on phospho medians")

ggsave("PCA_phospho.png", gp_pca, width = 7, height = 6)
ggsave("PCA_phospho.pdf", gp_pca, width = 7, height = 6)

um <- umap(pm_scaled, n_neighbors = 15, min_dist = 0.3)
um_dt <- data.table(UMAP1 = um[,1], UMAP2 = um[,2])
um_dt <- cbind(pm[, .(file_base, cluster, time_min, cytokine)], um_dt)
fwrite(um_dt, "UMAP_phospho_coordinates.csv")

gp_um <- ggplot(um_dt, aes(UMAP1, UMAP2, colour = cytokine)) +
  geom_point(alpha = 0.7, size = 0.8) +
  theme_bw() +
  labs(title = "UMAP on phospho medians")

ggsave("UMAP_phospho.png", gp_um, width = 7, height = 6)
ggsave("UMAP_phospho.pdf", gp_um, width = 7, height = 6)

###############################################################################
## LOGFC HEATMAPS
###############################################################################

dir.create("logFC_heatmaps", showWarnings = FALSE)

pm_mean <- pm_long[, .(value = mean(value, na.rm = TRUE)),
                   by = .(cytokine, cluster, time_min, marker)]

baseline <- pm_mean[time_min %in% c("0", 0), .(baseline = mean(value, na.rm = TRUE)),
                    by = .(cytokine, cluster, marker)]

pm_fc <- merge(pm_mean, baseline,
               by = c("cytokine","cluster","marker"),
               all.x = TRUE)
pm_fc[, logFC := value - baseline]

for (cy in sort(unique(pm_fc$cytokine))) {
  sub_cy <- pm_fc[cytokine == cy]
  for (tt in sort(unique(sub_cy$time_min))) {
    if (is.na(tt) || tt %in% c("0",0)) next
    sub_tt <- sub_cy[time_min == tt]
    if (nrow(sub_tt) == 0) next

    mat <- dcast(sub_tt,
                 cluster ~ marker,
                 value.var = "logFC",
                 fill = 0)
    rn <- mat$cluster
    mat <- as.matrix(mat[, -1])
    rownames(mat) <- rn

    fname <- paste0("logFC_heatmaps/logFC_", cy, "_t", tt)
    png(paste0(fname, ".png"), width = 1200, height = 900)
    pheatmap(mat, main = paste("logFC", cy, "time", tt),
             cluster_rows = TRUE, cluster_cols = TRUE)
    dev.off()
    pdf(paste0(fname, ".pdf"), width = 10, height = 7)
    pheatmap(mat, main = paste("logFC", cy, "time", tt),
             cluster_rows = TRUE, cluster_cols = TRUE)
    dev.off()
  }
}

###############################################################################
## TRAJECTORY ANALYSIS WITH GAM
###############################################################################

dir.create("trajectory_plots", showWarnings = FALSE)

auc_list <- list()

for (cy in sort(unique(pm_long$cytokine))) {
  for (mk in phos_markers) {
    sub <- pm_long[cytokine == cy & marker == mk & !is.na(time_min)]
    if (nrow(sub) < 10 || uniqueN(sub$time_min) < 3) next

    ## Fit GAM per cluster
    clusters_here <- sort(unique(sub$cluster))
    pred_all <- list()
    auc_cy_mk <- list()

    for (cl in clusters_here) {
      sub_cl <- sub[cluster == cl]
      if (uniqueN(sub_cl$time_min) < 3) next
      sub_cl$time_num <- as.numeric(sub_cl$time_min)
      try({
        fit <- gam(value ~ s(time_num, k = min(5, uniqueN(sub_cl$time_num))),
                   data = sub_cl)
        grid_t <- seq(min(sub_cl$time_num), max(sub_cl$time_num), length.out = 50)
        pred <- predict(fit, newdata = data.frame(time_num = grid_t))
        pred_all[[length(pred_all) + 1]] <- data.table(
          cytokine = cy, marker = mk, cluster = cl,
          time_num = grid_t, value = pred
        )
        ## trapezoidal AUC
        auc_val <- sum(diff(grid_t) * (head(pred, -1) + tail(pred, -1)) / 2)
        auc_cy_mk[[length(auc_cy_mk) + 1]] <- data.table(
          cytokine = cy, marker = mk, cluster = cl, AUC = auc_val
        )
      }, silent = TRUE)
    }

    if (length(pred_all) == 0) next

    pred_dt <- rbindlist(pred_all)
    auc_dt  <- rbindlist(auc_cy_mk)
    auc_list[[length(auc_list) + 1]] <- auc_dt

    p_traj <- ggplot(pred_dt,
                     aes(x = time_num, y = value,
                         colour = factor(cluster), group = cluster)) +
      geom_line() +
      theme_bw() +
      labs(title = paste("Trajectory:", cy, mk),
           x = "time_min", y = "smooth asinh value",
           colour = "cluster")

    fname_tr <- paste0("trajectory_plots/trajectory_", cy, "_", mk)
    ggsave(paste0(fname_tr, ".png"), p_traj, width = 8, height = 6)
    ggsave(paste0(fname_tr, ".pdf"), p_traj, width = 8, height = 6)
  }
}

if (length(auc_list)) {
  auc_all <- rbindlist(auc_list, use.names = TRUE, fill = TRUE)
  ## Aggregate to cluster-level activation score per cytokine
  act_rank <- auc_all[, .(activation_score = sum(abs(AUC), na.rm = TRUE)),
                      by = .(cytokine, cluster)]
  setorder(act_rank, cytokine, -activation_score)
  fwrite(act_rank, "cluster_activation_rankings.csv")
}

###############################################################################
## FUNCTIONAL CLUSTER ANNOTATIONS & CYTOKINE CONTRASTS
###############################################################################

## define marker groups for interpretation
stat3_markers  <- c("pSTAT3", "Er166Di", "Er170Di")
stat5_markers  <- c("pSTAT5", "Yb173Di")
tcr_markers    <- c("pZAP70","pSLP76","pPLCy","pNFATc1")
metab_markers  <- c("pAkt","pGSK3","pAMPK","pPGC1","pCREB")
stress_markers <- c("IkBa")

## use logFC table pm_fc to derive simple phenotypes
## focus on an early time (e.g., t2) where contrasts exist
early_fc <- pm_fc[time_min %in% c("2", 2)]

score_group <- function(df, markers) {
  m_int <- intersect(markers, df$marker)
  if (!length(m_int)) return(NA_real_)
  df[marker %in% m_int, mean(logFC, na.rm = TRUE)]
}

fun_list <- list()
for (cy in sort(unique(early_fc$cytokine))) {
  for (cl in sort(unique(early_fc$cluster))) {
    sub <- early_fc[cytokine == cy & cluster == cl]
    if (!nrow(sub)) next
    stat3 <- score_group(sub, stat3_markers)
    stat5 <- score_group(sub, stat5_markers)
    tcr   <- score_group(sub, tcr_markers)
    metab <- score_group(sub, metab_markers)
    stress<- score_group(sub, stress_markers)
    fun_list[[length(fun_list)+1]] <- data.table(
      cytokine = cy, cluster = cl,
      stat3_score = stat3,
      stat5_score = stat5,
      tcr_score   = tcr,
      metab_score = metab,
      stress_score= stress
    )
  }
}

fun_dt <- if (length(fun_list)) rbindlist(fun_list, use.names = TRUE, fill = TRUE) else
  data.table()

## combine with activation rankings to label super-responders
if (nrow(fun_dt) && exists("act_rank")) {
  fun_merged <- merge(fun_dt, act_rank,
                      by = c("cytokine","cluster"),
                      all.x = TRUE)
} else {
  fun_merged <- fun_dt
}

## assign dominant phenotype by max group score
if (nrow(fun_merged)) {
  fun_merged[, phenotype := {
    vals <- c(stat3_score, stat5_score, tcr_score, metab_score, stress_score)
    labs <- c("STAT3-dominant","STAT5-dominant","TCR-proximal",
              "Metabolic/CREB","Stress/low-responder")
    if (all(is.na(vals))) "Unclassified" else labs[which.max(ifelse(is.na(vals), -Inf, vals))]
  }, by = .(cytokine, cluster)]

  fwrite(fun_merged, "functional_cluster_annotations.csv")
}

## cytokine contrasts at logFC level: IL21-IL2, IL2-IL15, IL21-CD3, CD3-IL15
dir.create("cytokine_contrast_heatmaps", showWarnings = FALSE)

contrast_pairs <- list(
  c("IL21","IL2"),
  c("IL2","IL15"),
  c("IL21","CD3"),
  c("CD3","IL15")
)

for (pair in contrast_pairs) {
  cyA <- pair[1]; cyB <- pair[2]
  subA <- pm_fc[cytokine == cyA]
  subB <- pm_fc[cytokine == cyB]
  key_cols <- c("cluster","time_min","marker")
  setkeyv(subA, c("cytokine", key_cols))
  setkeyv(subB, c("cytokine", key_cols))
  merged <- merge(subA, subB,
                  by = key_cols,
                  suffixes = c("_A","_B"))
  if (!nrow(merged)) next
  merged[, diff := logFC_A - logFC_B]
  for (tt in sort(unique(merged$time_min))) {
    sub_tt <- merged[time_min == tt]
    if (!nrow(sub_tt) || tt %in% c("0",0) || is.na(tt)) next
    mat <- dcast(sub_tt, cluster ~ marker,
                 value.var = "diff", fill = 0)
    rn <- mat$cluster
    mat <- as.matrix(mat[, -1])
    rownames(mat) <- rn
    fname <- paste0("cytokine_contrast_heatmaps/contrast_",
                    cyA, "_vs_", cyB, "_t", tt)
    png(paste0(fname, ".png"), width = 1200, height = 900)
    pheatmap(mat,
             main = paste(cyA, "-", cyB, "time", tt),
             cluster_rows = TRUE, cluster_cols = TRUE)
    dev.off()
    pdf(paste0(fname, ".pdf"), width = 10, height = 7)
    pheatmap(mat,
             main = paste(cyA, "-", cyB, "time", tt),
             cluster_rows = TRUE, cluster_cols = TRUE)
    dev.off()
  }
}

###############################################################################
## FINAL REPORT SKELETON
###############################################################################

report_lines <- c(
  "# CyTOF Signaling Diagnostics Report",
  "",
  "This report summarizes integrity checks, diagnostics, and signaling phenotypes",
  "computed from the existing CSV outputs (no re-running of CytoNorm or FlowSOM).",
  "",
  "Key outputs:",
  "- `check_sample_table.txt`, `check_phos_medians.txt`",
  "- `diagnostic_plots/` (coverage, cluster counts, marker distributions, FlowSOM vs PhenoGraph)",
  "- `PCA_phospho.*`, `UMAP_phospho.*`",
  "- `logFC_heatmaps/` and `cytokine_contrast_heatmaps/`",
  "- `trajectory_plots/` and `cluster_activation_rankings.csv`",
  "- `functional_cluster_annotations.csv`",
  "- `FlowSOM_vs_Phenograph_contingency.csv`, `cluster_stability_report.txt`"
)

writeLines(report_lines, "FINAL_REPORT.md")

cat("Diagnostics and analyses completed.\n")
