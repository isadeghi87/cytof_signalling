#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

set.seed(1)

req_files <- c(
  "UMAP_phospho_coordinates.csv",
  "functional_cluster_annotations.csv",
  "cluster_cytokine_bias.csv",
  "cluster_signaling_summary.csv"
)
missing <- req_files[!file.exists(req_files)]
if (length(missing)) {
  stop("Missing required input files: ", paste(missing, collapse = ", "))
}

dir.create("summary_figures", showWarnings = FALSE)

um <- fread("UMAP_phospho_coordinates.csv")
fun <- fread("functional_cluster_annotations.csv")
bias <- fread("cluster_cytokine_bias.csv")
sum_dt <- fread("cluster_signaling_summary.csv")

###############################################################################
## a) UMAP by phenotype
###############################################################################

## join phenotype by cytokine+cluster
um[, cluster := as.integer(cluster)]
fun[, cluster := as.integer(cluster)]

um_pheno <- merge(um, fun[, .(cytokine, cluster, phenotype)],
                  by = c("cytokine","cluster"),
                  all.x = TRUE)

p1 <- ggplot(um_pheno,
             aes(UMAP1, UMAP2, colour = phenotype)) +
  geom_point(alpha = 0.7, size = 0.5) +
  theme_bw() +
  labs(title = "UMAP coloured by functional phenotype")

ggsave("summary_figures/UMAP_by_phenotype.png", p1, width = 7, height = 6)

###############################################################################
## b) UMAP by STAT3 cytokine preference
###############################################################################

bias[, cluster_id := as.integer(cluster_id)]
um_pref <- merge(um, bias[, .(cluster_id, best_cytokine_STAT3)],
                 by.x = "cluster", by.y = "cluster_id",
                 all.x = TRUE)

p2 <- ggplot(um_pref,
             aes(UMAP1, UMAP2, colour = best_cytokine_STAT3)) +
  geom_point(alpha = 0.7, size = 0.5) +
  theme_bw() +
  labs(title = "UMAP coloured by best STAT3 cytokine")

ggsave("summary_figures/UMAP_by_cytokine_preference.png", p2,
       width = 7, height = 6)

###############################################################################
## c) Activation bars by cytokine
###############################################################################

act <- fread("cluster_activation_rankings.csv")
act[, cluster := as.integer(cluster)]

act_merged <- merge(act,
                    sum_dt[, .(cluster_id, cytokine, phenotype)],
                    by.x = c("cluster","cytokine"),
                    by.y = c("cluster_id","cytokine"),
                    all.x = TRUE)

## mark top 3 clusters per cytokine
act_merged[, rank := rank(-activation_score, ties.method = "first"),
           by = cytokine]

p3 <- ggplot(act_merged,
             aes(x = factor(cluster), y = activation_score,
                 fill = phenotype)) +
  geom_col() +
  facet_wrap(~ cytokine, scales = "free_y") +
  theme_bw() +
  labs(title = "Activation (AUC) per cluster and cytokine",
       x = "FlowSOM cluster", y = "activation_score_total") +
  geom_text(data = act_merged[rank <= 3],
            aes(label = cluster), vjust = -0.3, size = 2)

ggsave("summary_figures/activation_bars_by_cytokine.png", p3,
       width = 10, height = 7)

###############################################################################
## d) Key trajectories for STAT3 and STAT5
###############################################################################

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) {
  stop("Could not find phos_medians_by_sample_cluster.csv for trajectories.")
}
pm <- fread(pm_path)

id_cols <- c("file_base","cluster","time_min","cytokine")
phos_cols <- setdiff(names(pm), id_cols)

pm_long <- melt(pm,
                id.vars = id_cols,
                measure.vars = phos_cols,
                variable.name = "marker",
                value.name = "value")

pm_mean <- pm_long[, .(value = mean(value, na.rm = TRUE)),
                   by = .(cytokine, cluster, time_min, marker)]

## compute logFC vs 0 for STAT3- and STAT5-linked channels
marker_stat3 <- "Er166Di"
marker_stat5 <- "Yb173Di"

base <- pm_mean[marker %in% c(marker_stat3, marker_stat5) &
                time_min %in% c("0",0),
                .(baseline = mean(value, na.rm = TRUE)),
                by = .(cytokine, cluster, marker)]

pm_fc <- merge(pm_mean[marker %in% c(marker_stat3, marker_stat5)],
               base,
               by = c("cytokine","cluster","marker"),
               all.x = TRUE)
pm_fc[, logFC := value - baseline]
pm_fc[, time_num := as.numeric(time_min)]

## identify top 3 STAT3 and STAT5 clusters per cytokine using AUC
auc_stat <- pm_fc[, {
  ord <- order(time_num)
  t_ord <- time_num[ord]
  v_ord <- logFC[ord]
  auc <- if (length(t_ord) > 1) {
    sum(diff(t_ord) * (head(v_ord, -1) + tail(v_ord, -1)) / 2)
  } else NA_real_
  .(AUC = auc)
}, by = .(cytokine, cluster, marker)]

top_sel <- list()
for (cy in sort(unique(auc_stat$cytokine))) {
  for (mk in c(marker_stat3, marker_stat5)) {
    sub <- auc_stat[cytokine == cy & marker == mk]
    sub <- sub[order(-AUC)]
    cl_top <- head(sub$cluster, 3)
    top_sel[[length(top_sel)+1]] <- data.table(
      cytokine = cy,
      marker = mk,
      cluster = cl_top
    )
  }
}
top_sel_dt <- rbindlist(top_sel, use.names = TRUE, fill = TRUE)

traj_dt <- merge(pm_fc, top_sel_dt,
                 by = c("cytokine","cluster","marker"),
                 all.y = TRUE)

## plot trajectories
p4 <- ggplot(traj_dt,
             aes(x = time_num, y = logFC,
                 colour = factor(cluster),
                 linetype = marker, group = interaction(cluster, marker))) +
  geom_line() +
  facet_wrap(~ cytokine, scales = "free_y") +
  theme_bw() +
  labs(title = "Key STAT3 / STAT5 trajectories",
       x = "time_min", y = "logFC vs 0",
       colour = "cluster")

ggsave("summary_figures/key_trajectories_STAT3_STAT5.png", p4,
       width = 10, height = 7)

cat("Summary figures written to summary_figures/.\n")
