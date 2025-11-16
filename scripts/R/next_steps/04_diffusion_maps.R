#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(1)

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
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

pm_wide <- dcast(pm_mean,
                 cytokine + cluster + time_min ~ marker,
                 value.var = "value",
                 fill = 0)

## try destiny, otherwise fall back to UMAP as proxy
have_destiny <- requireNamespace("destiny", quietly = TRUE)
if (!have_destiny && !requireNamespace("uwot", quietly = TRUE)) {
  stop("Need either destiny or uwot installed for diffusion maps.")
}

make_dm_plot <- function(cy) {
  sub <- pm_wide[cytokine == cy & !is.na(time_min)]
  if (nrow(sub) < 5) return(NULL)
  mat <- as.matrix(sub[, ..phos_cols])
  mat_scaled <- scale(mat)

  if (have_destiny) {
    dm <- destiny::DiffusionMap(mat_scaled)
    emb <- destiny::eigenvectors(dm)[, 1:2, drop = FALSE]
    df <- data.table(DC1 = emb[,1], DC2 = emb[,2])
  } else {
    emb <- uwot::umap(mat_scaled, n_neighbors = 10, min_dist = 0.1)
    df <- data.table(DC1 = emb[,1], DC2 = emb[,2])
  }

  df[, cluster := sub$cluster]
  df[, time_min := sub$time_min]

  p <- ggplot(df, aes(DC1, DC2,
                      colour = factor(time_min),
                      shape = factor(cluster))) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    labs(title = paste("Diffusion map (proxy) -", cy),
         colour = "time_min", shape = "cluster")

  out <- paste0("diffusion_map_", cy, ".png")
  ggsave(out, p, width = 7, height = 6)
}

for (cy in sort(unique(pm_wide$cytokine))) {
  if (cy %in% c("CD3","IL2","IL15","IL21")) {
    make_dm_plot(cy)
  }
}

cat("Diffusion map plots generated.\n")

