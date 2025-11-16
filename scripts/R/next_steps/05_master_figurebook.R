#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

set.seed(1)

## PCA and UMAP on phospho medians
pca_scores <- fread("PCA_phospho_coordinates.csv")
umap_scores <- fread("UMAP_phospho_coordinates.csv")

gp_pca <- ggplot(pca_scores, aes(PC1, PC2, colour = cytokine)) +
  geom_point(alpha = 0.7, size = 0.8) +
  theme_bw() +
  labs(title = "PCA on phospho medians")

gp_umap <- ggplot(umap_scores, aes(UMAP1, UMAP2, colour = cytokine)) +
  geom_point(alpha = 0.7, size = 0.8) +
  theme_bw() +
  labs(title = "UMAP on phospho medians")

## cluster composition per cytokine/time
pm <- fread("diff_phos_results/phos_medians_by_sample_cluster.csv")
cl_counts <- pm[, .N, by = .(cluster, cytokine, time_min)]
gp_counts <- ggplot(cl_counts,
                    aes(x = time_min, y = N, fill = factor(cluster))) +
  geom_col() +
  facet_wrap(~ cytokine, scales = "free_y") +
  theme_bw() +
  labs(title = "Cluster representation per cytokine × time",
       x = "time_min", y = "count")

## functional phenotypes
fun <- fread("cluster_functional_annotations.csv")
gp_fun <- ggplot(fun,
                 aes(x = cluster, y = cytokine, fill = phenotype)) +
  geom_tile() +
  theme_bw() +
  labs(title = "Functional cluster annotations") +
  scale_fill_brewer(palette = "Set2")

## open multi-page PDF
pdf("MASTER_FIGUREBOOK.pdf", width = 8, height = 6)

print(gp_pca)
print(gp_umap)
print(gp_counts)
print(gp_fun)

## logFC heatmaps summary: show one example per cytokine
for (cy in c("CD3","IL2","IL15","IL21")) {
  files <- list.files("logFC_heatmaps",
                      pattern = paste0("logFC_", cy, "_t"),
                      full.names = TRUE)
  files <- files[grepl("\\.png$", files)]
  if (!length(files)) next
  ## just indicate file names; actual heatmaps are separate PNGs
  df <- data.table(file = basename(files))
  p <- ggplot(df, aes(x = 1, y = file)) +
    geom_text(aes(label = file)) +
    theme_void() +
    labs(title = paste("logFC heatmaps for", cy))
  print(p)
}

dev.off()

cat("MASTER_FIGUREBOOK.pdf created.\n")

