#!/usr/bin/env Rscript
# t-SNE / UMAP + clustering for sampled events across samples

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
})

option_list <- list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv", help="Path to sample_table.csv"),
  make_option(c("--channels"), type="character", default="channels_to_adjust.tsv", help="Channels file"),
  make_option(c("--outdir"), type="character", default="qc_plots", help="Output dir for plots"),
  make_option(c("--n-samples"), type="integer", default=20, help="Max number of samples to include"),
  make_option(c("--n-events"), type="integer", default=1000, help="Events per sample to sample") ,
  make_option(c("--k"), type="integer", default=10, help="k for kmeans clustering")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!file.exists(opt$`sample-table`)) stop('sample_table.csv not found')
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

tab <- fread(opt$`sample-table`)
channels <- if (file.exists(opt$channels)) as.character(fread(opt$channels, header = FALSE)[[1]]) else stop('channels file not found')

if (!requireNamespace('flowCore', quietly = TRUE)) stop('flowCore required')
if (!requireNamespace('uwot', quietly = TRUE)) stop('uwot required for UMAP (install from CRAN)')
if (!requireNamespace('Rtsne', quietly = TRUE)) message('Rtsne not available; t-SNE will be skipped')

sel_samples <- head(tab$file, opt$`n-samples`)

all_list <- list()
for (s in sel_samples) {
  if (!file.exists(s)) { message('[WARN] sample missing: ', s); next }
  f <- flowCore::read.FCS(s, transformation = FALSE, truncate_max_range = FALSE)
  df <- as.data.frame(flowCore::exprs(f))
  present <- intersect(channels, colnames(df))
  if (length(present) == 0) { message('[WARN] no channels for ', s); next }
  n <- min(nrow(df), opt$`n-events`)
  idx <- if (nrow(df) > n) sample(seq_len(nrow(df)), n) else seq_len(nrow(df))
  sub <- df[idx, present, drop = FALSE]
  dt <- as.data.table(sub)
  dt[, sample := basename(s)]
  dt[, event_id := .I]
  all_list[[basename(s)]] <- dt
}

if (length(all_list) == 0) stop('No samples read for dimensionality reduction')
big <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

# select channels actually present
present_channels <- intersect(channels, colnames(big))
if (length(present_channels) == 0) stop('None of the requested channels are present in the data')
mat <- as.matrix(big[, ..present_channels, with = FALSE])

# PCA for speed/stability
mat[is.na(mat)] <- 0
pc <- prcomp(mat, center = TRUE, scale. = TRUE)
pcn <- pc$x[, 1:min(50, ncol(pc$x)), drop = FALSE]

# UMAP
um <- uwot::umap(pcn, n_components = 2, verbose = TRUE)
um_dt <- data.table(UMAP1 = um[,1], UMAP2 = um[,2], sample = big$sample)
um_dt[, cluster := as.factor(kmeans(pcn, centers = opt$k)$cluster)]

plt1 <- ggplot(um_dt, aes(x = UMAP1, y = UMAP2, color = cluster)) + geom_point(size=0.6, alpha=0.7) + theme_minimal() + labs(title='UMAP colored by cluster')
plt2 <- ggplot(um_dt, aes(x = UMAP1, y = UMAP2, color = sample)) + geom_point(size=0.6, alpha=0.7) + theme_minimal() + labs(title='UMAP colored by sample')

ggsave(file.path(opt$outdir, 'umap_clusters.png'), plt1, width = 8, height = 6, dpi = 150)
ggsave(file.path(opt$outdir, 'umap_by_sample.png'), plt2, width = 10, height = 6, dpi = 150)
message('[ok] UMAP plots saved to ', opt$outdir)

# t-SNE (optional)
if (requireNamespace('Rtsne', quietly = TRUE)) {
  # Rtsne cannot accept duplicate rows — deduplicate PCA rows then map back
  df_pcn <- as.data.frame(pcn)
  key_all <- do.call(paste, c(df_pcn, sep = "__"))
  dup <- duplicated(key_all)
  if (any(dup)) {
    uniq_df <- df_pcn[!dup, , drop = FALSE]
    message(sprintf('[INFO] Detected %d duplicate PCA rows — running Rtsne on %d unique rows', sum(dup), nrow(uniq_df)))
    ts <- Rtsne::Rtsne(as.matrix(uniq_df), perplexity = 30, verbose = TRUE)
    # map back
    key_uniq <- do.call(paste, c(as.data.frame(uniq_df), sep = "__"))
    map_idx <- match(key_all, key_uniq)
    ts_coords_full <- ts$Y[map_idx, , drop = FALSE]
  } else {
    ts <- Rtsne::Rtsne(pcn, perplexity = 30, verbose = TRUE)
    ts_coords_full <- ts$Y
  }

  ts_dt <- data.table(TSNE1 = ts_coords_full[,1], TSNE2 = ts_coords_full[,2], sample = big$sample)
  ts_dt[, cluster := as.factor(kmeans(pcn, centers = opt$k)$cluster)]
  pta <- ggplot(ts_dt, aes(x = TSNE1, y = TSNE2, color = cluster)) + geom_point(size=0.6, alpha=0.7) + theme_minimal() + labs(title='t-SNE colored by cluster')
  ggsave(file.path(opt$outdir, 'tsne_clusters.png'), pta, width = 8, height = 6, dpi = 150)
  message('[ok] t-SNE plot saved to ', opt$outdir)
} else message('[INFO] Rtsne not installed; skipped t-SNE')
## end of script
