#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(flowCore)
  library(ggplot2)
  library(Rtsne)
  library(uwot)
  library(FNN)
})

opt_list <- list(
  make_option("--norm-dir", type="character", default="normalized_fcs"),
  make_option("--sample-table", type="character", default="sample_table.csv"),
  make_option("--flowsom-rds", type="character", default="flowsom_results/flowsom_fitted.rds"),
  make_option("--cluster-labels", type="character", default="flowsom_results/cluster_labels_template.csv"),
  make_option("--outdir", type="character", default="tsne_umap_plots"),
  make_option("--cofactor", type="double", default=5),
  make_option("--n-events", type="integer", default=60000),
  make_option("--seed", type="integer", default=1)
)

opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
set.seed(opt$seed)

################################################################################
# Load FlowSOM and metadata
################################################################################

fs <- readRDS(opt$`flowsom-rds`)
meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) meta[, file_base := basename(file)]

meta_s <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_s$file_base))
idx <- file.exists(norm_files)
meta_s <- meta_s[idx]
norm_files <- norm_files[idx]

################################################################################
# Identify lineage channels
################################################################################

ff <- read.FCS(norm_files[1], transformation=FALSE)
par <- pData(parameters(ff))
cn <- par$name

lineage_ch <- intersect(
  c("Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
    "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di",
    "Cd113Di","Gd157Di"), cn
)

tf <- arcsinhTransform(a=0, b=1/opt$cofactor, c=0)
tlist <- transformList(lineage_ch, tf)

################################################################################
# Downsample
################################################################################

per_file <- round(opt$`n-events` / length(norm_files))
expr_list <- list()
meta_list <- list()

for (i in seq_along(norm_files)) {
  ff0 <- read.FCS(norm_files[i], transformation=FALSE)
  ffT <- transform(ff0, tlist)
  M <- exprs(ffT)[, lineage_ch, drop=FALSE]

  take <- min(nrow(M), per_file)
  idx2 <- sample(seq_len(nrow(M)), take)

  expr_list[[i]] <- M[idx2, ]
  meta_list[[i]] <- data.frame(
    file = meta_s$file_base[i],
    cytokine = if ("cytokine" %in% names(meta_s)) meta_s$cytokine[i] else meta_s$stim[i],
    time_min = meta_s$time_min[i]
  )[rep(1, take), ]
}

expr <- do.call(rbind, expr_list)
meta_dr <- do.call(rbind, meta_list)

################################################################################
# Map FlowSOM clusters
################################################################################

codes <- fs$map$codes
nn <- FNN::get.knnx(codes, expr, k=1)$nn.index
meta_dr$cluster <- fs$metaclustering[nn]

labs <- fread(opt$`cluster-labels`)
labs[, cluster := factor(cluster)]
meta_dr$cluster <- factor(meta_dr$cluster)
meta_dr <- merge(meta_dr, labs, by="cluster", all.x=TRUE)

################################################################################
# Run tSNE
################################################################################

ts <- Rtsne(expr, perplexity=30, check_duplicates=FALSE)
meta_dr$TSNE1 <- ts$Y[,1]
meta_dr$TSNE2 <- ts$Y[,2]

################################################################################
# Run UMAP
################################################################################

um <- umap(expr, n_neighbors=30, min_dist=0.3)
meta_dr$UMAP1 <- um[,1]
meta_dr$UMAP2 <- um[,2]

################################################################################
# Save plots
################################################################################

plot_save <- function(p, fname) {
  ggsave(file.path(opt$outdir, fname), p, width=7, height=6)
}

plot_save(
  ggplot(meta_dr, aes(TSNE1, TSNE2, col=cluster)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "tsne_by_cluster.png"
)

plot_save(
  ggplot(meta_dr, aes(UMAP1, UMAP2, col=cluster)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "umap_by_cluster.png"
)

plot_save(
  ggplot(meta_dr, aes(TSNE1, TSNE2, col=label)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "tsne_by_identity.png"
)

plot_save(
  ggplot(meta_dr, aes(UMAP1, UMAP2, col=label)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "umap_by_identity.png"
)

plot_save(
  ggplot(meta_dr, aes(TSNE1, TSNE2, col=factor(time_min))) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "tsne_by_time.png"
)

plot_save(
  ggplot(meta_dr, aes(UMAP1, UMAP2, col=factor(time_min))) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "umap_by_time.png"
)

plot_save(
  ggplot(meta_dr, aes(TSNE1, TSNE2, col=cytokine)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "tsne_by_cytokine.png"
)

plot_save(
  ggplot(meta_dr, aes(UMAP1, UMAP2, col=cytokine)) +
    geom_point(size=0.3, alpha=0.7) + theme_bw(),
  "umap_by_cytokine.png"
)

saveRDS(meta_dr, file.path(opt$outdir,"dimred_metadata.rds"))

cat("Done: tSNE + UMAP plots\n")
