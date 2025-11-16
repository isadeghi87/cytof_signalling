#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(data.table)
  library(dplyr)
  library(Rtsne)
  library(uwot)
  library(ggplot2)
  library(FNN)
})

opt_list <- list(
  make_option("--norm-dir", type="character", default="normalized_fcs"),
  make_option("--sample-table", type="character", default="sample_table.csv"),
  make_option("--flowsom-rds", type="character",
              default="flowsom_results/flowsom_fitted.rds"),
  make_option("--cluster-labels", type="character",
              default="flowsom_results/cluster_labels_template.csv"),
  make_option("--outdir", type="character", default="dimred_flowsom"),
  make_option("--n-events", type="integer", default=60000),
  make_option("--cofactor", type="double", default=5),
  make_option("--seed", type="integer", default=1)
)

opt <- parse_args(OptionParser(option_list=opt_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
# Load FlowSOM and metadata
###############################################################################

fs <- readRDS(opt$`flowsom-rds`)

meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) meta[, file_base := basename(file)]
}

meta_s <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_s$file_base))
idx <- file.exists(norm_files)
meta_s <- meta_s[idx]
norm_files <- norm_files[idx]

###############################################################################
# Identify lineage channels
###############################################################################

ff <- read.FCS(norm_files[1], transformation=FALSE)
par_tab <- pData(parameters(ff))
desc <- par_tab$desc
cn <- par_tab$name

lineage_markers <- c(
  "Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
  "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di","Cd113Di","Gd157Di"
)

lineage_ch <- intersect(lineage_markers, cn)

tf <- arcsinhTransform(a=0, b=1/opt$cofactor, c=0)
tlist <- transformList(lineage_ch, tf)

###############################################################################
# Read + downsample for DR
###############################################################################

expr_list <- list()
meta_list <- list()

per_file <- round(opt$`n-events` / length(norm_files))

for (i in seq_along(norm_files)) {

  ff0 <- read.FCS(norm_files[i], transformation=FALSE)
  ffT <- transform(ff0, tlist)
  m <- exprs(ffT)[, lineage_ch, drop=FALSE]

  n <- nrow(m)
  take <- min(n, per_file)
  idx0 <- sample(seq_len(n), take)

  expr_list[[i]] <- m[idx0, ]

  meta_list[[i]] <- data.frame(
    file = meta_s$file_base[i],
    time_min = meta_s$time_min[i],
    cytokine = if ("cytokine" %in% names(meta_s)) meta_s$cytokine[i] else meta_s$stim[i],
    stringsAsFactors=FALSE
  )[rep(1, take), ]
}

expr_mat <- do.call(rbind, expr_list)
meta_dr <- do.call(rbind, meta_list)

###############################################################################
# Map FlowSOM clusters
###############################################################################

codes <- fs$map$codes
nn <- FNN::get.knnx(codes, expr_mat, k=1)$nn.index
meta_dr$cluster <- fs$metaclustering[nn]

# Load biological labels
lab <- fread(opt$`cluster-labels`)
lab[, cluster := factor(cluster)]
meta_dr$cluster <- factor(meta_dr$cluster)
meta_dr <- left_join(meta_dr, lab, by=c("cluster"="cluster"))

###############################################################################
# tSNE
###############################################################################

ts <- Rtsne(expr_mat, perplexity=30, pca=TRUE, check_duplicates=FALSE)
meta_dr$TSNE1 <- ts$Y[,1]
meta_dr$TSNE2 <- ts$Y[,2]

###############################################################################
# UMAP
###############################################################################

um <- umap(expr_mat, n_neighbors=30, min_dist=0.3)
meta_dr$UMAP1 <- um[,1]
meta_dr$UMAP2 <- um[,2]

###############################################################################
# Plots
###############################################################################

ggsave(file.path(opt$outdir, "tsne_by_cluster.png"),
       ggplot(meta_dr, aes(TSNE1, TSNE2, color=cluster)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "tsne_by_identity.png"),
       ggplot(meta_dr, aes(TSNE1, TSNE2, color=label)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "tsne_by_time.png"),
       ggplot(meta_dr, aes(TSNE1, TSNE2, color=factor(time_min))) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "tsne_by_cytokine.png"),
       ggplot(meta_dr, aes(TSNE1, TSNE2, color=cytokine)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "umap_by_cluster.png"),
       ggplot(meta_dr, aes(UMAP1, UMAP2, color=cluster)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "umap_by_identity.png"),
       ggplot(meta_dr, aes(UMAP1, UMAP2, color=label)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "umap_by_time.png"),
       ggplot(meta_dr, aes(UMAP1, UMAP2, color=factor(time_min))) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

ggsave(file.path(opt$outdir, "umap_by_cytokine.png"),
       ggplot(meta_dr, aes(UMAP1, UMAP2, color=cytokine)) +
       geom_point(size=0.3, alpha=0.7) + theme_bw(), width=7, height=6)

saveRDS(meta_dr, file.path(opt$outdir, "dimred_metadata.rds"))

cat("\nDone DR.\n")
