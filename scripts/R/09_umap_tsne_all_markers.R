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
  make_option("--sample-table", type = "character", default = "sample_table.csv"),
  make_option("--norm-dir", type = "character", default = "normalized_fcs"),
  make_option("--flowsom-rds", type = "character",
              default = "flowsom_results/flowsom_fitted.rds"),
make_option("--cluster-labels", type = "character",
            default = "flowsom_results/cluster_labels_template.csv"),
make_option("--outdir", type = "character", default = "dimred_all_markers"),
make_option("--cofactor", type = "double", default = 5),
make_option("--per-file", type = "integer", default = -1,
            help = "Events per file (-1 = use all; beware memory)"),
  make_option("--pca-dims", type = "integer", default = 50,
              help = "PCA dimensions before t-SNE/UMAP (0 = no PCA)"),
  make_option("--perplexity", type = "double", default = 30),
  make_option("--umap-nn", type = "integer", default = 30),
  make_option("--umap-min-dist", type = "double", default = 0.3),
  make_option("--umap-threads", type = "integer", default = 0,
              help = "Threads for UMAP (0 = auto)"),
  make_option("--tsne-max-cells", type = "integer", default = 200000,
              help = "Max cells for t-SNE (downsample if exceeded; 0 = use all, may be very slow)"),
  make_option("--seed", type = "integer", default = 1)
)

opt <- parse_args(OptionParser(option_list = opt_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Load metadata and file paths (anchors + samples)
###############################################################################

meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) meta[, file_base := basename(file)] else
    stop("sample_table must contain file_base or file")
}
if (!"batch" %in% names(meta)) stop("sample_table missing batch")

meta[, batch := factor(batch)]
meta[, is_anchor := as.logical(is_anchor)]
if (!"time_min" %in% names(meta)) meta[, time_min := NA_integer_]
if (!"cytokine" %in% names(meta)) {
  if ("stim" %in% names(meta)) meta[, cytokine := stim] else meta[, cytokine := NA_character_]
}

norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta$file_base))
meta[, norm_file := norm_files]
meta <- meta[file.exists(norm_file)]
if (nrow(meta) == 0) stop("No normalized files found")

###############################################################################
# Determine marker set (all non-technical channels)
###############################################################################

ff0 <- read.FCS(meta$norm_file[1], transformation = FALSE, truncate_max_range = FALSE)
all_channels <- colnames(exprs(ff0))
remove_patterns <- c("Ir191", "Ir193", "Ce140", "Bead", "EQ", "DNA",
                     "Time", "Event_length", "length")
is_bad <- Reduce(`|`, lapply(remove_patterns, function(p) grepl(p, all_channels, ignore.case = TRUE)))
marker_cols <- all_channels[!is_bad]
if (!length(marker_cols)) stop("No marker channels remaining after filtering technical columns")
cat("Markers used (", length(marker_cols), "): ", paste(marker_cols, collapse = ", "), "\n", sep = "")

cofactor <- opt$cofactor

###############################################################################
# Collect expression matrix
###############################################################################

expr_list <- vector("list", nrow(meta))
meta_list <- vector("list", nrow(meta))
per_file <- opt$`per-file`
cat("Sampling per file:", if (per_file < 0) "all events" else per_file, "\n")

for (i in seq_len(nrow(meta))) {
  ff <- read.FCS(meta$norm_file[i], transformation = FALSE, truncate_max_range = FALSE,
                 which.lines = if (per_file > 0) seq_len(per_file) else NULL)
  mat_raw <- exprs(ff)[, marker_cols, drop = FALSE]

  # Enforce cap if requested (keeps first N rows for speed/memory)
  if (per_file > 0 && nrow(mat_raw) > per_file) {
    mat_raw <- mat_raw[seq_len(per_file), , drop = FALSE]
  }

  # arcsinh transform on the subset only
  mat <- asinh(mat_raw / cofactor)
  n_take <- nrow(mat)

  expr_list[[i]] <- mat
  meta_list[[i]] <- data.table(
    file_base = meta$file_base[i],
    batch = meta$batch[i],
    time_min = meta$time_min[i],
    cytokine = meta$cytokine[i],
    is_anchor = meta$is_anchor[i]
  )[rep(1, n_take)]
}

expr_mat <- do.call(rbind, expr_list)
meta_dr <- rbindlist(meta_list)
cat("Total cells used:", nrow(expr_mat), "\n")

###############################################################################
# Map FlowSOM clusters (if available)
###############################################################################

meta_dr[, cluster := NA_integer_]
meta_dr[, label := NA_character_]

if (file.exists(opt$`flowsom-rds`)) {
  fs <- readRDS(opt$`flowsom-rds`)
  lineage_markers <- c(
    "Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
    "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di",
    "Cd113Di","Gd157Di"
  )
  lineage_ch <- intersect(lineage_markers, marker_cols)
  if (length(lineage_ch) == ncol(fs$map$codes)) {
    nn <- FNN::get.knnx(fs$map$codes, expr_mat[, lineage_ch, drop = FALSE], k = 1)$nn.index
    meta_dr$cluster <- fs$metaclustering[nn]
    labs <- if (file.exists(opt$`cluster-labels`)) fread(opt$`cluster-labels`) else NULL
    if (!is.null(labs) && all(c("cluster", "label") %in% names(labs))) {
      labs[, cluster := as.integer(cluster)]
      meta_dr <- merge(meta_dr, labs, by = "cluster", all.x = TRUE)
      if ("label.x" %in% names(meta_dr) && "label.y" %in% names(meta_dr)) {
        meta_dr[, label := fifelse(!is.na(label.y), label.y, label.x)]
        meta_dr[, c("label.x", "label.y") := NULL]
      }
    }
  }
}

###############################################################################
# Dimensionality reduction on all markers
###############################################################################

X_full <- scale(expr_mat)

if (opt$`pca-dims` > 0 && ncol(X_full) > opt$`pca-dims`) {
  cat("Running PCA to", opt$`pca-dims`, "dims for DR input\n")
  pca <- prcomp(X_full, center = FALSE, scale. = FALSE, rank. = opt$`pca-dims`)
  X <- pca$x
} else {
  X <- X_full
}

um <- umap(X, n_neighbors = opt$`umap-nn`, min_dist = opt$`umap-min-dist`,
           n_threads = opt$`umap-threads`)
meta_dr$UMAP1 <- um[, 1]
meta_dr$UMAP2 <- um[, 2]

###############################################################################
# t-SNE (optionally downsample to keep runtime reasonable)
###############################################################################

ts_input <- X
ts_meta <- meta_dr
if (opt$`tsne-max-cells` > 0 && nrow(ts_input) > opt$`tsne-max-cells`) {
  set.seed(opt$seed)
  idx_ts <- sample.int(nrow(ts_input), opt$`tsne-max-cells`)
  ts_input <- ts_input[idx_ts, , drop = FALSE]
  ts_meta <- ts_meta[idx_ts, ]
  cat("t-SNE downsampled to", nrow(ts_input), "cells (limit =", opt$`tsne-max-cells`, ")\n")
} else {
  cat("t-SNE using all", nrow(ts_input), "cells\n")
}

ts <- Rtsne(ts_input, perplexity = opt$perplexity, check_duplicates = FALSE)
ts_meta$TSNE1 <- ts$Y[, 1]
ts_meta$TSNE2 <- ts$Y[, 2]

if (nrow(ts_meta) < nrow(meta_dr)) {
  # merge back for plotting (leave NA for non-tsne cells)
  meta_dr$TSNE1 <- NA_real_
  meta_dr$TSNE2 <- NA_real_
  meta_dr[seq_len(nrow(ts_meta)), c("TSNE1", "TSNE2") := ts_meta[, .(TSNE1, TSNE2)]]
} else {
  meta_dr$TSNE1 <- ts_meta$TSNE1
  meta_dr$TSNE2 <- ts_meta$TSNE2
}

fwrite(meta_dr, file.path(opt$outdir, "umap_tsne_all_markers.csv"))
saveRDS(meta_dr, file.path(opt$outdir, "umap_tsne_all_markers.rds"))

###############################################################################
# Plot helpers
###############################################################################

make_palette <- function(vals) {
  lv <- sort(unique(na.omit(vals)))
  if (!length(lv)) return(NULL)
  pal <- scales::hue_pal()(max(3, length(lv)))
  names(pal) <- as.character(lv)
  pal
}

plot_scatter <- function(df, x, y, color, title, fname, use_palette = FALSE) {
  if (use_palette) {
    df <- data.table::as.data.table(df)  # ensure copy
    df[[color]] <- factor(df[[color]])
  }
  p <- ggplot(df, aes_string(x = x, y = y, color = color)) +
    geom_point(size = 0.25, alpha = 0.5) +
    theme_bw() +
    labs(title = title, color = color)
  if (use_palette) {
    pal <- make_palette(df[[color]])
    if (!is.null(pal)) p <- p + scale_color_manual(values = pal)
  }
  ggsave(file.path(opt$outdir, fname), p, width = 8, height = 6)
}

plot_scatter(meta_dr, "UMAP1", "UMAP2", "batch", "UMAP by batch", "umap_by_batch.png")
plot_scatter(meta_dr, "UMAP1", "UMAP2", "time_min", "UMAP by time", "umap_by_time.png")
plot_scatter(meta_dr, "UMAP1", "UMAP2", "file_base", "UMAP by file", "umap_by_file.png")
plot_scatter(meta_dr, "UMAP1", "UMAP2", "cluster", "UMAP by cluster", "umap_by_cluster.png", use_palette = TRUE)
plot_scatter(meta_dr, "UMAP1", "UMAP2", "label", "UMAP by cell type/label", "umap_by_label.png", use_palette = TRUE)
plot_scatter(meta_dr, "UMAP1", "UMAP2", "cytokine", "UMAP by cytokine", "umap_by_cytokine.png")
plot_scatter(meta_dr, "UMAP1", "UMAP2", "is_anchor", "UMAP anchors vs samples", "umap_by_anchor.png")

plot_scatter(meta_dr, "TSNE1", "TSNE2", "batch", "t-SNE by batch", "tsne_by_batch.png")
plot_scatter(meta_dr, "TSNE1", "TSNE2", "time_min", "t-SNE by time", "tsne_by_time.png")
plot_scatter(meta_dr, "TSNE1", "TSNE2", "file_base", "t-SNE by file", "tsne_by_file.png")
plot_scatter(meta_dr, "TSNE1", "TSNE2", "cluster", "t-SNE by cluster", "tsne_by_cluster.png", use_palette = TRUE)
plot_scatter(meta_dr, "TSNE1", "TSNE2", "label", "t-SNE by cell type/label", "tsne_by_label.png", use_palette = TRUE)
plot_scatter(meta_dr, "TSNE1", "TSNE2", "cytokine", "t-SNE by cytokine", "tsne_by_cytokine.png")
plot_scatter(meta_dr, "TSNE1", "TSNE2", "is_anchor", "t-SNE anchors vs samples", "tsne_by_anchor.png")

cat("Saved UMAP/t-SNE embeddings and plots to", opt$outdir, "\n")
