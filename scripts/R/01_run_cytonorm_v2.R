#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(flowCore)
  library(CytoNorm)
  library(ggplot2)
  library(Rtsne)
  library(uwot)
})

option_list <- list(
  make_option("--sample-table", type = "character", default = "sample_table.csv",
              help = "Metadata table with columns: file_base, batch, is_anchor, sample_id, stim"),
  make_option("--files-list", type = "character", default = "files.txt",
              help = "Text file with one sample FCS filename per line (no anchors)"),
  make_option("--base-dir", type = "character", default = "batch",
              help = "Base directory that contains sample FCS files and 'anchors' subdir"),
  make_option("--channels", type = "character", default = NA,
              help = "Optional TSV with column 'channel' listing channels to normalize/check"),
  make_option("--outdir", type = "character", default = "normalized_fcs",
              help = "Output directory for normalized FCS"),
  make_option("--qcdir", type = "character", default = "qc_plots",
              help = "Output directory for QC plots"),
  make_option("--n-q", type = "integer", default = 99,
              help = "Number of quantiles for QuantileNorm.train (default 99)"),
  make_option("--n-events-qc", type = "integer", default = 30000,
              help = "Number of events sampled per condition for QC plots"),
  make_option("--seed", type = "integer", default = 1,
              help = "Random seed")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# define CyTOF transforms (asinh with cofactor 5 and its inverse)
cytof_cofactor <- 5
cytofTransform <- flowCore::arcsinhTransform(a = 0, b = 1 / cytof_cofactor, c = 0)
cytofTransform.reverse <- function(x) sinh(x) * cytof_cofactor

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$qcdir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== Step 1: Read metadata ===\n")

meta <- fread(opt$`sample-table`)

# Support both legacy (file_base) and current (file) sample_table schemas
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) {
    meta[, file_base := basename(file)]
  } else {
    stop("sample_table must contain either a 'file_base' or 'file' column.")
  }
}
required_cols <- c("file_base", "batch", "is_anchor")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("sample_table is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

meta[, batch := as.factor(batch)]
meta[, is_anchor := as.logical(is_anchor)]

cat("  - n rows in metadata:", nrow(meta), "\n")
cat("  - batches:", paste(levels(meta$batch), collapse = ", "), "\n")

# ------------------------------------------------------------
# Determine sample files from files.txt
# ------------------------------------------------------------
cat("\n=== Step 2: Determine sample files ===\n")

samples_dt <- meta[is_anchor == FALSE]

if (!is.null(opt$`files-list`) && !is.na(opt$`files-list`) && file.exists(opt$`files-list`)) {
  cat("  - Filtering samples using files list:", opt$`files-list`, "\n")
  files_vec <- readLines(opt$`files-list`, warn = FALSE)
  files_vec <- trimws(files_vec)
  files_vec <- files_vec[nchar(files_vec) > 0]
  files_base <- basename(files_vec)
  cat("  - n entries in files.txt:", length(files_base), "\n")
  samples_dt <- samples_dt[file_base %in% files_base]
}

cat("  - n sample files used for normalization:", nrow(samples_dt), "\n")

if (nrow(samples_dt) == 0) {
  stop("No non-anchor sample rows available for normalization (after optional files.txt filtering).")
}

# anchors: from sample_table; their FCS actually live in base-dir/anchors
anchors_dt <- meta[is_anchor == TRUE]

cat("  - n anchor files in metadata:", nrow(anchors_dt), "\n")

# Build full paths
anchors_dt[, file_full := file.path(opt$`base-dir`, "anchors", file_base)]
samples_dt[, file_full := file.path(opt$`base-dir`, file_base)]

# sanity check
chk_anchors <- anchors_dt[!file.exists(file_full)]
chk_samples <- samples_dt[!file.exists(file_full)]

if (nrow(chk_anchors) > 0) {
  cat("[WARN] Some anchor files not found on disk (showing first 5):\n")
  print(head(chk_anchors[, .(file_base, file_full)], 5))
}
if (nrow(chk_samples) > 0) {
  cat("[WARN] Some sample files not found on disk (showing first 5):\n")
  print(head(chk_samples[, .(file_base, file_full)], 5))
}

anchors_dt <- anchors_dt[file.exists(file_full)]
samples_dt <- samples_dt[file.exists(file_full)]

if (nrow(anchors_dt) == 0) stop("No existing anchor FCS files found.")
if (nrow(samples_dt) == 0) stop("No existing sample FCS files from files.txt found.")

# ------------------------------------------------------------
# Determine channels to normalize
# ------------------------------------------------------------
cat("\n=== Step 3: Determine channels to normalize ===\n")
### === Step 3: Determine channels to normalize (drop DNA / beads) ===

ff_example <- read.FCS(anchors_dt$file_full[1],
                       transformation = FALSE,
                       truncate_max_range = FALSE)
all_channels <- colnames(ff_example)

# things we NEVER want to normalize on
remove_patterns <- c("Ir191", "Ir193", "Ce140", "Bead", "EQ", "DNA",
                     "Time", "Event_length", "length")

is_bad <- Reduce(`|`, lapply(remove_patterns, function(p) {
  grepl(p, all_channels, ignore.case = TRUE)
}))

candidate_channels <- all_channels[!is_bad]

if (is.na(opt$channels)) {
  # Default: metal channels ending in "Di"
  norm_channels <- grep("Di$", candidate_channels, value = TRUE)
  cat("  - Using channels ending in 'Di' (technical removed); n =",
      length(norm_channels), "\n")
} else {
  ch_tab <- fread(opt$channels)
  if (!"channel" %in% colnames(ch_tab)) {
    stop("channels TSV must have a column named 'channel'")
  }
  ch_vec <- ch_tab$channel
  ch_vec <- ch_vec[ch_vec != "" & !is.na(ch_vec) & ch_vec != "channel"]
  norm_channels <- intersect(ch_vec, candidate_channels)
  cat("  - Using channels from file:", opt$channels, "\n")
  cat("  - n channels found in FCS:", length(norm_channels), "\n")
}

if (length(norm_channels) == 0) stop("No normalization channels found.")

# transforms (unchanged)
transformList         <- flowCore::transformList(norm_channels, cytofTransform)
transformList.reverse <- flowCore::transformList(norm_channels, cytofTransform.reverse)

# ------------------------------------------------------------
# Train QuantileNorm model on anchors
# ------------------------------------------------------------
cat("\n=== Step 4: Train QuantileNorm model on anchors ===\n")

model_qn <- QuantileNorm.train(
  files = anchors_dt$file_full,
  labels = anchors_dt$batch,
  channels = norm_channels,
  transformList = transformList,
  nQ = opt$`n-q`
)

cat("  - Model trained with", length(unique(anchors_dt$batch)), "batches\n")

# ------------------------------------------------------------
# Apply normalization to anchors + samples
# ------------------------------------------------------------
cat("\n=== Step 5: Normalize anchors + samples ===\n")

all_files <- c(anchors_dt$file_full, samples_dt$file_full)
all_labels <- c(as.character(anchors_dt$batch), as.character(samples_dt$batch))

QuantileNorm.normalize(
  model = model_qn,
  files = all_files,
  labels = all_labels,
  transformList = transformList,
  transformList.reverse = transformList.reverse,
  outputDir = opt$outdir,
  prefix = "Norm_",
  verbose = TRUE
)

cat("  - Normalized FCS files written to:", opt$outdir, "\n")

# ------------------------------------------------------------
# QC: density overlays raw vs normalized
# ------------------------------------------------------------
cat("\n=== Step 6: QC density overlays (raw vs norm) ===\n")

# Pick a small subset of files across batches for QC (up to 3 per batch)
qc_sample_files <- samples_dt[, .SD[sample(.N, min(3L, .N))], by = batch]

cat("  - Using", nrow(qc_sample_files), "sample files for density QC\n")

density_qc <- function(channel, n_events = opt$`n-events-qc`) {

  df_list <- list()

  for (i in seq_len(nrow(qc_sample_files))) {
    row <- qc_sample_files[i, ]
    raw_path  <- row$file_full
    norm_path <- file.path(opt$outdir, paste0("Norm_", basename(raw_path)))

    if (!file.exists(raw_path) || !file.exists(norm_path)) next

    ff_raw  <- read.FCS(raw_path, transformation = FALSE, truncate_max_range = FALSE)
    ff_norm <- read.FCS(norm_path, transformation = FALSE, truncate_max_range = FALSE)

    if (!(channel %in% colnames(ff_raw))) next

    # Apply cytof transform for plotting (asinh/5)
    ff_raw_t  <- flowCore::transform(ff_raw, transformList)
    ff_norm_t <- flowCore::transform(ff_norm, transformList)

    vals_raw  <- exprs(ff_raw_t)[, channel]
    vals_norm <- exprs(ff_norm_t)[, channel]

    # downsample
    n_take <- min(length(vals_raw), n_events)
    idx_r  <- sample(seq_along(vals_raw), n_take)
    idx_n  <- sample(seq_along(vals_norm), n_take)

    df_list[[length(df_list) + 1]] <-
      data.frame(value = c(vals_raw[idx_r], vals_norm[idx_n]),
                 source = rep(c("raw", "norm"), each = n_take),
                 file = basename(raw_path),
                 stringsAsFactors = FALSE)
  }

  if (length(df_list) == 0) return(NULL)
  df <- bind_rows(df_list)

  p <- ggplot(df, aes(x = value, colour = source)) +
    geom_density() +
    labs(
      title = paste0("Density overlay: ", channel, " (asinh/5)"),
      x = paste0(channel, " (asinh/5)"),
      y = "Density"
    ) +
    theme_bw(base_size = 14)

  ggsave(file.path(opt$qcdir, paste0("density_", channel, ".png")),
         p, width = 8, height = 5)
}

# Limit number of channels for QC to keep runtime reasonable
qc_channels <- norm_channels[1:min(8, length(norm_channels))]
cat("  - Plotting density overlays for channels:\n    ",
    paste(qc_channels, collapse = ", "), "\n")

invisible(lapply(qc_channels, density_qc))

# ------------------------------------------------------------
# QC: t-SNE & UMAP on normalized data
# ------------------------------------------------------------
cat("\n=== Step 7: Dimensionality reduction & clustering on normalized data ===\n")

# Use a subset of normalized sample files for DR (up to 3 per batch)
dr_files <- samples_dt[, .SD[sample(.N, min(3L, .N))], by = batch]

cat("  - Using", nrow(dr_files), "normalized sample files for TSNE/UMAP\n")

expr_list <- list()
meta_list <- list()

for (i in seq_len(nrow(dr_files))) {
  row <- dr_files[i, ]
  norm_path <- file.path(opt$outdir, paste0("Norm_", basename(row$file_full)))
  if (!file.exists(norm_path)) next

  ff <- read.FCS(norm_path, transformation = FALSE, truncate_max_range = FALSE)
  ff_t <- flowCore::transform(ff, transformList)

  mat <- exprs(ff_t)[, norm_channels, drop = FALSE]

  # downsample per file for DR
  n_take <- min(nrow(mat), opt$`n-events-qc`)
  idx <- sample(seq_len(nrow(mat)), n_take)

  expr_list[[length(expr_list) + 1]] <- mat[idx, ]
  meta_list[[length(meta_list) + 1]] <- data.frame(
    file      = basename(row$file_full),
    batch     = row$batch,
    time_min  = if ("time_min" %in% names(row)) row$time_min else NA,
    stringsAsFactors = FALSE
  )[rep(1, n_take), ]
}

if (length(expr_list) > 0) {
  expr_mat <- do.call(rbind, expr_list)
  meta_dr  <- do.call(rbind, meta_list)

  # simple k-means clustering on normalized expression
  k <- 12
  km <- kmeans(expr_mat, centers = k, nstart = 10)
  meta_dr$cluster <- factor(km$cluster)

  cat("  - Running t-SNE...\n")
  tsne_res <- Rtsne(expr_mat, check_duplicates = FALSE, pca = TRUE, perplexity = 30)
  meta_dr$TSNE1 <- tsne_res$Y[,1]
  meta_dr$TSNE2 <- tsne_res$Y[,2]

  cat("  - Running UMAP...\n")
  umap_res <- umap(expr_mat, n_neighbors = 30, min_dist = 0.3, metric = "euclidean")
  meta_dr$UMAP1 <- umap_res[,1]
  meta_dr$UMAP2 <- umap_res[,2]

  # t-SNE plot
  p_tsne <- ggplot(meta_dr, aes(x = TSNE1, y = TSNE2, colour = cluster)) +
    geom_point(size = 0.4, alpha = 0.8) +
    theme_bw(base_size = 14) +
    ggtitle("t-SNE colored by cluster")

  ggsave(file.path(opt$qcdir, "tsne_clusters.png"),
         p_tsne, width = 8, height = 6)

  # UMAP colored by sample
  p_umap_sample <- ggplot(meta_dr, aes(x = UMAP1, y = UMAP2, colour = file)) +
    geom_point(size = 0.4, alpha = 0.8) +
    theme_bw(base_size = 14) +
    ggtitle("UMAP colored by sample")

  ggsave(file.path(opt$qcdir, "umap_by_sample.png"),
         p_umap_sample, width = 10, height = 6)

  # UMAP colored by cluster
  p_umap_cluster <- ggplot(meta_dr, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 0.4, alpha = 0.8) +
    theme_bw(base_size = 14) +
    ggtitle("UMAP colored by cluster")

  ggsave(file.path(opt$qcdir, "umap_clusters.png"),
         p_umap_cluster, width = 8, height = 6)

  # t-SNE colored by batch
  p_tsne_batch <- ggplot(meta_dr, aes(x = TSNE1, y = TSNE2, colour = batch)) +
    geom_point(size = 0.4, alpha = 0.8) +
    theme_bw(base_size = 14) +
    ggtitle("t-SNE colored by batch")

  ggsave(file.path(opt$qcdir, "tsne_by_batch.png"),
         p_tsne_batch, width = 8, height = 6)

  # UMAP colored by batch
  p_umap_batch <- ggplot(meta_dr, aes(x = UMAP1, y = UMAP2, colour = batch)) +
    geom_point(size = 0.4, alpha = 0.8) +
    theme_bw(base_size = 14) +
    ggtitle("UMAP colored by batch")

  ggsave(file.path(opt$qcdir, "umap_by_batch.png"),
         p_umap_batch, width = 8, height = 6)

  # t-SNE and UMAP colored by time point (if available)
  if ("time_min" %in% colnames(meta_dr)) {
    meta_dr$time_factor <- factor(meta_dr$time_min)

    p_tsne_time <- ggplot(meta_dr, aes(x = TSNE1, y = TSNE2, colour = time_factor)) +
      geom_point(size = 0.4, alpha = 0.8) +
      theme_bw(base_size = 14) +
      ggtitle("t-SNE colored by time point")

    ggsave(file.path(opt$qcdir, "tsne_by_timepoint.png"),
           p_tsne_time, width = 8, height = 6)

    p_umap_time <- ggplot(meta_dr, aes(x = UMAP1, y = UMAP2, colour = time_factor)) +
      geom_point(size = 0.4, alpha = 0.8) +
      theme_bw(base_size = 14) +
      ggtitle("UMAP colored by time point")

    ggsave(file.path(opt$qcdir, "umap_by_timepoint.png"),
           p_umap_time, width = 8, height = 6)
  }

  cat("  - Saved TSNE/UMAP plots in:", opt$qcdir, "\n")
} else {
  cat("[WARN] No data for DR – check normalized files.\n")
}



cat("\n=== DONE ===\n")
