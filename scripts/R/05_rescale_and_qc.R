#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(purrr); library(stringr)
  library(flowCore);  library(matrixStats); library(stats)
  library(uwot);      library(cluster);    library(ggplot2); library(pheatmap)
})

# ---------- Inputs ----------
sample_csv   <- "sample_table.csv"              # from your prepare_metadata.R
norm_dir_in  <- "normalized_fcs"                # post-CytoNorm FCS
out_dir_fcs  <- "normalized_fcs_rescaled"       # output (after channel-specific re-scaling)
qc_dir       <- "qc"
cofactor     <- 5
set.seed(1)

dir.create(out_dir_fcs, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

tab <- fread(sample_csv)

# Resolve normalized file path for each input (handles CytoNorm 'Norm__' and 'Norm_Batch_*' cases)
norm_files <- list.files(norm_dir_in, full.names = TRUE)
resolve_norm_for <- function(bn) {
  # direct match on basename
  cand <- norm_files[basename(norm_files) == bn]
  if (length(cand)) return(cand[1])
  # contains-match
  cand <- norm_files[grepl(bn, basename(norm_files), fixed = TRUE)]
  if (length(cand)) return(cand[which.max(file.info(cand)$size)])
  # relaxed token match
  toks <- unlist(strsplit(tolower(bn), "[^a-z0-9]+")); toks <- toks[toks != ""]
  if (length(toks)) {
    ok <- vapply(basename(norm_files), function(x){ any(toks %in% unlist(strsplit(tolower(x), "[^a-z0-9]+"))) }, logical(1))
    cand <- norm_files[ok]
    if (length(cand)) return(cand[which.max(file.info(cand)$size)])
  }
  NA_character_
}

tab[, norm_file := vapply(basename(file), resolve_norm_for, character(1))]
missing <- tab[is.na(norm_file) | !file.exists(norm_file)]
if (nrow(missing)) {
  stop(sprintf("%d normalized files missing. e.g.: %s", nrow(missing), paste(head(basename(missing$file), 3), collapse=", ")))
}

lineage <- if (file.exists("effective_lineage_markers.txt"))
  readLines("effective_lineage_markers.txt") else character(0)
phospho <- if (file.exists("effective_phospho_markers.txt"))
  readLines("effective_phospho_markers.txt") else character(0)
markers_use <- unique(c(lineage, phospho))

# Helper to read a channel vector (already arcsinh in your pipeline)
read_ch <- function(f, ch){
  M <- exprs(read.FCS(f, transformation = FALSE))
  if (!ch %in% colnames(M)) stop("Missing channel: ", ch, " in ", f)
  M[, ch]
}

norm_path <- function(x) x
out_path  <- function(x) file.path(out_dir_fcs, basename(x))

# ---------- 1) Diagnostics (pre-rescaling) ----------
message("Computing per-file medians & 99th percentiles (pre-rescale)...")

# Compute per-file statistics for all channels we might use
file_stats <- function(files, chans){
  res <- lapply(files, function(f){
    M <- exprs(read.FCS(f, transformation=FALSE))
    keep <- intersect(chans, colnames(M))
    X <- M[, keep, drop=FALSE]
    q99 <- matrixStats::colQuantiles(X, probs = 0.99, na.rm=TRUE, drop=FALSE)[, 1]
    data.frame(file = basename(f),
               channel = keep,
               median = colMedians(X, na.rm=TRUE),
               p99    = q99,
               stringsAsFactors = FALSE)
  })
  bind_rows(res)
}

pre_stats <- file_stats(tab$norm_file, markers_use)
# normalize 'file' key to original basename so joins with post_stats work
file_map <- data.table(file_pre = basename(tab$norm_file), file_target = basename(tab$file))
pre_stats$file <- file_map$file_target[match(pre_stats$file, file_map$file_pre)]
fwrite(pre_stats, file.path(qc_dir, "pre_file_stats.csv"))

# ANOVA (batch effect on medians) per channel
message("Running per-channel ANOVA on medians (pre-rescale)...")
anova_pre <- pre_stats %>%
  left_join(tab[, .(file = basename(file), batch)], by="file") %>%
  group_by(channel) %>%
  summarise(
    p_value = tryCatch({
      fit <- aov(median ~ factor(batch))
      summary(fit)[[1]][["Pr(>F)"]][1]
    }, error=function(e) NA_real_),
    n_batches = n_distinct(batch)
  ) %>%
  mutate(p_adj = p.adjust(p_value, method="fdr"))
fwrite(anova_pre, file.path(qc_dir, "anova_pre.csv"))

# 99th-percentile range across files per channel
p99_pre <- pre_stats %>%
  group_by(channel) %>%
  summarise(p99_range = max(p99, na.rm=TRUE) - min(p99, na.rm=TRUE))
fwrite(p99_pre, file.path(qc_dir, "p99_range_pre.csv"))

# ---------- 2) Select channels needing re-scale ----------
# Criteria: significant batch effect OR large p99 spread among *phospho* markers
bad_ch <- anova_pre %>% dplyr::filter(!is.na(p_adj) & p_adj < 0.05) %>% pull(channel)
spread_ch <- p99_pre %>% arrange(desc(p99_range)) %>% dplyr::filter(p99_range > quantile(p99_range, 0.8, na.rm=TRUE)) %>% pull(channel)
cand_rescale <- intersect(unique(c(bad_ch, spread_ch)), phospho)  # phosphos only

message("Channels flagged for re-scaling: ", paste(cand_rescale, collapse=", "))

# ---------- 3) Learn per-batch linear maps from UNSTIM anchors ----------
is_unstim <- !is.na(tab$stim) & grepl("^unstim$", tab$stim, ignore.case=TRUE)
anchors   <- tab[is_anchor == TRUE & is_unstim == TRUE]

# Reference (global) anchor percentiles per channel
ref_quants <- function(ch, probs=c(0.05, 0.95)){
  vals <- unlist(lapply(anchors$norm_file, function(f) read_ch(f, ch)))
  stats::quantile(vals, probs=probs, na.rm=TRUE, names=FALSE)
}

# Batch-specific anchor percentiles per channel
batch_quants <- function(batch_id, ch, probs=c(0.05,0.95)){
  fs <- anchors[batch == batch_id]$norm_file
  vals <- unlist(lapply(fs, function(f) read_ch(f, ch)))
  stats::quantile(vals, probs=probs, na.rm=TRUE, names=FALSE)
}

# Compute linear parameters a,b so that a*x + b maps batch 5th/95th to reference 5th/95th
fit_ab <- function(q_src, q_ref){
  # Solve: a*q_src_low + b = q_ref_low; a*q_src_hi + b = q_ref_hi
  a <- (q_ref[2] - q_ref[1]) / (q_src[2] - q_src[1] + 1e-8)
  b <- q_ref[1] - a * q_src[1]
  c(a=a, b=b)
}

# Build mapping table for all flagged channels × batches
map_tbl <- data.table()
if (nrow(anchors) == 0) {
  warning("No unstim anchors found in sample_table; skipping rescaling maps.")
} else if (length(cand_rescale) == 0) {
  message("No phospho channels flagged; skipping rescaling maps.")
} else {
  map_tbl <- rbindlist(lapply(cand_rescale, function(ch){
    q_ref <- ref_quants(ch)
    rbindlist(lapply(unique(tab$batch), function(b){
      q_src <- batch_quants(b, ch)
      pars  <- fit_ab(q_src, q_ref)
      data.table(batch=b, channel=ch, a=pars["a"], b=pars["b"],
                 src_q05=q_src[1], src_q95=q_src[2], ref_q05=q_ref[1], ref_q95=q_ref[2])
    }))
  }))
  if (nrow(map_tbl)) fwrite(map_tbl, file.path(qc_dir, "rescale_maps.csv"))
}

# ---------- 4) Apply re-scaling to all files (only flagged channels) ----------
rescale_apply <- function(f_in, f_out, maps_rowset){
  ff <- read.FCS(f_in, transformation=FALSE)
  M  <- exprs(ff)
  for (j in seq_len(nrow(maps_rowset))) {
    ch <- maps_rowset$channel[j]
    if (!ch %in% colnames(M)) next
    a <- maps_rowset$a[j]; b <- maps_rowset$b[j]
    M[, ch] <- a * M[, ch] + b
  }
  exprs(ff) <- M
  suppressWarnings(write.FCS(ff, filename=f_out))
}

if (nrow(map_tbl)) {
  message("Applying per-channel linear re-scaling...")
  for (i in seq_len(nrow(tab))) {
    f_in  <- tab$norm_file[i]
    f_out <- out_path(tab$file[i])
    maps  <- map_tbl[batch == tab$batch[i]]
    if (nrow(maps)) rescale_apply(f_in, f_out, maps) else file.copy(f_in, f_out, overwrite=TRUE)
    cat("[ok] ->", f_out, "\n")
  }
} else {
  if (length(cand_rescale) > 0 && nrow(anchors) == 0) {
    message("Rescaling skipped (no unstim anchors); copying inputs to ", out_dir_fcs)
  } else {
    message("No channels flagged; copying inputs to ", out_dir_fcs)
  }
  file.copy(tab$norm_file, out_path(tab$file), overwrite=TRUE)
}

# ---------- 5) Diagnostics (post-rescaling) ----------
message("Recomputing diagnostics after re-scaling...")

post_stats <- file_stats(out_path(tab$file), markers_use)
fwrite(post_stats, file.path(qc_dir, "post_file_stats.csv"))

anova_post <- post_stats %>%
  left_join(tab[, .(file = basename(file), batch)], by="file") %>%
  group_by(channel) %>%
  summarise(
    p_value = tryCatch({
      fit <- aov(median ~ factor(batch))
      summary(fit)[[1]][["Pr(>F)"]][1]
    }, error=function(e) NA_real_),
    n_batches = n_distinct(batch)
  ) %>%
  mutate(p_adj = p.adjust(p_value, method="fdr"))
fwrite(anova_post, file.path(qc_dir, "anova_post.csv"))

p99_post <- post_stats %>%
  group_by(channel) %>%
  summarise(p99_range = max(p99, na.rm=TRUE) - min(p99, na.rm=TRUE))
fwrite(p99_post, file.path(qc_dir, "p99_range_post.csv"))

# Quick comparison heatmap (delta medians by file/channel)
preM  <- pre_stats %>% select(file, channel, median)  %>% mutate(type="pre")
postM <- post_stats %>% select(file, channel, median) %>% mutate(type="post")
delta <- postM %>% inner_join(preM, by=c("file","channel"), suffix=c("_post","_pre")) %>%
  mutate(delta = median_post - median_pre)
DM <- reshape2::acast(delta, file ~ channel, value.var="delta")
heatmap_path <- file.path(qc_dir, "median_delta_heatmap.pdf")
finite_vals <- DM[is.finite(DM)]
if (length(finite_vals) == 0 || (max(finite_vals) - min(finite_vals)) == 0) {
  warning("Median delta matrix is constant; skipping heatmap at ", heatmap_path)
} else {
  pdf(heatmap_path, width=11, height=7)
  pheatmap::pheatmap(t(DM), main="Median (post - pre) per channel/file")
  dev.off()
}

# ---------- 6) UMAP + Silhouette by sample ----------
# Sample per file to keep it light
message("Computing UMAP and silhouette by sample...")
max_cells_per_file <- 5000
max_cells_total <- 80000

sample_cells <- function(f, feats){
  M <- exprs(read.FCS(f, transformation=FALSE))
  feats <- intersect(feats, colnames(M))
  X <- M[, feats, drop=FALSE]
  n <- nrow(X)
  if (n > max_cells_per_file) X <- X[sample.int(n, max_cells_per_file), , drop=FALSE]
  X
}

umap_feats <- unique(c(lineage, phospho))
X_list <- lapply(out_path(tab$file), sample_cells, feats=umap_feats)
rows_per_file <- vapply(X_list, nrow, integer(1))
X <- do.call(rbind, X_list)
lab_sample <- rep(basename(tab$file), times = rows_per_file)
if (nrow(X) > max_cells_total) {
  idx <- sample.int(nrow(X), max_cells_total)
  X <- X[idx, , drop=FALSE]
  lab_sample <- lab_sample[idx]
}
# standardize features
X <- scale(X)

emb <- uwot::umap(X, n_neighbors = 15, min_dist = 0.2, n_components = 2, metric = "euclidean", verbose = TRUE)
colnames(emb) <- c("UMAP1","UMAP2")

# Silhouette by sample (labels = sample/file)
# --- UMAP already computed into 'emb' and labels 'lab_sample' ---

sil_subset <- 5000L  # <= 5k is usually safe
if (nrow(emb) > sil_subset) {
  set.seed(1)
  keep <- sample.int(nrow(emb), sil_subset)
  emb_s  <- emb[keep, , drop=FALSE]
  labs_s <- lab_sample[keep]
} else {
  emb_s  <- emb
  labs_s <- lab_sample
}

# compute distances ONLY on the subset
dmat <- dist(emb_s)
sil  <- cluster::silhouette(as.integer(factor(labs_s)), dmat)
sil_df <- data.frame(sample = labs_s, sil_width = sil[,3])
sil_summary <- dplyr::summarise(dplyr::group_by(sil_df, sample),
                                mean_sil = mean(sil_width, na.rm=TRUE))

data.table::fwrite(sil_df,      file.path(qc_dir, "silhouette_per_cell.csv"))
data.table::fwrite(sil_summary, file.path(qc_dir, "silhouette_by_sample.csv"))
cat("Mean silhouette on subset: ",
    round(mean(sil_summary$mean_sil, na.rm=TRUE), 3), "\n")


# Plot UMAP colored by sample (optional)
p <- data.frame(emb, sample = lab_sample) %>%
  ggplot(aes(UMAP1, UMAP2, color = sample)) + geom_point(size=0.2, alpha=0.6) +
  theme_minimal(base_size=9) + guides(color=FALSE) + ggtitle("UMAP colored by sample (rescaled)")
ggsave(file.path(qc_dir, "umap_by_sample_rescaled.png"), p, width=7, height=5, dpi=200)


# ---------- Nice UMAP + clean legends (kNN = 15) ----------
if (!requireNamespace("irlba", quietly=TRUE)) {
  stop("Please install.packages('irlba') for memory-efficient PCA")
}

# PCA -> UMAP (kNN=15) with robust fallback for small matrices
Xs <- scale(X)  # X already sampled & arcsinh
n_pcs <- max(2, min(30, ncol(Xs) - 1, nrow(Xs) - 1))
ok_umap2 <- TRUE
emb <- NULL
tryCatch({
  pc <- irlba::prcomp_irlba(Xs, n = n_pcs, center = FALSE, scale. = FALSE)
  emb <- uwot::umap(
    pc$x,
    n_neighbors = 15,    # as requested
    min_dist    = 0.35,  # a bit smoother clusters; tweak 0.3–0.5
    n_components= 2,
    metric      = "euclidean",
    verbose     = TRUE
  )
  colnames(emb) <- c("UMAP1","UMAP2")
}, error=function(e){
  warning("Secondary PCA+UMAP failed: ", conditionMessage(e)); ok_umap2 <<- FALSE
})

# Build plotting dataframe with metadata
df <- data.frame(if (is.null(emb)) matrix(NA_real_, nrow=0, ncol=2, dimnames=list(NULL, c("UMAP1","UMAP2"))) else emb,
                 sample = lab_sample[seq_len(ifelse(is.null(emb), 0, nrow(emb)))], stringsAsFactors = FALSE)
meta_map <- tab[, .(sample = basename(file),
                    batch  = factor(batch, levels = sort(unique(batch))),
                    stim   = factor(ifelse(is.na(stim) | stim=="", "unknown", tolower(stim))))]
if ("stim_class" %in% names(tab)) {
  meta_map[, stim_class := tolower(tab$stim_class)]
} else {
  meta_map[, stim_class := NA_character_]
}
# Derive a simplified stim class for plotting (first token; anchors/US -> unstim)
meta_map[, stim_simple := stim]
meta_map[!is.na(stim_simple) & grepl("_", stim_simple), stim_simple := sub("_.*$", "", stim_simple)]
meta_map[, stim_simple := factor(ifelse(stim_class == "unstim", "unstim", as.character(stim_simple)))]
df <- dplyr::left_join(df, meta_map, by = "sample")

# A fast point layer that degrades gracefully if pkg missing
point_layer <- function() {
  if (requireNamespace("ggpointdensity", quietly=TRUE)) {
    ggpointdensity::geom_pointdensity(size=0.2) + scale_color_viridis_c()
  } else if (requireNamespace("scattermore", quietly=TRUE)) {
    scattermore::geom_scattermore(pointsize = 1.2, alpha = 0.6)
  } else {
    geom_point(size = 0.2, alpha = 0.6)
  }
}

# THEME + legend helpers
theme_umap <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(4,  "mm")
  )

# 1) UMAP by BATCH  (compact legend)
p_batch <- ggplot(df, aes(UMAP1, UMAP2, color = batch)) +
  point_layer() +
  guides(color = guide_legend(
    title = "Batch",
    override.aes = list(size = 3, alpha = 1),
    nrow = 1
  )) +
  labs(title = "UMAP colored by batch (rescaled)") +
  theme_umap + theme(legend.position = "top")
ggsave(file.path(qc_dir, "umap_by_batch_rescaled.png"),
       p_batch, width = 7.5, height = 6, dpi = 300)

# 2) UMAP by STIM  (compact legend; consistent colors)
stim_pal <- c("cd3"="#D95F02","cd328"="#1B9E77","il15"="#7570B3",
              "il2"="#66A61E","unknown"="#999999","unstim"="#E7298A")
stim_pal <- stim_pal[names(stim_pal) %in% levels(df$stim_simple)]

if (ok_umap2 && nrow(df) > 0) {
  p_stim <- ggplot(df, aes(UMAP1, UMAP2, color = stim_simple)) +
    point_layer() +
    scale_color_manual(values = stim_pal) +
    guides(color = guide_legend(
      title = "Stim",
      override.aes = list(size = 3, alpha = 1),
      nrow = 1
    )) +
    labs(title = "UMAP colored by stimulation (rescaled)") +
    theme_umap + theme(legend.position = "top")
  ggsave(file.path(qc_dir, "umap_by_stim_rescaled.png"),
         p_stim, width = 7.5, height = 6, dpi = 300)
}

# 3) UMAP by SAMPLE (legend OFF – too many categories)
#    If you want a legend, save it as a separate file with top N samples only.
if (ok_umap2 && nrow(df) > 0) {
  p_sample <- ggplot(df, aes(UMAP1, UMAP2, color = sample)) +
    point_layer() +
    guides(color = "none") +
    labs(title = "UMAP colored by sample (rescaled)") +
    theme_umap
  ggsave(file.path(qc_dir, "umap_by_sample_rescaled.png"),
         p_sample, width = 8, height = 6, dpi = 300)
}

# OPTIONAL: save a compact legend for top-N samples (by cell count)
topN <- 12
sample_counts <- df |>
  dplyr::count(sample, name="n") |>
  dplyr::arrange(desc(n)) |>
  dplyr::slice(1:topN)
df_top <- df[df$sample %in% sample_counts$sample, ]
if (ok_umap2 && nrow(df_top) > 0) {
  p_sample_legend <- ggplot(df_top, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.4, alpha = 0.7) +
    theme_umap +
    guides(color = guide_legend(
      title = paste0("Top ", topN, " samples"),
      override.aes = list(size = 3, alpha = 1),
      ncol = 2
    )) +
    labs(title = "UMAP by top samples (legend only)")
  ggsave(file.path(qc_dir, "umap_by_sample_topN_with_legend.png"),
         p_sample_legend, width = 7.5, height = 6.5, dpi = 300)
}




# Summary to console
cat("\n=== SUMMARY ===\n")
cat("Flagged & rescaled channels: ", ifelse(length(cand_rescale), paste(cand_rescale, collapse=", "), "none"), "\n")
cat("Median batch ANOVA (channels with FDR < 0.05): PRE=", sum(anova_pre$p_adj < 0.05, na.rm=TRUE),
    " POST=", sum(anova_post$p_adj < 0.05, na.rm=TRUE), "\n")
cat("Silhouette by sample (mean across samples) ~ ",
    round(mean(sil_summary$mean_sil, na.rm=TRUE), 3), " (closer to 0 is better)\n")
cat("Outputs:\n - Rescaled FCS: ", out_dir_fcs,
    "\n - ANOVA CSVs: ", file.path(qc_dir, "anova_pre.csv"), " & ", file.path(qc_dir, "anova_post.csv"),
    "\n - P99 ranges: ", file.path(qc_dir, "p99_range_pre.csv"), " & ", file.path(qc_dir, "p99_range_post.csv"),
    "\n - Silhouette: ", file.path(qc_dir, "silhouette_by_sample.csv"), "\n")
