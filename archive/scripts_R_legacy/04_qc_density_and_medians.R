#!/usr/bin/env Rscript

## NOTE (archived):
## Older QC script for density overlays and medians; superseded
## by 04_qc_plots.R and 05_plot_normalization_efficiency.R in
## the current pipeline.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table); library(ggplot2)
  library(flowCore);   library(pheatmap)
})

option_list <- list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv"),
  make_option(c("--normalized-dir"), type="character", default="normalized_fcs"),
  make_option(c("--outdir"), type="character", default="qc"),
  make_option(c("--phospho-list"), type="character", default="effective_phospho_markers.txt"),
  make_option(c("--cofactor"), type="double", default=5)
)
opt <- parse_args(OptionParser(option_list = option_list))

tab <- data.table::fread(opt$`sample-table`)
outdir <- opt$outdir; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
normdir <- opt$`normalized-dir`

if (file.exists(opt$`phospho-list`)) {
  phospho <- readLines(opt$`phospho-list`)
} else {
  stop("Phospho list not found: ", opt$`phospho-list`)
}
cofactor <- opt$cofactor

# resolve a file path from sample table: accept absolute, or look under 'batch/', 'normalized_fcs/', or cwd
resolve_file <- function(fp){
  if (is.na(fp) || fp == "") return(NA_character_)
  if (file.exists(fp)) return(normalizePath(fp))
  bn <- basename(fp)
  cand <- c(file.path("batch", bn), file.path(normdir, bn), bn)
  for (c in cand) if (file.exists(c)) return(normalizePath(c))
  return(NA_character_)
}

read_ch <- function(f, ch){
  realf <- resolve_file(f)
  if (is.na(realf)) return(numeric(0))
  fcs <- tryCatch(read.FCS(realf, transformation=FALSE), error=function(e) NULL)
  if (is.null(fcs)) return(numeric(0))
  mat <- exprs(fcs)
  if (!(ch %in% colnames(mat))) return(numeric(0))
  as.numeric(mat[, ch])
}

# Density overlays per phospho channel (anchors highlighted)
for (ch in phospho) {
  df <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
    f <- tab$file[i]
    v <- read_ch(f, ch); if (length(v) > 20000) v <- sample(v, 20000)
    data.table(intensity = asinh(v/cofactor), file = basename(f), type = "raw", anchor = tab$is_anchor[i])
  }))
  df2 <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
    f <- file.path(normdir, basename(tab$file[i]))
    v <- read_ch(f, ch); if (length(v) > 20000) v <- sample(v, 20000)
    data.table(intensity = asinh(v/cofactor), file = basename(f), type = "norm", anchor = tab$is_anchor[i])
  }))
  D <- rbind(df, df2)
  p <- ggplot(D, aes(intensity, group=file, colour=type, alpha=anchor)) +
    geom_density(size=0.6) +
    scale_alpha_manual(values=c(`FALSE`=0.5, `TRUE`=1.0)) +
    theme_minimal(base_size=9) + labs(title=paste("Density:", ch), y="Density", x="asinh/5")
  ggsave(file.path(outdir, paste0("density_", ch, ".pdf")), p, width=4, height=3)
}

# Median shift heatmap (post - pre) over phosphos
Mpre  <- sapply(phospho, function(ch) sapply(tab$file, function(f) median(read_ch(f, ch), na.rm=TRUE)))
Mpost <- sapply(phospho, function(ch) sapply(file.path(normdir, basename(tab$file)), function(f) median(read_ch(f, ch), na.rm=TRUE)))

# sapply returns a matrix with rows = samples and cols = phospho channels
rownames(Mpre) <- tab$file
rownames(Mpost) <- tab$file
colnames(Mpre) <- colnames(Mpost) <- phospho

# compute post - pre (rows: samples, cols: phospho), then transpose so rows=phospho for the heatmap
shift <- Mpost - Mpre
pdf(file.path(outdir, "median_shifts_phospho.pdf"), width=10, height=7)
pheatmap::pheatmap(t(shift), cluster_rows=TRUE, cluster_cols=TRUE, main="Median(post - pre): phospho only", fontsize=8)
dev.off()
