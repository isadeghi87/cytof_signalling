#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(data.table)
  library(Rphenograph)
  library(igraph)
  library(ggplot2)
})

opt_list <- list(
  make_option("--norm-dir", type="character", default="normalized_fcs"),
  make_option("--sample-table", type="character", default="sample_table.csv"),
  make_option("--outdir", type="character", default="phenograph_results"),
  make_option("--cofactor", type="double", default=5),
  make_option("--n-events", type="integer", default=50000)
)

opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
# Load metadata and files
###############################################################################

meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) meta[, file_base := basename(file)]

meta_s <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_s$file_base))
idx <- file.exists(norm_files)
meta_s <- meta_s[idx]
norm_files <- norm_files[idx]

###############################################################################
# Identify lineage panel
###############################################################################

ff <- read.FCS(norm_files[1], transformation=FALSE)
cn <- pData(parameters(ff))$name

lineage_ch <- intersect(
  c("Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
    "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di",
    "Cd113Di","Gd157Di"), cn)

tf <- arcsinhTransform(a=0,b=1/opt$cofactor,c=0)
tlist <- transformList(lineage_ch, tf)

###############################################################################
# Downsample
###############################################################################

expr_list <- list()
per_file <- round(opt$`n-events` / length(norm_files))

for (i in seq_along(norm_files)) {
  ff0 <- read.FCS(norm_files[i], transformation=FALSE)
  ffT <- transform(ff0, tlist)
  m <- exprs(ffT)[, lineage_ch]

  take <- min(nrow(m), per_file)
  expr_list[[i]] <- m[sample(seq_len(nrow(m)), take), ]
}

mat <- do.call(rbind, expr_list)

###############################################################################
# Run PhenoGraph
###############################################################################

ph <- tryCatch({
  Rphenograph(mat, k=30)
}, error = function(e) {
  message("[WARN] Rphenograph failed: ", e$message, " – falling back to k-means clustering.")
  NULL
})

if (is.null(ph)) {
  set.seed(1)
  km <- kmeans(mat, centers = 20)
  membership_vec <- factor(km$cluster)
} else {
  membership_vec <- factor(igraph::membership(ph[[2]]))
}

dt <- data.table(mat)
dt$cluster <- membership_vec

colnames(dt)[1:2] <- c("X1","X2")

ggplot(dt, aes(X1, X2, col=cluster)) +
  geom_point(size=0.3, alpha=0.7) +
  theme_bw()

ggsave(file.path(opt$outdir,"phenograph_scatter.png"), width=8,height=6)

fwrite(data.table(cluster=membership_vec), file.path(opt$outdir,"phenograph_labels.csv"))

cat("Done.\n")
