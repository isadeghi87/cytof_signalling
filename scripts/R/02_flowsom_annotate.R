#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(FlowSOM)
  library(data.table)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})

###############################################################################
# Options
###############################################################################

opt_list <- list(
  make_option("--norm-dir", type="character", default="normalized_fcs"),
  make_option("--sample-table", type="character", default="sample_table.csv"),
  make_option("--outdir", type="character", default="flowsom_results"),
  make_option("--n-events", type="integer", default=150000),
  make_option("--metak", type="integer", default=20),
  make_option("--cofactor", type="double", default=5),
  make_option("--seed", type="integer", default=1)
)

opt <- parse_args(OptionParser(option_list=opt_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
# Step 1: load metadata + normalized files
###############################################################################

cat("\n=== Load metadata ===\n")
meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) meta[, file_base := basename(file)] else
    stop("sample_table must contain file_base or file")
}

meta_s <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_s$file_base))
idx <- file.exists(norm_files)
meta_s <- meta_s[idx]
norm_files <- norm_files[idx]

if (length(norm_files)==0) stop("No normalized files detected")

###############################################################################
# Step 2: identify lineage + phospho markers
###############################################################################

cat("\n=== Identify marker channels ===\n")

ff <- read.FCS(norm_files[1], transformation=FALSE)

par_tab <- pData(parameters(ff))
desc <- par_tab$desc
cn   <- par_tab$name

match_channel <- function(x) {
  i <- grep(x, desc, ignore.case=TRUE)
  if (length(i)==0) i <- grep(x, cn, ignore.case=TRUE)
  if (length(i)==0) return(NA)
  cn[i[1]]
}

# LINEAGE channels recommended from your experiment
lineage_markers <- c(
  "Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
  "Gd155Di","Nd150Di","Yb174Di","Tm169Di",
  "Er168Di","Cd113Di","Gd157Di"
)

# PHOSPHO markers
phos_markers <- c(
  "Er166Di","Er170Di","Eu153Di","Gd156Di","Gd158Di","Ho165Di",
  "Lu175Di","Nd143Di","Nd145Di","Nd146Di","Nd148Di",
  "Sm152Di","Sm154Di","Yb172Di","Yb173Di"
)

lineage_ch <- intersect(lineage_markers, cn)
phos_ch    <- intersect(phos_markers, cn)

cat("Lineage channels:", paste(lineage_ch, collapse=","), "\n")
cat("Phospho channels:", paste(phos_ch, collapse=","), "\n")

tf <- arcsinhTransform(a=0, b=1/opt$cofactor, c=0)
tlist <- transformList(c(lineage_ch, phos_ch), tf)

###############################################################################
# Step 3: downsample and build FlowSOM
###############################################################################

cat("\n=== Build FlowSOM ===\n")

expr_list <- lapply(norm_files, function(f){
  ff0 <- read.FCS(f, transformation=FALSE)
  ffT <- transform(ff0, tlist)
  exprs(ffT)[, lineage_ch, drop=FALSE]
})

total_ev <- sum(sapply(expr_list, nrow))
cat("Total events:", total_ev, "\n")

if (total_ev > opt$`n-events`) {
  frac <- opt$`n-events` / total_ev
  expr_sub <- do.call(rbind, lapply(expr_list, function(m){
    keep <- max(2000, round(nrow(m)*frac))
    keep <- min(nrow(m), keep)
    m[sample(seq_len(nrow(m)), keep), , drop=FALSE]
  }))
} else {
  expr_sub <- do.call(rbind, expr_list)
}

cat("Events for FlowSOM:", nrow(expr_sub), "\n")

fs <- FlowSOM::ReadInput(expr_sub, transform=FALSE, scale=FALSE)
fs <- FlowSOM::BuildSOM(fs, colsToUse=seq_along(lineage_ch), xdim=10, ydim=10)
fs <- FlowSOM::BuildMST(fs)

metaCl <- metaClustering_consensus(fs$map$codes, k=opt$metak)
fs$metaclustering <- metaCl

saveRDS(fs, file.path(opt$outdir, "flowsom_fitted.rds"))

###############################################################################
# Plot MST Tree
###############################################################################

png(file.path(opt$outdir, "flowsom_MST_tree.png"), width=1800, height=1200)
PlotStars(fs, backgroundValues=fs$metaclustering)
dev.off()

###############################################################################
# Step 4: marker medians per metacluster
###############################################################################

cat("\n=== Compute per-metacluster medians ===\n")

cl_ids <- fs$metaclustering[fs$map$mapping[,1]]
dt <- as.data.table(expr_sub)
dt[, cluster := cl_ids]

med_lineage <- dt[, lapply(.SD, median), by=cluster, .SDcols=lineage_ch]
setorder(med_lineage, cluster)

mat <- as.matrix(med_lineage[, ..lineage_ch])
rownames(mat) <- paste0("MC", med_lineage$cluster)

pheatmap(
  mat,
  scale="row",
  clustering_distance_rows="euclidean",
  clustering_distance_cols="correlation",
  filename=file.path(opt$outdir,"metacluster_lineage_heatmap.png"),
  width=10, height=8
)

###############################################################################
# Step 5: phospho heatmap (optional)
###############################################################################

if (length(phos_ch)) {
  phos_vals <- lapply(norm_files, function(f){
    ff0 <- read.FCS(f, transformation=FALSE)
    ffT <- transform(ff0, tlist)
    exprs(ffT)[, phos_ch, drop=FALSE]
  })

  phos_mat <- do.call(rbind, phos_vals)

  if (nrow(phos_mat) == nrow(dt)) {
    cl_ids2 <- FNN::get.knnx(fs$map$codes, dt[, lineage_ch, with=FALSE], k=1)$nn.index
    dt_phos <- as.data.table(phos_mat)
    dt_phos[, cluster := fs$metaclustering[cl_ids2]]

    med_phos <- dt_phos[, lapply(.SD, median), by=cluster, .SDcols=phos_ch]
    setorder(med_phos, cluster)

    pheatmap(
      as.matrix(med_phos[, ..phos_ch]),
      scale="row",
      filename=file.path(opt$outdir, "metacluster_phospho_heatmap.png"),
      width=10, height=8
    )
  } else {
    message("[WARN] Skipping phospho metacluster heatmap: row mismatch between lineage (",
            nrow(dt), ") and phospho (", nrow(phos_mat), ") events.")
  }
}

###############################################################################
# Save cluster label template
###############################################################################

annot <- data.table(
  cluster = med_lineage$cluster,
  label = paste0("Cluster_", med_lineage$cluster)
)
fwrite(annot, file.path(opt$outdir,"cluster_labels_template.csv"))

cat("\nDone FlowSOM annotation.\n")
