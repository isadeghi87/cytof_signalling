#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(data.table)
  library(FNN)
  library(dplyr)
})

opt_list <- list(
  make_option("--norm-dir", type="character", default="normalized_fcs"),
  make_option("--sample-table", type="character", default="sample_table.csv"),
  make_option("--flowsom-rds", type="character",
              default="flowsom_results/flowsom_fitted.rds"),
  make_option("--outdir", type="character", default="phospho_medians"),
  make_option("--cofactor", type="double", default=5),
  make_option("--seed", type="integer", default=1)
)

opt <- parse_args(OptionParser(option_list=opt_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
# Metadata
###############################################################################

meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) {
    meta[, file_base := basename(file)]
  } else {
    stop("sample_table must contain 'file_base' or 'file'.")
  }
}

if (!"is_anchor" %in% names(meta)) {
  stop("sample_table must contain 'is_anchor' logical column.")
}

if (!"time_min" %in% names(meta)) {
  stop("sample_table must contain 'time_min' column.")
}

if (!"cytokine" %in% names(meta)) {
  if ("stim" %in% names(meta)) {
    stim_to_cytokine <- function(stim_vec) {
      v <- tolower(as.character(stim_vec))
      res <- ifelse(grepl("il21", v), "IL21",
             ifelse(grepl("il15", v), "IL15",
             ifelse(grepl("il2",  v), "IL2",
             ifelse(grepl("pac", v), "PAC",
             ifelse(grepl("cd3", v) & !grepl("il", v), "CD3",
             ifelse(grepl("unstim", v) | v == "" | is.na(v),
                    "UNSTIM", "OTHER"))))))
      res
    }
    meta[, cytokine := stim_to_cytokine(stim)]
  } else {
    stop("sample_table needs 'cytokine' column, or a 'stim' column to derive it.")
  }
}

meta_s <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_s$file_base))
idx <- file.exists(norm_files)
meta_s <- meta_s[idx]
norm_files <- norm_files[idx]

###############################################################################
# Load FlowSOM
###############################################################################

fs <- readRDS(opt$`flowsom-rds`)
codes <- fs$map$codes

###############################################################################
# Identify markers
###############################################################################

ff <- read.FCS(norm_files[1], transformation=FALSE)
par <- pData(parameters(ff))
cn <- par$name

lineage_ch <- intersect(
  c("Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
    "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di",
    "Cd113Di","Gd157Di"), cn )

phos_ch <- intersect(
  c("Er166Di","Er170Di","Eu153Di","Gd156Di","Gd158Di",
    "Ho165Di","Lu175Di","Nd143Di","Nd145Di","Nd146Di",
    "Nd148Di","Sm152Di","Sm154Di","Yb172Di","Yb173Di"),
  cn)

if (!length(lineage_ch)) {
  stop("No lineage channels found in normalized FCS; cannot compute medians.")
}
if (!length(phos_ch)) {
  warning("No phospho channels found in normalized FCS; output will be empty.")
}

tf <- arcsinhTransform(a=0, b=1/opt$cofactor, c=0)
tlist <- transformList(c(lineage_ch, phos_ch), tf)

###############################################################################
# Compute medians: sample × cluster × time × cytokine
###############################################################################

res <- list()

for (i in seq_along(norm_files)) {
  nf <- norm_files[i]

  ff0 <- read.FCS(nf, transformation=FALSE, truncate_max_range=FALSE)
  ffT <- transform(ff0, tlist)
  M <- exprs(ffT)

  if (!all(lineage_ch %in% colnames(M)) || !all(phos_ch %in% colnames(M))) {
    warning("Skipping file without complete lineage/phospho channels: ", nf)
    next
  }

  # FlowSOM cluster mapping
  nn <- FNN::get.knnx(codes, M[, lineage_ch], k=1)$nn.index
  cl <- fs$metaclustering[nn]

  dt <- as.data.table(M[, phos_ch, drop=FALSE])
  dt[, cluster := cl]

  med <- dt[, lapply(.SD, median), by=cluster]
  med[, file_base := meta_s$file_base[i]]
  med[, time_min  := meta_s$time_min[i]]
  med[, cytokine  := meta_s$cytokine[i]]

  res[[length(res)+1]] <- med
}

if (!length(res)) {
  warning("No medians could be computed; check input files and channels.")
  out <- data.table(
    cluster = integer(),
    file_base = character(),
    time_min = numeric(),
    cytokine = character()
  )
} else {
  out <- rbindlist(res, fill=TRUE, use.names=TRUE)
}

fwrite(out, file.path(opt$outdir,"phos_medians_by_sample_cluster.csv"))
cat("Saved medians.\n")
