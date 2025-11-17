#!/usr/bin/env Rscript

## NOTE (archived):
## Legacy end-to-end differential phospho-signaling script.
## The current pipeline instead uses 05_compute_phospho_medians.R
## followed by 06_diff_phospho_analysis.R and downstream helpers.

suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(data.table)
  library(dplyr)
  library(limma)
})

option_list <- list(
  make_option("--norm-dir", type = "character", default = "normalized_fcs_full"),
  make_option("--sample-table", type = "character", default = "sample_table.csv"),
  make_option("--flowsom-rds", type = "character",
              default = "flowsom_results/flowsom_fitted.rds"),
  make_option("--outdir", type = "character", default = "diff_phos_results"),
  make_option("--cofactor", type = "double", default = 5),
  make_option("--seed", type = "integer", default = 1)
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cofactor <- opt$cofactor

# === Step 1: metadata ===
meta <- fread(opt$`sample-table`)
if (!"file_base" %in% names(meta)) {
  if ("file" %in% names(meta)) meta[, file_base := basename(file)] else
    stop("sample_table must contain 'file_base' or 'file'.")
}
if (!all(c("time_min") %in% names(meta))) {
  stop("sample_table needs 'time_min' column (minutes).")
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
             ifelse(grepl("unstim", v) | v == "" | is.na(v), "UNSTIM", "OTHER"))))))
      res
    }
    meta[, cytokine := stim_to_cytokine(stim)]
  } else {
    stop("sample_table needs 'cytokine' column (IL2, IL15, IL21, PAC, etc.), or a 'stim' column to derive it from.")
  }
}

meta_samp <- meta[is_anchor == FALSE]
norm_files <- file.path(opt$`norm-dir`, paste0("Norm_", meta_samp$file_base))
exist_idx <- file.exists(norm_files)
meta_samp <- meta_samp[exist_idx]
norm_files <- norm_files[exist_idx]

# === Step 2: load FlowSOM & define phospho channels ===
fsom <- readRDS(opt$`flowsom-rds`)
codes <- fsom$map$codes

param <- pData(parameters(read.FCS(norm_files[1])))
desc  <- param$desc; cn <- param$name

marker_match <- function(pattern) {
  ix <- grep(pattern, desc, ignore.case = TRUE)
  if (length(ix) == 0) ix <- grep(pattern, cn, ignore.case = TRUE)
  if (length(ix) == 0) return(NA)
  cn[ix[1]]
}

# lineage channels (must match those used for FlowSOM)
lineage_markers <- c("Cd112Di","Cd111Di","Cd116Di","Sm147Di","Tb159Di",
                     "Gd155Di","Nd150Di","Yb174Di","Tm169Di","Er168Di",
                     "Cd113Di","Gd157Di")
lineage_channels <- intersect(lineage_markers, cn)

phos_names <- c("pZAP70", "pSLP76", "pAkt", "pNFATc1", "pSTAT3",
                "pPLCy", "pCREB", "pGSK3", "pAMPK",
                "pMK2", "pPGC1", "IkBa", "pMEK", "pS6", "pSTAT5")

phos_patterns <- c("pZAP", "pSLP", "pAkt", "pNFAT", "pSTAT3",
                   "pPLC", "pCREB", "pGSK3", "pAMPK",
                   "pMK2", "pPGC1", "IkB", "pMEK", "pS6", "pSTAT5")

phos_channels <- sapply(phos_patterns, marker_match)
keep <- !is.na(phos_channels)
phos_channels <- phos_channels[keep]
phos_names    <- phos_names[keep]

cat("Phospho channels used:\n")
print(data.frame(marker = phos_names, channel = phos_channels))

tf <- flowCore::arcsinhTransform(a = 0, b = 1 / cofactor, c = 0)
tlist <- flowCore::transformList(c(lineage_channels, phos_channels), tf)

# === Step 3: aggregate medians per sample × cluster × time × cytokine ===
agg_list <- list()

for (i in seq_along(norm_files)) {
  nf <- norm_files[i]
  fb <- meta_samp$file_base[i]
  time_min <- meta_samp$time_min[i]
  cytok    <- meta_samp$cytokine[i]

  ff  <- read.FCS(nf, transformation = FALSE, truncate_max_range = FALSE)
  ffT <- flowCore::transform(ff, tlist)
  expr_mat <- exprs(ffT)
  if (!all(lineage_channels %in% colnames(expr_mat))) next
  if (!all(phos_channels %in% colnames(expr_mat))) next

  m_lineage <- expr_mat[, lineage_channels, drop = FALSE]
  m_phos    <- expr_mat[, phos_channels,    drop = FALSE]

  # map events to closest FlowSOM code (metacluster)
  nn <- FNN::get.knnx(data = codes, query = m_lineage, k = 1)$nn.index
  mc <- fsom$metaclustering[nn]

  dt <- as.data.table(m_phos)
  dt[, cluster := mc]

  med <- dt[, lapply(.SD, median), by = cluster, .SDcols = phos_channels]
  med[, file_base := fb]
  med[, time_min  := time_min]
  med[, cytokine  := cytok]

  agg_list[[length(agg_list) + 1]] <- med
}

agg <- rbindlist(agg_list, use.names = TRUE, fill = TRUE)
setcolorder(agg, c("file_base", "cluster", "time_min", "cytokine", phos_channels))

# rename columns to human-readable phospho names
setnames(agg, old = phos_channels, new = phos_names)

fwrite(agg, file.path(opt$outdir, "phos_medians_by_sample_cluster.csv"))
cat("Wrote aggregated medians to phos_medians_by_sample_cluster.csv\n")

# === Step 4: limma differential analysis ===
results_all <- list()

clusters <- sort(unique(agg$cluster))
cyts     <- sort(unique(agg$cytokine))

for (cl in clusters) {
  for (cy in cyts) {
    sub <- agg[cluster == cl & cytokine == cy]
    if (nrow(sub) < 4) next

    time_fac <- factor(sub$time_min, levels = sort(unique(sub$time_min)))
    if (!any(levels(time_fac) %in% c("0", 0))) {
      # use earliest time as baseline
      baseline <- levels(time_fac)[1]
    } else {
      baseline <- "0"
    }

    design <- model.matrix(~ 0 + time_fac)
    colnames(design) <- paste0("t", levels(time_fac))

    Y <- t(as.matrix(sub[, ..phos_names]))
    fit <- lmFit(Y, design)

    # contrasts: each time vs baseline
    lev <- levels(time_fac)
    lev_other <- setdiff(lev, baseline)
    if (length(lev_other) == 0) next

    contr_str <- paste(
      sprintf("t%s - t%s", lev_other, baseline),
      collapse = ";"
    )
    contrast_mat <- makeContrasts(contrasts = contr_str, levels = design)

    fit2 <- contrasts.fit(fit, contrast_mat)
    fit2 <- eBayes(fit2)

    for (j in seq_len(ncol(contrast_mat))) {
      cname <- colnames(contrast_mat)[j]
      tt <- topTable(fit2, coef = j, number = Inf, sort.by = "none")
      tt$cluster  <- cl
      tt$cytokine <- cy
      tt$contrast <- cname
      tt$marker   <- rownames(tt)
      results_all[[length(results_all) + 1]] <- tt
    }
  }
}

if (length(results_all) > 0) {
  res <- rbindlist(results_all, use.names = TRUE, fill = TRUE)
  fwrite(res, file.path(opt$outdir, "diff_phos_limma_results.csv"))
  cat("Wrote limma results to diff_phos_limma_results.csv\n")
} else {
  cat("No contrasts could be computed (check data).\n")
}
