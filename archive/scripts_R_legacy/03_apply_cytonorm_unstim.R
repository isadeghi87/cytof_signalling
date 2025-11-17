#!/usr/bin/env Rscript

## NOTE (archived):
## Legacy helper to apply CytoNorm for unstim-only workflows.
## Not used in the current end-to-end pipeline; kept for reference.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(flowCore)
  library(CytoNorm)
  library(stringr)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv"),
  make_option(c("--model"),        type="character", default="models/cytonorm_unstim_phospho.rds"),
  make_option(c("--outdir"),       type="character", default="normalized_fcs"),
  make_option(c("--cofactor"),     type="double",    default=5)
)))

# Avoid truncation warnings when reading FCS
options("flowCore.truncate.max.range" = FALSE)

tab <- fread(opt$`sample-table`)
need <- c("file","batch","is_anchor","stim")
stopifnot(all(need %in% names(tab)))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# --- Load model & effective marker lists ---
model   <- readRDS(opt$model)
phospho <- if (file.exists("effective_phospho_markers.txt")) readLines("effective_phospho_markers.txt") else character(0)
lineage <- if (file.exists("effective_lineage_markers.txt")) readLines("effective_lineage_markers.txt") else character(0)
if (!length(phospho)) stop("Missing effective_phospho_markers.txt (run training first).")
if (!length(lineage)) stop("Missing effective_lineage_markers.txt (run training first).")

# --- Build transform lists (forward + reverse) ---
# arcsinhTransform is defined as: asinh(a + b*x) + c
# Use a=0, b=1/cofactor, c=0 (default), and supply the analytic inverse
a_tr <- 0
b_tr <- 1/opt$cofactor
c_tr <- 0
tfun   <- flowCore::arcsinhTransform(a = a_tr, b = b_tr, c = c_tr)
tf     <- flowCore::transformList(from = c(lineage, phospho), tfun = tfun)
# Important: provide a plain function so transformList stores it as a callable
inv_fun <- function(x) (sinh(x - c_tr) - a_tr) / b_tr
tf_rev  <- flowCore::transformList(from = c(lineage, phospho), tfun = inv_fun)
# --- Helper: try to recover batch when NA ---
# 1) from filename
recover_batch_from_filename <- function(p) {
  b <- basename(p)
  m <- stringr::str_match(b, "Batch_([0-9]+)")[,2]
  ifelse(is.na(m), NA_character_, m)
}

# 2) from metadata columns in sample_table (anything containing 'batch' except the main 'batch')
meta_batch_cols <- setdiff(grep("batch", names(tab), ignore.case = TRUE, value = TRUE), "batch")
recover_batch_from_meta <- function(i) {
  if (!length(meta_batch_cols)) return(NA_character_)
  vals <- as.character(unlist(tab[i, ..meta_batch_cols]))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return(NA_character_)
  # take the first that looks like a number/id
  cand <- vals[grepl("^[A-Za-z0-9._-]+$", vals)]
  ifelse(length(cand), cand[1], vals[1])
}

# 3) from FCS $DATE: map to nearest anchor day and inherit that anchor's batch
anchors <- tab[is_anchor==TRUE]
get_fcs_date <- function(f) {
  # returns Date or NA
  k <- tryCatch(keyword(read.FCS(f, transformation = FALSE)), error=function(e) list())
  d <- k[["$DATE"]]
  if (is.null(d) || is.na(d)) return(NA)
  # $DATE can be dd/mm/yyyy or mm/dd/yyyy depending on software; be lenient
  d <- gsub("[^0-9/.-]", "", d)
  for (fmt in c("%d/%m/%Y","%m/%d/%Y","%Y-%m-%d","%d-%m-%Y","%m-%d-%Y")) {
    res <- try(as.Date(d, format=fmt), silent=TRUE)
    if (!inherits(res, "try-error") && !is.na(res)) return(res)
  }
  # last resort; still guard against errors and return NA instead
  res <- suppressWarnings(try(as.Date(d), silent=TRUE))
  if (!inherits(res, "try-error") && !is.na(res)) return(res)
  return(NA)
}
anchor_dates <- NULL
if (nrow(anchors)) {
  ad <- lapply(anchors$file, get_fcs_date)
  anchor_dates <- data.table(file=anchors$file, batch=anchors$batch, date=as.Date(unlist(ad)))
  anchor_dates <- anchor_dates[!is.na(date)]
}

recover_batch_from_fcsdate <- function(f) {
  if (is.null(anchor_dates) || !nrow(anchor_dates)) return(NA_character_)
  d <- get_fcs_date(f)
  if (is.na(d)) return(NA_character_)
  # find closest anchor date
  anchor_dates[, dist := abs(as.numeric(d - date))]
  anchor_dates[which.min(dist)]$batch
}

# Apply recovery for NA batches
if (any(is.na(tab$batch))) {
  for (i in which(is.na(tab$batch))) {
    b <- recover_batch_from_filename(tab$file[i])
    if (is.na(b)) b <- recover_batch_from_meta(i)
    if (is.na(b)) b <- recover_batch_from_fcsdate(tab$file[i])
    tab$batch[i] <- b
  }
}

# Still NA? skip with warning
if (any(is.na(tab$batch))) {
  bad <- tab[is.na(batch), file]
  warning("These files have batch=NA and will be skipped:\n - ", paste(basename(bad), collapse="\n - "))
  tab <- tab[!is.na(batch)]
}
stopifnot(nrow(tab) > 0)

# --- Get training label levels from model (several strategies) ---
train_levels <- NULL
try_list <- list(
  function(m) try(names(m@models), silent=TRUE),
  function(m) try(levels(m@transformation$label), silent=TRUE),
  function(m) try(rownames(m@transformation$weights), silent=TRUE)
)
for (f in try_list) {
  x <- f(model)
  if (!inherits(x, "try-error") && length(x)) { train_levels <- as.character(x); break }
}
if (is.null(train_levels) || !length(train_levels)) {
  warning("[WARN] Could not infer training label levels from model; using batches found in sample table.")
  train_levels <- sort(unique(tab$batch))
}

tab[, label := factor(batch, levels = train_levels)]

# skip files from unseen batches
to_skip <- tab[is.na(label)]
if (nrow(to_skip)) {
  warning("Skipping files from unseen batches: ",
          paste(unique(to_skip$batch), collapse=", "),
          "\nAffected files:\n - ", paste(basename(to_skip$file), collapse="\n - "))
  tab <- tab[!is.na(label)]
}
stopifnot(nrow(tab) > 0)

# --- Apply CytoNorm (normalize PHOSPHO only; consistent with training) ---
for (i in seq_len(nrow(tab))) {
  f_in  <- tab$file[i]
  lbl   <- tab$label[i]
  f_out <- file.path(opt$outdir, basename(f_in))
  CytoNorm.normalize(
    model                 = model,
    files                 = f_in,
    labels                = lbl,
    channels              = phospho,
    transformList         = tf,
    transformList.reverse = tf_rev,
    outputDir             = opt$outdir,
    outputSuffix          = "",
    write                 = TRUE,
    verbose               = FALSE
  )
  cat("[ok] normalized ->", f_out, "\n")
}

cat("\n[done] Wrote normalized FCS to: ", opt$outdir, "\n", sep="")
