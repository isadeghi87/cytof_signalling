#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table); library(readr)
  library(flowCore);   library(FlowSOM)
  library(CytoNorm)
})

# ---------- CLI ----------
opt <- parse_args(OptionParser(option_list = list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv"),
  make_option(c("--out-model"),     type="character", default="models/cytonorm_unstim_phospho.rds"),
  make_option(c("--phospho"),       type="character", default=NULL, help="Path to phospho marker list (.txt or .tsv)"),
  make_option(c("--lineage"),       type="character", default=NULL, help="Path to lineage marker list (.txt or .tsv)"),
  make_option(c("--cofactor"),      type="double",    default=5),
  make_option(c("--allow-fallback"),action="store_true", default=FALSE,
              help="If no UNSTIM anchors found, fall back to ALL anchors for training")
)))

sample_table <- opt$`sample-table`
out_model    <- opt$`out-model`
cofactor     <- opt$cofactor
dir.create(dirname(out_model), recursive = TRUE, showWarnings = FALSE)

# Sanity: provide clear error if heavy deps missing
required_pkgs <- c("flowCore", "FlowSOM", "CytoNorm")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) {
  stop(sprintf("Missing required R packages: %s. Please install them (e.g. scripts/R/00_install_pkgs.R).",
               paste(missing_pkgs, collapse=", ")))
}

# ---------- helpers ----------
read_marker_list <- function(path){
  if (is.null(path) || !file.exists(path)) return(character(0))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("txt","lst")) {
    x <- readLines(path, warn=FALSE)
    x <- trimws(x); x <- x[nchar(x)>0]
    return(unique(x))
  } else if (ext %in% c("tsv","csv")) {
    df <- tryCatch(suppressMessages(readr::read_tsv(path, show_col_types = FALSE)),
                   error=function(e) NULL)
    if (!is.null(df)) {
      # prefer 'channel' column if present, else take first column
      if ("channel" %in% names(df)) {
        x <- df$channel
      } else {
        x <- df[[1]]
      }
      x <- trimws(as.character(x)); x <- x[nchar(x)>0]
      return(unique(x))
    } else {
      # fall back to plain lines
      x <- readLines(path, warn=FALSE)
      x <- trimws(x); x <- x[nchar(x)>0]
      return(unique(x))
    }
  } else {
    x <- readLines(path, warn=FALSE)
    x <- trimws(x); x <- x[nchar(x)>0]
    return(unique(x))
  }
}

log_vec <- function(lbl, v){
  if (length(v)) message(sprintf("[%s] %d: %s", lbl, length(v), paste(head(v, 50), collapse=", ")))
  else message(sprintf("[%s] 0", lbl))
}

# ---------- load table ----------
tab <- data.table::fread(sample_table)
need <- c("file","batch","is_anchor")
stopifnot(all(need %in% names(tab)))

# Use stim_class if available; else fall back to exact 'stim' tokens
if ("stim_class" %in% names(tab)) {
  is_unstim <- tolower(tab$stim_class) == "unstim"
} else {
  is_unstim <- !is.na(tab$stim) & grepl("^unstim$", tab$stim, ignore.case = TRUE)
}

# ---------- pick training anchors: UNSTIM only ----------
train <- tab[is_anchor == TRUE & is_unstim == TRUE]
if (nrow(train) == 0) {
  if (!opt$`allow-fallback`) {
    stop("No anchors explicitly labeled as 'unstim' were found. ",
         "Fix stim parsing in prepare_metadata.R or pass --allow-fallback to use all anchors.")
  } else {
    message("[WARN] No explicit UNSTIM anchors; using ALL anchors for training due to --allow-fallback.")
    train <- tab[is_anchor == TRUE]
  }
}
# batches must be covered
unique_batches <- unique(tab$batch[!is.na(tab$batch)])
missing_batches <- setdiff(unique_batches, unique(train$batch))
if (length(missing_batches)) {
  stop("These batches lack a training anchor: ", paste(missing_batches, collapse=", "))
}
message("[OK] Training anchors per batch:")
print(train[, .N, by=batch][order(batch)])

# ---------- markers ----------
# load phospho list (priority: user-supplied -> channels_phospho*.tsv -> channels_to_adjust.tsv)
phospho <- character(0)
cands <- c(opt$phospho, "channels_phospho3.tsv", "channels_phospho.tsv", "phospho_markers.txt", "channels_to_adjust.tsv")
for (p in cands) {
  if (length(phospho)==0 && !is.null(p) && file.exists(p)) phospho <- read_marker_list(p)
}
if (length(phospho)==0) stop("No phospho marker list found. Provide --phospho or ensure channels_phospho*.tsv exists.")

# read panel from first training FCS
panel <- colnames(read.FCS(train$file[1], transformation = FALSE))

phospho <- intersect(unique(phospho), panel)
if (length(phospho) < 3) stop("Too few phospho markers after intersection with panel. Got: ", paste(phospho, collapse=", "))

# lineage
if (!is.null(opt$lineage) && file.exists(opt$lineage)) {
  lineage <- read_marker_list(opt$lineage)
  lineage <- intersect(unique(lineage), panel)
} else {
  # auto-derive lineage: everything not phospho minus obvious instrument/aux channels
  forbidden <- c("Time","Event_length","Center","Offset","Residual","BCKG190Di","ArAr80Di","Xe131Di","Xe134Di","I127Di","Ir191Di","Ir193Di")
  cand <- setdiff(panel, phospho)
  cand <- cand[!cand %in% forbidden]
  lineage <- head(cand, 12)
}
if (length(lineage) < 6) stop("Too few lineage markers (", length(lineage), "). Provide --lineage with a proper list.")

log_vec("LINEAGE (for FlowSOM)", lineage)
log_vec("PHOSPHO (to normalize)", phospho)

# ---------- screen phosphos (drop unstable) ----------
screen_phosphos <- function(files, chans, cof=5) {
  meds <- list(); iqrs <- list()
  for (f in files) {
    M <- exprs(read.FCS(f, transformation=FALSE))
    M <- asinh(M / cof)
    keep <- intersect(colnames(M), chans)
    meds[[basename(f)]] <- apply(M[, keep, drop=FALSE], 2, median, na.rm=TRUE)
    iqrs[[basename(f)]] <- apply(M[, keep, drop=FALSE], 2, IQR,    na.rm=TRUE)
  }
  med_df <- do.call(rbind, lapply(meds, \(v) data.frame(t(v), check.names=FALSE)))
  iqr_df <- do.call(rbind, lapply(iqrs, \(v) data.frame(t(v), check.names=FALSE)))
  med_sd <- apply(med_df, 2, sd, na.rm=TRUE)
  iqr_mu <- colMeans(iqr_df, na.rm=TRUE)
  keep <- names(which(iqr_mu > 0.02 & med_sd < 0.35))  # thresholds
  list(keep = keep, drop = setdiff(chans, keep))
}
scr <- screen_phosphos(train$file, phospho, cof = cofactor)
if (length(scr$drop)) message("[Drop] Noisy phospho channels removed: ", paste(scr$drop, collapse=", "))
phospho <- scr$keep
if (length(phospho) < 3) stop("After screening, <3 phospho markers remain. Loosen thresholds or review panel.")

# ---------- train ----------
train$label <- factor(train$batch)
fs_params <- list(xdim=10, ydim=10, nClus=50, scale=FALSE, nCells=20000, colsToUse=lineage)

set.seed(1)
model <- CytoNorm.train(
  files          = train$file,
  labels         = train$label,
  channels       = phospho,                                # normalize PHOSPHO only
  transformList  = flowCore::transformList(from = c(lineage, phospho),
                                           tfun = flowCore::arcsinhTransform(a=0, b=1/cofactor)),
  FlowSOM.params = fs_params,                              # cluster on LINEAGE
  normParams     = list(nQ=101, goal="mean", limitQuantile=c(0.01, 0.99)),
  seed           = 1
)

saveRDS(model, out_model)
writeLines(lineage, "effective_lineage_markers.txt")
writeLines(phospho, "effective_phospho_markers.txt")
cat("[ok] Trained model saved ->", out_model, "\n")
