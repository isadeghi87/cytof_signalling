#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(stringr)
})

# Check for heavy runtime packages and provide a clear error if missing. Installing
# packages requires internet access; if running in an air-gapped environment, ask
# the user to install the binaries manually or run `scripts/R/00_install_pkgs.R`.
required_pkgs <- c("flowCore", "FlowSOM", "CytoNorm")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) {
  stop(sprintf("Missing required R packages: %s.\nPlease install them (e.g. run 'Rscript scripts/R/00_install_pkgs.R' with internet access) or install from a local repository.", paste(missing_pkgs, collapse=", ")))
}
suppressPackageStartupMessages({
  library(flowCore); library(FlowSOM); library(CytoNorm)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("--sample-table"), type="character"),
  make_option(c("--channels-tsv"), type="character", default="channels_to_adjust.tsv"),
  make_option(c("--cofactor"), type="double", default=5),
  make_option(c("--n-clusters"), type="integer", default=10),
  make_option(c("--out-model"), type="character", default="models/cytonorm_model.rds")
)))
dir.create(dirname(opt$`out-model`), recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(opt$`sample-table`), file.exists(opt$`channels-tsv`))
tab <- data.table::fread(opt$`sample-table`)
if (!"is_anchor" %in% names(tab)) stop("'is_anchor' column missing in sample-table; run 01_prepare_metadata.R")
anchors <- tab[is_anchor==TRUE]
if (nrow(anchors) == 0) stop("No anchor files flagged in sample-table.")
channels <- unique(readr::read_tsv(opt$`channels-tsv`, col_names=TRUE)$channel); channels <- channels[nchar(channels)>0]
read_ff <- function(p){ suppressWarnings(read.FCS(p, transformation=FALSE, alter.names=TRUE)) }
panel <- colnames(read_ff(anchors$file[1]))
miss <- setdiff(channels, panel)
if (length(miss)) { warning("Missing channels in anchor FCS: ", paste(miss, collapse=", ")); channels <- setdiff(channels, miss) }
if (length(channels) < 3) stop("Too few valid channels after filtering.")
labels <- anchors$batch; if (any(is.na(labels))) labels <- factor(anchors$plate); if (any(is.na(labels))) labels <- factor(basename(anchors$file))
set.seed(1L)
model <- CytoNorm::CytoNorm.train(
  files = anchors$file,
  labels = labels,
  channels = channels,
  transformList = flowCore::transformList(from=channels, tfun=flowCore::arcsinhTransform(a=0, b=1/opt$cofactor)),
  # nCells controls how many cells FlowSOM samples per file; set a reasonable
  # default (10k) to ensure clustering is stable while keeping runtime bounded.
  FlowSOM.params = list(xdim=7, ydim=7, nClus=opt$`n-clusters`, scale=FALSE, nCells=10000),
  # Use the newer CytoNorm API: normParams is a list controlling quantile
  # normalization; CytoNorm::CytoNorm.train expects e.g. list(nQ=99).
  normParams = list(nQ=101),
  seed = 1
)
saveRDS(model, file=opt$`out-model`)
cat("[ok] Trained CytoNorm model ->", opt$`out-model`, "\n")
