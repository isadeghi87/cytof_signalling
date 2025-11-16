#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(flowCore); library(CytoNorm); library(stringr)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("--sample-table"), type="character"),
  make_option(c("--model"), type="character"),
  make_option(c("--outdir"), type="character", default="normalized_fcs"),
  make_option(c("--cofactor"), type="double", default=5)
)))
stopifnot(file.exists(opt$`sample-table`), file.exists(opt$model))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
tab <- data.table::fread(opt$`sample-table`)
channels <- unique(readr::read_tsv("channels_to_adjust.tsv", col_names=TRUE)$channel)
channels <- channels[nchar(channels)>0]
model <- readRDS(opt$model)
for (i in seq_len(nrow(tab))) {
  f <- tab$file[i]; lbl <- tab$batch[i]; if (is.na(lbl)) lbl <- tab$plate[i]
  out_file <- file.path(opt$outdir, basename(f))
  tryCatch({
    CytoNorm::CytoNorm.normalize(
      model = model,
      files = f,
      labels = lbl,
      channels = channels,
      transformList = flowCore::transformList(from=channels, tfun=flowCore::arcsinhTransform(a=0, b=1/opt$cofactor)),
      outputDir = opt$outdir,
      outputSuffix = "",
      write = TRUE, verbose = FALSE
    )
    if (!file.exists(out_file)) {
      # CytoNorm may write many files per input (e.g. suffixed with _fsomN)
      # and with a 'Norm__' prefix that encodes the original path. Find
      # candidates whose basename contains the original basename and pick
      # the largest one (heuristic: likely the full normalized FCS), then
      # rename it to the expected out_file.
      all_cand <- list.files(opt$outdir, full.names=TRUE)
      cand <- all_cand[vapply(all_cand, function(x) grepl(basename(f), basename(x), fixed=TRUE), logical(1))]
      if (length(cand)) {
        sizes <- file.info(cand)$size
        # pick the largest candidate
        sel <- cand[which.max(sizes)]
        if (file.exists(sel)) {
          file.rename(sel, out_file)
        }
      }
    }
    cat("[ok] normalized ->", out_file, "\n")
  }, error=function(e){ warning("[skip] ", f, " : ", e$message) })
}
