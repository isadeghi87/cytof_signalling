#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(flowCore); library(ggplot2)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("--sample-table"), type="character"),
  make_option(c("--normalized-dir"), type="character", default="normalized_fcs"),
  make_option(c("--out"), type="character", default="summary")
)))
dir.create(opt$out, recursive=TRUE, showWarnings=FALSE)
tab <- data.table::fread(opt$`sample-table`)
channels <- unique(readr::read_tsv("channels_to_adjust.tsv", col_names=TRUE)$channel); channels <- channels[nchar(channels)>0]
file_medians <- function(f){
  ff <- tryCatch(read.FCS(f, transformation=FALSE, alter.names=TRUE), error=function(e) NULL)
  if (is.null(ff)) return(NULL)
  M <- exprs(ff)
  med <- sapply(channels, function(ch) if (ch %in% colnames(M)) median(M[,ch], na.rm=TRUE) else NA_real_)
  data.table(file=basename(f), channel=channels, median=as.numeric(med))
}
pre_dt <- rbindlist(lapply(tab$file, file_medians), fill=TRUE)
post_dt <- rbindlist(lapply(file.path(opt$`normalized-dir`, basename(tab$file)), file_medians), fill=TRUE)
pre_dt[, status:="pre"]; post_dt[, status:="post"]
all_dt <- rbind(pre_dt, post_dt, fill=TRUE)
fwrite(all_dt, file.path(opt$out, "channel_medians.tsv"), sep="\t")
rle_dt <- all_dt[, .(value = median - median(median, na.rm=TRUE)), by=.(status, channel, file)]
p <- ggplot(rle_dt, aes(channel, value, fill=status)) + geom_boxplot(outlier.size=0.5) + coord_flip() +
  theme_bw() + labs(title="RLE-style centered medians across files", y="relative median")
ggsave(file.path(opt$out, "RLE_boxplot.png"), p, width=8, height=10, dpi=200)
cat("[ok] Summary written to ", opt$out, "\n")
