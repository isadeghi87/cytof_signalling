#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(flowCore); library(ggplot2); library(scales); library(cowplot); library(RColorBrewer)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option(c("--sample-table"), type="character"),
  make_option(c("--normalized-dir"), type="character", default="normalized_fcs"),
  make_option(c("--outdir"), type="character", default="qc_plots"),
  make_option(c("--cofactor"), type="double", default=5),
  make_option(c("--channels"), type="character", default="channels_to_adjust.tsv",
              help="Path to channels list (TSV with column 'channel' or single-column file)")
)))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
tab <- data.table::fread(opt$`sample-table`)
ch_df <- suppressMessages(try(readr::read_tsv(opt$channels, show_col_types = FALSE), silent = TRUE))
if (inherits(ch_df, "try-error")) {
  ch <- readLines(opt$channels, warn = FALSE)
} else {
  if ("channel" %in% names(ch_df)) ch <- ch_df$channel else ch <- ch_df[[1]]
}
channels <- unique(trimws(as.character(ch)))
channels <- channels[nchar(channels) > 0]
sample_rows <- tab[sample(.N, min(.N, 12))]
dens_df <- function(f){
  ff <- tryCatch(read.FCS(f, transformation=FALSE, alter.names=TRUE), error=function(e) NULL)
  if (is.null(ff)) return(NULL)
  M <- exprs(ff)
  rbindlist(lapply(channels, function(ch){
    if (!(ch %in% colnames(M))) return(NULL)
    data.table(channel=ch, value=asinh(M[,ch]/opt$cofactor), file=basename(f))
  }), fill=TRUE)
}

# sample more files to stabilize densities (up to 24)
ns <- min(nrow(tab), 24)
set.seed(1)
sample_rows <- tab[sample(.N, ns)]

pre <- rbindlist(lapply(sample_rows$file, dens_df), fill=TRUE); if (nrow(pre)==0) stop("No pre-normalization data read.")
post <- rbindlist(lapply(file.path(opt$`normalized-dir`, basename(sample_rows$file)), dens_df), fill=TRUE); if (nrow(post)==0) stop("No post-normalization data read.")
pre[, status:="pre"]; post[, status:="post"]; df <- rbind(pre, post, fill=TRUE)

for (ch in unique(df$channel)) {
  d <- df[channel==ch]
  # robust x-limits from pooled 1% and 99% quantiles
  q <- quantile(d$value, probs=c(0.01, 0.99), na.rm=TRUE)
  med <- d[, .(median = median(value, na.rm=TRUE)), by = status]
  p <- ggplot(d, aes(value, fill=status, color=status)) +
    geom_density(adjust=1.1, alpha=0.35, size=0.7, na.rm=TRUE) +
    geom_vline(data=med, aes(xintercept=median, color=status), linetype="dashed", linewidth=0.6, show.legend=FALSE) +
    coord_cartesian(xlim = c(q[1], q[2])) +
    theme_bw(base_size = 12) +
    labs(title=paste0("Density (asinh/", opt$cofactor, "): ", ch), x="asinh intensity", y="Density") +
    scale_fill_manual(values=c("pre"="grey65","post"="#1b9e77")) +
    scale_color_manual(values=c("pre"="grey40","post"="#1b9e77")) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))
  ggsave(file.path(opt$outdir, paste0("density_", gsub("[^A-Za-z0-9]+","_", ch), ".png")), p, width=8, height=5, dpi=300)
}
cat("[ok] QC plots saved in ", opt$outdir, "\n")
