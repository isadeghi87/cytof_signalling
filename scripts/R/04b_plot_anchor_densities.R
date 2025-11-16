#!/usr/bin/env Rscript
# Plot density overlays for anchor samples per channel

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv", help="Path to sample_table.csv"),
  make_option(c("--anchors"), type="character", default="templates/anchors.txt", help="Anchors list (one per line) or 'auto' to use sample_table is_anchor column"),
  make_option(c("--channels"), type="character", default="channels_to_adjust.tsv", help="Channels file (one channel per line or tab file)") ,
  make_option(c("--outdir"), type="character", default="qc_plots", help="Output dir for plots"),
  make_option(c("--n-events"), type="integer", default=5000, help="Max events to sample per anchor")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!file.exists(opt$`sample-table`)) stop("sample_table.csv not found: ", opt$`sample-table`)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

tab <- fread(opt$`sample-table`)

anchors_list <- NULL
if (!is.null(opt$anchors) && file.exists(opt$anchors)) {
  anchors_list <- fread(opt$anchors, header = FALSE)$V1
} else if ("is_anchor" %in% names(tab)) {
  anchors_list <- tab[file == file & (is_anchor==TRUE | is_anchor==1), file]
} else {
  stop("No anchors found: pass --anchors or ensure 'is_anchor' column in sample_table.csv")
}

channels <- NULL
if (file.exists(opt$channels)) {
  # try as single-column or tab-separated
  ch <- fread(opt$channels, header = FALSE)
  channels <- as.character(ch[[1]])
} else stop("Channels file not found: ", opt$channels)

# helper to read FCS safely
read_fcs_exprs <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace('flowCore', quietly = TRUE)) stop('flowCore required to read FCS')
  f <- try(flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE), silent = TRUE)
  if (inherits(f, 'try-error')) return(NULL)
  as.data.frame(flowCore::exprs(f))
}

libraryNamespaceOK <- function(pkg) requireNamespace(pkg, quietly = TRUE)
if (!libraryNamespaceOK('flowCore')) stop('Please install flowCore to run anchor density plotting')

plot_dt <- list()
for (anc in anchors_list) {
  # resolve path: sample_table may contain full path or basename; try several fallbacks
  candidate <- NULL
  if (file.exists(anc)) candidate <- anc
  else {
    # try matching basename in sample_table
    b <- basename(anc)
    cand_row <- tab[basename(file) == b]
    if (nrow(cand_row) > 0 && file.exists(cand_row$file[1])) candidate <- cand_row$file[1]
  }
  if (is.null(candidate)) {
    message('[WARN] anchor not found on disk: ', anc)
    next
  }
  df <- read_fcs_exprs(candidate)
  if (is.null(df)) { message('[WARN] unable to read: ', candidate); next }
  # pick channels present
  present <- intersect(channels, colnames(df))
  if (length(present) == 0) { message('[WARN] no channels present in ', candidate); next }
  # sample rows
  set.seed(1)
  n <- min(nrow(df), opt$`n-events`)
  idx <- if (nrow(df) > n) sample(seq_len(nrow(df)), n) else seq_len(nrow(df))
  sub <- df[idx, present, drop = FALSE]
  long <- as.data.table(sub)
  long[, sample := basename(candidate)]
  long[, rowid := .I]
  long_melt <- melt(long, id.vars = c('sample','rowid'), measure.vars = present, variable.name = 'channel', value.name = 'value')
  plot_dt[[basename(candidate)]] <- long_melt
}

if (length(plot_dt) == 0) stop('No anchor data available to plot')

all_dt <- rbindlist(plot_dt)

plt <- ggplot(all_dt, aes(x = value, color = sample)) +
  geom_density(size=0.6, alpha=0.8) +
  facet_wrap(~channel, scales = 'free') +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(title = 'Anchor density overlays', x = 'Intensity', y = 'Density')

outf <- file.path(opt$outdir, 'anchor_densities.png')
ggsave(outf, plt, width = 14, height = 10, dpi = 150)
message('[ok] Anchor density plot saved to ', outf)

## end of script
