#!/usr/bin/env Rscript
# Compare pre- and post-normalization densities by sample groups (P1..P4)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(flowCore)
  library(ggplot2)
  library(gridExtra)
  library(stringr)
  library(scales)
})

option_list <- list(
  make_option(c("--sample-table"), type="character", default="sample_table.csv"),
  make_option(c("--pre-dir"), type="character", default="batch", help="Directory with original FCS files"),
  make_option(c("--post-dir"), type="character", default="normalized_fcs", help="Directory with normalized FCS files (primary)"),
  make_option(c("--post-dir2"), type="character", default=NULL, help="Optional second post dir (e.g. K=30) to compare side-by-side"),
  make_option(c("--channels"), type="character", default="channels_to_adjust.tsv"),
  make_option(c("--out"), type="character", default="summary/normalization_efficiency.pdf"),
  make_option(c("--cofactor"), type="double", default=0, help="Arcsinh cofactor (0 = no transform)") ,
  make_option(c("--n-events"), type="integer", default=3000),
  make_option(c("--mode"), type="character", default="aggregated", help="aggregated|per-sample overlay mode")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$`sample-table`)) stop('sample_table.csv missing')
tab <- fread(opt$`sample-table`)

if (!file.exists(opt$channels)) stop('channels file missing')
ch_df <- suppressMessages(try(readr::read_tsv(opt$channels, show_col_types = FALSE), silent = TRUE))
if (inherits(ch_df, "try-error")) {
  channels <- as.character(fread(opt$channels, header=FALSE)[[1]])
} else {
  channels <- if ("channel" %in% names(ch_df)) as.character(ch_df$channel) else as.character(ch_df[[1]])
}
channels_to_plot <- head(channels, 8)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# helper to resolve pre/post paths for a sample row
resolve_paths <- function(basename_or_path) {
  # if path exists, use it; else try to find in pre-dir or post-dir
  b <- basename(basename_or_path)
  pre <- if (file.exists(basename_or_path)) basename_or_path else file.path(opt$`pre-dir`, b)
  post <- file.path(opt$`post-dir`, b)
  list(pre = if (file.exists(pre)) pre else NA_character_, post = if (file.exists(post)) post else NA_character_)
}

# determine group P1..P4 from filename
get_group <- function(fname) {
  b <- basename(fname)
  m <- str_match(b, "_P(\\d)")
  if (!is.na(m[1,2])) return(paste0('P', m[1,2]))
  # fallback: look for P1_ or P1
  m2 <- str_match(b, "P(\\d)")
  if (!is.na(m2[1,2])) return(paste0('P', m2[1,2]))
  return('other')
}

read_sample_subset <- function(path, channels, n, cofactor = 0) {
  if (!file.exists(path)) return(NULL)
  f <- try(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE), silent = TRUE)
  if (inherits(f, 'try-error')) return(NULL)
  expr <- as.data.frame(exprs(f))
  present <- intersect(channels, colnames(expr))
  if (length(present) == 0) return(NULL)
  n0 <- min(nrow(expr), n)
  idx <- if (nrow(expr) > n0) sample(seq_len(nrow(expr)), n0) else seq_len(nrow(expr))
  dt <- as.data.table(expr[idx, present, drop = FALSE])
  # optional arcsinh transform
  if (!is.null(cofactor) && cofactor > 0) {
    for (c in present) dt[[c]] <- asinh(dt[[c]]/cofactor)
  }
  dt
}

plots <- list()

groups <- unique(sapply(tab$file, get_group))
groups <- intersect(groups, c('P1','P2','P3','P4','other'))

pdf(opt$out, width = 11, height = 8.5)
for (g in groups) {
  rows <- tab[sapply(file, function(x) get_group(x) == g)]
  if (nrow(rows) == 0) next
  # for each channel, make side-by-side pre/post overlay
  for (ch in channels_to_plot) {
    pre_list <- list(); post_list <- list(); post2_list <- list(); labels <- c()
    for (i in seq_len(nrow(rows))) {
      r <- rows$file[i]
      rp <- resolve_paths(r)
      labels <- c(labels, basename(r))
      pre_dt <- read_sample_subset(rp$pre, ch, opt$`n-events`, opt$cofactor)
      post_dt <- read_sample_subset(rp$post, ch, opt$`n-events`, opt$cofactor)
      post2_dt <- NULL
      if (!is.null(opt$`post-dir2`) && nzchar(opt$`post-dir2`)) {
        rp2 <- file.path(opt$`post-dir2`, basename(r))
        post2_dt <- read_sample_subset(rp2, ch, opt$`n-events`, opt$cofactor)
      }
      pre_list[[i]] <- if (!is.null(pre_dt)) data.table(value = pre_dt[[ch]], sample = basename(r), which = 'pre') else NULL
      post_list[[i]] <- if (!is.null(post_dt)) data.table(value = post_dt[[ch]], sample = basename(r), which = 'post') else NULL
      post2_list[[i]] <- if (!is.null(post2_dt)) data.table(value = post2_dt[[ch]], sample = basename(r), which = 'post2') else NULL
    }
    pre_dt_all <- rbindlist(pre_list, use.names = TRUE, fill = TRUE)
    post_dt_all <- rbindlist(post_list, use.names = TRUE, fill = TRUE)
    post2_dt_all <- rbindlist(post2_list, use.names = TRUE, fill = TRUE)
    if (nrow(pre_dt_all) == 0 && nrow(post_dt_all) == 0) next

    if (tolower(opt$mode) == 'aggregated') {
      # Overlay pooled pre vs post with medians and shared x-limits (1–99% quantiles)
      pre_dt_all$status <- 'pre'; post_dt_all$status <- 'post'
      agg <- rbind(pre_dt_all[,.(value)], post_dt_all[,.(value)])
      q <- quantile(agg$value, probs = c(0.01, 0.99), na.rm = TRUE)
      med <- rbind(
        data.table(status='pre',  median = median(pre_dt_all$value,  na.rm=TRUE)),
        data.table(status='post', median = median(post_dt_all$value, na.rm=TRUE))
      )
      p_overlay <- ggplot(rbind(pre_dt_all[,.(value,status='pre')], post_dt_all[,.(value,status='post')]),
                          aes(x=value, fill=status, color=status)) +
        geom_density(alpha=0.35, size=0.7, adjust=1.1) +
        geom_vline(data=med, aes(xintercept=median, color=status), linetype='dashed', linewidth=0.6, show.legend=FALSE) +
        coord_cartesian(xlim=c(q[1], q[2])) +
        scale_fill_manual(values=c(pre='grey65', post='#1b9e77')) +
        scale_color_manual(values=c(pre='grey40', post='#1b9e77')) +
        theme_minimal(base_size=12) + theme(legend.position='top') +
        labs(title=sprintf('%s - %s (asinh/%g) pooled', g, ch, opt$cofactor), x=ch, y='Density')
      grid.arrange(p_overlay, ncol = 1)
    } else {
      # Original per-sample side-by-side panels with consistent colors
      sample_names <- unique(c(as.character(pre_dt_all$sample), as.character(post_dt_all$sample), as.character(post2_dt_all$sample)))
      sample_names <- sample_names[!is.na(sample_names)]
      cols <- if(length(sample_names)>0) setNames(hue_pal()(length(sample_names)), sample_names) else NULL
      p1 <- ggplot(pre_dt_all, aes(x = value, color = sample, fill = sample)) + geom_density(alpha = 0.25) + theme_minimal() + labs(title = sprintf('%s - %s - pre', g, ch))
      p2 <- ggplot(post_dt_all, aes(x = value, color = sample, fill = sample)) + geom_density(alpha = 0.25) + theme_minimal() + labs(title = sprintf('%s - %s - post', g, ch))
      if(!is.null(cols)) { p1 <- p1 + scale_color_manual(values=cols) + scale_fill_manual(values=cols) + theme(legend.position='top', legend.title=element_blank()) }
      if(!is.null(cols)) { p2 <- p2 + scale_color_manual(values=cols) + scale_fill_manual(values=cols) + theme(legend.position='top', legend.title=element_blank()) }
      if (nrow(post2_dt_all) > 0) {
        p3 <- ggplot(post2_dt_all, aes(x = value, color = sample, fill = sample)) + geom_density(alpha = 0.25) + theme_minimal() + labs(title = sprintf('%s - %s - post2', g, ch))
        if(!is.null(cols)) { p3 <- p3 + scale_color_manual(values=cols) + scale_fill_manual(values=cols) + theme(legend.position='top', legend.title=element_blank()) }
        grid.arrange(p1, p2, p3, ncol = 3)
      } else {
        grid.arrange(p1, p2, ncol = 2)
      }
    }
  }
}
dev.off()
message('[ok] normalization efficiency plots written to ', opt$out)
