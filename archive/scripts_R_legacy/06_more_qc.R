#!/usr/bin/env Rscript

## NOTE (archived):
## Additional QC routines developed during method exploration.
## Not part of the streamlined pipeline; kept for completeness.
## More QC: pre/post medians, CVs, event counts, heatmaps, density overlays, UMAP
suppressPackageStartupMessages({
  library(optparse)
  library(flowCore)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(uwot)
})

option_list <- list(
  make_option(c("--raw-dir"), type="character", default="batch/", help="directory with raw FCS files"),
  make_option(c("--norm-dir"), type="character", default="normalized_fcs/", help="directory with normalized FCS files"),
  make_option(c("--channels"), type="character", default="channels_to_adjust.tsv", help="channels to check (one per line)") ,
  make_option(c("--outdir"), type="character", default="summary/more_qc/", help="output directory for QC results"),
  make_option(c("--subsample"), type="integer", default=3000, help="per-file subsample for UMAP/density")
)

opt <- parse_args(OptionParser(option_list=option_list))

# map hyphen-option names to clean variable names
raw_dir <- opt[["raw-dir"]]
norm_dir <- opt[["norm-dir"]]
outdir <- opt$outdir

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

read_channels <- function(chfile){
  if(!file.exists(chfile)) stop("channels file not found: ", chfile)
  ch <- scan(chfile, what=character(), sep="\n", quiet=TRUE)
  ch <- trimws(ch)
  ch[ch!=""]
}

channels <- read_channels(opt$channels)

list_fcs <- function(d){
  if(!dir.exists(d)) return(character())
  fs <- list.files(d, pattern="\\.fcs$", full.names=TRUE, recursive=FALSE)
  fs
}

raw_files <- list_fcs(raw_dir)
norm_files <- list_fcs(norm_dir)

basename_key <- function(path) {
  bn <- basename(path)
  # strip common patterns like Batch_..._ prefix to allow matching
  bn
}

match_norm <- function(raw_path){
  bn <- basename_key(raw_path)
  # try exact basename match
  cand <- norm_files[basename(norm_files) == bn]
  if(length(cand)) return(cand[1])
  # try matching by sample id token (first or last few tokens)
  token <- sub("\\.fcs$", "", bn)
  cand <- norm_files[grepl(token, basename(norm_files), fixed=TRUE)]
  if(length(cand)) return(cand[1])
  return(NA_character_)
}

read_medians <- function(fpath, channels){
  if(is.na(fpath) || !file.exists(fpath)) return(NULL)
  f <- tryCatch(read.FCS(fpath, transformation=FALSE, truncate_max_range=FALSE), error=function(e) NULL)
  if(is.null(f)) return(NULL)
  expr <- exprs(f)
  present <- intersect(channels, colnames(expr))
  if(length(present)==0) return(NULL)
  med <- apply(expr[, present, drop=FALSE], 2, median, na.rm=TRUE)
  ev <- nrow(expr)
  list(medians=med, events=ev)
}

pre_list <- list()
post_list <- list()
rows <- list()

for(r in raw_files){
  key <- basename_key(r)
  res_pre <- read_medians(r, channels)
  res_post <- read_medians(match_norm(r), channels)
  rows[[key]] <- list(raw=r, norm=match_norm(r), pre=res_pre, post=res_post)
}

# build summary tables
build_table <- function(rows, which="pre"){
  df_list <- list()
  for(k in names(rows)){
    item <- rows[[k]][[which]]
    if(is.null(item)) next
    med <- item$medians
    ev <- item$events
    r <- data.table(file=k, events=ev)
    for(ch in names(med)) r[[ch]] <- med[[ch]]
    df_list[[k]] <- r
  }
  if(length(df_list)==0) return(NULL)
  dt <- rbindlist(df_list, fill=TRUE)
  setcolorder(dt, c("file","events", setdiff(names(dt), c("file","events"))))
  dt
}

dt_pre <- build_table(rows, "pre")
dt_post <- build_table(rows, "post")

fwrite(dt_pre, file.path(opt$outdir, "medians_pre.tsv"), sep="\t")
fwrite(dt_post, file.path(opt$outdir, "medians_post.tsv"), sep="\t")

if(!is.null(dt_pre) && !is.null(dt_post)){
  # align columns
  common_ch <- intersect(colnames(dt_pre), colnames(dt_post))
  common_ch <- setdiff(common_ch, c("file","events"))
  if(length(common_ch)>0){
    # compute delta (post - pre)
    setkey(dt_pre, file); setkey(dt_post, file)
    merged <- merge(dt_pre, dt_post, by="file", suffixes=c("_pre","_post"), all=TRUE)
    delta <- data.table(file=merged$file)
    for(ch in common_ch){
      delta[[ch]] <- merged[[paste0(ch,"_post")]] - merged[[paste0(ch,"_pre")]]
    }
    fwrite(delta, file.path(opt$outdir, "medians_delta.tsv"), sep="\t")

    # heatmap of median shifts (channels x samples)
    mat <- as.matrix(t(delta[, ..common_ch]))
    colnames(mat) <- delta$file
    # handle NA/Inf: replace non-finite with 0 but warn
    if(!all(is.finite(mat))){
      warning("Non-finite values found in median delta matrix; replacing with 0 for heatmap")
      mat[!is.finite(mat)] <- 0
    }
    pdf(file.path(outdir, "median_shifts_heatmap.pdf"), width=10, height=6)
    pheatmap(mat, main="Median (post - pre)", fontsize=8)
    dev.off()

    # per-channel boxplot of shifts
    mlong <- melt(delta, id.vars="file", variable.name="channel", value.name="delta")
    png(file.path(opt$outdir, "median_shifts_boxplot.png"), width=1400, height=700)
    ggplot(mlong, aes(x=channel, y=delta)) + geom_boxplot() + theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Per-channel median shifts (post - pre)")
    dev.off()
  }
}

# per-sample event counts pre/post
if(!is.null(dt_pre)){
  ev_pre <- dt_pre[, .(file, events)]
  fwrite(ev_pre, file.path(opt$outdir, "events_pre.tsv"), sep="\t")
}
if(!is.null(dt_post)){
  ev_post <- dt_post[, .(file, events)]
  fwrite(ev_post, file.path(opt$outdir, "events_post.tsv"), sep="\t")
}

# density overlays for top channels by variance
if(!is.null(dt_pre)){
  chan_cols <- setdiff(colnames(dt_pre), c("file","events"))
  if(length(chan_cols)>0){
    # choose top 6 channels by median absolute deviation across files (pre)
    mad_vals <- sapply(chan_cols, function(ch) mad(dt_pre[[ch]], na.rm=TRUE))
    nchoose <- min(6, length(mad_vals))
    if(nchoose > 0){
      topch <- names(sort(mad_vals, decreasing=TRUE))[seq_len(nchoose)]
    } else {
      topch <- character()
    }

    # create densities subdir for per-channel PNGs
    dir.create(file.path(outdir, "densities"), recursive=TRUE, showWarnings=FALSE)

    # use cairo_pdf for more robust multi-page PDF creation
    pdf_path <- file.path(outdir, "density_overlays_top_channels.pdf")
    cairo_opened <- FALSE
    try({
      cairo_pdf(pdf_path, width=10, height=8)
      cairo_opened <- TRUE
      for(ch in topch){
      # aggregate a subsample of events per-file for raw and normalized, then plot combined density
      dt_vals <- list()
      for(r in raw_files){
        # read raw
        fraw <- tryCatch(read.FCS(r, transformation=FALSE, truncate_max_range=FALSE), error=function(e) NULL)
        if(!is.null(fraw) && ch %in% colnames(exprs(fraw))){
          e <- exprs(fraw)[, ch]
          if(length(e) > opt$subsample){
            e <- sample(e, opt$subsample)
          }
          dt_vals[[length(dt_vals)+1]] <- data.table(value = as.numeric(e), source = "raw", sample = basename(r))
        }
        # read normalized counterpart if exists
        rn <- match_norm(r)
        if(!is.na(rn) && file.exists(rn)){
          fn <- tryCatch(read.FCS(rn, transformation=FALSE, truncate_max_range=FALSE), error=function(e) NULL)
          if(!is.null(fn) && ch %in% colnames(exprs(fn))){
            e2 <- exprs(fn)[, ch]
            if(length(e2) > opt$subsample){
              e2 <- sample(e2, opt$subsample)
            }
            dt_vals[[length(dt_vals)+1]] <- data.table(value = as.numeric(e2), source = "norm", sample = basename(rn))
          }
        }
      }
      if(length(dt_vals)==0){
        # write a simple placeholder PNG and a blank PDF page
        png(file.path(outdir, "densities", paste0(gsub("[^A-Za-z0-9_\\-]", "_", ch), "_empty.png")), width=800, height=500)
        plot.new(); title(main=paste("No data for channel:", ch))
        dev.off()
        plot.new(); title(main=paste("No data for channel:", ch))
        next
      }
      dt_all <- rbindlist(dt_vals, use.names=TRUE, fill=TRUE)
      # cap total points to avoid huge plotting: sample down if >100000
      if(nrow(dt_all) > 100000){
        dt_all <- dt_all[sample(.N, 100000)]
      }

      # transform values (arcsinh with cofactor 5), compute medians per source
      cofactor <- 5
      dt_all[, value_asinh := asinh(value / cofactor)]
      medians_src <- dt_all[, .(median = median(value_asinh, na.rm=TRUE)), by=source]

      # determine sensible x-limits from 1st-99th percentiles after transform
      q <- quantile(dt_all$value_asinh, probs=c(0.01, 0.99), na.rm=TRUE)
      xlim_low <- as.numeric(q[1]); xlim_high <- as.numeric(q[2])

      p <- ggplot(dt_all, aes(x = value_asinh, color = source, fill = source)) +
        geom_density(alpha=0.25, linewidth=0.6) +
        theme_minimal(base_size=12) +
        ggtitle(paste0("Density overlay: ", ch, " (asinh/", cofactor, ")")) +
        xlab(paste0(ch, " (asinh/", cofactor, ")")) + ylab("Density") +
        xlim(c(xlim_low, xlim_high)) +
        theme(axis.text.x = element_text(angle=0, hjust=0.5),
              legend.position = "top", legend.title = element_blank())

      # add dashed vertical median lines per source
      p <- p + geom_vline(data=medians_src, aes(xintercept=median, color=source), linetype="dashed", linewidth=0.6)
      # save PNG high-res
      safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", ch)
      png_file <- file.path(outdir, "densities", paste0(safe_name, ".png"))
      tryCatch({
        ggsave(filename=png_file, plot=p, width=10, height=6, dpi=200)
      }, error=function(e){
        warning("Failed to ggsave PNG for ", ch, ": ", e$message)
      })
      # print to the open PDF device
      print(p)
    }
    }, silent=FALSE)
    if(cairo_opened) dev.off()
  }
}

cat("More QC run completed. Outputs written to:", opt$outdir, "\n")
