#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt_list <- list(
  make_option("--medians", type="character",
              default="diff_phos_results/phos_medians_by_sample_cluster.csv"),
  make_option("--outdir", type="character", default="riverplots")
)

opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
# Load medians and restrict to stimulated cytokines
###############################################################################

dt <- fread(opt$medians)

# keep only stimulated conditions and rows with valid time
stim_keep <- c("CD3","IL2","IL15","IL21","PAC")
dt <- dt[cytokine %in% stim_keep]
dt <- dt[!is.na(time_min)]

phos_cols <- setdiff(colnames(dt), c("file_base","cluster","time_min","cytokine"))

if (!length(phos_cols) || nrow(dt) == 0) {
  message("No suitable phospho medians for time‑course plots; exiting.")
  quit(save = "no")
}

###############################################################################
# Time‑course plots per cytokine (cluster × marker vs time)
###############################################################################

cyts <- sort(unique(dt$cytokine))

for (cy in cyts) {
  sub <- dt[cytokine == cy]

  # require at least 2 time points and 2 clusters
  if (uniqueN(sub$time_min) < 2 || uniqueN(sub$cluster) < 2) {
    message("[WARN] Skipping cytokine ", cy,
            " for time‑course plot (insufficient time points or clusters).")
    next
  }

  long <- melt(
    sub,
    id.vars = c("file_base","cluster","time_min","cytokine"),
    measure.vars = phos_cols,
    variable.name = "marker",
    value.name = "value"
  )

  # aggregate over files: median per cluster × time × marker
  long_agg <- long[, .(value = median(value, na.rm = TRUE)),
                   by = .(cluster, time_min, cytokine, marker)]

  p <- ggplot(long_agg,
              aes(x = time_min, y = value,
                  group = cluster, colour = factor(cluster))) +
    geom_line(alpha = 0.7) +
    geom_point(size = 0.8) +
    facet_wrap(~ marker, scales = "free_y") +
    theme_bw() +
    labs(title = paste("Phospho time‑courses:", cy),
         x = "time (min)",
         y = "asinh median",
         colour = "FlowSOM cluster")

  ggsave(file.path(opt$outdir, paste0("timecourses_", cy, ".png")),
         p, width = 10, height = 6)
}

cat("Time‑course phospho plots written to: ", opt$outdir, "\n")
