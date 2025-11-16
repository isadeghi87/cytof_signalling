#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(1)

if (!file.exists("cluster_signaling_summary.csv")) {
  stop("cluster_signaling_summary.csv not found; run 01_merge_functional_and_lineage.R first.")
}

sum_dt <- fread("cluster_signaling_summary.csv")

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) {
  stop("Could not find phos_medians_by_sample_cluster.csv in diff_phos_results/ or phospho_medians/.")
}
pm <- fread(pm_path)

id_cols <- c("file_base","cluster","time_min","cytokine")
phos_cols <- setdiff(names(pm), id_cols)

pm_long <- melt(pm,
                id.vars = id_cols,
                measure.vars = phos_cols,
                variable.name = "marker",
                value.name = "value")

pm_mean <- pm_long[, .(value = mean(value, na.rm = TRUE)),
                   by = .(cytokine, cluster, time_min, marker)]

base <- pm_mean[time_min %in% c("0",0),
                .(baseline = mean(value, na.rm = TRUE)),
                by = .(cytokine, cluster, marker)]

pm_fc <- merge(pm_mean, base,
               by = c("cytokine","cluster","marker"),
               all.x = TRUE)
pm_fc[, logFC := value - baseline]

stat3_markers  <- c("pSTAT3","Er166Di","Er170Di")
stat5_markers  <- c("pSTAT5","Yb173Di")
tcr_markers    <- c("pZAP70","pSLP76","pPLCy","pNFATc1","pS6","pERK")

score_group_total <- function(dt, markers) {
  m_int <- intersect(markers, dt$marker)
  if (!length(m_int)) return(NA_real_)
  dt[marker %in% m_int,
     sum(abs(logFC), na.rm = TRUE)]
}

clusters <- sort(unique(pm_fc$cluster))
cytokines <- sort(unique(pm_fc$cytokine))

res_bias <- list()

for (cl in clusters) {
  for (grp_name in c("STAT3","STAT5","TCR")) {
    scores <- numeric(length(cytokines))
    names(scores) <- cytokines
    for (cy in cytokines) {
      sub <- pm_fc[cluster == cl & cytokine == cy]
      if (!nrow(sub)) {
        scores[cy] <- NA_real_
        next
      }
      markers <- switch(grp_name,
                        STAT3 = stat3_markers,
                        STAT5 = stat5_markers,
                        TCR   = tcr_markers)
      scores[cy] <- score_group_total(sub, markers)
    }
    # choose best cytokine
    if (all(is.na(scores))) {
      best <- NA_character_
      note <- "flat"
    } else {
      ord <- order(scores, decreasing = TRUE, na.last = NA)
      nonzero <- scores[ord]
      if (length(nonzero) == 0 || max(nonzero, na.rm = TRUE) < 1e-3) {
        best <- NA_character_
        note <- "flat"
      } else {
        best <- names(nonzero)[1]
        # simple description
        if (sum(nonzero > 1e-3, na.rm = TRUE) == 1) {
          note <- paste0(best, " only")
        } else {
          note <- paste(names(nonzero), collapse = " > ")
        }
      }
    }
    res_bias[[length(res_bias)+1]] <- data.table(
      cluster_id = cl,
      pathway = grp_name,
      best_cytokine = best,
      notes = note
    )
  }
}

bias_dt <- rbindlist(res_bias, use.names = TRUE, fill = TRUE)

## reshape to wide per cluster
bias_wide <- dcast(bias_dt,
                   cluster_id ~ pathway,
                   value.var = "best_cytokine")
setnames(bias_wide,
         c("STAT3","STAT5","TCR"),
         c("best_cytokine_STAT3","best_cytokine_STAT5","best_cytokine_TCR"))

## combine notes for STAT3 group as main notes
notes_stat3 <- bias_dt[pathway == "STAT3",
                       .(notes = notes[1]), by = cluster_id]

out <- merge(bias_wide, notes_stat3, by = "cluster_id", all.x = TRUE)

## add manual_label from cluster_signaling_summary (unique per cluster)
lab_dt <- unique(sum_dt[, .(cluster_id, manual_label)])
out <- merge(out, lab_dt, by = "cluster_id", all.x = TRUE)

setcolorder(out,
            c("cluster_id","manual_label",
              "best_cytokine_STAT3","best_cytokine_STAT5","best_cytokine_TCR",
              "notes"))

fwrite(out, "cluster_cytokine_bias.csv")

cat("cluster_cytokine_bias.csv written.\n")

