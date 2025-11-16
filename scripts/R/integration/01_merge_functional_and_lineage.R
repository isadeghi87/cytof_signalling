#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(1)

## Required inputs
req_files <- c(
  "functional_cluster_annotations.csv",
  "cluster_activation_rankings.csv",
  "flowsom_results/cluster_labels_template.csv",
  "flowsom_results/metacluster_medians_asinh.csv",
  "lineage_channels.tsv"
)

missing <- req_files[!file.exists(req_files)]
if (length(missing)) {
  stop("Missing required input files: ", paste(missing, collapse = ", "))
}

fun <- fread("functional_cluster_annotations.csv")
act <- fread("cluster_activation_rankings.csv")
cl_lab <- fread("flowsom_results/cluster_labels_template.csv")
mc_med <- fread("flowsom_results/metacluster_medians_asinh.csv")
lineage_ch <- fread("lineage_channels.tsv")$channel

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) {
  stop("Could not find phos_medians_by_sample_cluster.csv in diff_phos_results/ or phospho_medians/.")
}
pm <- fread(pm_path)

## identify phospho markers
id_cols <- c("file_base","cluster","time_min","cytokine")
phos_cols <- setdiff(names(pm), id_cols)

## compute logFC vs 0 per cytokine/cluster/time/marker
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

## marker groups per phenotype
stat3_markers  <- c("pSTAT3","Er166Di","Er170Di")
stat5_markers  <- c("pSTAT5","Yb173Di")
tcr_markers    <- c("pZAP70","pSLP76","pPLCy","pNFATc1","pS6","pERK")
metab_markers  <- c("pAkt","pGSK3","pAMPK","pPGC1","pCREB")
stress_markers <- c("IkBa")

group_for_pheno <- function(ph) {
  if (grepl("STAT3", ph)) stat3_markers else
  if (grepl("STAT5", ph)) stat5_markers else
  if (grepl("TCR", ph))   tcr_markers   else
  if (grepl("Metabolic", ph)) metab_markers else
  if (grepl("Stress", ph)) stress_markers else
    phos_cols
}

## lineage highlights per metacluster (top 3 lineage markers by median)
lineage_markers <- intersect(lineage_ch, colnames(mc_med))
lin_list <- list()
for (cl in mc_med$cluster) {
  row <- mc_med[cluster == cl]
  vals <- as.numeric(row[, ..lineage_markers])
  names(vals) <- lineage_markers
  vals <- sort(vals, decreasing = TRUE)
  top3 <- head(names(vals), 3)
  lin_list[[length(lin_list)+1]] <- data.table(
    cluster = cl,
    lineage_top_markers = paste(top3, collapse = ",")
  )
}
lin_dt <- rbindlist(lin_list)

## merge activation totals
setnames(act, "activation_score", "activation_score_total")

## compute activation_score_early and dominant_markers & best_timepoint
res_rows <- list()

for (i in seq_len(nrow(fun))) {
  row <- fun[i]
  cy  <- row$cytokine
  cl  <- row$cluster
  ph  <- row$phenotype
  grp_markers <- group_for_pheno(ph)

  sub_fc <- pm_fc[cytokine == cy & cluster == cl & marker %in% grp_markers]
  # early activation (2,5 min)
  early <- sub_fc[time_min %in% c(2,5),
                  .(activation_score_early = mean(abs(logFC), na.rm = TRUE))]
  if (!nrow(early)) early[, activation_score_early := NA_real_]

  # dominant markers: top 3 by mean |logFC| in early window (fall back to all times)
  dom_tbl <- sub_fc[time_min %in% c(2,5)]
  if (!nrow(dom_tbl)) dom_tbl <- sub_fc
  if (nrow(dom_tbl)) {
    dom_summary <- dom_tbl[, .(m = mean(abs(logFC), na.rm = TRUE)), by = marker]
    dom_summary <- dom_summary[order(-m)]
    dom_markers <- paste(head(dom_summary$marker, 3), collapse = ",")
  } else {
    dom_markers <- NA_character_
  }

  # best timepoint: max mean logFC across group markers
  bt_tbl <- sub_fc[, .(mean_logFC = mean(logFC, na.rm = TRUE)),
                   by = time_min][order(-mean_logFC)]
  best_tp <- if (nrow(bt_tbl)) bt_tbl$time_min[1] else NA

  # activation_score_total from act
  act_row <- act[cytokine == cy & cluster == cl]
  act_tot <- if (nrow(act_row)) act_row$activation_score_total[1] else NA_real_

  res_rows[[length(res_rows)+1]] <- data.table(
    cytokine = cy,
    cluster_id = cl,
    phenotype = ph,
    activation_score_early = early$activation_score_early[1],
    activation_score_total = act_tot,
    dominant_markers = dom_markers,
    best_timepoint = best_tp
  )
}

sig_dt <- rbindlist(res_rows, use.names = TRUE, fill = TRUE)

## add manual labels and metacluster info
setnames(cl_lab, c("cluster","label"), c("cluster_id","manual_label"))
sig_dt <- merge(sig_dt, cl_lab, by = "cluster_id", all.x = TRUE)

## metacluster_id (here equal to FlowSOM metacluster id)
sig_dt[, metacluster_id := cluster_id]

## add lineage_top_markers
sig_dt <- merge(sig_dt, lin_dt,
                by.x = "metacluster_id", by.y = "cluster",
                all.x = TRUE)

setcolorder(sig_dt,
            c("cytokine","cluster_id","metacluster_id","manual_label",
              "phenotype","activation_score_early","activation_score_total",
              "dominant_markers","best_timepoint","lineage_top_markers"))

fwrite(sig_dt, "cluster_signaling_summary.csv")

cat("cluster_signaling_summary.csv written.\n")

