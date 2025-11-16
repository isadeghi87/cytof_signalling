#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

set.seed(1)

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
pm <- fread(pm_path)

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

## helper to average across selected times
avg_times <- function(dt, tvec) {
  dt[time_min %in% tvec,
     .(score = mean(logFC, na.rm = TRUE)),
     by = .(cytokine, cluster)]
}

early_scores <- avg_times(pm_fc, c(2,5))
setnames(early_scores, "score", "early_activation")

late_scores  <- avg_times(pm_fc, c(30,60))
setnames(late_scores, "score", "late_activation")

## AUC across all time points per marker
pm_fc[ , time_num := as.numeric(time_min)]
auc_list <- pm_fc[!is.na(time_num),
                  {
                    ord <- order(time_num)
                    t_ord <- time_num[ord]
                    v_ord <- logFC[ord]
                    auc <- if (length(t_ord) > 1) {
                      sum(diff(t_ord) * (head(v_ord, -1) + tail(v_ord, -1)) / 2)
                    } else NA_real_
                    .(AUC = auc)
                  },
                  by = .(cytokine, cluster, marker)]

auc_tot <- auc_list[, .(AUC_total = sum(abs(AUC), na.rm = TRUE)),
                    by = .(cytokine, cluster)]

## merge summaries
fs <- merge(early_scores, late_scores,
            by = c("cytokine","cluster"), all = TRUE)
fs <- merge(fs, auc_tot, by = c("cytokine","cluster"), all = TRUE)

## marker groups for phenotypes
## Use channel names actually present in phos_medians_by_sample_cluster.csv
## STAT3 and STAT5 channels are Er166Di and Yb173Di in this panel.
## For a simple TCR-proximal proxy we use Er170Di, Eu153Di and Yb172Di,
## which correspond to early signalling / proximal markers in this panel.
tcr_markers   <- c("Er170Di","Eu153Di","Yb172Di")
stat5_markers <- c("Yb173Di")
stat3_markers <- c("Er166Di")

score_group <- function(dt, markers) {
  m_int <- intersect(markers, dt$marker)
  if (!length(m_int)) return(NA_real_)
  dt[marker %in% m_int, mean(logFC, na.rm = TRUE)]
}

## compute group scores at peak (e.g. 10 min)
peak_fc <- pm_fc[time_min %in% c(10)]

phen_list <- list()
for (cy in sort(unique(peak_fc$cytokine))) {
  for (cl in sort(unique(peak_fc$cluster))) {
    sub <- peak_fc[cytokine == cy & cluster == cl]
    if (!nrow(sub)) next
    stat5 <- score_group(sub, stat5_markers)
    stat3 <- score_group(sub, stat3_markers)
    tcr   <- score_group(sub, tcr_markers)
    phen_list[[length(phen_list)+1]] <- data.table(
      cytokine = cy, cluster = cl,
      STAT5_peak = stat5,
      STAT3_peak = stat3,
      TCR_peak   = tcr
    )
  }
}

phen_dt <- if (length(phen_list)) rbindlist(phen_list, use.names = TRUE, fill = TRUE) else
  data.table()

fs <- merge(fs, phen_dt, by = c("cytokine","cluster"), all = TRUE)

## phenotype categories
if (nrow(fs)) {
  fs[, phenotype := {
    s5  <- STAT5_peak
    s3  <- STAT3_peak
    tcr <- TCR_peak
    if (is.na(s5) & is.na(s3) & is.na(tcr)) {
      "Weak/non-responder"
    } else if (!is.na(s5) & !is.na(s3) & s5 > 0 & s3 > 0) {
      "Dual-STAT"
    } else if (!is.na(s5) & (is.na(s3) | s5 > s3) & s5 > 0) {
      "STAT5-biased"
    } else if (!is.na(s3) & (is.na(s5) | s3 > s5) & s3 > 0) {
      "STAT3-biased"
    } else if (!is.na(tcr) & tcr > 0) {
      "TCR-proximal"
    } else {
      "Weak/non-responder"
    }
  }, by = .(cytokine, cluster)]
}

fwrite(fs, "cluster_functional_annotations.csv")

## functional signatures heatmap (clusters × features)
if (nrow(fs)) {
  mat_fs <- dcast(fs,
                  cluster ~ cytokine,
                  value.var = "AUC_total",
                  fill = 0)
  rn <- mat_fs$cluster
  mat <- as.matrix(mat_fs[, -1])
  rownames(mat) <- rn
  png("functional_signatures.png", width = 1200, height = 800)
  pheatmap(mat, main = "Cluster functional signatures (AUC total)")
  dev.off()
}

## cytokine × cluster comparison radar-like plot
## summarise peak stats
rad <- fs[, .(STAT3_peak = STAT3_peak,
              STAT5_peak = STAT5_peak,
              TCR_peak   = TCR_peak),
          by = .(cytokine, cluster)]

rad_long <- melt(rad,
                 id.vars = c("cytokine","cluster"),
                 variable.name = "feature",
                 value.name = "value")

gp <- ggplot(rad_long,
             aes(x = feature, y = value,
                 group = interaction(cytokine, cluster),
                 colour = cytokine)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ cluster, ncol = 4) +
  coord_polar() +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 14),
    axis.text  = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) +
  labs(title = "Cytokine × cluster comparison",
       y = "peak logFC", x = "feature")

ggsave("cytokine_cluster_comparison.png", gp,
       width = 18, height = 12, dpi = 300)

cat("Functional signatures computed.\n")
