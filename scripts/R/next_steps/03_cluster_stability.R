#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

pl <- fread("phenograph_results/phenograph_labels.csv")
dim_meta <- readRDS("dimred_flowsom/dimred_metadata.rds")
dim_dt <- as.data.table(dim_meta)
dim_dt[, flowsom_cluster := as.integer(as.character(cluster))]

## if lengths match, compute full contingency and ARI/NMI
if (nrow(pl) == nrow(dim_dt)) {
  pl[, phenograph_cluster := cluster]
  pl[, cluster := NULL]
  full <- cbind(dim_dt, pl)
  tab <- full[, .N, by = .(flowsom_cluster, phenograph_cluster)]
  fwrite(tab, "FlowSOM_Phenograph_contingency.csv")

  fs <- full$flowsom_cluster
  ph <- full$phenograph_cluster

  ari_simple <- function(x, y) {
    tab <- table(x, y)
    n <- sum(tab)
    sum_comb <- sum(choose(tab, 2))
    sum_row  <- sum(choose(rowSums(tab), 2))
    sum_col  <- sum(choose(colSums(tab), 2))
    expected <- sum_row * sum_col / choose(n, 2)
    max_val  <- (sum_row + sum_col) / 2
    (sum_comb - expected) / (max_val - expected)
  }

  ari_val <- ari_simple(fs, ph)

  ## simple NMI
  nmi_simple <- function(x, y) {
    tab <- table(x, y)
    px <- rowSums(tab) / sum(tab)
    py <- colSums(tab) / sum(tab)
    pxy <- tab / sum(tab)
    nz <- pxy > 0
    mi <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]])))
    hx <- -sum(px * log(px))
    hy <- -sum(py * log(py))
    mi / sqrt(hx * hy)
  }
  nmi_val <- nmi_simple(fs, ph)

  rep <- list(
    n_events = nrow(full),
    ari = ari_val,
    nmi = nmi_val
  )
  writeLines(jsonlite::toJSON(rep, auto_unbox = TRUE, pretty = TRUE),
             "cluster_stability_report.txt")
} else {
  ## length mismatch – cannot compute true overlap
  fs_counts <- dim_dt[, .N, by = flowsom_cluster]
  ph_counts <- pl[, .N, by = cluster]
  fs_counts[, label_type := "FlowSOM"]
  setnames(fs_counts, "flowsom_cluster", "cluster")
  ph_counts[, label_type := "PhenoGraph"]
  all_counts <- rbind(fs_counts, ph_counts, fill = TRUE)
  fwrite(all_counts, "FlowSOM_Phenograph_contingency.csv")

  rep <- list(
    message = "Row counts of FlowSOM-labelled sample (dimred_metadata) and phenograph_labels differ; cannot compute ARI/NMI reliably.",
    n_flowsom = nrow(dim_dt),
    n_phenograph = nrow(pl),
    suggestion = "Regenerate phenograph_labels with the same event subset used for dimred_flowsom, or export per-event FlowSOM labels for the PhenoGraph subset."
  )
  writeLines(jsonlite::toJSON(rep, auto_unbox = TRUE, pretty = TRUE),
             "cluster_stability_report.txt")
}

cat("Cluster stability analysis completed.\n")

