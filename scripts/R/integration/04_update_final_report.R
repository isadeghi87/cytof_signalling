#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(1)

## Required inputs
req_files <- c(
  "cluster_signaling_summary.csv",
  "cluster_cytokine_bias.csv",
  "cluster_activation_rankings.csv"
)
missing <- req_files[!file.exists(req_files)]
if (length(missing)) {
  stop("Missing required input files: ", paste(missing, collapse = ", "))
}

sig <- fread("cluster_signaling_summary.csv")
bias <- fread("cluster_cytokine_bias.csv")
act  <- fread("cluster_activation_rankings.csv")

lines <- c(
  "# Final Report v3",
  "",
  "## 1. Overview",
  "",
  "This report summarises functional signaling phenotypes, cytokine preferences,",
  "and compact summary figures derived from cluster-level phospho medians.",
  "",
  "Key new outputs:",
  "- `cluster_signaling_summary.csv`",
  "- `cluster_cytokine_bias.csv`",
  "- `summary_figures/UMAP_by_phenotype.png`",
  "- `summary_figures/UMAP_by_cytokine_preference.png`",
  "- `summary_figures/activation_bars_by_cytokine.png`",
  "- `summary_figures/key_trajectories_STAT3_STAT5.png`",
  ""
)

## Phenotype distribution
phen_tab <- sig[, .N, by = phenotype][order(-N)]

lines <- c(lines,
  "## 2. Functional Phenotypes per Cytokine",
  "",
  "### Phenotype counts (all cytokines combined)",
  ""
)

lines <- c(lines,
  paste0("- ", phen_tab$phenotype, ": ", phen_tab$N, " cluster–cytokine pairs"))

## brief per-cytokine summary: most frequent phenotypes
cy_phen <- sig[, .N, by = .(cytokine, phenotype)][order(cytokine, -N)]
lines <- c(lines,
  "",
  "### Dominant phenotypes per cytokine",
  ""
)
for (cy in unique(cy_phen$cytokine)) {
  sub <- cy_phen[cytokine == cy]
  top <- head(sub, 3)
  desc <- paste(paste0(top$phenotype, " (", top$N, ")"), collapse = "; ")
  lines <- c(lines,
             paste0("- ", cy, ": ", desc))
}

## Cytokine bias notes
lines <- c(lines,
  "",
  "## 3. Cytokine Preferences per Cluster",
  "",
  "Based on integrated logFC for STAT3/STAT5/TCR marker groups,",
  "each cluster has a preferred cytokine for each pathway.",
  "",
  "Examples (first 10 clusters):",
  ""
)

head_bias <- head(bias[order(cluster_id)], 10)
for (i in seq_len(nrow(head_bias))) {
  r <- head_bias[i]
  lines <- c(lines,
    sprintf("- Cluster %s (%s): STAT3 → %s; STAT5 → %s; TCR → %s; notes: %s",
            r$cluster_id,
            ifelse(is.na(r$manual_label), "unlabelled", r$manual_label),
            r$best_cytokine_STAT3,
            r$best_cytokine_STAT5,
            r$best_cytokine_TCR,
            r$notes))
}

## Activation highlights
lines <- c(lines,
  "",
  "## 4. Most Responsive Clusters",
  "",
  "Using `cluster_activation_rankings.csv`, the following clusters show highest",
  "integrated activation (AUC) per cytokine:",
  ""
)

act_top <- act[, .SD[order(-activation_score)][1:3], by = cytokine]
for (cy in unique(act_top$cytokine)) {
  sub <- act_top[cytokine == cy]
  desc <- paste(paste0("cluster ", sub$cluster, " (", round(sub$activation_score,1), ")"),
                collapse = "; ")
  lines <- c(lines, paste0("- ", cy, ": ", desc))
}

## Figure references
lines <- c(lines,
  "",
  "## 5. Figure References",
  "",
  "- `summary_figures/UMAP_by_phenotype.png`: UMAP coloured by functional phenotype.",
  "- `summary_figures/UMAP_by_cytokine_preference.png`: UMAP coloured by best STAT3 cytokine.",
  "- `summary_figures/activation_bars_by_cytokine.png`: AUC per cluster×cytokine, annotated with top responders.",
  "- `summary_figures/key_trajectories_STAT3_STAT5.png`: STAT3/STAT5 time courses for top responder clusters.",
  ""
)

## Limitations
lines <- c(lines,
  "## 6. Limitations",
  "",
  "- FlowSOM vs PhenoGraph stability is only partially assessed:",
  "  the available PhenoGraph labels are not perfectly row-aligned with the",
  "  FlowSOM-labelled event set, so ARI/NMI metrics are approximate at best.",
  "- All signaling summaries are based on cluster medians, not raw single-cell",
  "  trajectories; sub-cluster heterogeneity or rare subpopulations may be",
  "  smoothed out.",
  "- Functional phenotypes and cytokine preferences are derived from a limited",
  "  marker panel and early time-point contrasts; additional markers or time",
  "  points could refine these assignments.",
  "- No batch-specific or subject-level random effects are modeled; all",
  "  statistics are aggregate across samples.",
  ""
)

writeLines(lines, "FINAL_REPORT_v3.md")

cat("FINAL_REPORT_v3.md written.\n")

