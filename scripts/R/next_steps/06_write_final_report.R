#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

## basic summaries
integrity <- if (file.exists("integrity_summary.txt"))
  readLines("integrity_summary.txt") else character()

fun <- if (file.exists("cluster_functional_annotations.csv"))
  fread("cluster_functional_annotations.csv") else NULL

lim_full <- if (file.exists("limma_full_results.csv"))
  fread("limma_full_results.csv") else NULL

act <- if (file.exists("cluster_activation_rankings.csv"))
  fread("cluster_activation_rankings.csv") else NULL

lines <- c(
  "# Final Report v2",
  "",
  "## 1. QC and Metadata Integrity",
  "",
  "- Sample and cytokine/time coverage are summarized in `integrity_summary.txt`.",
  "- Diagnostic coverage plots are in `diagnostic_plots/`.",
  "",
  "### Integrity summary (excerpt)",
  ""
)

lines <- c(lines, paste0("- ", integrity))

lines <- c(lines,
  "",
  "## 2. Differential Signaling (limma)",
  "",
  if (!is.null(lim_full)) {
    cy <- paste(sort(unique(lim_full$cytokine)), collapse = ", ")
    c(
      paste0("- Full limma results: `limma_full_results.csv` (cytokines: ", cy, ")."),
      "- Time contrasts: 5, 10, 15, 30, 60 minutes vs 0 per cytokine.",
      "- Cytokine contrasts: IL2–CD3, IL15–CD3, IL21–CD3 at matched times."
    )
  } else {
    "- limma_full_results.csv not found."
  },
  "",
  "## 3. Functional Cluster Phenotypes",
  ""
)

if (!is.null(fun) && nrow(fun)) {
  phen_tab <- fun[, .N, by = phenotype]
  lines <- c(lines,
             "- Cluster functional annotations stored in `cluster_functional_annotations.csv`.",
             "- Functional signature heatmap: `functional_signatures.png`.",
             "",
             "### Phenotype counts",
             "",
             paste0("- ", apply(phen_tab, 1,
                                function(r) paste0(r["phenotype"], ": ", r["N"], " entries"))))
}

lines <- c(lines,
  "",
  "## 4. Cytokine Comparisons",
  "",
  "- Direct cytokine contrast heatmaps are in `cytokine_contrast_heatmaps/`.",
  "- `cytokine_cluster_comparison.png` summarises STAT3/STAT5/TCR peaks per cluster.",
  "",
  "## 5. Trajectories and Diffusion Maps",
  "",
  "- GAM-based trajectories and AUC rankings: `trajectory_plots/` and `cluster_activation_rankings.csv`.",
  "- Diffusion-map-like embeddings per cytokine: `diffusion_map_CD3.png`, `diffusion_map_IL2.png`, `diffusion_map_IL15.png`, `diffusion_map_IL21.png`.",
  "",
  "## 6. Cluster Stability (FlowSOM vs PhenoGraph)",
  "",
  "- Contingency information: `FlowSOM_Phenograph_contingency.csv`.",
  "- Stability report JSON: `cluster_stability_report.txt`.",
  "",
  "## 7. Summary Figure Book",
  "",
  "- Multi-page figure compilation: `MASTER_FIGUREBOOK.pdf`.",
  "",
  "## 8. Key Biological Highlights (to be filled by user)",
  "",
  "- Review STAT3/STAT5-biased clusters and link to lineage identities from `cluster_labels_template.csv`.",
  "- Inspect cytokine-specific contrasts to identify preferential responders.",
  "- Use activation rankings to prioritise clusters for downstream validation."
)

writeLines(lines, "FINAL_REPORT_v2.md")

cat("FINAL_REPORT_v2.md written.\n")

