#!/usr/bin/env Rscript
options(warn=1)
message("[setup] Installing packages...")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
pkgs_bioc <- c("flowCore","FlowSOM","CytoNorm","PeacoQC")
pkgs_cran <- c("data.table","readr","stringr","ggplot2","RColorBrewer","scales","cowplot","optparse")
for (p in pkgs_bioc) {
  if (!suppressWarnings(requireNamespace(p, quietly=TRUE))) {
    BiocManager::install(p, ask=FALSE, update=FALSE)
  }
}
for (p in pkgs_cran) {
  if (!suppressWarnings(requireNamespace(p, quietly=TRUE))) {
    install.packages(p, dependencies=TRUE, repos="https://cloud.r-project.org")
  }
}
message("[setup] Done.")
