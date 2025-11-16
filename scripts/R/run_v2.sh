#!/bin/bash

set -euo pipefail

echo "-----------------------------------------------------"
echo "CyTOF PIPELINE STARTED"
echo "-----------------------------------------------------"

LOGDIR="logs"
mkdir -p ${LOGDIR}

# Activate your R env if needed
# source ~/.bashrc
# conda activate r-omics

# 1. CYTONORM NORMALIZATION
echo "[1] Normalization with CytoNorm"
Rscript 01_run_cytonorm.R \
  --sample-table sample_table.csv \
  --files-list files.txt \
  --base-dir batch \
  --channels channels_to_adjust.tsv \
  --outdir normalized_fcs \
  --qcdir qc_plots \
  > ${LOGDIR}/01_run_cytonorm.log 2>&1

#echo "[1] Completed."


# 2. FLOWSOM + METACLUSTERS
echo "[2] FlowSOM clustering + MST + heatmaps"
Rscript scripts/R/02_flowsom_annotate.R \
  --norm-dir normalized_fcs \
  --sample-table sample_table.csv \
  --outdir flowsom_results \
  --metak 20 \
  > ${LOGDIR}/02_flowsom_annotate.log 2>&1

echo "[2] Completed."


# 3. DIMENSIONALITY REDUCTION WITH FLOWSOM CLUSTERS
echo "[3] t-SNE + UMAP using FlowSOM clusters"
Rscript scripts/R/03_dimred_flowsom.R \
  --norm-dir normalized_fcs \
  --sample-table sample_table.csv \
  --flowsom-rds flowsom_results/flowsom_fitted.rds \
  --cluster-labels flowsom_results/cluster_labels_template.csv \
  --outdir dimred_flowsom \
  > ${LOGDIR}/03_dimred_flowsom.log 2>&1

echo "[3] Completed."


# 4. DEDICATED ADVANCED t-SNE + UMAP SET (SCRIPT 4)
echo "[4] Extra dimensionality reduction plots"
Rscript scripts/R/04_tns_umap_plots.R \
  --norm-dir normalized_fcs \
  --sample-table sample_table.csv \
  --flowsom-rds flowsom_results/flowsom_fitted.rds \
  --cluster-labels flowsom_results/cluster_labels_template.csv \
  --outdir tsne_umap_plots \
  > ${LOGDIR}/04_tsne_umap_plots.log 2>&1

echo "[4] Completed."


# 5. PHOSPHO MEDIANS
echo "[5] Compute phospho medians per cluster × sample"
Rscript scripts/R/05_compute_phospho_medians.R \
  --norm-dir normalized_fcs \
  --sample-table sample_table.csv \
  --flowsom-rds flowsom_results/flowsom_fitted.rds \
  --outdir phospho_medians \
  > ${LOGDIR}/05_compute_phospho_medians.log 2>&1

echo "[5] Completed."


# 6. DIFFERENTIAL SIGNALING
echo "[6] Differential phospho-signaling (limma)"
Rscript scripts/R/06_diff_phospho_analysis.R \
  --medians phospho_medians/phos_medians_by_sample_cluster.csv \
  --outdir diff_phos_results \
  > ${LOGDIR}/06_diff_phospho_analysis.log 2>&1

echo "[6] Completed."


# 7. PHENOGRAPH ALTERNATIVE CLUSTERING
echo "[7] PhenoGraph clustering"
Rscript scripts/R/07_phenograph_clustering.R \
  --norm-dir normalized_fcs \
  --sample-table sample_table.csv \
  --outdir phenograph_results \
  > ${LOGDIR}/07_phenograph_clustering.log 2>&1

echo "[7] Completed."


# 8. RIVERPLOTS
echo "[8] Riverplots for signaling trajectories"
Rscript scripts/R/08_riverplots_signaling.R \
  --medians diff_phos_results/phos_medians_by_sample_cluster.csv \
  --outdir riverplots \
  > ${LOGDIR}/08_riverplots_signaling.log 2>&1

echo "[8] Completed."


echo "-----------------------------------------------------"
echo "CyTOF PIPELINE COMPLETED SUCCESSFULLY"
echo "All logs saved in ./logs/"
echo "-----------------------------------------------------"
