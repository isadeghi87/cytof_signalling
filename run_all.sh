#!/usr/bin/env bash
set -euo pipefail
FILES=${1:-/home/isadeghi/projects/cyTof/cyTof/batch/files.txt}
ANCHORS=${2:-/home/isadeghi/projects/cyTof/cyTof/batch/anchors/anchors.txt}
CHANNELS=${3:-/home/isadeghi/projects/cyTof/cyTof/batch/ChannelsToAdjust.txt}
META=${4:-/home/isadeghi/projects/cyTof/cyTof/batch/OMIQ_metadata-20250919_LS_pTCR_cHCV.csv}
BASEDIR=${5:-/home/isadeghi/projects/cyTof/cyTof/batch}
#Rscript scripts/R/00_install_pkgs.R
Rscript scripts/R/01_prepare_metadata.R --files-list "$FILES" --anchors-list "$ANCHORS" --channels "$CHANNELS" --metadata "$META" --base-dir "$BASEDIR" --out sample_table.csv
Rscript scripts/R/02_train_cytonorm.R --sample-table sample_table.csv --out-model models/cytonorm_model.rds
Rscript scripts/R/03_apply_cytonorm.R --sample-table sample_table.csv --model models/cytonorm_model.rds --outdir normalized_fcs
## Consolidate CytoNorm outputs to canonical basenames (pick largest candidate)
Rscript scripts/R/03b_consolidate_normalized.R sample_table.csv normalized_fcs || true
# archive leftover temporary Norm__ files (keep backups for traceability)
mkdir -p normalized_fcs/backups_norm || true
shopt -s nullglob || true
for f in normalized_fcs/Norm__*; do mv "$f" normalized_fcs/backups_norm/ || true; done

# QC and extra plots
Rscript scripts/R/04b_plot_anchor_densities.R --sample-table sample_table.csv --anchors auto --channels channels_to_adjust.tsv --outdir qc_plots || true
Rscript scripts/R/04c_tsne_umap_clusters.R --sample-table sample_table.csv --channels channels_to_adjust.tsv --outdir qc_plots --n-samples 20 --n-events 1000 --k 12 || true
Rscript scripts/R/05_plot_normalization_efficiency.R --sample-table sample_table.csv --pre-dir batch --post-dir normalized_fcs --channels channels_to_adjust.tsv --out summary/normalization_efficiency.pdf --n-events 2000 --cofactor 5 || true
Rscript scripts/R/04_qc_plots.R --sample-table sample_table.csv --normalized-dir normalized_fcs --outdir qc_plots
Rscript scripts/R/05_export_summary.R --sample-table sample_table.csv --normalized-dir normalized_fcs --out summary
echo "[ok] Finished. See normalized_fcs/, qc_plots/, summary/"
