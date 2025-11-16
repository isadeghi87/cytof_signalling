# CyTOF Anchor‑Based Normalization (CytoNorm) — Pro Pack

This pack ingests your file lists, anchor references, channel list to adjust, and metadata,
then performs **anchor‑based batch normalization** with **CytoNorm** + **FlowSOM**, optional
**PeacoQC** cleaning, and exports normalized FCS files with extensive QC plots.

## Inputs
- `files.txt` — list of FCS files. Supports either full paths (one per line) **or** `ls -l` style rows.
- `anchors/anchors.txt` — list of anchor FCS paths (one per line). Ideally one anchor per batch.
- `ChannelsToAdjust.txt` — channels to normalize (one per line), e.g. metal channel names.
- `metadata.csv` — sample annotations: at least `sample_id` and `group` (batch optional).

## Quick start
```bash
# 0) install packages (one-time)
Rscript scripts/R/00_install_pkgs.R

# 1) prepare metadata table
Rscript scripts/R/01_prepare_metadata.R   --files-list /home/isadeghi/projects/cyTof/cyTof/batch/files.txt   --anchors-list /home/isadeghi/projects/cyTof/cyTof/batch/anchors/anchors.txt   --channels /home/isadeghi/projects/cyTof/cyTof/batch/ChannelsToAdjust.txt   --metadata /home/isadeghi/projects/cyTof/cyTof/batch/OMIQ_metadata-20250919_LS_pTCR_cHCV.csv   --base-dir /home/isadeghi/projects/cyTof/cyTof/batch   --out sample_table.csv

# 2) train CytoNorm model on anchors
Rscript scripts/R/02_train_cytonorm.R   --sample-table sample_table.csv   --out-model models/cytonorm_model.rds

# 3) apply normalization to all files
Rscript scripts/R/03_apply_cytonorm.R   --sample-table sample_table.csv   --model models/cytonorm_model.rds   --outdir normalized_fcs

# 4) QC plots (before vs after)
Rscript scripts/R/04_qc_plots.R   --sample-table sample_table.csv   --normalized-dir normalized_fcs   --outdir qc_plots

# 5) export summary per file/channel
Rscript scripts/R/05_export_summary.R   --sample-table sample_table.csv   --normalized-dir normalized_fcs   --out summary
```

All outputs go to `normalized_fcs/`, `qc_plots/`, and `summary/` by default.
