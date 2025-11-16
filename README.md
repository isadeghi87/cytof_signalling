# γc Cytokine Signalling CyTOF Pipeline

This repository contains an end‑to‑end CyTOF analysis pipeline for human T cells stimulated with γc cytokines (CD3, IL‑2, IL‑15, IL‑21).  
The code starts from raw FCS files, performs CytoNorm‑based batch normalisation, FlowSOM clustering and dimensionality reduction, computes cluster‑level phospho‑signalling, runs limma differential analyses, derives functional STAT3/STAT5/TCR modules, performs automated lineage annotation, and assembles summary figures and a manuscript.

The code and figures are intended to be fully reproducible once the raw FCS and metadata are available locally (these are **not** included in the repository and should not be pushed to GitHub).

---

## 1. Repository Structure

- `batch/` – raw FCS files and anchor files (local only; do not publish).
  - `batch/anchors/` – anchor FCS files for CytoNorm.
  - `OMIQ_metadata-*.csv` – experimental metadata exported from OMIQ.
- `normalized_fcs/`, `normalized_fcs_rescaled/` – CytoNorm‑normalised FCS files (**large; keep local**).
- `scripts/R/` – all R scripts for normalisation, clustering, QC, and downstream analysis.
  - `00_install_pkgs.R` – installs all required R packages.
  - `01_prepare_metadata.R` – builds `sample_table.csv` from file lists, anchors, and metadata.
  - `01_run_cytonorm.R` – main CytoNorm normalisation entrypoint.
  - `02_train_cytonorm*.R`, `03_apply_cytonorm*.R` – lower‑level CytoNorm steps.
  - `03_dimred_flowsom.R`, `04_tns_umap_plots.R`, `04c_tsne_umap_clusters.R` – t‑SNE / UMAP plots.
  - `05_compute_phospho_medians.R`, `06_diff_phospho_analysis.R`, `diff_phos.R` – phospho medians and limma.
  - `integration/` – automated lineage mapping and annotation.
  - `next_steps/` – functional module scoring and higher‑level summaries.
  - `run_v2.sh` – full downstream analysis pipeline.
- `run_all.sh` – one‑shot script to go from raw FCS + metadata → CytoNorm normalised FCS + QC.
- `flowsom_results/` – FlowSOM clustering, MST tree and metacluster medians.
- `dimred_flowsom/`, `tsne_umap_plots/` – dimension‑reduction embeddings.
- `phospho_medians/`, `diff_phos_results/` – cluster‑level medians and limma differential signalling.
- `lineage_annotation/` – channel mapping and automated lineage/differentiation labels.
- `summary/`, `summary_figures/` – QC summaries, module scores, cytokine bias plots.
- `diagnostic_plots/`, `qc_plots/`, `qc_plots_debug/` – additional QC and diagnostics.
- `figures/` – final multi‑panel figures (`Figure1.png`–`Figure6.png`) used in the manuscript.
- `MANUSCRIPT_science_immunology.*` – manuscript source (`.md`/`.tex`) and compiled outputs (`.pdf`, `.docx`).
- `FINAL_REPORT_v3.*` – internal summary report.

---

## 2. Requirements & Installation

**Software**

- R ≥ 4.1 (tested on recent R 4.x).
- Recommended: Unix‑like environment (Linux/macOS) with `bash` and `make`.

**R packages**

Key packages include: `data.table`, `dplyr`, `flowCore`, `CytoNorm`, `FlowSOM`, `Rtsne`, `uwot`, `ggplot2`, `limma`, and others.

Install all dependencies once using:

```bash
Rscript scripts/R/00_install_pkgs.R
```

If you use conda or a module system, activate your R environment before running the pipeline (see comments inside `scripts/R/run_v2.sh`).

---

## 3. Input Data and Metadata

The pipeline assumes the following inputs **on your local machine**:

- Raw FCS files under `batch/`:
  - Sample FCS: `batch/Batch_X_*.fcs`
  - Anchor FCS: `batch/anchors/Batch_X_Anchor-*.fcs`
- A metadata CSV exported from OMIQ (or similar), e.g.:
  - `batch/OMIQ_metadata-20250919_LS_pTCR_cHCV.csv`
- A text file listing the non‑anchor sample FCS files:
  - `files.txt` – one FCS filename per line.
- A list of anchor FCS files:
  - `batch/anchors/anchors.txt` – one FCS filename per line (ideally at least one per batch).
- A TSV listing channels to normalise and inspect:
  - `channels_to_adjust.tsv` – column `channel` with metal channel names to normalise.

These inputs are used to construct `sample_table.csv`, which becomes the central design table for all downstream analyses.

> **Important:** Raw FCS and per‑event data should generally **not** be pushed to a public GitHub repository. Keep `batch/`, `normalized_fcs/`, and `normalized_fcs_*` local or on secure storage.

---

## 4. Normalisation and QC: `run_all.sh`

To go from raw FCS + metadata to normalised FCS and QC summaries:

```bash
./run_all.sh \
  /path/to/files.txt \
  /path/to/anchors/anchors.txt \
  /path/to/ChannelsToAdjust.txt \
  /path/to/OMIQ_metadata.csv \
  /path/to/batch_dir
```

With defaults (as currently configured in the project) you can simply run:

```bash
./run_all.sh
```

This script performs:

1. **Metadata assembly** – `scripts/R/01_prepare_metadata.R`  
   Creates `sample_table.csv` with `file_base`, `batch`, `is_anchor`, `sample_id`, `stimulus`, `time_min`, etc.
2. **CytoNorm training & application** – `scripts/R/02_train_cytonorm.R`, `03_apply_cytonorm.R`  
   Trains a quantile normalisation model on anchors and applies it to all samples, writing to `normalized_fcs/`.
3. **Consolidation** – `scripts/R/03b_consolidate_normalized.R`  
   Collapses multiple CytoNorm outputs per sample to a canonical filename.
4. **QC & diagnostics** – `scripts/R/04b_plot_anchor_densities.R`, `04c_tsne_umap_clusters.R`, `05_plot_normalization_efficiency.R`, `04_qc_plots.R`, `05_export_summary.R`  
   Produces:
   - Marker density plots and anchor overlays (`qc_plots/`, `summary/more_qc/`).
   - t‑SNE/UMAP coloured by batch, cluster, and time (`qc_plots/`).
   - Normalisation efficiency plots (`summary/normalization_efficiency.pdf`, `summary/normalization_efficiency_lineage.pdf`).
   - Channel‑level medians and event counts (`summary/channel_medians.tsv`, `summary/events_post.tsv`).

Outputs from this stage:

- `sample_table.csv` – master metadata table.  
- `normalized_fcs/` – normalised FCS per sample.  
- `qc_plots/`, `diagnostic_plots/`, `summary/` – QC plots and summaries.

You can use `run_all.sh` on new datasets by updating `files.txt`, anchor lists, channel lists and metadata.

---

## 5. Full Downstream Pipeline: `scripts/R/run_v2.sh`

Once `sample_table.csv` and `normalized_fcs/` are available, run the full analysis pipeline:

```bash
cd scripts/R
./run_v2.sh
```

This script orchestrates the following steps (logging to `logs/`):

1. **CytoNorm normalisation (optional repeat)** – `01_run_cytonorm.R`  
   Confirms/updates normalised FCS and QC plots.
2. **FlowSOM clustering** – `02_flowsom_annotate.R`  
   - Trains a FlowSOM model on normalised data.  
   - Derives ~20 metaclusters.  
   - Generates MST trees and heatmaps of metacluster medians:  
     - `flowsom_results/flowsom_MST_tree.png`  
     - `flowsom_results/metacluster_heatmap.png`  
     - `flowsom_results/metacluster_lineage_heatmap.png`
3. **Dimensionality reduction** – `03_dimred_flowsom.R`, `04_tns_umap_plots.R`  
   - t‑SNE and UMAP using FlowSOM clusters as features.  
   - Outputs in `dimred_flowsom/` and `tsne_umap_plots/` coloured by cluster, cytokine, time, lineage, etc.
4. **Cluster‑level phospho medians** – `05_compute_phospho_medians.R`  
   - Computes arcsinh‑transformed medians per cluster × sample for phospho markers.  
   - Writes `phospho_medians/phos_medians_by_sample_cluster.csv`.
5. **Differential phospho‑signalling (limma)** – `06_diff_phospho_analysis.R`, `diff_phos.R`  
   - Fits limma models for cytokine/time contrasts per marker.  
   - Outputs `diff_phos_results/diff_phospho_limma_results.csv` and logFC heatmaps (`logFC_heatmaps/`).
6. **Alternative clustering (PhenoGraph)** – `07_phenograph_clustering.R`  
   - Optional validation of FlowSOM results via PhenoGraph.  
   - Outputs labels and scatter plots in `phenograph_results/`.
7. **Signalling trajectories (riverplots)** – `08_riverplots_signaling.R`  
   - Builds cluster‑level timecourses and riverplots of STAT3/STAT5/TCR signalling.  
   - Outputs `riverplots/` and `trajectory_plots/`.

Additional scripts under `scripts/R/integration/` and `scripts/R/next_steps/` perform:

- Automated lineage annotation (`lineage_annotation/cluster_lineage_annotation_final.csv`).  
- Functional module scoring (STAT3, STAT5, TCR‑proximal, etc.).  
- Cytokine‑bias summary per cluster (`cluster_cytokine_bias.csv`).  
- High‑level summaries (`cluster_signaling_summary.csv`, `cluster_functional_annotations.csv`).

---

## 6. Key Outputs and How They Relate to the Manuscript

After the full pipeline, the most important outputs are:

**Tables**

- `sample_table.csv` – design table used by all scripts.  
- `flowsom_results/metacluster_medians_asinh.csv` – cluster × marker medians.  
- `phospho_medians/phos_medians_by_sample_cluster.csv` – per‑sample × cluster phospho medians.  
- `diff_phos_results/diff_phospho_limma_results.csv` – limma differential signalling results.  
- `functional_cluster_annotations.csv` – STAT3/STAT5/TCR dominance per cluster × cytokine.  
- `cluster_signaling_summary.csv` – AUCs, dominant markers, best timepoints per cluster.  
- `cluster_cytokine_bias.csv` – best cytokine (IL‑2/IL‑15/IL‑21/CD3) per pathway and cluster.  
- `lineage_annotation/cluster_lineage_annotation_final.csv` – lineage and differentiation labels for each metacluster.

**Figures**

- Normalisation & QC:  
  - `summary/normalization_efficiency.pdf`, `summary/normalization_efficiency_lineage.pdf`  
  - `summary/more_qc/density_overlays_top_channels.png`, `qc_plots/anchor_densities.png`
- Clustering & embeddings:  
  - `flowsom_results/flowsom_MST_tree.png`, `flowsom_results/metacluster_heatmap.png`  
  - `dimred_flowsom/umap_by_cluster.png`, `tsne_umap_plots/umap_by_cytokine.png`, etc.
- Functional modules & cytokine bias:  
  - `functional_signatures.png`, `cytokine_cluster_comparison.png`, `activation_bars_by_cytokine.png`
- Lineage‑resolved signalling:  
  - `lineage_annotation/lineage_heatmap_scaled.png`  
  - `lineage_annotation/UMAP_by_auto_lineage.png`  
  - `summary_figures/UMAP_by_phenotype.png`, `summary_figures/UMAP_by_cytokine_preference.png`
- Trajectories and diffusion maps:  
  - `summary_figures/key_trajectories_STAT3_STAT5.png`  
  - `diffusion_map_CD3.png`, `diffusion_map_IL2.png`, `diffusion_map_IL21.png`

**Final assembled multi‑panel figures**

- `figures/Figure1.png` – experimental design, coverage, and normalisation QC.  
- `figures/Figure2.png` – FlowSOM clustering, MST, and QC.  
- `figures/Figure3.png` – PCA/UMAP of phospho signalling and functional modules.  
- `figures/Figure4.png` – time‑resolved limma results and cytokine bias.  
- `figures/Figure5.png` – lineage annotation and lineage × phenotype × cytokine preference.  
- `figures/Figure6.png` – IL‑21‑biased dual‑STAT trajectories and diffusion maps.

These correspond directly to Figures 1–6 referenced in `MANUSCRIPT_science_immunology.md` and the compiled `.pdf`.

To regenerate the manuscript after updating the text, knit `MANUSCRIPT_science_immunology.md` (e.g. via RStudio or `rmarkdown::render`), or rebuild the LaTeX version `MANUSCRIPT_science_immunology.tex`.

---

## 7. Adapting the Pipeline to New Datasets

To reuse this pipeline on another CyTOF dataset:

1. Copy the repository to a new project directory.  
2. Populate `batch/` with your new raw FCS files and anchors (under `batch/anchors/`).  
3. Export a metadata CSV from your acquisition/analysis software and place it under `batch/`.  
4. Create/update:
   - `files.txt` – list of sample FCS files (no anchors).  
   - `batch/anchors/anchors.txt` – list of anchor FCS files.  
   - `channels_to_adjust.tsv` – channels to normalise.  
5. Run `./run_all.sh` to generate `sample_table.csv`, `normalized_fcs/`, and QC plots.  
6. Run `scripts/R/run_v2.sh` to generate clusters, medians, differential signalling, lineage annotations, and summary figures.  
7. Update the manuscript text and figure legends if your design differs (number of donors, cytokines, timepoints, etc.).

If your panel uses different channels for STAT3/STAT5/TCR, update the marker definitions in the relevant R scripts under `scripts/R/next_steps/` and `scripts/R/integration/`.

---

## 8. Data, Privacy, and .gitignore Recommendations

For a public GitHub repository you will usually want to:

- **Keep raw and normalised FCS files private.**  
  Do not commit `batch/`, `normalized_fcs/`, `normalized_fcs_rescaled/`, or other per‑cell data directories.
- Consider ignoring large intermediate outputs that can be regenerated:
  - `qc_plots/`, `qc_plots_debug/`, `diagnostic_plots/`  
  - `dimred_flowsom/`, `tsne_umap_plots/`, `phenograph_results/`, `riverplots/`, `trajectory_plots/`

A minimal `.gitignore` could include:

```gitignore
batch/
normalized_fcs*/
qc_plots*/
diagnostic_plots/
dimred_flowsom/
tsne_umap_plots/
phenograph_results/
riverplots*/
trajectory_plots/
logs/
*.Rproj
.Rhistory
.RData
```

You may choose to keep key summary tables and final figures under version control so that others can inspect the main results without raw data.

---

## 9. Citation and Contact

If you use this pipeline or derived figures/tables in a publication, please cite the associated γc cytokine signalling manuscript once it is published, and acknowledge the use of the CytoNorm/FlowSOM‑based pipeline.

For questions about the code or analysis, please contact the corresponding author or maintainer of this repository.

This repository currently has no explicit open‑source license; please do not redistribute or reuse the code for commercial purposes without permission.
