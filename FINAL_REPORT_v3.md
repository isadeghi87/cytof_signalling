# Final Report v3

## 1. Overview

This report summarises functional signaling phenotypes, cytokine preferences,
and compact summary figures derived from cluster-level phospho medians.

Key new outputs:
- `cluster_signaling_summary.csv`
- `cluster_cytokine_bias.csv`
- `summary_figures/UMAP_by_phenotype.png`
- `summary_figures/UMAP_by_cytokine_preference.png`
- `summary_figures/activation_bars_by_cytokine.png`
- `summary_figures/key_trajectories_STAT3_STAT5.png`

## 2. Functional Phenotypes per Cytokine

### Phenotype counts (all cytokines combined)

- STAT5-dominant: 40 cluster–cytokine pairs
- STAT3-dominant: 40 cluster–cytokine pairs

### Dominant phenotypes per cytokine

- CD3: STAT5-dominant (10); STAT3-dominant (10)
- IL15: STAT5-dominant (15); STAT3-dominant (5)
- IL2: STAT3-dominant (15); STAT5-dominant (5)
- IL21: STAT5-dominant (10); STAT3-dominant (10)

## 3. Cytokine Preferences per Cluster

Based on integrated logFC for STAT3/STAT5/TCR marker groups,
each cluster has a preferred cytokine for each pathway.

Examples (first 10 clusters):

- Cluster 1 (Cluster_1): STAT3 → IL15; STAT5 → IL21; TCR → NA; notes: IL15 > IL2 > CD3 > IL21 > UNSTIM
- Cluster 2 (Cluster_2): STAT3 → IL15; STAT5 → IL2; TCR → NA; notes: IL15 > IL2 > IL21 > CD3 > UNSTIM
- Cluster 3 (Cluster_3): STAT3 → IL2; STAT5 → IL21; TCR → NA; notes: IL2 > IL15 > IL21 > CD3 > UNSTIM
- Cluster 4 (Cluster_4): STAT3 → IL2; STAT5 → IL2; TCR → NA; notes: IL2 > IL21 > IL15 > CD3 > UNSTIM
- Cluster 5 (Cluster_5): STAT3 → IL2; STAT5 → IL2; TCR → NA; notes: IL2 > CD3 > IL15 > IL21 > UNSTIM
- Cluster 6 (Cluster_6): STAT3 → IL15; STAT5 → IL15; TCR → NA; notes: IL15 > IL2 > IL21 > CD3 > UNSTIM
- Cluster 7 (Cluster_7): STAT3 → IL15; STAT5 → IL21; TCR → NA; notes: IL15 > CD3 > IL2 > IL21 > UNSTIM
- Cluster 8 (Cluster_8): STAT3 → IL15; STAT5 → IL21; TCR → NA; notes: IL15 > IL2 > IL21 > CD3 > UNSTIM
- Cluster 9 (Cluster_9): STAT3 → IL2; STAT5 → IL21; TCR → NA; notes: IL2 > IL15 > IL21 > CD3 > UNSTIM
- Cluster 10 (Cluster_10): STAT3 → IL15; STAT5 → IL15; TCR → NA; notes: IL15 > CD3 > IL2 > IL21 > UNSTIM

## 4. Most Responsive Clusters

Using `cluster_activation_rankings.csv`, the following clusters show highest
integrated activation (AUC) per cytokine:

- CD3: cluster 18 (1096.3); cluster 4 (1056); cluster 11 (1007.3)
- IL15: cluster 4 (1316.8); cluster 18 (1295.7); cluster 11 (1166.3)
- IL2: cluster 18 (1194.5); cluster 4 (1191.9); cluster 11 (1064.1)
- IL21: cluster 4 (1298); cluster 18 (1170); cluster 5 (1109.3)

## 5. Figure References

- `summary_figures/UMAP_by_phenotype.png`: UMAP coloured by functional phenotype.
- `summary_figures/UMAP_by_cytokine_preference.png`: UMAP coloured by best STAT3 cytokine.
- `summary_figures/activation_bars_by_cytokine.png`: AUC per cluster×cytokine, annotated with top responders.
- `summary_figures/key_trajectories_STAT3_STAT5.png`: STAT3/STAT5 time courses for top responder clusters.

## 6. Limitations

- FlowSOM vs PhenoGraph stability is only partially assessed:
  the available PhenoGraph labels are not perfectly row-aligned with the
  FlowSOM-labelled event set, so ARI/NMI metrics are approximate at best.
- All signaling summaries are based on cluster medians, not raw single-cell
  trajectories; sub-cluster heterogeneity or rare subpopulations may be
  smoothed out.
- Functional phenotypes and cytokine preferences are derived from a limited
  marker panel and early time-point contrasts; additional markers or time
  points could refine these assignments.
- No batch-specific or subject-level random effects are modeled; all
  statistics are aggregate across samples.

