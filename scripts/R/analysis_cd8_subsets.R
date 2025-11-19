#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(glue)
})

# Helper ------------------------------------------------------------------
scale_safe <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  sd_x <- sd(x, na.rm = TRUE)
  if (is.na(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

state_from_z <- function(z, pos = 0.5, neg = -0.5) {
  dplyr::case_when(
    is.na(z) ~ NA_character_,
    z >= pos  ~ "high",
    z <= neg  ~ "low",
    TRUE      ~ "dim"
  )
}

assign_core_subset <- function(df) {
  auto_naive <- rep(FALSE, nrow(df))
  if ("auto_label" %in% names(df)) {
    auto_lab <- ifelse(is.na(df$auto_label), "", df$auto_label)
    auto_naive <- grepl("CD8_T_naive_like", auto_lab)
  }

  dplyr::case_when(
    auto_naive ~ "Naive",
    df$naive_flag ~ "Naive",
    df$non_naive_flag & df$CCR7_state == "high" &
      df$CD45RA_state %in% c("low", "dim") &
      df$CD27_state %in% c("high", "dim") ~ "Tcm",
    df$non_naive_flag & df$CCR7_state == "low" &
      df$CD45RA_state %in% c("low", "dim") &
      df$CD27_state %in% c("high", "dim") ~ "Tem_CD27pos",
    df$non_naive_flag & df$CCR7_state == "low" &
      df$CD45RA_state %in% c("low", "dim") &
      df$CD27_state == "low" ~ "Tem_CD27neg",
    df$non_naive_flag & df$CCR7_state == "low" &
      df$CD45RA_state == "high" &
      df$CD27_state == "low" ~ "Temra_CD27neg",
    TRUE ~ "Unresolved"
  )
}

assign_exhaustion_subset <- function(df, core_subset) {
  dplyr::case_when(
    df$non_naive_flag &
      df$pd1_exhaustion == "PD1_hi" &
      (df$cd57_state == "high" | df$tbet_state == "high") ~ "Tex1",
    df$non_naive_flag &
      df$pd1_exhaustion %in% c("PD1_hi", "PD1_mid") &
      df$cd57_state %in% c("low", "dim") &
      df$tbet_state %in% c("low", "dim") &
      df$tig_it_state %in% c("high", "dim") ~ "Tpex1",
    TRUE ~ core_subset
  )
}

# Paths --------------------------------------------------------------------
lineage_file   <- "lineage_annotation/cluster_lineage_annotation_final.csv"
flowsom_file   <- "flowsom_results/metacluster_medians_asinh.csv"
phos_file      <- "phospho_medians/phos_medians_by_sample_cluster.csv"
fun_file       <- "functional_cluster_annotations.csv"
signaling_file <- "cluster_signaling_summary.csv"
umap_file      <- "UMAP_phospho_coordinates.csv"

stopifnot(file.exists(lineage_file))
stopifnot(file.exists(flowsom_file))
stopifnot(file.exists(phos_file))
stopifnot(file.exists(fun_file))
stopifnot(file.exists(signaling_file))

# Output directory ---------------------------------------------------------
out_dir <- "CD8_subset_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data ----------------------------------------------------------------
lineage <- readr::read_csv(lineage_file, show_col_types = FALSE) %>%
  mutate(cluster = as.integer(cluster))
flowsom <- readr::read_csv(flowsom_file, show_col_types = FALSE) %>%
  mutate(cluster = as.integer(cluster))
phos <- readr::read_csv(phos_file, show_col_types = FALSE) %>%
  mutate(
    cluster = as.integer(cluster),
    time_min = suppressWarnings(as.numeric(time_min)),
    cytokine = toupper(cytokine)
  )
fun_tab <- readr::read_csv(fun_file, show_col_types = FALSE) %>%
  mutate(cluster = as.integer(cluster))
signaling <- readr::read_csv(signaling_file, show_col_types = FALSE) %>%
  mutate(cluster_id = as.integer(cluster_id))

phos_marker_cols <- setdiff(colnames(phos), c("cluster", "file_base", "time_min", "cytokine"))

# Aggregate phospho markers of interest (CCR7, CD27, PD1, TIGIT, CD57, Tbet)
state_markers <- c("Ho165Di", "Nd143Di", "Yb172Di", "Nd146Di", "Nd148Di", "Lu175Di")
missing_markers <- setdiff(state_markers, colnames(phos))
if (length(missing_markers) > 0) {
  warning("Markers not found in phos_medians: ", paste(missing_markers, collapse = ", "))
}

phos_summary <- phos %>%
  group_by(cluster) %>%
  summarise(across(all_of(intersect(state_markers, colnames(phos))),
                   ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
            .groups = "drop")

cd8_clusters <- lineage %>%
  filter(main_lineage_group == "CD8_lineage") %>%
  select(cluster, lineage_confidence, naive_score, memory_score,
         activation_score, auto_label)

cd8_markers <- cd8_clusters %>%
  left_join(flowsom %>% select(cluster, Yb174Di), by = "cluster") %>%
  left_join(phos_summary, by = "cluster") %>%
  rename(
    CD45RA = Yb174Di,
    CCR7   = Ho165Di,
    CD27   = Nd143Di,
    PD1    = Yb172Di,
    TIGIT  = Nd146Di,
    CD57   = Nd148Di,
    Tbet   = Lu175Di
  )

# Compute z-scores ---------------------------------------------------------
z_cols <- c("CD45RA", "CCR7", "CD27", "PD1", "TIGIT", "CD57", "Tbet")
cd8_markers <- cd8_markers %>%
  mutate(across(all_of(intersect(z_cols, colnames(cd8_markers))), scale_safe, .names = "{.col}_z"))

# Gating logic -------------------------------------------------------------
cd8_gated <- cd8_markers %>%
  mutate(
    CCR7_state   = state_from_z(CCR7_z),
    CD45RA_state = state_from_z(CD45RA_z),
    CD27_state   = state_from_z(CD27_z),
    PD1_state    = state_from_z(PD1_z, pos = 0.75, neg = -0.75),
    TIGIT_state  = state_from_z(TIGIT_z),
    CD57_state   = state_from_z(CD57_z),
    Tbet_state   = state_from_z(Tbet_z),
    naive_flag    = CCR7_state == "high" & CD45RA_state == "high" & CD27_state %in% c("high", "dim"),
    non_naive_flag = !naive_flag,
    pd1_exhaustion = case_when(
      is.na(PD1_state)     ~ "PD1_NA",
      PD1_state == "high"  ~ "PD1_hi",
      PD1_state == "low"   ~ "PD1_lo",
      TRUE                 ~ "PD1_mid"
    ),
    tig_it_state  = TIGIT_state,
    cd57_state    = CD57_state,
    tbet_state    = Tbet_state
  ) %>%
  mutate(
    core_subset = assign_core_subset(
      pick(naive_flag, non_naive_flag, CCR7_state,
           CD45RA_state, CD27_state, auto_label)
    ),
    subset_call = assign_exhaustion_subset(
      pick(pd1_exhaustion, cd57_state, tbet_state, tig_it_state, non_naive_flag),
      core_subset
    )
  )

# Merge functional annotations --------------------------------------------
cd8_functional <- fun_tab %>%
  inner_join(cd8_gated %>% select(cluster, subset_call), by = "cluster")

subset_functional_summary <- cd8_functional %>%
  count(subset_call, cytokine, phenotype, name = "n") %>%
  group_by(subset_call, cytokine) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

readr::write_csv(subset_functional_summary,
                 file.path(out_dir, "cd8_subset_cytokine_functional_summary.csv"))

# Prepare final assignment table ------------------------------------------
assignment_cols <- c(
  "cluster", "lineage_confidence", "core_subset", "subset_call",
  "naive_flag", "non_naive_flag", "pd1_exhaustion", "tig_it_state",
  "cd57_state", "tbet_state", paste0(z_cols, c("", "_z"))
)
assignment_cols <- assignment_cols[assignment_cols %in% colnames(cd8_gated)]

subset_order <- c("Naive", "Tcm", "Tem_CD27pos", "Tem_CD27neg",
                  "Temra_CD27neg", "Tpex1", "Tex1", "Unresolved")

cd8_assignments <- cd8_gated %>%
  select(all_of(assignment_cols)) %>%
  arrange(match(subset_call, subset_order), cluster)

readr::write_csv(cd8_assignments, file.path(out_dir, "cd8_subset_assignments.csv"))

# Heatmap of marker states -------------------------------------------------
heatmap_df <- cd8_gated %>%
  select(cluster, subset_call, ends_with("_z")) %>%
  pivot_longer(cols = ends_with("_z"), names_to = "marker", values_to = "z") %>%
  mutate(marker = str_remove(marker, "_z"))

cluster_levels <- cd8_assignments$cluster
marker_levels <- c("CCR7", "CD45RA", "CD27", "PD1", "TIGIT", "CD57", "Tbet")

heatmap_plot <- ggplot(heatmap_df,
                       aes(x = factor(marker, levels = marker_levels),
                           y = factor(cluster, levels = cluster_levels),
                           fill = z)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-2.5, 2.5), oob = scales::squish) +
  facet_grid(subset_call ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Marker", y = "Cluster", fill = "z-score",
       title = "CD8 metacluster marker z-scores",
       subtitle = "Thresholds: |z| >= 0.5") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(filename = file.path(out_dir, "cd8_marker_heatmap.png"),
                plot = heatmap_plot, width = 6, height = 6.5, dpi = 300)

# UMAP overlay -------------------------------------------------------------
if (file.exists(umap_file)) {
  umap_tab <- readr::read_csv(umap_file, show_col_types = FALSE) %>%
    mutate(cluster = as.integer(cluster)) %>%
    left_join(cd8_assignments %>% select(cluster, subset_call), by = "cluster") %>%
    mutate(subset_plot = coalesce(subset_call, "Other lineages"))

  subset_palette <- c(
    "Naive" = "#1b9e77",
    "Tcm" = "#7570b3",
    "Tem_CD27pos" = "#d95f02",
    "Tem_CD27neg" = "#e7298a",
    "Temra_CD27neg" = "#66a61e",
    "Tpex1" = "#e6ab02",
    "Tex1" = "#a6761d",
    "Unresolved" = "#666666",
    "Other lineages" = "#bdbdbd"
  )

  umap_plot <- ggplot(umap_tab, aes(x = UMAP1, y = UMAP2, color = subset_plot)) +
    geom_point(size = 0.2, alpha = 0.5) +
    theme_classic(base_size = 11) +
    scale_color_manual(values = subset_palette, drop = FALSE) +
    labs(title = "UMAP colored by CD8 subset calls", color = "Subset")

  ggplot2::ggsave(filename = file.path(out_dir, "UMAP_cd8_subsets.png"),
                  plot = umap_plot, width = 6.5, height = 5.5, dpi = 300)
}

# Trajectory summaries -----------------------------------------------------
subset_long <- phos %>%
  filter(cluster %in% cd8_assignments$cluster) %>%
  left_join(cd8_assignments %>% select(cluster, subset_call), by = "cluster") %>%
  filter(!is.na(subset_call)) %>%
  pivot_longer(cols = all_of(phos_marker_cols), names_to = "marker", values_to = "signal")

subset_sample_marker <- subset_long %>%
  group_by(file_base, cytokine, time_min, subset_call, marker) %>%
  summarise(median_signal = mean(signal, na.rm = TRUE), .groups = "drop")

trajectory_marker <- subset_sample_marker %>%
  group_by(cytokine, subset_call, marker, time_min) %>%
  summarise(
    mean_signal = mean(median_signal, na.rm = TRUE),
    sd_signal = sd(median_signal, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

trajectory_overall <- subset_sample_marker %>%
  group_by(file_base, cytokine, time_min, subset_call) %>%
  summarise(mean_phospho = mean(median_signal, na.rm = TRUE), .groups = "drop") %>%
  group_by(cytokine, subset_call, time_min) %>%
  summarise(
    mean_signal = mean(mean_phospho, na.rm = TRUE),
    sd_signal = sd(mean_phospho, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

readr::write_csv(trajectory_marker, file.path(out_dir, "cd8_subset_marker_trajectories.csv"))
readr::write_csv(trajectory_overall, file.path(out_dir, "cd8_subset_overall_trajectory.csv"))

# Plots: global trajectories ----------------------------------------------
overall_plot <- ggplot(trajectory_overall %>% filter(!is.na(time_min)),
                       aes(x = time_min, y = mean_signal,
                           color = subset_call)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~ cytokine, scales = "free_y") +
  scale_color_manual(values = c(
    "Naive" = "#1b9e77",
    "Tcm" = "#7570b3",
    "Tem_CD27pos" = "#d95f02",
    "Tem_CD27neg" = "#e7298a",
    "Temra_CD27neg" = "#66a61e",
    "Tpex1" = "#e6ab02",
    "Tex1" = "#a6761d",
    "Unresolved" = "#666666"
  ), drop = FALSE) +
  labs(title = "Overall phospho trajectories by subset",
       x = "Time (min)", y = "Mean phospho signal (all markers)",
       color = "Subset") +
  theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(out_dir, "overall_subset_trajectories.png"),
                overall_plot, width = 10, height = 6.5, dpi = 300)

marker_plot <- ggplot(trajectory_marker %>% filter(!is.na(time_min)),
                      aes(x = time_min, y = mean_signal,
                          color = subset_call)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(marker ~ cytokine, scales = "free_y") +
  scale_color_manual(values = c(
    "Naive" = "#1b9e77",
    "Tcm" = "#7570b3",
    "Tem_CD27pos" = "#d95f02",
    "Tem_CD27neg" = "#e7298a",
    "Temra_CD27neg" = "#66a61e",
    "Tpex1" = "#e6ab02",
    "Tex1" = "#a6761d",
    "Unresolved" = "#666666"
  ), drop = FALSE) +
  labs(title = "Marker-level trajectories across CD8 subsets",
       x = "Time (min)", y = "Mean phospho signal") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 7),
        legend.position = "bottom")

ggplot2::ggsave(file.path(out_dir, "marker_subset_trajectories.png"),
                marker_plot, width = 12, height = 8.5, dpi = 300)

# Condition-specific comparisons ------------------------------------------
requested_conditions <- c(
  "CD3",
  "CD3+CD28",
  "CD3+CD28+IL2",
  "CD3+CD28+IL21",
  "CD3+CD28+IL15",
  "CD3+CD28+PAC IL21"
)

condition_map <- c(
  "CD3" = "CD3",
  "CD3+CD28" = "CD3",
  "CD3+CD28+IL2" = "IL2",
  "CD3+CD28+IL21" = "IL21",
  "CD3+CD28+IL15" = "IL15",
  "CD3+CD28+PAC IL21" = "IL21"
)

missing_conditions <- character()

for (cond in requested_conditions) {
  mapped <- condition_map[[cond]]
  cond_df <- trajectory_overall %>%
    filter(cytokine == mapped, !is.na(time_min))

  if (nrow(cond_df) == 0) {
    missing_conditions <- c(missing_conditions, cond)
    next
  }

  cond_plot <- ggplot(cond_df,
                      aes(x = time_min, y = mean_signal,
                          color = subset_call)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c(
      "Naive" = "#1b9e77",
      "Tcm" = "#7570b3",
      "Tem_CD27pos" = "#d95f02",
      "Tem_CD27neg" = "#e7298a",
      "Temra_CD27neg" = "#66a61e",
      "Tpex1" = "#e6ab02",
      "Tex1" = "#a6761d",
      "Unresolved" = "#666666"
    ), drop = FALSE) +
    labs(
      title = glue("Subset trajectories: {cond} (mapped to '{mapped}')"),
      x = "Time (min)",
      y = "Mean phospho signal",
      color = "Subset"
    ) +
    theme_minimal(base_size = 11)

  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cond)
  ggplot2::ggsave(file.path(out_dir, glue("trajectory_{safe_name}.png")),
                  cond_plot, width = 7, height = 5, dpi = 300)
}

# Summary markdown --------------------------------------------------------
subset_counts <- cd8_assignments %>%
  count(subset_call, name = "n_clusters") %>%
  arrange(match(subset_call, subset_order))

pd1_summary <- cd8_assignments %>%
  count(subset_call, pd1_exhaustion, name = "n") %>%
  group_by(subset_call) %>%
  mutate(freq = scales::percent(n / sum(n), accuracy = 0.1)) %>%
  ungroup()

trajectory_note <- if (length(missing_conditions) > 0) {
  paste0("- Requested conditions without data: ",
         paste(missing_conditions, collapse = ", "),
         ". Check raw metadata if these stimuli exist under different labels.")
} else {
  "- All requested condition comparisons were generated."
}

summary_lines <- c(
  "# CD8 subset gating and trajectory summary",
  "",
  glue("Total CD8 metaclusters: {nrow(cd8_assignments)}"),
  glue("Lineage confidence (median): {round(stats::median(cd8_gated$lineage_confidence, na.rm = TRUE), 3)}"),
  "",
  "## Subset counts",
  subset_counts %>% mutate(line = glue("- {subset_call}: {n_clusters} clusters")) %>% pull(line),
  "",
  "## PD1 distribution per subset",
  pd1_summary %>% mutate(line = glue("- {subset_call} / {pd1_exhaustion}: {n} clusters ({freq})")) %>% pull(line),
  "",
  "## Trajectory notes",
  trajectory_note,
  "- Tex/Tpex assignments approximate PD1/TIGIT/CD57/Tbet states because TCF1/CD127 are unavailable in this panel.",
  "- Condition labels were mapped to available cytokine metadata (CD3, IL2, IL15, IL21, UNSTIM)."
)

readr::write_lines(summary_lines, file.path(out_dir, "cd8_subset_summary.md"))

message("CD8 subset analysis complete. Outputs written to ", out_dir)
