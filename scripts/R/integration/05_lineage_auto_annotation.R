#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

message("1) Loading data")

metacluster_file <- "flowsom_results/metacluster_medians_asinh.csv"
template_file    <- "flowsom_results/cluster_labels_template.csv"
mapping_file     <- "lineage_annotation/lineage_mapping.csv"

stopifnot(file.exists(metacluster_file))

meta <- readr::read_csv(metacluster_file, show_col_types = FALSE)

if (!"cluster" %in% colnames(meta)) {
  stop("metacluster_medians_asinh.csv must have a column named 'cluster'")
}

clusters <- meta$cluster
marker_cols <- setdiff(colnames(meta), "cluster")

## load explicit channel→antigen mapping
if (!file.exists(mapping_file)) {
  stop("lineage_mapping.csv not found at ", mapping_file)
}

map_tab <- readr::read_csv(mapping_file, show_col_types = FALSE)
if (!all(c("channel_name","antigen","group") %in% colnames(map_tab))) {
  stop("lineage_mapping.csv must contain columns: channel_name, antigen, group")
}

## keep only channels present in metacluster medians
map_tab <- map_tab %>%
  filter(channel_name %in% marker_cols)

if (nrow(map_tab) == 0) {
  stop("No channels in lineage_mapping.csv match metacluster_medians_asinh.csv; cannot score lineages.")
}

missing_ch <- setdiff(readr::read_csv(mapping_file, show_col_types = FALSE)$channel_name,
                      map_tab$channel_name)
if (length(missing_ch) > 0) {
  warning("Channels in lineage_mapping.csv not present in metacluster_medians_asinh.csv will be ignored: ",
          paste(unique(missing_ch), collapse = ", "))
}

lineage_markers <- unique(map_tab$channel_name)

lin_mat <- meta %>%
  select(cluster, all_of(lineage_markers)) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

message("Detected ", ncol(lin_mat), " lineage markers for ", nrow(lin_mat), " clusters")

## 2) Scale markers across clusters
message("2) Scaling lineage markers (z score)")

scale_rows <- function(m) {
  t(scale(t(m)))
}

lin_scaled <- scale_rows(lin_mat)
lin_scaled[is.na(lin_scaled)] <- 0

## 3) Heatmap of lineage structure
message("3) Producing lineage heatmap")

dir.create("lineage_annotation", showWarnings = FALSE)

pheatmap::pheatmap(
  lin_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  filename = "lineage_annotation/lineage_heatmap_scaled.png",
  width = 10,
  height = 8
)

pheatmap::pheatmap(
  lin_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  filename = "lineage_annotation/lineage_heatmap_scaled.pdf",
  width = 10,
  height = 8
)

## 4) Define marker groups using mapping

all_markers <- colnames(lin_scaled)

group_to_channels <- split(map_tab$channel_name, map_tab$group)

## ensure we have explicit columns for expected groups, even if empty
all_groups <- unique(c(
  names(group_to_channels),
  "CD4_lineage","CD8_lineage","NK_lineage","Treg_lineage","gd_lineage",
  "Naive_markers","Memory_markers","Activation","Tfh_like"
))

safe_mean <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

## 5) Compute group scores per cluster

message("4) Computing group scores for each cluster")

scores <- map_dfr(
  seq_len(nrow(lin_scaled)),
  function(i) {
    z <- lin_scaled[i, ]
    vals <- map_dbl(all_groups, function(g) {
      ch <- group_to_channels[[g]]
      if (is.null(ch)) return(NA_real_)
      idx <- match(ch, all_markers)
      idx <- idx[!is.na(idx)]
      safe_mean(z[idx])
    })
    tibble(
      cluster = rownames(lin_scaled)[i],
      !!!setNames(as.list(vals), all_groups)
    )
  }
)

## 6) Build automatic labels

message("5) Assigning automatic lineage labels")

## 5) Derive lineage and differentiation labels directly from scores

anno_extra <- scores %>%
  mutate(
    main_lineage_group = dplyr::case_when(
      !is.na(CD4_lineage) | !is.na(CD8_lineage) ~
        if_else(CD4_lineage >= CD8_lineage, "CD4_lineage", "CD8_lineage"),
      TRUE ~ "unknown_lineage"
    ),
    lineage_confidence = dplyr::case_when(
      main_lineage_group %in% c("CD4_lineage","CD8_lineage") ~
        abs(CD4_lineage - CD8_lineage),
      TRUE ~ NA_real_
    ),
    naive_score  = Naive_markers,
    memory_score = Memory_markers,
    activation_score = Activation,
    tfh_score    = Tfh_like,
    diff_state = dplyr::case_when(
      !is.na(naive_score) & !is.na(memory_score) &
        naive_score - memory_score > 0.3 ~ "naive_like",
      !is.na(naive_score) & !is.na(memory_score) &
        memory_score - naive_score > 0.3 ~ "memory_effector_like",
      !is.na(naive_score) & !is.na(memory_score) ~ "intermediate",
      TRUE ~ "unknown_diff"
    ),
    lineage_label = dplyr::case_when(
      main_lineage_group == "CD4_lineage"  ~ "CD4_T",
      main_lineage_group == "CD8_lineage"  ~ "CD8_T",
      main_lineage_group == "NK_lineage"   ~ "NK_like",
      main_lineage_group == "Treg_lineage" ~ "Treg_like",
      main_lineage_group == "gd_lineage"   ~ "gamma_delta_T_like",
      TRUE                                 ~ "unknown_lineage"
    ),
    auto_label = paste(lineage_label, diff_state, sep = "_")
  )

## 7) Merge top markers and optional manual labels

message("6) Extracting top markers per cluster and merging manual template if available")

top_marker_tbl <- map_dfr(
  seq_len(nrow(lin_scaled)),
  function(i) {
    z <- lin_scaled[i, ]
    ord <- order(z, decreasing = TRUE)
    top_idx <- ord[seq_len(min(10, length(ord)))]
    tibble(
      cluster = rownames(lin_scaled)[i],
      top_markers = paste(all_markers[top_idx], collapse = ";")
    )
  }
)

if (file.exists(template_file)) {
  template <- readr::read_csv(template_file, show_col_types = FALSE)
  if (!"cluster" %in% colnames(template)) {
    warning("cluster_labels_template.csv has no 'cluster' column. Ignoring manual labels.")
    template <- tibble(cluster = clusters, manual_label = NA_character_)
  } else {
    if (!"manual_label" %in% colnames(template)) {
      template$manual_label <- NA_character_
    }
    template <- template %>%
      mutate(cluster = as.character(cluster))
  }
} else {
  template <- tibble(cluster = clusters, manual_label = NA_character_)
}

final_anno <- scores %>%
  left_join(
    anno_extra %>%
      select(cluster, main_lineage_group, lineage_confidence,
             naive_score, memory_score, activation_score, tfh_score,
             auto_label),
    by = "cluster"
  ) %>%
  left_join(top_marker_tbl, by = "cluster") %>%
  left_join(template %>% select(cluster, manual_label), by = "cluster") %>%
  arrange(as.integer(cluster))

## 9) Add pathway phenotype per cluster using functional_cluster_annotations
fun_file <- "functional_cluster_annotations.csv"
if (file.exists(fun_file)) {
  fun_tab <- readr::read_csv(fun_file, show_col_types = FALSE)
  if (all(c("cluster","phenotype") %in% colnames(fun_tab))) {
    fun_tab <- fun_tab %>%
      mutate(cluster = as.character(cluster)) %>%
      group_by(cluster, phenotype) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(cluster) %>%
      slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
      ungroup()

    map_pathway <- function(ph) {
      case_when(
        grepl("STAT3", ph, ignore.case = TRUE) ~ "STAT3",
        grepl("STAT5", ph, ignore.case = TRUE) ~ "STAT5",
        grepl("TCR",   ph, ignore.case = TRUE) ~ "TCR",
        grepl("Metabolic", ph, ignore.case = TRUE) ~ "Metabolic",
        grepl("Stress", ph, ignore.case = TRUE) ~ "Stress",
        TRUE ~ "Mixed/other"
      )
    }

    fun_tab <- fun_tab %>%
      mutate(pathway_phenotype = map_pathway(phenotype))

    final_anno <- final_anno %>%
      left_join(fun_tab %>% select(cluster, pathway_phenotype),
                by = "cluster")
  } else {
    warning("functional_cluster_annotations.csv found but lacks 'cluster' or 'phenotype'; skipping pathway phenotype assignment.")
  }
} else {
  warning("functional_cluster_annotations.csv not found; skipping pathway phenotype assignment.")
}

## 8) Save final table

out_csv <- "lineage_annotation/cluster_lineage_annotation_final.csv"
readr::write_csv(final_anno, out_csv)
message("Saved automatic lineage annotation to: ", out_csv)

## 11) Write a simple per-cluster summary in markdown
summary_lines <- c(
  "# Cluster lineage summary",
  "",
  "This file summarises automatic lineage assignments and pathway phenotypes per FlowSOM cluster.",
  ""
)

lineage_counts <- final_anno %>%
  count(main_lineage_group, name = "n") %>%
  arrange(desc(n))

summary_lines <- c(
  summary_lines,
  "## Lineage group counts",
  ""
)
summary_lines <- c(
  summary_lines,
  purrr::map_chr(seq_len(nrow(lineage_counts)), function(i) {
    row <- lineage_counts[i, ]
    paste0("- ", row$main_lineage_group, ": ", row$n)
  }),
  ""
)

if ("pathway_phenotype" %in% colnames(final_anno)) {
  path_counts <- final_anno %>%
    count(pathway_phenotype, name = "n") %>%
    arrange(desc(n))
  summary_lines <- c(
    summary_lines,
    "## Pathway phenotype counts",
    ""
  )
  summary_lines <- c(
    summary_lines,
    purrr::map_chr(seq_len(nrow(path_counts)), function(i) {
      row <- path_counts[i, ]
      paste0("- ", row$pathway_phenotype, ": ", row$n)
    }),
    ""
  )
}

readr::write_lines(summary_lines, "lineage_annotation/cluster_lineage_summary.md")
message("Wrote lineage_annotation/cluster_lineage_summary.md")

## 10) Optional simple UMAP/TSNE overlay if coordinates are available

umap_file <- "UMAP_phospho_coordinates.csv"
if (file.exists(umap_file)) {
  message("9) Found UMAP_phospho_coordinates.csv, generating UMAP_by_auto_lineage.png")
  umap <- readr::read_csv(umap_file, show_col_types = FALSE)
  # expect at least: cluster, cytokine, UMAP1, UMAP2
  if (all(c("cluster", "UMAP1", "UMAP2") %in% colnames(umap))) {
    umap2 <- umap %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(final_anno %>% select(cluster, auto_label), by = "cluster")

    p <- ggplot(umap2, aes(x = UMAP1, y = UMAP2, color = auto_label)) +
      geom_point(size = 0.3, alpha = 0.6) +
      theme_classic() +
      theme(legend.position = "right") +
      labs(title = "UMAP colored by automatic lineage annotation")

    ggsave("lineage_annotation/UMAP_by_auto_lineage.png", p,
           width = 7, height = 6, dpi = 300)

    ## combined lineage + pathway phenotype overlay if available
    if ("pathway_phenotype" %in% colnames(final_anno)) {
      umap3 <- umap2 %>%
        left_join(final_anno %>% select(cluster, pathway_phenotype),
                  by = "cluster") %>%
        mutate(lineage_pathway = interaction(auto_label, pathway_phenotype,
                                             drop = TRUE, lex.order = TRUE))

      p2 <- ggplot(umap3, aes(x = UMAP1, y = UMAP2, color = lineage_pathway)) +
        geom_point(size = 0.3, alpha = 0.6) +
        theme_classic() +
        theme(legend.position = "right") +
        labs(title = "UMAP colored by lineage and pathway phenotype")

      ggsave("lineage_annotation/UMAP_by_lineage_and_phenotype.png", p2,
             width = 7, height = 6, dpi = 300)
    }
  } else {
    warning("UMAP_phospho_coordinates.csv exists but lacks cluster/UMAP1/UMAP2 columns, skipping UMAP overlay")
  }
} else {
  message("UMAP_phospho_coordinates.csv not found, skipping UMAP overlay")
}

tsne_file <- "TSNE_phospho_coordinates.csv"
if (file.exists(tsne_file)) {
  message("10) Found TSNE_phospho_coordinates.csv, generating TSNE overlays")
  tsne <- readr::read_csv(tsne_file, show_col_types = FALSE)
  if (all(c("cluster", "TSNE1", "TSNE2") %in% colnames(tsne))) {
    tsne2 <- tsne %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(final_anno %>% select(cluster, auto_label), by = "cluster")

    p_ts <- ggplot(tsne2, aes(x = TSNE1, y = TSNE2, color = auto_label)) +
      geom_point(size = 0.3, alpha = 0.6) +
      theme_classic() +
      theme(legend.position = "right") +
      labs(title = "t-SNE colored by automatic lineage annotation")

    ggsave("lineage_annotation/TSNE_by_auto_lineage.png", p_ts,
           width = 7, height = 6, dpi = 300)

    if ("pathway_phenotype" %in% colnames(final_anno)) {
      tsne3 <- tsne2 %>%
        left_join(final_anno %>% select(cluster, pathway_phenotype),
                  by = "cluster") %>%
        mutate(lineage_pathway = interaction(auto_label, pathway_phenotype,
                                             drop = TRUE, lex.order = TRUE))

      p_ts2 <- ggplot(tsne3, aes(x = TSNE1, y = TSNE2, color = lineage_pathway)) +
        geom_point(size = 0.3, alpha = 0.6) +
        theme_classic() +
        theme(legend.position = "right") +
        labs(title = "t-SNE colored by lineage and pathway phenotype")

      ggsave("lineage_annotation/TSNE_by_lineage_and_phenotype.png", p_ts2,
             width = 7, height = 6, dpi = 300)
    }
  } else {
    warning("TSNE_phospho_coordinates.csv exists but lacks cluster/TSNE1/TSNE2 columns, skipping TSNE overlay")
  }
} else {
  message("TSNE_phospho_coordinates.csv not found, skipping TSNE overlay")
}

message("Done.")
