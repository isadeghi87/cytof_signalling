#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

set.seed(1)

###############################################################################
## 1. LOAD AND VALIDATE STRUCTURE
###############################################################################

st <- fread("sample_table.csv")
if (!"file_base" %in% names(st)) {
  if ("file" %in% names(st)) {
    st[, file_base := basename(file)]
  } else {
    stop("sample_table.csv must contain 'file' or 'file_base'.")
  }
}

if (!"cytokine" %in% names(st)) {
  if (!"stim" %in% names(st)) {
    stop("sample_table.csv must contain 'cytokine' or 'stim'.")
  }
  stim_to_cytokine <- function(stim_vec) {
    v <- tolower(as.character(stim_vec))
    res <- ifelse(grepl("il21", v), "IL21",
           ifelse(grepl("il15", v), "IL15",
           ifelse(grepl("il2",  v), "IL2",
           ifelse(grepl("pac", v), "PAC",
           ifelse(grepl("cd3", v) & !grepl("il", v), "CD3",
           ifelse(grepl("unstim", v) | v == "" | is.na(v),
                  "UNSTIM", "OTHER"))))))
    res
  }
  st[, cytokine := stim_to_cytokine(stim)]
}

pm_path <- "diff_phos_results/phos_medians_by_sample_cluster.csv"
if (!file.exists(pm_path)) pm_path <- "phospho_medians/phos_medians_by_sample_cluster.csv"
pm <- fread(pm_path)

mc <- fread("flowsom_results/metacluster_medians_asinh.csv")
cl_tmpl <- fread("flowsom_results/cluster_labels_template.csv")
pl <- fread("phenograph_results/phenograph_labels.csv")

int_lines <- character()
int_lines <- c(int_lines,
  sprintf("sample_table: n=%d rows, unique file_base=%d",
          nrow(st), length(unique(st$file_base))),
  sprintf("sample_table: cytokines=%s",
          paste(sort(unique(st$cytokine)), collapse=", ")),
  sprintf("sample_table: time_min=%s",
          paste(sort(unique(st$time_min)), collapse=", ")))

## check is_anchor
is_anchor_tab <- capture.output(print(table(st$is_anchor, useNA="ifany")))
int_lines <- c(int_lines,
  sprintf("sample_table: is_anchor table: %s", paste(is_anchor_tab, collapse=" ")))

## phos medians uniqueness
key_pm <- c("file_base","cluster","cytokine","time_min")
dup_pm <- pm[duplicated(pm[, ..key_pm])]
int_lines <- c(int_lines,
  sprintf("phos_medians: n=%d, duplicated key rows=%d",
          nrow(pm), nrow(dup_pm)))

## phenograph alignment (length only)
int_lines <- c(int_lines,
  sprintf("phenograph_labels rows=%d", nrow(pl)))

## metaclusters
int_lines <- c(int_lines,
  sprintf("metacluster_medians: n_clusters=%d (expected 20)",
          length(unique(mc$cluster))))

## cluster_labels_template check
int_lines <- c(int_lines,
  sprintf("cluster_labels_template: n_clusters=%d",
          length(unique(cl_tmpl$cluster))))

writeLines(int_lines, "integrity_summary.txt")

###############################################################################
## 2. FULL LIMMA: TIME CONTRASTS VS 0
###############################################################################

id_cols <- c("file_base","cluster","time_min","cytokine")
phos_cols <- setdiff(names(pm), id_cols)

res_all <- list()

clusters <- sort(unique(pm$cluster))
cyt_list <- sort(unique(pm$cytokine))

for (cl in clusters) {
  for (cy in cyt_list) {
    sub <- pm[cluster == cl & cytokine == cy]
    if (nrow(sub) < 2) next

    tlev <- sort(unique(sub$time_min))
    if (!0 %in% tlev && !"0" %in% tlev) next
    base <- if (0 %in% tlev) 0 else 0

    time_fac <- factor(sub$time_min, levels = tlev)
    design <- model.matrix(~ 0 + time_fac)
    colnames(design) <- paste0("t", tlev)

    Y <- t(as.matrix(sub[, ..phos_cols]))
    fit <- lmFit(Y, design)

    other_t <- setdiff(tlev, base)
    if (!length(other_t)) next

    contr_str <- paste(sprintf("t%s - t%s", other_t, base), collapse = ";")
    cm <- makeContrasts(contrasts = contr_str, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm))

    n_coef <- ncol(cm)
    for (j in seq_len(n_coef)) {
      tt <- tryCatch(
        topTable(fit2, coef = j, number = Inf, sort.by = "none"),
        error = function(e) NULL
      )
      if (is.null(tt)) next
      tt$cluster  <- cl
      tt$cytokine <- cy
      tt$contrast <- paste0("t", other_t[j], "-t", base)
      tt$time_min <- other_t[j]
      tt$marker   <- rownames(tt)
      res_all[[length(res_all)+1]] <- as.data.table(tt)
    }
  }
}

###############################################################################
## 3. LIMMA: CYTOKINE CONTRASTS AT MATCHED TIME POINTS
###############################################################################

pair_list <- list(
  c("IL2","CD3"),
  c("IL15","CD3"),
  c("IL21","CD3")
)

for (cl in clusters) {
  for (pair in pair_list) {
    cyA <- pair[1]; cyB <- pair[2]
    sub <- pm[cluster == cl & cytokine %in% c(cyA, cyB)]
    if (nrow(sub) < 2) next

    times_here <- sort(unique(sub$time_min))
    for (tt in times_here) {
      if (is.na(tt)) next
      sub_tt <- sub[time_min == tt]
      if (length(unique(sub_tt$cytokine)) < 2) next

      design <- model.matrix(~ 0 + factor(cytokine, levels = c(cyA, cyB)),
                             data = sub_tt)
      colnames(design) <- c(cyA, cyB)

      Y <- t(as.matrix(sub_tt[, ..phos_cols]))
      fit <- lmFit(Y, design)
      cm <- makeContrasts(contrasts = paste0(cyA, "-", cyB),
                          levels = design)
      fit2 <- eBayes(contrasts.fit(fit, cm))

      tt_tab <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
      tt_tab$cluster  <- cl
      tt_tab$cytokine <- paste0(cyA, "_vs_", cyB)
      tt_tab$contrast <- paste0(cyA, "-", cyB, "_t", tt)
      tt_tab$time_min <- tt
      tt_tab$marker   <- rownames(tt_tab)
      res_all[[length(res_all)+1]] <- as.data.table(tt_tab)
    }
  }
}

###############################################################################
## 4. SAVE FULL LIMMA RESULTS
###############################################################################

if (length(res_all)) {
  lim_full <- rbindlist(res_all, use.names = TRUE, fill = TRUE)
  setcolorder(lim_full,
              c("cluster","cytokine","time_min","contrast","marker",
                "logFC","AveExpr","t","P.Value","adj.P.Val","B"))
  fwrite(lim_full, "limma_full_results.csv")
} else {
  warning("No limma results could be computed.")
}

cat("Full limma modeling completed.\n")
