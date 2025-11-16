#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(limma)
})

opt_list <- list(
  make_option("--medians", type="character",
              default="phospho_medians/phos_medians_by_sample_cluster.csv"),
  make_option("--outdir", type="character", default="diff_phos_results")
)

opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

dt <- fread(opt$medians)

# ensure expected identifier columns exist
needed_id <- c("file_base","cluster","time_min","cytokine")
missing_id <- setdiff(needed_id, colnames(dt))
if (length(missing_id)) {
  stop("Medians table is missing required columns: ",
       paste(missing_id, collapse=", "))
}

# write a copy of medians into the diff_phos_results folder for step 8
setcolorder(dt, c("file_base","cluster","time_min","cytokine",
                  setdiff(colnames(dt), c("file_base","cluster","time_min","cytokine"))))
fwrite(dt, file.path(opt$outdir, "phos_medians_by_sample_cluster.csv"))

phos_cols <- setdiff(colnames(dt), c("file_base","cluster","time_min","cytokine"))

res_all <- list()

clusters <- sort(unique(dt$cluster))
cyt_list <- sort(unique(dt$cytokine))

for (cl in clusters) {
  for (cy in cyt_list) {

    sub <- dt[cluster == cl & cytokine == cy]

    # need enough samples per condition
    if (nrow(sub) < 3) next

    tlev <- sort(unique(sub$time_min))
    tlev <- as.character(tlev)

    # require at least 2 timepoints
    if (length(tlev) < 2) next

    base <- if ("0" %in% tlev) "0" else tlev[1]

    design <- model.matrix(~ 0 + factor(time_min, levels=tlev), data=sub)
    colnames(design) <- paste0("t",tlev)

    Y <- t(as.matrix(sub[, ..phos_cols]))

    fit <- lmFit(Y, design)

    others <- setdiff(tlev, base)
    if (length(others)==0) next

    contr <- paste(sprintf("t%s - t%s", others, base), collapse=";")
    cm <- makeContrasts(contrasts=contr, levels=design)
    fit2 <- eBayes(contrasts.fit(fit, cm))

    for (j in seq_along(others)) {
      tt <- tryCatch(
        topTable(fit2, coef=j, number=Inf, sort.by="none"),
        error = function(e) {
          message("Skipping contrast for cluster=", cl,
                  " cytokine=", cy,
                  " coef index=", j,
                  " due to error: ", e$message)
          NULL
        }
      )
      if (is.null(tt)) next
      tt$cluster <- cl
      tt$cytokine <- cy
      tt$contrast <- paste0("t",others[j],"-t",base)
      tt$marker <- rownames(tt)
      res_all[[length(res_all)+1]] <- tt
    }

  }
}

if (length(res_all) == 0) {
  message("No valid contrasts could be computed; writing empty limma results.")
  out <- data.table(
    logFC = numeric(),
    AveExpr = numeric(),
    t = numeric(),
    P.Value = numeric(),
    adj.P.Val = numeric(),
    B = numeric(),
    cluster = integer(),
    cytokine = character(),
    contrast = character(),
    marker = character()
  )
} else {
  out <- rbindlist(res_all, fill=TRUE)
}

fwrite(out, file.path(opt$outdir,"diff_phospho_limma_results.csv"))

cat("Done.\n")
