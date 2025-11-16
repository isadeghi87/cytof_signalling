#!/usr/bin/env Rscript
# Consolidate CytoNorm outputs in normalized_fcs/ by mapping each sample in
# sample_table.csv to exactly one normalized FCS. For each sample, find files
# under normalized_fcs/ whose basename contains the original sample basename;
# select the largest candidate (by file size) and rename it to
# normalized_fcs/<original_basename>.fcs. Log successes and failures.

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
sample_table <- ifelse(length(args) >= 1, args[1], "sample_table.csv")
norm_dir <- ifelse(length(args) >= 2, args[2], "normalized_fcs")

if (!file.exists(sample_table)) stop("sample_table.csv not found: ", sample_table)
if (!dir.exists(norm_dir)) stop("normalized_fcs directory not found: ", norm_dir)

tab <- fread(sample_table)
if (!"file" %in% names(tab)) stop("sample_table.csv must contain a 'file' column with input paths")

all_norm <- list.files(norm_dir, full.names = TRUE)

log_msgs <- character()
failures <- 0L
renamed <- 0L
diag <- list()

# helper: normalize a string for tolerant matching
norm_str <- function(x) {
  x <- tolower(x)
  # replace common separators with underscore
  x <- gsub("[[:punct:]\\s]+", "_", x)
  # collapse underscores
  x <- gsub("_+", "_", x)
  x
}

# tokens: alphanumeric tokens from basename
tokens_of <- function(x) {
  toks <- unlist(strsplit(x, "[^a-zA-Z0-9]+"))
  toks[toks != ""]
}

for (i in seq_len(nrow(tab))) {
  orig <- tab$file[i]
  b <- basename(orig)
  b_norm <- norm_str(b)
  # find candidates containing the original basename
  cand <- all_norm[str_detect(tolower(basename(all_norm)), fixed(tolower(b)))]
  note <- NA_character_
  # if no direct contains-match, try tolerant strategies
  if (length(cand) == 0) {
    # token-based match: require at least one token match
    toks <- tokens_of(b_norm)
    if (length(toks) > 0) {
      cand_tokens <- sapply(basename(all_norm), function(x) any(toks %in% tokens_of(norm_str(x))))
      cand2 <- all_norm[which(cand_tokens)]
      if (length(cand2) > 0) {
        cand <- cand2
        note <- "token-match"
      }
    }
  }
  # relaxed substring match: try splitting sample basename into pieces and match any piece
  if (length(cand) == 0) {
    pieces <- unlist(strsplit(b_norm, "_+"))
    pieces <- pieces[pieces != ""]
    if (length(pieces) > 0) {
      cand_sub <- all_norm[sapply(basename(all_norm), function(x) any(sapply(pieces, function(p) grepl(p, norm_str(x), fixed = TRUE))))]
      if (length(cand_sub) > 0) {
        cand <- cand_sub
        note <- ifelse(is.na(note), "piece-substring", paste(note, "piece-substring", sep=","))
      }
    }
  }
  # final fallback: try relaxed glob-like: replace non-alnum in b with '.*' and grep
  if (length(cand) == 0) {
    pat <- gsub("[^a-zA-Z0-9]+", ".*", b_norm)
    cand_glob <- all_norm[grepl(pat, sapply(basename(all_norm), norm_str))]
    if (length(cand_glob) > 0) {
      cand <- cand_glob
      note <- ifelse(is.na(note), "regex-glob", paste(note, "regex-glob", sep=","))
    }
  }
  if (length(cand) == 0) {
    # record diagnostic row with no candidates
    msg <- sprintf("[WARN] No normalized candidate found for %s", b)
    message(msg)
    log_msgs <- c(log_msgs, msg)
    failures <- failures + 1L
    diag[[length(diag) + 1L]] <- list(sample_file = orig,
                                       sample_basename = b,
                                       candidates = paste(basename(all_norm), collapse = ";"),
                                       candidate_sizes = paste(file.info(all_norm)$size, collapse = ";"),
                                       chosen = NA_character_,
                                       outcome = "no_candidate",
                                       note = NA_character_)
    next
  }
  sizes <- file.info(cand)$size
  sel <- cand[which.max(sizes)]
  dest <- file.path(norm_dir, b)
  if (file.exists(dest)) {
    msg <- sprintf("[SKIP] Destination exists, skipping: %s", dest)
    message(msg)
    log_msgs <- c(log_msgs, msg)
    next
  }
  ok <- file.rename(from = sel, to = dest)
  if (isTRUE(ok)) {
    msg <- sprintf("[OK] Renamed %s -> %s", basename(sel), basename(dest))
    message(msg)
    log_msgs <- c(log_msgs, msg)
    renamed <- renamed + 1L
    # remove sel from all_norm so we don't pick it again
    all_norm <- setdiff(all_norm, sel)
  } else {
    msg <- sprintf("[ERR] Failed to rename %s -> %s", sel, dest)
    message(msg)
    log_msgs <- c(log_msgs, msg)
    failures <- failures + 1L
    diag[[length(diag) + 1L]] <- list(sample_file = orig,
                                       sample_basename = b,
                                       candidates = paste(basename(cand), collapse = ";"),
                                       candidate_sizes = paste(sizes, collapse = ";"),
                                       chosen = basename(sel),
                                       outcome = "rename_failed",
                                       note = note)
  }
}

summary_msg <- sprintf("Consolidation finished: renamed=%d, failures=%d", renamed, failures)
message(summary_msg)
invisible({
  # write diagnostics to CSV for later review
  if (length(diag) > 0) {
    dt <- rbindlist(lapply(diag, as.data.table), fill = TRUE)
    fwrite(dt, file = file.path(norm_dir, "consolidation_diagnostics.csv"))
    message(sprintf("[ok] Wrote consolidation diagnostics to %s", file.path(norm_dir, "consolidation_diagnostics.csv")))
  }
  list(log = log_msgs, renamed = renamed, failures = failures)
})
