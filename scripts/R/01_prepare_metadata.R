#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tools)
})

# ---------- CLI ----------
opt <- OptionParser(option_list = list(
  make_option("--files-list", type = "character", help = "Path to files.txt (non-anchor FCS basenames, one per line)."),
  make_option("--anchors-list", type = "character", help = "Path to anchors.txt (anchor FCS basenames, one per line)."),
  make_option("--base-dir", type = "character", help = "Base directory containing batch FCS files (and anchors/ subdir for anchors)."),
  make_option("--channels", type = "character", default = "channels_to_adjust.tsv", help = "Channels TSV to copy/normalize [default: %default]"),
  make_option("--metadata", type = "character", help = "OMIQ metadata CSV (must contain a 'Filename' column or similar)."),
  make_option("--out", type = "character", default = "sample_table.csv", help = "Output CSV for sample table [default: %default]"),
  make_option("--normalized-dir", type = "character", default = "normalized_fcs", help = "Output dir for canonical symlinks/copies [default: %default]"),
  make_option("--copy-files", action = "store_true", default = FALSE, help = "Copy files instead of creating symlinks.")
))
args <- parse_args(opt)

must <- c("files-list","anchors-list","base-dir","metadata")
miss <- must[!nzchar(unlist(args[must]))]
if (length(miss)) {
  stop("Missing required arguments: ", paste(miss, collapse=", "))
}

# ---------- Helpers ----------
read_list <- function(p) {
  x <- fread(p, header = FALSE, sep = "\n", col.names = "entry", data.table = TRUE, encoding = "UTF-8")
  x[, entry := trimws(entry)]
  x <- x[nchar(entry) > 0]
  x
}

safe_int <- function(x) {
  suppressWarnings(as.integer(x))
}

pad_min <- function(m) {
  ifelse(is.na(m), NA_character_, sprintf("%03dmin", as.integer(m)))
}

stim_token <- function(stim) {
  s <- tolower(trimws(stim))
  s <- gsub("\\s+", "", s)
  s <- gsub("\\+", "_", s)
  s <- gsub("^unstim.*$", "unstim", s)
  s
}

# Canonicalize detailed stimulation into tokens like: cd328_il2_il21_pac
# Accepts values from OMIQ (e.g., "CD3/28+IL21 PAC") or filenames.
canon_stimulation <- function(x) {
  if (is.null(x)) return(NA_character_)
  v <- as.character(x)
  out <- v
  for (i in seq_along(v)) {
    s <- tolower(trimws(v[i]))
    if (is.na(s) || !nzchar(s)) { out[i] <- NA_character_; next }
    s <- gsub("[()]+", " ", s)
    # Early exit for explicit unstim
    if (grepl("\\bunstim\\b", s)) { out[i] <- "unstim"; next }
    toks <- character()
    if (grepl("cd3[\\./ ]?28", s) || grepl("\\bcd328\\b", s)) toks <- c(toks, "cd328")
    else if (grepl("\\bcd3\\b", s)) toks <- c(toks, "cd3")
    if (grepl("(?<![A-Za-z0-9])(?:il\\s*2|il2)(?![0-9A-Za-z])", s, perl = TRUE))  toks <- c(toks, "il2")
    if (grepl("(?<![A-Za-z0-9])(?:il\\s*15|il15)(?![0-9A-Za-z])", s, perl = TRUE)) toks <- c(toks, "il15")
    if (grepl("(?<![A-Za-z0-9])(?:il\\s*21|il21)(?![0-9A-Za-z])", s, perl = TRUE)) toks <- c(toks, "il21")
    if (grepl("\\bpac\\b", s))   toks <- c(toks, "pac")
    out[i] <- if (length(toks) == 0) NA_character_ else paste(toks, collapse = "_")
  }
  out
}

make_stem <- function(is_anchor, batch, sample_id, plate, stim, time_min, anchor_k = NA_integer_) {
  b   <- as.character(batch)
  sid <- as.character(sample_id)
  pl  <- ifelse(is.na(plate) | plate=="", NA, as.integer(plate))
  tm  <- pad_min(time_min)
  st  <- stim_token(stim)

  if (isTRUE(is_anchor)) {
    ak <- ifelse(is.na(anchor_k), "", paste0("Anchor-", anchor_k, "_"))
    return(paste0("Norm_Batch-", b, "_", ak, "Sample-", sid, ".fcs"))
  } else {
    pl_s <- ifelse(is.na(pl), "", paste0("_P", pl))
    st_s <- ifelse(is.na(st) | st=="", "", paste0("_", st))
    tm_s <- ifelse(is.na(tm), "", paste0("_", tm))
    return(paste0("Norm_Batch-", b, "_Sample-", sid, pl_s, st_s, tm_s, ".fcs"))
  }
}

sanitize_fallback <- function(bn) {
  x <- bn
  x <- gsub("\\s+", "_", x)
  x <- gsub("[']", "", x)
  x
}

# Parse from filename (robust to your patterns)
parse_from_basename <- function(bn, is_anchor = FALSE) {
  # Initialize
  batch <- sample_id <- plate <- time_min <- stim <- NA

  # Batch & Sample
  m <- regexec("Batch[_-]([0-9]+).*?Sample[_-]?([0-9]+)", bn, ignore.case = TRUE)
  mt <- regmatches(bn, m)
  if (length(mt[[1]]) >= 3) {
    batch     <- safe_int(mt[[1]][2])
    sample_id <- safe_int(mt[[1]][3])
  }

  # Plate
  m <- regexec("P([0-9]+)", bn, ignore.case = TRUE)
  mt <- regmatches(bn, m)
  if (length(mt[[1]]) >= 2) plate <- safe_int(mt[[1]][2])

  # Time (min)
  m <- regexec("([0-9]+)\\s*min", bn, ignore.case = TRUE)
  mt <- regmatches(bn, m)
  if (length(mt[[1]]) >= 2) time_min <- safe_int(mt[[1]][2])

  # Unstim patterns (US-1, unstim-1, _unstim, -US-2, etc.)
  if (grepl("\\bunstim\\b", bn, ignore.case = TRUE) ||
      grepl("\\bUS-?[0-9]+\\b", bn, ignore.case = TRUE)) {
    stim <- "unstim"
    time_min <- NA_integer_
  }

  # CD3 vs CD3.28
  has_cd328 <- grepl("cd3[\\._]?28", bn, ignore.case = TRUE)
  has_cd3   <- grepl("cd3(?![\\._]?28)", bn, perl = TRUE, ignore.case = TRUE) | grepl("cd3_", bn, ignore.case = TRUE)
  # Cytokines & PAC
  has_il2  <- grepl("(?i)\\bIL\\s*2\\b|IL2", bn)
  has_il15 <- grepl("(?i)\\bIL\\s*15\\b|IL15", bn)
  has_il21 <- grepl("(?i)\\bIL\\s*21\\b|IL21", bn)
  has_pac  <- grepl("(?i)\\bPAC\\b|_PAC[_\\s]", paste0("_", bn, "_"))

  # Build stim
  if (is.na(stim)) {
    toks <- character()
    if (has_cd328) toks <- c(toks, "cd328")
    else if (has_cd3) toks <- c(toks, "cd3")
    if (has_il2)  toks <- c(toks, "il2")
    if (has_il15) toks <- c(toks, "il15")
    if (has_il21) toks <- c(toks, "il21")
    if (has_pac)  toks <- c(toks, "pac")
    if (length(toks) == 0) {
      # Leave NA, will try OMIQ rescue later; otherwise set to 'unstim' when we see US
      stim <- NA_character_
    } else {
      stim <- paste(toks, collapse = "+")
    }
  }

  # Anchors: force unstim & time NA
  if (is_anchor) {
    stim <- "unstim"
    time_min <- NA_integer_
  }

  data.table(batch = batch, sample_id = sample_id, plate = plate, time_min = time_min, stim = stim)
}

# Merge preferring OMIQ metadata values when present and non-empty
prefer <- function(a, b) {
  # prefer b when b is not NA / not empty
  ifelse(!is.na(b) & b != "", b, a)
}

# ---------- Input ----------
base_dir      <- normalizePath(args[["base-dir"]], mustWork = TRUE)
anchors_dir   <- file.path(base_dir, "anchors")

files_dt   <- read_list(args[["files-list"]])
anchors_dt <- read_list(args[["anchors-list"]])

# Resolve full paths
files_dt[, is_anchor := FALSE]
files_dt[, file_base := entry]
files_dt[, file := file.path(base_dir, file_base)]
files_dt[, exists := file.exists(file)]

anchors_dt[, is_anchor := TRUE]
anchors_dt[, file_base := entry]
anchors_dt[, file := file.path(anchors_dir, file_base)]
anchors_dt[, exists := file.exists(file)]

if (any(!files_dt$exists)) {
  cat("[WARN] These entries in files.txt could not be resolved:\n")
  for (x in files_dt[exists == FALSE, file_base]) cat(" - ", x, "\n", sep = "")
}
if (any(!anchors_dt$exists)) {
  cat("[WARN] These entries in anchors.txt could not be resolved:\n")
  for (x in anchors_dt[exists == FALSE, file_base]) cat(" - ", x, "\n", sep = "")
}

all_dt <- rbindlist(list(files_dt, anchors_dt), use.names = TRUE, fill = TRUE)
all_dt <- all_dt[exists == TRUE]

# ---------- Parse from basename ----------
all_dt[, bn := basename(file)]
parsed <- rbindlist(lapply(seq_len(nrow(all_dt)), function(i) {
  parse_from_basename(all_dt$bn[i], is_anchor = all_dt$is_anchor[i])
}), use.names = TRUE)

all_dt <- cbind(all_dt, parsed)

# ---------- Read OMIQ metadata ----------
meta <- tryCatch({
  fread(args[["metadata"]])
}, error = function(e) {
  stop("Failed to read metadata CSV: ", e$message)
})

## ---------- Normalize filenames to enable robust matching ----------
## OMIQ CSV uses filenames like:
##  - 20250917_pTCR_cHCV_Sample1_HBUF8157_US-1.fcs
##  - 20250917_pTCR_cHCV_SampleX_..._Anchor-1.fcs (sometimes prefixed by 'export')
## Local disk files can be prefixed with 'Batch_<n>_' and anchors can have
## 'Anchor-<k>_' at the start. We derive a core key by:
##   1) lowercasing and taking basename
##   2) removing leading 'batch_<n>_'
##   3) removing leading 'export'
##   4) removing any 'anchor-<k>' token regardless of position
## After (4) the remaining core should match between disk and OMIQ.
normalize_core_key <- function(x) {
  b <- tolower(basename(x))
  # drop leading Batch_<n>_
  b <- sub("^batch_([0-9]+)_", "", b, ignore.case = TRUE)
  # drop leading export
  b <- sub("^export", "", b, ignore.case = TRUE)
  # remove any anchor token anywhere (prefix/suffix); allow '-', '_' or space before number
  b <- gsub("anchor[- _]?[0-9]+", "", b, ignore.case = TRUE)
  # collapse duplicate underscores introduced by removals
  b <- gsub("__+", "_", b)
  # trim stray underscores left by removal
  b <- gsub("^_+|_+$", "", b)
  # remove underscore immediately before extension after anchor removal
  b <- sub("_+\\.fcs$", ".fcs", b, ignore.case = TRUE)
  b
}

# Find filename-like column in OMIQ CSV
fn_cols <- names(meta)[grepl("filename", names(meta), ignore.case = TRUE)]
if (length(fn_cols) == 0) {
  # Fallback: guess first string column
  txt_cols <- names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
  if (length(txt_cols) == 0) stop("No suitable 'Filename' column found in metadata.")
  fn_col <- txt_cols[1]
} else {
  fn_col <- fn_cols[1]
}

meta[, filename_basename := basename(get(fn_col))]
setnames(meta, fn_col, "OMIQ_Filename_Original")
meta[, key_core := normalize_core_key(filename_basename)]

# Normalize prospective columns from OMIQ (robust naming)
col_map <- list(
  batch     = c("Batch","batch"),
  sample_id = c("Sample","SampleID","Sample_Id","sample_id","Sample_ID","SampleNumber"),
  plate     = c("Plate","plate","WellPlate","P"),
  time_min  = c("Time","time","Time_min","TimeMin","Minutes","minute","min","timepoint","Timepoint"),
  stim      = c("Stim","stim","Stimulation","Condition","cond")
)

pick_col <- function(dt, candidates) {
  for (c in candidates) {
    if (c %in% names(dt)) return(c)
  }
  return(NA_character_)
}

meta_cols <- lapply(col_map, pick_col, dt = meta)
# Cast/clean OMIQ fields if present
if (!is.na(meta_cols$batch))     meta[, OMIQ_batch     := safe_int(get(meta_cols$batch))]
if (!is.na(meta_cols$sample_id)) meta[, OMIQ_sample_id := safe_int(get(meta_cols$sample_id))]
if (!is.na(meta_cols$plate))     meta[, OMIQ_plate     := safe_int(get(meta_cols$plate))]
if (!is.na(meta_cols$time_min))  {
  raw_time <- meta[[meta_cols$time_min]]
  if (is.numeric(raw_time)) {
    meta[, OMIQ_time_min := safe_int(raw_time)]
  } else {
    # extract first integer from strings like "0min", "60 min", or "x"
    v <- as.character(raw_time)
    digits <- suppressWarnings(as.integer(gsub("^.*?([0-9]+).*$", "\\1", v)))
    # if no digits, set NA
    digits[!grepl("[0-9]+", v)] <- NA_integer_
    meta[, OMIQ_time_min := digits]
  }
} else {
  meta[, OMIQ_time_min := NA_integer_]
}
if (!is.na(meta_cols$stim)) {
  meta[, OMIQ_stim_detail := canon_stimulation(meta[[meta_cols$stim]])]
} else {
  meta[, OMIQ_stim_detail := NA_character_]
}

# ---------- Ensure OMIQ_* columns exist even if missing in CSV ----------
need_cols <- c("OMIQ_batch","OMIQ_sample_id","OMIQ_plate","OMIQ_time_min","OMIQ_stim_detail")
for (nm in need_cols) {
  if (!nm %in% names(meta)) meta[, (nm) := NA]
}

# Type-safety (so merging + prefer() behave consistently)
if ("OMIQ_batch"     %in% names(meta))     meta[, OMIQ_batch     := safe_int(OMIQ_batch)]
if ("OMIQ_sample_id" %in% names(meta))     meta[, OMIQ_sample_id := safe_int(OMIQ_sample_id)]
if ("OMIQ_plate"     %in% names(meta))     meta[, OMIQ_plate     := safe_int(OMIQ_plate)]
if ("OMIQ_time_min"  %in% names(meta))     meta[, OMIQ_time_min  := safe_int(OMIQ_time_min)]
if ("OMIQ_stim_detail" %in% names(meta))   meta[, OMIQ_stim_detail := canon_stimulation(OMIQ_stim_detail)]

# ---------- Derive robust join key for local files ----------
setDT(all_dt)
setDT(meta)
all_dt[, key_core := normalize_core_key(bn)]
meta_small <- meta[, .(key_core, OMIQ_batch, OMIQ_sample_id, OMIQ_plate, OMIQ_time_min, OMIQ_stim_detail)]
# Deduplicate in case multiple OMIQ rows map to same core (e.g. 'export' variants)
meta_small <- unique(meta_small, by = "key_core")
all_dt <- merge(all_dt, meta_small, by = "key_core", all.x = TRUE)

# ---------- Prefer OMIQ where available ----------
all_dt[, batch     := prefer(batch,     OMIQ_batch)]
all_dt[, sample_id := prefer(sample_id, OMIQ_sample_id)]
all_dt[, plate     := prefer(plate,     OMIQ_plate)]
all_dt[, time_min  := prefer(time_min,  OMIQ_time_min)]

# Build detailed stimulation from filename fallback if OMIQ missing
all_dt[, stim_detail_fn := canon_stimulation(bn)]
all_dt[, stim_detail    := prefer(stim_detail_fn, OMIQ_stim_detail)]
# anchors are always unstim; fix any remaining NAs due to metadata typos
all_dt[is_anchor == TRUE & (is.na(stim_detail) | stim_detail == ""), stim_detail := "unstim"]

# If still missing stim, last-resort: set unstim for anchor or explicit US/unstim; else leave NA
# Class column: 'unstim' for anchors or files explicitly marked unstim; else 'stim'
all_dt[, stim_class := ifelse(is_anchor | (!is.na(stim_detail) & stim_detail == "unstim"), "unstim", "stim")]

# ---------- Build sample table ----------
sample_dt <- all_dt[, .(
  file = normalizePath(file),
  is_anchor,
  batch = safe_int(batch),
  sample_id = safe_int(sample_id),
  plate = safe_int(plate),
  time_min = safe_int(time_min),
  stim = stim_detail,                 # detailed stimulation tokens
  stim_class = stim_class             # broad class used by downstream scripts
)]

# Column order exactly as requested:
setcolorder(sample_dt, c("file","is_anchor","batch","sample_id","plate","time_min","stim","stim_class"))

# ---------- Channels TSV ----------
channels_out <- args[["channels"]]
if (file.exists(args[["channels"]])) {
  # keep as-is; ensure TSV delimiter is tab
  # (we don't modify — just confirm it exists)
  message("[ok] Using channels file: ", channels_out)
} else {
  message("[WARN] Channels file not found, creating placeholder at: ", channels_out)
  fwrite(data.table(channel = character(), action = character()), file = channels_out, sep = "\t")
}

# ---------- Write sample table ----------
fwrite(sample_dt, file = args[["out"]])
message("[ok] Wrote ", args[["out"]], " with ", nrow(sample_dt), " rows")

# ---------- Normalized outputs (symlinks or copies) ----------
out_dir <- args[["normalized-dir"]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# assign anchor index within each (batch, sample_id)
sample_dt[, anchor_k := NA_integer_]
sample_dt[is_anchor == TRUE, anchor_k := seq_len(.N), by = .(batch, sample_id)]

sample_dt[, out_stem := ifelse(
  is_anchor,
  make_stem(TRUE, batch, sample_id, plate, "unstim", NA, anchor_k),
  make_stem(FALSE, batch, sample_id, plate, stim, time_min)
)]
sample_dt[, out_path := file.path(out_dir, out_stem)]

# collision guard
dups <- sample_dt[duplicated(out_path) | duplicated(out_path, fromLast = TRUE)]
if (nrow(dups)) {
  message("[WARN] Name collisions detected for ", nrow(dups), " rows. Appending _repX.")
  sample_dt[, rep_idx := seq_len(.N), by = out_path]
  sample_dt[rep_idx > 1, out_path := sub("\\.fcs$", paste0("_rep", rep_idx, ".fcs"), out_path)]
}

# create links/copies
use_symlink <- !isTRUE(args[["copy-files"]])
for (i in seq_len(nrow(sample_dt))) {
  src <- sample_dt$file[i]
  dst <- sample_dt$out_path[i]
  if (!file.exists(src)) {
    warning("Missing source: ", src)
    next
  }
  if (file.exists(dst)) next
  if (use_symlink) {
    ok <- tryCatch({ file.symlink(src, dst); TRUE }, error = function(e) FALSE)
    if (!ok) invisible(file.copy(src, dst, overwrite = TRUE))
  } else {
    invisible(file.copy(src, dst, overwrite = TRUE))
  }
}

# ---------- Quick diagnostics ----------
cat("\n=== Quick diagnostics ===\n")
tb_batch <- sample_dt[!is.na(batch), .N, by = batch][order(batch)]
if (nrow(tb_batch)) {
  cat("\nBATCH\n\n")
  print(setNames(as.integer(tb_batch$N), tb_batch$batch))
} else {
  cat("\nBATCH\n\n<all NA>\n")
}

tb_stim <- sample_dt[, .N, by = stim][order(is.na(stim), stim)]
cat("\nSTIM\n\n")
if (nrow(tb_stim)) {
  print(setNames(as.integer(tb_stim$N), ifelse(is.na(tb_stim$stim), "<NA>", tb_stim$stim)))
} else {
  cat("<none>\n")
}

na_rows <- sample_dt[is.na(stim) | is.na(batch) | is.na(sample_id)]
if (nrow(na_rows)) {
  cat("\n[WARN] Rows with NA in stim/batch/sample_id (showing up to 20):\n")
  print(na_rows[1:min(20, .N), .(file_base = basename(file), batch, sample_id, stim)])
}

# ensure filenames have no spaces/quotes (safety pass)
old <- list.files(out_dir, full.names = TRUE)
bad <- old[grepl(" ", basename(old), fixed = TRUE) | grepl("'", basename(old), fixed = TRUE)]
if (length(bad)) {
  for (f in bad) {
    bn <- basename(f)
    new_bn <- sanitize_fallback(bn)
    new_f <- file.path(dirname(f), new_bn)
    if (!file.exists(new_f)) file.rename(f, new_f)
  }
  message("[ok] Renamed ", length(bad), " pre-existing normalized files to remove spaces/quotes.")
}

# Final echo
resolved_n  <- nrow(sample_dt)
anchors_n   <- sample_dt[is_anchor == TRUE, .N]
nonanchor_n <- sample_dt[is_anchor == FALSE, .N]
message("\n[ok] Resolved ", resolved_n, " files (", anchors_n, " anchors, ", nonanchor_n, " others)")
message("[ok] Normalized outputs in: ", out_dir)
