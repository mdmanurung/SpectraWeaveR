# Purpose: per-dataset config for benchmark 05. One copy of every script serves 3D, 24H and EXV; the
#          dataset comes from the required env var SW_DATASET (3d, 24h, exv). There is no default.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") near the top of a script,
#          after the library() calls. SW_DATASET_DIR (set by the sbatch wrappers) points at a frozen copy
#          of this folder; unset, it is the repo copy.
# Inputs:  <dataset dir>/<ds>.json         config values, each with a reason in "reasons"
#          <dataset dir>/<ds>_runs.csv     run registry: key, path, producer, approved_on
#          <dataset dir>/<ds>_samples.csv  sample, donor, condition, well, acq_batch, stem [00_samples.R]
# Outputs: none; defines DS (the config list) and the helpers below in the global environment.
#
# Helpers
#   cfg("name")        config value; stops when the field is missing or null (never a silent default)
#   data_path(...)     file.path(<data_dir>, ...)
#   run_in("key")      path of a registered run folder; stops when the key is unknown or the folder is missing
#   run_out("key")     new run folder experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-<key>[SW_RUN_SUFFIX]; stops if it exists;
#                      does not register it (registering is the user's approval step)
#   has_run("key")     TRUE when the key is registered (the folder may still be missing: run_in() checks that)
#   samples_tbl()      the sample table
#   cut_of(arm, sample, channel, kind, variant)   positive or negative cut of one cell; see 97e_common_gs.R.
#                      A rule of type "both" needs variant = one of its variants' names (e.g. "fixed15").

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

.sw05_dir <- Sys.getenv(
  "SW_DATASET_DIR",
  "experiments/benchmarks/05_bmv-unmix-compare/datasets"
)
.sw05_ds <- Sys.getenv("SW_DATASET")
if (!.sw05_ds %in% c("3d", "24h", "exv")) {
  stop(
    "SW_DATASET must be one of 3d, 24h, exv; got '",
    .sw05_ds,
    "'",
    call. = FALSE
  )
}
DS <- jsonlite::read_json(
  file.path(.sw05_dir, paste0(.sw05_ds, ".json")),
  simplifyVector = TRUE
)
stopifnot("config names itself" = identical(DS$dataset, .sw05_ds))

cfg <- function(name) {
  v <- DS[[name]]
  if (is.null(v)) {
    stop(
      "config '",
      name,
      "' is missing or null for dataset ",
      .sw05_ds,
      " (",
      file.path(.sw05_dir, paste0(.sw05_ds, ".json")),
      "); an impactful value needs the user's decision",
      call. = FALSE
    )
  }
  v
}

data_path <- function(...) file.path(cfg("data_dir"), ...)

.sw05_runs <- read_csv(
  file.path(.sw05_dir, paste0(.sw05_ds, "_runs.csv")),
  show_col_types = FALSE
)

has_run <- function(key) key %in% .sw05_runs$key

run_in <- function(key) {
  hit <- .sw05_runs |>
    dplyr::filter(.data$key == .env$key)
  if (nrow(hit) != 1) {
    stop(
      "run key '",
      key,
      "' is not registered once in ",
      .sw05_ds,
      "_runs.csv (found ",
      nrow(hit),
      ")",
      call. = FALSE
    )
  }
  if (!dir.exists(hit$path)) {
    stop("registered run folder is missing: ", hit$path, call. = FALSE)
  }
  hit$path
}

run_out <- function(key) {
  day <- Sys.getenv("SW_RUN_DATE")
  if (!grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", day)) {
    stop("SW_RUN_DATE must be YYYY-MM-DD", call. = FALSE)
  }
  out <- file.path(
    "experiments/runs",
    paste0(day, "_bmv", cfg("run_tag"), "-", key, Sys.getenv("SW_RUN_SUFFIX"))
  )
  if (dir.exists(out)) {
    stop("run dir exists: ", out, call. = FALSE)
  }
  out
}

samples_tbl <- function() {
  read_csv(
    file.path(.sw05_dir, paste0(.sw05_ds, "_samples.csv")),
    show_col_types = FALSE
  )
}

# Filled in by the 97e step (cuts.csv). A fixed rule is answered without the table.
# A rule of type "both" (24H negatives) has one named variant per cut; cuts.csv holds one column per
# variant (pos_cut, neg_cut_<variant>) and a call without variant stops rather than return NA.
cut_of <- function(
  arm,
  sample,
  channel,
  kind = c("pos", "neg"),
  variant = NULL
) {
  kind <- match.arg(kind)
  rule <- cfg(paste0(kind, "_rule"))
  col <- paste0(kind, "_cut")
  if (identical(rule$type, "both")) {
    if (is.null(variant) || !variant %in% rule$variants$name) {
      stop(
        "cut_of(kind = '",
        kind,
        "') needs variant = one of: ",
        paste(rule$variants$name, collapse = ", "),
        call. = FALSE
      )
    }
    v <- rule$variants[rule$variants$name == variant, ]
    if (identical(v$type, "fixed")) {
      return(v$value)
    }
    col <- paste0(col, "_", variant)
  } else if (identical(rule$type, "fixed")) {
    return(rule$value)
  }
  cuts <- read_csv(
    file.path(run_in("common-cuts"), "cuts.csv"),
    show_col_types = FALSE
  )
  hit <- cuts |>
    dplyr::filter(
      .data$arm == .env$arm,
      .data$sample == .env$sample,
      .data$channel == .env$channel
    )
  stopifnot("one cut per arm, sample and channel" = nrow(hit) == 1)
  if (is.na(hit[[col]])) {
    stop(
      "cuts.csv has NA in ",
      col,
      " for ",
      arm,
      " ",
      sample,
      " ",
      channel,
      call. = FALSE
    )
  }
  hit[[col]]
}

# Gate labels of the configured label run (cfg("label_key")), with the columns asked for.
# 3D labels (97d) sit on clean3 rows and carry raw_row; pooled AS 1.8.8 labels (99a) sit on the uncleaned rows of the
# arm files, so the row is the raw row and raw_row is made from it.
read_labels <- function(cols, key = cfg("label_key")) {
  f <- file.path(run_in(key), "labels.csv.gz")
  have <- names(readr::read_csv(f, n_max = 0, show_col_types = FALSE))
  pooled <- !"raw_row" %in% have
  want <- unique(if (pooled) sub("^raw_row$", "row", cols) else cols)
  types <- c(sample = "c", row = "d", raw_row = "d", FSC_A = "d", SSC_A = "d")
  spec <- stats::setNames(
    ifelse(want %in% names(types), types[want], "l"),
    want
  )
  d <- readr::read_csv(f, col_types = do.call(readr::cols_only, as.list(spec)))
  if (pooled && "raw_row" %in% cols) {
    d <- dplyr::rename(d, raw_row = row)
  }
  d
}
