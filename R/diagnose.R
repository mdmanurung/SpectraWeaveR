# R/diagnose.R
# One call that runs the file-level and batch-level diagnostics and stacks their
# numbers into a single table. It only calls the modules; every setting stays a
# required argument of the module it belongs to.

#' Run the technical diagnostics and collect their numbers in one table
#'
#' @description
#' Run any of eight diagnostic modules and stack what they return into one tidy
#' table (`module`, `subject`, `check`, `item`, `value`, `text`). A module runs
#' only when you give it a settings list; pass `NULL` to skip it. Nothing is
#' marked pass or fail. Use [sw_diagnosis_top()] to see the files or batches
#' with the largest value per check.
#'
#' | Module | Runs | Needs |
#' |---|---|---|
#' | `header` | [sw_fcs_header_audit()] | FCS paths on disk |
#' | `signal` | [sw_fcs_signal_audit()] | `data` |
#' | `anomaly` | [sw_anomaly_features()] + [sw_anomaly_flag_outlier()] (+ [sw_anomaly_flag_novelty()]) | `data` |
#' | `batch` | [sw_batch_scores()] | `data` and batch labels |
#' | `cluster` | [sw_cluster_scores()] | `data` |
#' | `relativise` | [sw_relativise_diagnose()] | a fitted relativisation model |
#' | `otd` | [sw_otd_scores()] | `data` |
#' | `negative` | [sw_negative_audit()] | `data` |
#'
#' @details
#' **No defaults anywhere.** Each settings list holds the arguments of its
#' module, named as in that module's help page, and each is checked by that
#' module. A missing argument stops with the module's own message. The
#' `anomaly` list also takes the forest arguments `ntrees`, `sample_size`,
#' `ndim` and `seed`; its `cutoff` is not needed because only the score is
#' kept. Add `novelty_reference` and `novelty_quantile` to also compute novelty
#' log-densities.
#'
#' **Modules are independent.** A fault one module misses may show in another.
#' Run all of them when you can.
#'
#' @param data A named list of per-file numeric matrices (event data on the
#'   analysis scale), or `NULL` when only `header` is used.
#' @param header `NULL` or a list with `paths`, `keywords`, `reference`.
#' @param signal `NULL` or a list with `channels`, `range_max`, `range_min`,
#'   `time_channel`, `time_bin`.
#' @param anomaly `NULL` or a list with `channels`, `features`,
#'   `quantile_step`, `n_events`, `seed`, `scale`, `aggregate_from`,
#'   `agg_size`, `reference` (when needed), `ntrees`, `sample_size`, `ndim`.
#' @param batch `NULL` or a list with `batch`, `markers`, `n_events`, `seed`,
#'   `order`.
#' @param cluster `NULL` or a list with `markers`, `k`, `n_events`, `seed`,
#'   `nstart`, `min_cluster_events`.
#' @param relativise `NULL` or a list with `model` (from
#'   [sw_relativise_fit()]), `batch`, `order`. Subjects are units (samples or
#'   batches); it does not use `data`.
#' @param otd `NULL` or a list with `markers`, `n_events`, `reference_events`,
#'   `blur`, `seed`, `max_iter`, `tol`.
#' @param negative `NULL` or a list with `channels`, `n_events`, `seed`,
#'   `hyper_k`, `pairs`.
#'
#' @return A tibble of class `sw_diagnosis` with columns `module`, `subject`
#'   (file, batch or batch pair), `check`, `item` (channel or marker, else
#'   `NA`), `value` (numeric) and `text` (extra detail). Attribute `modules`
#'   lists the modules that ran.
#'
#' @seealso [sw_diagnosis_top()]
#'
#' @examples
#' \dontrun{
#' dx <- sw_diagnose(
#'   data = files,
#'   header = NULL,
#'   signal = list(channels = chans, range_max = 6, range_min = NULL,
#'                 time_channel = NULL, time_bin = NULL),
#'   anomaly = NULL,
#'   batch = list(batch = batch, markers = chans, n_events = 500, seed = 1, order = NULL),
#'   cluster = NULL,
#'   relativise = NULL,
#'   otd = NULL,
#'   negative = NULL
#' )
#' sw_diagnosis_top(dx, n = 3)
#' }
#'
#' @export
sw_diagnose <- function(
  data,
  header,
  signal,
  anomaly,
  batch,
  cluster,
  relativise,
  otd,
  negative
) {
  .sw_stop_missing(c(
    data = missing(data),
    header = missing(header),
    signal = missing(signal),
    anomaly = missing(anomaly),
    batch = missing(batch),
    cluster = missing(cluster),
    relativise = missing(relativise),
    otd = missing(otd),
    negative = missing(negative)
  ))
  modules <- c(
    "header",
    "signal",
    "anomaly",
    "batch",
    "cluster",
    "relativise",
    "otd",
    "negative"
  )
  for (nm in modules) {
    v <- get(nm)
    if (!is.null(v) && !is.list(v)) {
      stop("`", nm, "` must be NULL or a list of settings.", call. = FALSE)
    }
  }
  if (all(vapply(modules, \(nm) is.null(get(nm)), logical(1)))) {
    stop("Give at least one module a settings list.", call. = FALSE)
  }
  if (
    is.null(data) &&
      (!is.null(signal) ||
        !is.null(anomaly) ||
        !is.null(batch) ||
        !is.null(cluster) ||
        !is.null(otd) ||
        !is.null(negative))
  ) {
    stop(
      "`data` is needed for the signal, anomaly, batch, cluster, otd and negative modules.",
      call. = FALSE
    )
  }

  parts <- list()
  ran <- character(0)

  if (!is.null(header)) {
    h <- do.call(sw_fcs_header_audit, header)
    parts$header <- tibble::tibble(
      module = "header",
      subject = h$file,
      check = paste(h$kind, h$key),
      item = h$channel,
      value = as.numeric(h$status != "match"),
      text = paste0(h$status, ": ", h$value, " vs ", h$reference_value)
    )
    ran <- c(ran, "header")
  }

  if (!is.null(signal)) {
    s <- do.call(sw_fcs_signal_audit, c(list(data = data), signal))
    parts$signal <- tibble::tibble(
      module = "signal",
      subject = s$file,
      check = s$check,
      item = s$channel,
      value = s$value,
      text = s$detail
    )
    ran <- c(ran, "signal")
  }

  if (!is.null(anomaly)) {
    forest <- c("ntrees", "sample_size", "ndim", "seed")
    novelty <- c("novelty_reference", "novelty_quantile")
    # `seed` feeds both the event draw and the forest, so it stays in both.
    feat_args <- anomaly[setdiff(
      names(anomaly),
      c(setdiff(forest, "seed"), novelty)
    )]
    feats <- do.call(sw_anomaly_features, c(list(data = data), feat_args))
    # The cutoff only sets a flag; the score is what this table keeps.
    fl <- do.call(
      sw_anomaly_flag_outlier,
      c(
        list(features = feats, cutoff = 0.5),
        anomaly[intersect(names(anomaly), forest)]
      )
    )
    parts$anomaly <- tibble::tibble(
      module = "anomaly",
      subject = fl$file,
      check = "outlier_score",
      item = NA_character_,
      value = fl$score,
      text = NA_character_
    )
    if (!is.null(anomaly$novelty_reference)) {
      nv <- sw_anomaly_flag_novelty(
        feats,
        reference = anomaly$novelty_reference,
        quantile = anomaly$novelty_quantile
      )
      parts$novelty <- tibble::tibble(
        module = "anomaly",
        subject = nv$file,
        check = "novelty_log_density",
        item = NA_character_,
        value = nv$log_density,
        text = nv$role
      )
    }
    ran <- c(ran, "anomaly")
  }

  if (!is.null(batch)) {
    b <- do.call(sw_batch_scores, c(list(data = data), batch))
    parts$batch <- tibble::tibble(
      module = "batch",
      subject = ifelse(is.na(b$batch), "all batches", b$batch),
      check = b$check,
      item = b$marker,
      value = b$value,
      text = paste0("n_files=", b$n_files)
    )
    ran <- c(ran, "batch")
  }

  if (!is.null(cluster)) {
    cl <- do.call(sw_cluster_scores, c(list(data = data), cluster))
    parts$cluster <- tibble::tibble(
      module = "cluster",
      subject = cl$file,
      check = cl$check,
      item = cl$marker,
      value = cl$value,
      text = paste0("cluster ", cl$cluster, "; n_events=", cl$n_events)
    )
    ran <- c(ran, "cluster")
  }

  if (!is.null(relativise)) {
    r <- do.call(sw_relativise_diagnose, relativise)
    parts$relativise <- tibble::tibble(
      module = "relativise",
      subject = ifelse(is.na(r$unit), "all units", r$unit),
      check = r$check,
      item = r$channel,
      value = r$value,
      text = NA_character_
    )
    ran <- c(ran, "relativise")
  }

  if (!is.null(otd)) {
    o <- do.call(sw_otd_scores, c(list(data = data), otd))
    parts$otd <- tibble::tibble(
      module = "otd",
      subject = o$file,
      check = "otd",
      item = NA_character_,
      value = o$otd,
      text = paste0("converged=", o$converged)
    )
    ran <- c(ran, "otd")
  }

  if (!is.null(negative)) {
    ng <- do.call(sw_negative_audit, c(list(data = data), negative))
    parts$negative <- tibble::tibble(
      module = "negative",
      subject = ng$file,
      check = ng$check,
      item = ng$channel,
      value = ng$value,
      text = ng$detail
    )
    ran <- c(ran, "negative")
  }

  out <- dplyr::bind_rows(parts)
  attr(out, "modules") <- ran
  class(out) <- c("sw_diagnosis", class(out))
  out
}

# Checks whose sign only gives a direction: the size of the deviation is what
# ranks. Novelty log-density is the other way round: low means unusual.
.sw_signed_checks <- c(
  "median_z",
  "sd_z",
  "freq_z",
  "time_trend",
  "order_rho",
  "gain_z",
  "gain_order_rho",
  "neg_median_z",
  "neg_sd_z",
  "neg_corr_z"
)
.sw_diag_key <- function(check, value) {
  # A gain ratio departs from 1 either way, so its log size ranks.
  ifelse(
    check == "novelty_log_density",
    -value,
    ifelse(
      check %in% .sw_signed_checks,
      abs(value),
      ifelse(check == "gain_batch_ratio", abs(log(value)), value)
    )
  )
}

#' Largest values per check in a diagnosis
#'
#' For each module and check, list the `n` subjects (files, batches or
#' markers) with the largest value. No cutoff is applied. Signed checks
#' (`median_z`, `sd_z`, `freq_z`, `time_trend`, `order_rho`) are ranked by the
#' size of the value, ignoring its sign. For `novelty_log_density` the smallest
#' values are the unusual ones, so that check is ranked from the lowest value.
#'
#' @param diagnosis An `sw_diagnosis` result.
#' @param n Whole number, rows kept per check.
#'
#' @return A tibble with the same columns as the diagnosis plus `rank`.
#'
#' @examples
#' \dontrun{
#' sw_diagnosis_top(dx, n = 3)
#' }
#'
#' @export
sw_diagnosis_top <- function(diagnosis, n) {
  if (missing(n)) {
    stop("`n` is required and has no default.", call. = FALSE)
  }
  if (!inherits(diagnosis, "sw_diagnosis")) {
    stop("`diagnosis` must come from sw_diagnose().", call. = FALSE)
  }
  .sw_check_count(n, "n")
  d <- diagnosis[!is.na(diagnosis$value), ]
  d$.key <- .sw_diag_key(d$check, d$value)
  # Group by module, check and item (NA item = file-level check).
  d$.group <- paste(
    d$module,
    d$check,
    ifelse(is.na(d$item), "", d$item),
    sep = "|"
  )
  d <- d[order(d$.group, -d$.key), ]
  d$rank <- stats::ave(seq_len(nrow(d)), d$.group, FUN = seq_along)
  d <- d[d$rank <= n, ]
  d$.group <- NULL
  tibble::as_tibble(d[setdiff(names(d), ".key")])
}
