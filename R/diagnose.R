# R/diagnose.R
# One call that runs the file-level and batch-level diagnostics and stacks their
# numbers into a single table. It only calls the modules; every setting stays a
# required argument of the module it belongs to.

#' Run the technical diagnostics and collect their numbers in one table
#'
#' @description
#' Run any of fourteen diagnostic modules and stack what they return into one tidy
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
#' | `gate` | [sw_gate_qc()] | `data` and `sw_gate` objects |
#' | `reference` | [sw_ref_consistency()] | `data` and a manifest with reference ids |
#' | `batch_pop` | [sw_batch_scores()] with `population` | `data`, batch labels and population labels |
#' | `flag_iqr` | [sw_flag_iqr()] | event-level `cells` (not `data`) |
#' | `batch_otd` | [sw_batch_otd()] | `data` and batch labels |
#' | `batch_mixing` | [sw_batch_mixing()] | `data` and batch labels |
#'
#' @details
#' **No defaults anywhere.** Each settings list holds the arguments of its
#' module, named as in that module's help page, and each is checked by that
#' module. A missing argument stops with the module's own message. The
#' `anomaly` list also takes the forest arguments `ntrees`, `sample_size`,
#' `ndim` and `seed`; its `cutoff` is not needed because only the score is
#' kept. Add `novelty_reference` and `novelty_quantile` to also compute novelty
#' log-densities. When `features` includes `"emd"`, the module also returns
#' `emd_z`: per channel, the robust z over files of each file's EMD to the
#' pooled sample. It finds a marker split into two modes better than `sd_z`.
#' With a `reference_group` in the `anomaly` list, the pooled sample and the z
#' are both taken within each file's group.
#'
#' **Modules are independent.** A fault one module misses may show in another.
#' Run all of them when you can.
#'
#' **Import tables are report-only.** The tables from the import step
#' (integrity, signatures, compensation, batch dates, dye map, bind) have no
#' check, batch and value layout, so this function does not stack them. Read
#' them from the import report.
#'
#' @param data A named list of per-file numeric matrices (event data on the
#'   analysis scale), or `NULL` when only `header` is used.
#' @param header `NULL` or a list with `paths`, `keywords`, `reference`.
#' @param signal `NULL` or a list with `channels`, `range_max`, `range_min`,
#'   `time_channel`, `time_bin`, `area_height`, `area_height_quantile`,
#'   `reference_group`, `reference`, `summary` (`NULL` or a list; see
#'   [sw_fcs_signal_audit()]).
#' @param anomaly `NULL` or a list with `channels`, `features`,
#'   `quantile_step`, `n_events`, `seed`, `scale`, `aggregate_from`,
#'   `agg_size`, `reference` (when needed), `reference_group`, `ntrees`,
#'   `sample_size`, `ndim`.
#' @param batch `NULL` or a list with `batch`, `markers`, `n_events`, `seed`,
#'   `order`, `reference_group`. A `population` entry is an error: this module
#'   stacks all-event rows only, so call [sw_batch_scores()] for population rows.
#' @param cluster `NULL` or a list with `markers`, `k`, `n_events`, `seed`,
#'   `nstart`, `min_cluster_events`, `reference`, `reference_group`.
#' @param relativise `NULL` or a list with `model` (from
#'   [sw_relativise_fit()]), `batch`, `order`. Subjects are units (samples or
#'   batches); it does not use `data`.
#' @param otd `NULL` or a list with `markers`, `n_events`, `reference_events`,
#'   `blur`, `seed`, `max_iter`, `tol`.
#' @param negative `NULL` or a list with `channels`, `n_events`, `seed`,
#'   `hyper_k`, `pairs`, `spread_quantile`, `min_separation`.
#' @param gate `NULL` or a list with `gates`, `reference`, `edge_width`. The
#'   `item` column is the gate id, with the channel after a `/` for
#'   per-channel checks.
#' @param reference `NULL` or a list with `manifest`, `markers`, `n_events`,
#'   `seed`, `order`, `min_batches_per_reference`, `min_events`, `roles`,
#'   `transform`, `cofactor`. Subjects are batches or batch pairs; `item` is the
#'   reference id and the marker, so each reference and marker is ranked on its
#'   own. `count_ratio` rows have no marker.
#' @param batch_pop `NULL` or a list with the arguments of
#'   [sw_batch_scores()], including a non-`NULL` `population` list and
#'   `min_events` and `group_summary`. Only rows inside a population are
#'   stacked (the all-event rows are the `batch` module's), plus the
#'   `batch_order_rho` rows. `item` is the marker and the population, or the
#'   population alone for `pop_share` and `share_rest_z`.
#' @param flag_iqr `NULL` or a list with the arguments of [sw_flag_iqr()]
#'   (`cells`, `markers`, `thresholds`, `iqr_multiplier`, `min_cells_per_side`,
#'   `quantile_type`, `min_effect`, and optionally `batch_col`, `control_col`,
#'   `sample_col`). The stacked value is `effect`, the distance beyond the
#'   Tukey fence; the fence flag itself is not stacked. `item` is the control,
#'   metric and marker.
#' @param batch_otd `NULL` or a list with `batch`, `markers`, `n_events`,
#'   `blur`, `seed`, `max_iter`, `tol`, `mode`. The check is `otd_` plus the
#'   mode.
#' @param batch_mixing `NULL` or a list with `batch`, `markers`, `k`,
#'   `perplexity`, `alpha`, `n_events`, `seed`, `min_events`, `population`.
#'   `item` is the population (`NA` without one).
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
#'                 time_channel = NULL, time_bin = NULL,
#'                 area_height = NULL, area_height_quantile = NULL,
#'                 reference_group = NULL, reference = NULL, summary = NULL),
#'   anomaly = NULL,
#'   batch = list(batch = batch, markers = chans, n_events = 500, seed = 1, order = NULL,
#'                reference_group = NULL),
#'   cluster = NULL,
#'   relativise = NULL,
#'   otd = NULL,
#'   negative = NULL,
#'   gate = NULL,
#'   reference = NULL,
#'   batch_pop = NULL,
#'   flag_iqr = NULL,
#'   batch_otd = NULL,
#'   batch_mixing = NULL
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
  negative,
  gate,
  reference,
  batch_pop,
  flag_iqr,
  batch_otd,
  batch_mixing
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
    negative = missing(negative),
    gate = missing(gate),
    reference = missing(reference),
    batch_pop = missing(batch_pop),
    flag_iqr = missing(flag_iqr),
    batch_otd = missing(batch_otd),
    batch_mixing = missing(batch_mixing)
  ))
  modules <- c(
    "header",
    "signal",
    "anomaly",
    "batch",
    "cluster",
    "relativise",
    "otd",
    "negative",
    "gate",
    "reference",
    "batch_pop",
    "flag_iqr",
    "batch_otd",
    "batch_mixing"
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
        !is.null(negative) ||
        !is.null(gate) ||
        !is.null(reference) ||
        !is.null(batch_pop) ||
        !is.null(batch_otd) ||
        !is.null(batch_mixing))
  ) {
    stop(
      "`data` is needed for the signal, anomaly, batch, cluster, otd, negative, gate, ",
      "reference, batch_pop, batch_otd and batch_mixing modules (not header or flag_iqr).",
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
    # Per-channel EMD to the pooled sample, as a robust z over files. Finds a
    # split into two modes better than sd_z (real files, ±0.5 SD: 0.70 vs 0.54).
    emd_cols <- grep("_emd$", names(feats), value = TRUE)
    if (length(emd_cols) > 0) {
      # With reference_group, the EMD is already measured against the file's own
      # group pool, and the z is taken within the group.
      emd_group <- .sw_check_reference_group(
        anomaly$reference_group,
        feats$file
      )
      emd_z <- lapply(emd_cols, function(col) {
        z <- .sw_robust_z(feats[[col]], emd_group)
        tibble::tibble(
          module = "anomaly",
          subject = feats$file,
          check = "emd_z",
          item = sub("_emd$", "", col),
          value = z,
          text = NA_character_
        )
      })
      parts$emd <- dplyr::bind_rows(emd_z)
    }
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
    # Population rows are not stacked here: this module runs the all-event scores.
    if (!is.null(batch[["population"]])) {
      stop(
        "`batch$population` is not supported: sw_diagnose() stacks all-event ",
        "batch rows only. Call sw_batch_scores() directly for population rows.",
        call. = FALSE
      )
    }
    b <- do.call(
      sw_batch_scores,
      c(list(data = data), batch, list(population = NULL))
    )
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
      # One ranking per marker and cluster: pooling a marker's clusters into
      # one ranking lost a shifted rare cluster on real files (AUROC 0.88 vs
      # 0.95 on the cluster's own row).
      item = ifelse(
        is.na(cl$marker),
        paste0("cluster ", cl$cluster),
        paste0(cl$marker, " / cluster ", cl$cluster)
      ),
      value = cl$value,
      text = paste0("n_events=", cl$n_events)
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

  if (!is.null(gate)) {
    gq <- do.call(sw_gate_qc, c(list(data = data), gate))
    parts$gate <- tibble::tibble(
      module = "gate",
      subject = gq$sample,
      check = gq$check,
      item = ifelse(
        is.na(gq$channel),
        gq$gate,
        paste0(gq$gate, "/", gq$channel)
      ),
      value = gq$value,
      text = NA_character_
    )
    ran <- c(ran, "gate")
  }

  if (!is.null(reference)) {
    rc <- do.call(sw_ref_consistency, c(list(data = data), reference))
    parts$reference <- tibble::tibble(
      module = "reference",
      subject = ifelse(is.na(rc$batch), "all batches", rc$batch),
      check = rc$check,
      # One ranking per reference and marker (count_ratio has no marker).
      item = ifelse(
        is.na(rc$marker),
        rc$reference_id,
        paste0(rc$reference_id, " / ", rc$marker)
      ),
      value = rc$value,
      text = ifelse(
        rc$status == "assessed",
        paste0("n_files=", rc$n_files),
        paste0("not_assessed: ", rc$reason, "; n_files=", rc$n_files)
      )
    )
    ran <- c(ran, "reference")
  }

  if (!is.null(batch_pop)) {
    if (is.null(batch_pop[["population"]])) {
      stop(
        "`batch_pop$population` is needed; without it the scores equal the ",
        "`batch` module.",
        call. = FALSE
      )
    }
    bp <- do.call(sw_batch_scores, c(list(data = data), batch_pop))
    # The all-event rows belong to the `batch` module; batch_order_rho is new.
    bp <- bp[!is.na(bp$population) | bp$check == "batch_order_rho", ]
    parts$batch_pop <- tibble::tibble(
      module = "batch_pop",
      subject = ifelse(is.na(bp$batch), "all batches", bp$batch),
      check = bp$check,
      # One ranking per marker and population.
      item = dplyr::case_when(
        is.na(bp$population) ~ bp$marker,
        is.na(bp$marker) ~ bp$population,
        .default = paste0(bp$marker, " / ", bp$population)
      ),
      value = bp$value,
      # An assessed row can still have no value (a robust z needs 3 batches and
      # a non-zero spread); say so, so no NA is silent.
      text = ifelse(
        bp$status == "assessed",
        paste0(
          ifelse(is.na(bp$value), "no value (assessed row); ", ""),
          "n_files=",
          bp$n_files
        ),
        paste0(
          "not_assessed: n_files=",
          bp$n_files,
          ifelse(is.na(bp$file), "", paste0("; file=", bp$file))
        )
      )
    )
    ran <- c(ran, "batch_pop")
  }

  if (!is.null(flag_iqr)) {
    fl_iqr <- do.call(sw_flag_iqr, flag_iqr)
    parts$flag_iqr <- tibble::tibble(
      module = "flag_iqr",
      subject = fl_iqr$batch,
      check = "iqr_effect",
      # One ranking per control, metric and marker; the fence flag is not kept.
      item = paste(fl_iqr$control, fl_iqr$metric, fl_iqr$marker, sep = " / "),
      value = fl_iqr$effect,
      text = paste0("status=", fl_iqr$status, "; n_cells=", fl_iqr$n_cells)
    )
    ran <- c(ran, "flag_iqr")
  }

  if (!is.null(batch_otd)) {
    bo <- do.call(sw_batch_otd, c(list(data = data), batch_otd))
    parts$batch_otd <- tibble::tibble(
      module = "batch_otd",
      subject = bo$batch,
      check = paste0("otd_", bo$mode),
      item = NA_character_,
      value = bo$otd,
      text = paste0("converged=", bo$converged, "; n_files=", bo$n_files)
    )
    ran <- c(ran, "batch_otd")
  }

  if (!is.null(batch_mixing)) {
    bm <- do.call(sw_batch_mixing, c(list(data = data), batch_mixing))
    parts$batch_mixing <- tibble::tibble(
      module = "batch_mixing",
      subject = ifelse(is.na(bm$batch), "all batches", bm$batch),
      check = bm$check,
      item = bm$population,
      value = bm$value,
      text = paste0(
        "n_events=",
        bm$n_events,
        "; n_batches=",
        bm$n_batches,
        "; status=",
        bm$status
      )
    )
    ran <- c(ran, "batch_mixing")
  }

  out <- dplyr::bind_rows(parts)
  attr(out, "modules") <- ran
  class(out) <- c("sw_diagnosis", class(out))
  out
}

# Checks whose sign only gives a direction: the size of the deviation is what
# ranks. Some checks are the other way round: low means unusual.
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
  "neg_corr",
  "neg_corr_z",
  "corr_pair_z",
  "global_sd_z",
  "spill_corr",
  "spill_corr_z",
  "kept_share_z",
  "kept_share_log_ratio",
  "centroid_shift",
  # A split blurs two populations, so the bimodality coefficient drops (G38).
  "bimodality_z",
  # Batch-level checks that depart either way (DR33).
  "median_shift_z",
  "share_rest_z",
  "batch_order_rho"
)
# Not here on purpose (G38): `spread_z` and `rest_emd_z` (also in the reference
# and batch_pop modules) rank by their raw value,
# because spillover spread only adds width and an EMD only grows with a fault; a
# negative z means a file or batch is closer to normal than the others.
# Novelty log-density: low is unusual. GMM separation: a small value means the
# channel did not split into a negative and a positive, so it ranks first.
# LISI: a low inverse Simpson index means the batches do not mix, so it ranks
# first. `kbet_rejection` and `otd_*` grow with a fault and rank by raw value.
.sw_low_first_checks <- c("novelty_log_density", "separation", "lisi_batch")
.sw_diag_key <- function(check, value) {
  key <- value
  signed <- check %in% .sw_signed_checks
  key[signed] <- abs(value[signed])
  low <- check %in% .sw_low_first_checks
  key[low] <- -value[low]
  # A ratio departs from 1 either way, so its log size ranks: a narrowed batch
  # (0.5) ranks as high as a widened one (2). Only these rows are logged, so
  # negative values elsewhere raise no warning.
  ratio <- check %in% c("gain_batch_ratio", "spread_ratio", "count_ratio")
  key[ratio] <- abs(log(value[ratio]))
  key
}

#' Largest values per check in a diagnosis
#'
#' For each module and check, list the `n` subjects (files, batches or
#' markers) with the largest value. No cutoff is applied. Signed checks
#' (`median_z`, `sd_z`, `freq_z`, `time_trend`, `order_rho`, `gain_z`,
#' `gain_order_rho`, `neg_median_z`, `neg_sd_z`, `neg_corr`, `neg_corr_z`,
#' `corr_pair_z`, `global_sd_z`, `spill_corr`, `spill_corr_z`, `kept_share_z`,
#' `kept_share_log_ratio`, `centroid_shift`, `bimodality_z`,
#' `median_shift_z`, `share_rest_z`, `batch_order_rho`) are
#' ranked by the size of the value, ignoring its sign. `gain_batch_ratio`,
#' `spread_ratio` and `count_ratio` are ranked by `abs(log(value))`, so 0.5 and
#' 2 rank equal. For `novelty_log_density`, `separation` and `lisi_batch` the
#' smallest values are the unusual ones, so those checks rank from the lowest
#' value. `rest_emd_z`, `kbet_rejection`, `iqr_effect` and `otd_*` rank by raw
#' value, highest first. Level values (`offset`, `gain`, `neg_median`, `neg_sd`, `neg_share`)
#' rank by raw value; read their `_z` rows to rank deviations.
#'
#' @param diagnosis An `sw_diagnosis` result.
#' @param n Whole number, rows kept per check.
#'
#' Rows with no value (`NA`, a thin file or an unjudged check, reason in `text`)
#' are not ranked and not dropped: they follow the ranked rows with `rank` `NA`.
#'
#' @return A tibble with the same columns as the diagnosis plus `rank`. Rows
#'   with `NA` value come last, with `rank` `NA`.
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
  # Rows with no value (not assessed, reason in `text`) are listed after the
  # ranked rows with rank NA, never dropped.
  not_assessed <- diagnosis[is.na(diagnosis$value), ]
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
  d <- tibble::as_tibble(d[setdiff(names(d), ".key")])
  if (nrow(not_assessed) > 0) {
    not_assessed <- tibble::as_tibble(not_assessed)
    not_assessed$rank <- NA_integer_
    d <- dplyr::bind_rows(d, not_assessed)
  }
  d
}
