# R/fcs_signal_audit.R
# Per-file technical checks on event data: saturation, dead channels, NA/Inf,
# negative share, event-rate stability and channel-pair correlation change.
# Written from scratch; it only computes summaries and sets no cutoffs.

#' Audit event data of each file for technical problems
#'
#' @description
#' Compute a set of plain per-file numbers that expose common technical faults.
#' The result is a long table, one row per file, check and channel. It holds
#' numbers only: no check is marked pass or fail, because a cutoff changes
#' which files look bad and is yours to choose.
#'
#' Checks (the `check` column):
#'
#' * `n_events`: events in the file (file level).
#' * `nonfinite_share`: share of `NA`, `NaN` or `Inf` values per channel.
#' * `zero_variance`: 1 when a channel is constant (a dead or disconnected
#'   detector), else 0, per channel.
#' * `negative_share`: share of values below 0 per channel.
#' * `saturation_share`: share of values at or above `range_max` per channel.
#'   Only when `range_max` is given.
#' * `floor_share`: share of values at or below `range_min` per channel. Only
#'   when `range_min` is given.
#' * `rate_cv`: coefficient of variation of event counts per time bin (file
#'   level). A clog or a pressure change makes the event rate uneven. Only when
#'   `time_channel` is given.
#' * `max_gap`: largest gap between consecutive event times, as a multiple of
#'   the median gap (file level). Only when `time_channel` is given.
#' * `median_z`, `sd_z`: robust z-score of the file's median and of its SD in
#'   a channel against all files, `(value - median) / (1.4826 * MAD)` over
#'   files. A shift of one marker shows in `median_z`; a stretch or a split
#'   into two modes widens the SD and shows in `sd_z`. The reference is the
#'   cohort you pass, batch effects included, so a whole batch that differs
#'   shifts the reference too. `NA` with fewer than 3 files.
#' * `global_sd_z`: per file, the median over channels of
#'   `log(file SD / cohort median SD)`, as a robust z over files (file level).
#'   A fault that widens or narrows every channel a little (a gain or focus
#'   change) shows here even when no single channel's `sd_z` stands out. It is
#'   a check of its own, not a summary of the other rows. `NA` with fewer than
#'   3 files.
#' * `time_trend`: Spearman correlation of a channel with time, per channel. A
#'   steady drift of a marker inside one acquisition (laser, fluidics or
#'   detector drift) gives a value far from 0. Only when `time_channel` is
#'   given.
#' * `corr_shift`: largest absolute difference between the file's Spearman
#'   correlation of a channel pair and the median of that pair over all files
#'   (file level, `detail` names the pair). A wrong spillover or unmixing
#'   matrix changes channel-pair correlations.
#' * `corr_pair_z`: the same correlations, one row per file and channel pair
#'   (`channel` is `"A ~ B"`): robust z of the file's correlation against all
#'   files. Rank these rows per pair; `corr_shift` keeps only the largest change
#'   per file, which hides a change in one pair among many (on 36 real files
#'   with 40 channels, a 0.3 spillover scored AUROC 0.66 with `corr_shift`).
#'   `NA` with fewer than 3 files.
#'
#' @details
#' **Scale is yours.** Saturation and floor limits must be on the same scale as
#' `data`. `negative_share` is only meaningful on a scale that keeps negative
#' values (raw or arcsinh).
#'
#' **No defaults.** `range_max`, `range_min`, `time_channel` and `time_bin`
#' each change what is reported, so each must be passed. Pass `NULL` for a
#' check you do not want. `corr_shift` uses every row of each file; subsample
#' first for large files.
#'
#' @param data A named list. Each element is a numeric matrix or data frame
#'   with one column per channel. Names are the file ids. At least two files.
#' @param channels Character vector of channel names to audit (not the time
#'   channel).
#' @param range_max `NULL`, or a single number or named numeric vector of the
#'   instrument maximum per channel.
#' @param range_min `NULL`, or a single number or named numeric vector of the
#'   lowest value counted as a floor event.
#' @param time_channel `NULL`, or the name of the time column in each file.
#' @param time_bin Number, width of a time bin in the units of the time
#'   channel. Needed when `time_channel` is given; otherwise pass `NULL`.
#'
#' @return A tibble of class `sw_fcs_signal_audit` with columns `file`,
#'   `check`, `channel` (`NA` for file-level checks), `value` and `detail`
#'   (text; the channel pair for `corr_shift`, else `NA`).
#'
#' @seealso [sw_anomaly_features()], [sw_fcs_header_audit()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:4, function(i) {
#'   matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:4)
#' sw_fcs_signal_audit(
#'   files,
#'   channels = c("A", "B"),
#'   range_max = NULL, range_min = NULL, time_channel = NULL, time_bin = NULL
#' )
#'
#' @export
sw_fcs_signal_audit <- function(
  data,
  channels,
  range_max,
  range_min,
  time_channel,
  time_bin
) {
  .sw_stop_missing(c(
    data = missing(data),
    channels = missing(channels),
    range_max = missing(range_max),
    range_min = missing(range_min),
    time_channel = missing(time_channel),
    time_bin = missing(time_bin)
  ))
  if (
    !is.list(data) ||
      is.null(names(data)) ||
      anyDuplicated(names(data)) ||
      length(data) < 2
  ) {
    stop("`data` must be a named list of at least two files.", call. = FALSE)
  }
  if (!is.character(channels) || length(channels) < 2) {
    stop("`channels` must name at least two channels.", call. = FALSE)
  }
  lim_max <- .sw_audit_limit(range_max, channels, "range_max")
  lim_min <- .sw_audit_limit(range_min, channels, "range_min")
  if (!is.null(time_channel)) {
    if (!is.character(time_channel) || length(time_channel) != 1) {
      stop("`time_channel` must be one column name or NULL.", call. = FALSE)
    }
    if (
      !is.numeric(time_bin) ||
        length(time_bin) != 1 ||
        is.na(time_bin) ||
        time_bin <= 0
    ) {
      stop(
        "`time_bin` must be a positive number when `time_channel` is given.",
        call. = FALSE
      )
    }
  }

  files <- sort(names(data))
  mats <- lapply(files, function(f) {
    m <- as.matrix(data[[f]])
    need <- c(channels, time_channel)
    if (!all(need %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' lacks column(s): ",
        paste(setdiff(need, colnames(m)), collapse = ", "),
        call. = FALSE
      )
    }
    m
  })
  names(mats) <- files

  rows <- list()
  add <- function(file, check, channel, value, detail = NA_character_) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      file = file,
      check = check,
      channel = channel,
      value = as.numeric(value), # as.numeric drops names carried by some inputs
      detail = detail
    )
  }
  # Per file and channel median and SD, for the cohort z-scores below.
  med_mat <- sd_mat <- matrix(
    NA_real_,
    length(files),
    length(channels),
    dimnames = list(files, channels)
  )
  for (f in files) {
    m <- mats[[f]]
    add(f, "n_events", NA_character_, nrow(m))
    for (ch in channels) {
      v <- m[, ch]
      finite <- is.finite(v)
      add(f, "nonfinite_share", ch, mean(!finite))
      vf <- v[finite]
      add(
        f,
        "zero_variance",
        ch,
        as.numeric(length(vf) > 0 && stats::var(vf) == 0)
      )
      add(f, "negative_share", ch, if (length(vf)) mean(vf < 0) else NA_real_)
      if (length(vf) > 1) {
        med_mat[f, ch] <- stats::median(vf)
        sd_mat[f, ch] <- stats::sd(vf)
      }
      if (!is.null(lim_max)) {
        add(
          f,
          "saturation_share",
          ch,
          if (length(vf)) mean(vf >= lim_max[[ch]]) else NA_real_
        )
      }
      if (!is.null(lim_min)) {
        add(
          f,
          "floor_share",
          ch,
          if (length(vf)) mean(vf <= lim_min[[ch]]) else NA_real_
        )
      }
    }
    if (!is.null(time_channel)) {
      tt <- sort(m[, time_channel])
      tt <- tt[is.finite(tt)]
      # Events per time bin; bins with no events inside the run count as zero.
      bin_id <- as.integer(floor((tt - min(tt)) / time_bin)) + 1L
      full <- tabulate(bin_id, nbins = max(bin_id))
      add(
        f,
        "rate_cv",
        NA_character_,
        if (length(full) > 1) stats::sd(full) / mean(full) else NA_real_
      )
      gaps <- diff(tt)
      med <- stats::median(gaps)
      add(
        f,
        "max_gap",
        NA_character_,
        if (length(gaps) && med > 0) max(gaps) / med else NA_real_
      )
      # Drift inside the file: rank correlation of each channel with time.
      for (ch in channels) {
        add(
          f,
          "time_trend",
          ch,
          as.numeric(suppressWarnings(stats::cor(
            m[, ch],
            m[, time_channel],
            method = "spearman",
            use = "pairwise.complete.obs"
          )))
        )
      }
    }
  }

  # Robust z of each file's median and SD against the other files, per channel:
  # (value - median over files) / (1.4826 * MAD over files). NA with fewer than
  # 3 files or when the spread over files is 0.
  robust_z <- function(x) {
    # stats::mad() scales by 1.4826 by default; constant = 1 keeps the raw MAD.
    s <- 1.4826 * stats::mad(x, constant = 1, na.rm = TRUE)
    if (sum(!is.na(x)) < 3 || is.na(s) || s == 0) {
      return(rep(NA_real_, length(x)))
    }
    (x - stats::median(x, na.rm = TRUE)) / s
  }
  for (ch in channels) {
    zm <- robust_z(med_mat[, ch])
    zs <- robust_z(sd_mat[, ch])
    for (j in seq_along(files)) {
      add(files[j], "median_z", ch, zm[j])
      add(files[j], "sd_z", ch, zs[j])
    }
  }
  # File-wide spread change (user decision 2026-11-10): per file, the median
  # over channels of log(file SD / cohort median SD), as a robust z over files.
  # A fault that widens every channel a little shows here, not in any one row.
  cohort_sd <- apply(sd_mat, 2, stats::median, na.rm = TRUE)
  log_ratio <- log(sweep(sd_mat, 2, cohort_sd, "/"))
  log_ratio[!is.finite(log_ratio)] <- NA_real_
  global <- apply(log_ratio, 1, stats::median, na.rm = TRUE)
  global[is.nan(global)] <- NA_real_
  zg <- robust_z(global)
  for (j in seq_along(files)) {
    add(files[j], "global_sd_z", NA_character_, zg[j])
  }

  # Channel-pair Spearman correlation against the cohort median per pair.
  cors <- lapply(mats, function(m) {
    suppressWarnings(stats::cor(
      m[, channels],
      method = "spearman",
      use = "pairwise.complete.obs"
    ))
  })
  idx <- which(upper.tri(cors[[1]]), arr.ind = TRUE)
  pair_name <- paste(channels[idx[, 1]], channels[idx[, 2]], sep = " ~ ")
  pair_vals <- vapply(cors, function(cm) cm[idx], numeric(nrow(idx)))
  pair_vals <- matrix(pair_vals, nrow = nrow(idx))
  center <- apply(pair_vals, 1, stats::median, na.rm = TRUE)
  for (j in seq_along(files)) {
    d <- abs(pair_vals[, j] - center)
    if (all(is.na(d))) {
      add(files[j], "corr_shift", NA_character_, NA_real_)
    } else {
      k <- which.max(d)
      add(files[j], "corr_shift", NA_character_, d[k], pair_name[k])
    }
  }

  # The same pairs one row each: robust z of the file's correlation against the
  # other files, per pair. Built in one step (pairs x files can be large).
  pair_z <- t(apply(pair_vals, 1, robust_z))
  if (nrow(idx) == 1) {
    pair_z <- matrix(pair_z, nrow = 1)
  }
  pair_rows <- tibble::tibble(
    file = rep(files, each = nrow(idx)),
    check = "corr_pair_z",
    channel = rep(pair_name, times = length(files)),
    value = as.numeric(pair_z),
    detail = NA_character_
  )

  out <- dplyr::bind_rows(c(rows, list(pair_rows)))
  class(out) <- c("sw_fcs_signal_audit", class(out))
  out
}

# NULL, one number, or a named vector covering every channel -> named vector.
.sw_audit_limit <- function(x, channels, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || anyNA(x)) {
    stop("`", name, "` must be NULL or numeric without NA.", call. = FALSE)
  }
  if (length(x) == 1 && is.null(names(x))) {
    return(stats::setNames(rep(x, length(channels)), channels))
  }
  if (!all(channels %in% names(x))) {
    stop("`", name, "` must name every channel.", call. = FALSE)
  }
  x[channels]
}
