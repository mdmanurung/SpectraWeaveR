#' @title Negative Shape and Marker Screen
#'
#' @description
#' Compare the negative events of several unmixing setups against a reference
#' setup, channel by channel, on paired cells (plan v2 steps 11 and 12,
#' user grill 2026-11-11).
#'
#' * [sw_unmix_negative_shape()] scores, per setup, sample, population and
#'   channel, the spread of the negatives (Nguyen rSD, 84th minus 50th
#'   percentile of events below the cut), its ratio to the reference, the share
#'   of events squeezed towards zero and the share above the cut.
#' * [sw_unmix_marker_screen()] flags a setup, population and channel when the
#'   squeezed share rises in every donor and the donor-median rise reaches
#'   `min_delta`.
#'
#' @details
#' Negatives are each setup's own events below the cut (grill R9). The rSD is
#' computed only when at least `min_negatives` events lie below the cut, in the
#' setup and in the reference (grill R10); otherwise the metrics that need it
#' are `NA`. There are no built-in thresholds: `cut`, `squeeze_frac`,
#' `min_negatives`, `min_delta` and `min_donors` are required (ROADMAP rule 2).
#'
#' The method ports the benchmark scripts `helper75_marker-scores.R` and
#' `86b_apply_flag_rule_v2.R` (experiments/benchmarks/05_bmv-unmix-compare).
#' One difference from 86b: too few donors gives `NA`, never `FALSE`.
#'
#' @name unmix_screen
#' @keywords internal
NULL

# ---- Internal helpers --------------------------------------------------------

# Nguyen rSD: 84th minus 50th percentile.
.sw_nguyen_rsd <- function(v) {
  unname(diff(stats::quantile(v, c(0.5, 0.84), names = FALSE)))
}

# One matrix for a setup and sample: a matrix as given, or read from an FCS path.
.sw_screen_matrix <- function(x, setup, sample) {
  if (is.matrix(x)) {
    return(x)
  }
  if (is.character(x) && length(x) == 1L) {
    if (!requireNamespace("flowCore", quietly = TRUE)) {
      stop("Reading FCS paths needs the flowCore package.", call. = FALSE)
    }
    return(flowCore::exprs(flowCore::read.FCS(
      x,
      transformation = FALSE,
      truncate_max_range = FALSE
    )))
  }
  stop(
    sprintf(
      "x[['%s']][['%s']] must be a matrix or one FCS path.",
      setup,
      sample
    ),
    call. = FALSE
  )
}

# Metrics for one channel on one population's rows; ref_rsd is the reference's rSD.
.sw_screen_metrics <- function(
  lin,
  cofactor,
  cut,
  ref_rsd,
  squeeze_frac,
  min_negatives
) {
  below <- asinh(lin / cofactor) < cut
  n_neg <- sum(below)
  rsd_neg <- if (n_neg >= min_negatives) {
    .sw_nguyen_rsd(lin[below])
  } else {
    NA_real_
  }
  # A reference spread of 0 or NA cannot scale the squeeze window or the ratio.
  ref_ok <- !is.na(ref_rsd) && ref_rsd > 0
  tibble::tibble(
    n_pop = length(lin),
    n_neg = n_neg,
    rsd_neg = rsd_neg,
    ref_rsd = ref_rsd,
    rsd_ratio = if (ref_ok) rsd_neg / ref_rsd else NA_real_,
    squeezed_share = if (ref_ok) {
      mean(abs(lin) < squeeze_frac * ref_rsd)
    } else {
      NA_real_
    },
    pos_share = mean(asinh(lin / cofactor) > cut)
  )
}

# ---- Negative shape -----------------------------------------------------------

#' Score the negatives of each setup against a reference setup
#'
#' @param x Named list by setup, each a named list by sample of unmixed
#'   matrices (events x channels, linear scale) or FCS paths. Rows are paired:
#'   row `i` is the same cell in every setup of a sample, checked on
#'   `key_cols`.
#' @param reference Name of the reference setup in `x`.
#' @param populations Named list by sample, each a named list by population
#'   of integer row indices into the paired rows.
#' @param channels Character vector of channels to score.
#' @param cofactor arcsinh cofactor.
#' @param cut Cut on the arcsinh scale: events below it are negatives.
#' @param squeeze_frac An event counts as squeezed when its absolute linear
#'   value is below `squeeze_frac` times the reference rSD.
#' @param min_negatives Fewest below-cut events needed to compute an rSD, in
#'   the setup and in the reference.
#' @param key_cols Columns that identify a cell and must match across setups.
#'
#' @return An `sw_negative_shape`: a list with `metrics` (one row per
#'   `setup`, `sample`, `pop`, `channel`: `n_pop`, `n_neg`, `rsd_neg`,
#'   `ref_rsd`, `rsd_ratio`, `squeezed_share`, `pos_share`; the reference's own
#'   rows included), `reference`, `cut`, `cofactor`, `squeeze_frac`,
#'   `min_negatives` and `provenance`.
#'
#' @examples
#' set.seed(1)
#' mk <- function(sd) cbind(Time = 1:200, TNF = rnorm(200, sd = sd))
#' x <- list(ref = list(d1 = mk(300)), alt = list(d1 = mk(150)))
#' x$alt$d1[, "Time"] <- x$ref$d1[, "Time"]
#' shape <- sw_unmix_negative_shape(
#'   x,
#'   reference = "ref",
#'   populations = list(d1 = list(T = 1:200)),
#'   channels = "TNF",
#'   cofactor = 6000,
#'   cut = 1,
#'   squeeze_frac = 0.25,
#'   min_negatives = 50,
#'   key_cols = "Time"
#' )
#' shape$metrics
#'
#' @export
sw_unmix_negative_shape <- function(
  x,
  reference,
  populations,
  channels,
  cofactor,
  cut,
  squeeze_frac,
  min_negatives,
  key_cols
) {
  .sw_stop_missing(c(
    x = missing(x),
    reference = missing(reference),
    populations = missing(populations),
    channels = missing(channels),
    cofactor = missing(cofactor),
    cut = missing(cut),
    squeeze_frac = missing(squeeze_frac),
    min_negatives = missing(min_negatives),
    key_cols = missing(key_cols)
  ))
  if (!is.list(x) || is.null(names(x)) || anyDuplicated(names(x))) {
    stop("'x' must be a list with unique setup names.", call. = FALSE)
  }
  if (
    !is.character(reference) ||
      length(reference) != 1L ||
      !reference %in% names(x)
  ) {
    stop("'reference' must name one setup in 'x'.", call. = FALSE)
  }
  for (arg in c("cofactor", "cut", "squeeze_frac", "min_negatives")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be one number.", arg), call. = FALSE)
    }
  }
  samples <- names(populations)
  if (is.null(samples) || anyDuplicated(samples)) {
    stop("'populations' must be a list named by sample.", call. = FALSE)
  }
  missing_samples <- unlist(lapply(names(x), \(s) {
    absent <- setdiff(samples, names(x[[s]]))
    if (length(absent)) paste0(s, "/", absent)
  }))
  if (length(missing_samples) > 0L) {
    stop(
      sprintf(
        "'x' lacks sample(s): %s.",
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  metrics <- list()
  for (s in samples) {
    mats <- lapply(names(x), \(setup) {
      .sw_screen_matrix(x[[setup]][[s]], setup, s)
    })
    names(mats) <- names(x)
    ref <- mats[[reference]]
    absent_cols <- setdiff(c(key_cols, channels), colnames(ref))
    if (length(absent_cols) > 0L) {
      stop(
        sprintf(
          "Sample '%s', setup '%s': column(s) not found: %s.",
          s,
          reference,
          paste(absent_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    # Rows must pair across setups: same row count and identical key columns.
    for (setup in names(mats)) {
      m <- mats[[setup]]
      paired <- nrow(m) == nrow(ref) &&
        all(c(key_cols, channels) %in% colnames(m)) &&
        isTRUE(all.equal(
          m[, key_cols, drop = FALSE],
          ref[, key_cols, drop = FALSE],
          check.attributes = FALSE
        ))
      if (!paired) {
        stop(
          sprintf(
            "Sample '%s': setup '%s' rows do not pair with the reference on key_cols.",
            s,
            setup
          ),
          call. = FALSE
        )
      }
    }
    pops <- populations[[s]]
    if (is.null(names(pops)) || anyDuplicated(names(pops))) {
      stop(
        sprintf("populations[['%s']] must be named by population.", s),
        call. = FALSE
      )
    }
    for (p in names(pops)) {
      rows <- pops[[p]]
      if (
        !is.numeric(rows) ||
          any(rows < 1 | rows > nrow(ref) | rows != round(rows))
      ) {
        stop(
          sprintf(
            "populations[['%s']][['%s']] must be row indices of the sample.",
            s,
            p
          ),
          call. = FALSE
        )
      }
      for (ch in channels) {
        # Reference spread on its own below-cut events (the squeeze window and ratio scale).
        ref_lin <- ref[rows, ch]
        ref_below <- asinh(ref_lin / cofactor) < cut
        ref_rsd <- if (sum(ref_below) >= min_negatives) {
          .sw_nguyen_rsd(ref_lin[ref_below])
        } else {
          NA_real_
        }
        for (setup in names(mats)) {
          metrics[[length(metrics) + 1L]] <- .sw_screen_metrics(
            mats[[setup]][rows, ch],
            cofactor,
            cut,
            ref_rsd,
            squeeze_frac,
            min_negatives
          ) |>
            dplyr::mutate(
              setup = .env$setup,
              sample = .env$s,
              pop = .env$p,
              channel = .env$ch,
              .before = 1
            )
        }
      }
    }
  }

  structure(
    list(
      metrics = dplyr::bind_rows(metrics),
      reference = reference,
      cut = cut,
      cofactor = cofactor,
      squeeze_frac = squeeze_frac,
      min_negatives = min_negatives,
      provenance = .sw_spike_provenance(reference = reference)
    ),
    class = "sw_negative_shape"
  )
}

#' @export
print.sw_negative_shape <- function(x, ...) {
  m <- x$metrics
  cat(
    "<sw_negative_shape> reference:",
    x$reference,
    "\n  setups:",
    paste(unique(m$setup), collapse = ", "),
    "| samples:",
    length(unique(m$sample)),
    "| populations:",
    paste(unique(m$pop), collapse = ", "),
    "| channels:",
    length(unique(m$channel)),
    "\n  cut:",
    x$cut,
    "arcsinh | squeeze_frac:",
    x$squeeze_frac,
    "| min_negatives:",
    x$min_negatives,
    "| rows with NA rSD:",
    sum(is.na(m$rsd_neg)),
    "\n"
  )
  invisible(x)
}

# ---- Marker screen ------------------------------------------------------------

#' Flag setups that squeeze a channel's negatives in every donor
#'
#' Per setup (reference excluded), population and channel: the change in
#' squeezed share versus the reference in each donor, its donor median, and
#' the flag. The flag is `TRUE` when the change is above 0 in every donor and
#' the donor median is at least `min_delta`; `FALSE` otherwise; `NA` when
#' fewer than `min_donors` donors have a value.
#'
#' @param shape An `sw_negative_shape` from [sw_unmix_negative_shape()].
#' @param min_delta Smallest donor-median rise in squeezed share that flags.
#' @param min_donors Fewest donors with a value needed for a flag.
#' @param channel_populations Optional named list, channel = populations in
#'   which that channel is screened (for example a marker only within its
#'   parent). Channels not named are screened in every population.
#'
#' @return An `sw_marker_screen`: a list with `flags` (one row per `setup`,
#'   `pop`, `channel`: `n_donors`, `median_d_squeezed`, `median_rsd_ratio`,
#'   `all_donors_up`, `flag`), `per_donor` (the per-donor changes), `rule`
#'   and `reference`.
#'
#' @examples
#' set.seed(1)
#' mk <- function(sd) cbind(Time = 1:200, TNF = rnorm(200, sd = sd))
#' x <- list(ref = list(d1 = mk(300), d2 = mk(300)), alt = list(d1 = mk(60), d2 = mk(60)))
#' for (d in c("d1", "d2")) x$alt[[d]][, "Time"] <- x$ref[[d]][, "Time"]
#' shape <- sw_unmix_negative_shape(
#'   x, "ref", list(d1 = list(T = 1:200), d2 = list(T = 1:200)), "TNF",
#'   cofactor = 6000, cut = 1, squeeze_frac = 0.25, min_negatives = 50, key_cols = "Time"
#' )
#' sw_unmix_marker_screen(shape, min_delta = 0.1, min_donors = 2)$flags
#'
#' @export
sw_unmix_marker_screen <- function(
  shape,
  min_delta,
  min_donors,
  channel_populations = NULL
) {
  .sw_stop_missing(c(
    shape = missing(shape),
    min_delta = missing(min_delta),
    min_donors = missing(min_donors)
  ))
  if (!inherits(shape, "sw_negative_shape")) {
    stop("'shape' must be an sw_negative_shape.", call. = FALSE)
  }
  for (arg in c("min_delta", "min_donors")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be one number.", arg), call. = FALSE)
    }
  }
  m <- shape$metrics
  if (!is.null(channel_populations)) {
    if (!is.list(channel_populations) || is.null(names(channel_populations))) {
      stop(
        "'channel_populations' must be a list named by channel.",
        call. = FALSE
      )
    }
    keep <- vapply(
      seq_len(nrow(m)),
      \(i) {
        ch <- m$channel[i]
        !ch %in% names(channel_populations) ||
          m$pop[i] %in% channel_populations[[ch]]
      },
      logical(1)
    )
    m <- m[keep, ]
  }

  # Each setup's donor rows against the same donor's reference row.
  ref <- m[
    m$setup == shape$reference,
    c("sample", "pop", "channel", "squeezed_share")
  ]
  names(ref)[4] <- "ref_squeezed"
  per_donor <- m[m$setup != shape$reference, ] |>
    dplyr::left_join(
      ref,
      by = c("sample", "pop", "channel"),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    dplyr::mutate(d_squeezed = .data$squeezed_share - .data$ref_squeezed)

  flags <- per_donor |>
    dplyr::summarise(
      n_donors = sum(!is.na(.data$d_squeezed)),
      median_d_squeezed = if (any(!is.na(.data$d_squeezed))) {
        stats::median(.data$d_squeezed, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_rsd_ratio = if (any(!is.na(.data$rsd_ratio))) {
        stats::median(.data$rsd_ratio, na.rm = TRUE)
      } else {
        NA_real_
      },
      all_donors_up = all(.data$d_squeezed[!is.na(.data$d_squeezed)] > 0),
      .by = c("setup", "pop", "channel")
    ) |>
    dplyr::mutate(
      flag = dplyr::if_else(
        .data$n_donors < min_donors,
        NA,
        .data$all_donors_up & .data$median_d_squeezed >= min_delta
      )
    )

  structure(
    list(
      flags = flags,
      per_donor = per_donor,
      rule = list(
        min_delta = min_delta,
        min_donors = min_donors,
        text = sprintf(
          "squeezed share rises (> 0) in every donor and the donor median rise is >= %s; NA under %s donors",
          format(min_delta),
          format(min_donors)
        )
      ),
      reference = shape$reference
    ),
    class = "sw_marker_screen"
  )
}

#' @export
print.sw_marker_screen <- function(x, ...) {
  f <- x$flags
  cat(
    "<sw_marker_screen> reference:",
    x$reference,
    "\n  rule:",
    x$rule$text,
    "\n  flagged:",
    sum(f$flag, na.rm = TRUE),
    "of",
    nrow(f),
    "| NA:",
    sum(is.na(f$flag)),
    "\n"
  )
  invisible(x)
}
