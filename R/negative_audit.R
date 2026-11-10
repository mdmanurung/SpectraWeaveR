# R/negative_audit.R
# Checks on the negative population of each channel: centring, width, events far
# below it (hyper-negatives), and correlation of two channels inside their joint
# negatives (an angled negative). The symptoms come from O. Burton's Colibri
# Cytometry posts (identifying-and-fixing-errors-in-flow-data,
# dealing-with-multiple-autofluorescences, panel-design part 1); the code is our
# own. The negative is the lower component of a two-component Gaussian mixture
# per channel (mclust, one dimension, unequal variances: model "V").

#' Audit the negative population of each channel
#'
#' @description
#' For each file and channel, fit a two-component Gaussian mixture with
#' unequal variances to the channel's values (mclust model `"V"`, the
#' one-dimensional form of `"VVV"`). The component with the lower mean is the
#' negative population; each event goes to the component with the higher
#' posterior. Then report numbers only:
#'
#' * `neg_median`, `neg_sd`: median and robust SD (`1.4826 * MAD`) of the
#'   negative events. A skewed unmixing moves the median off its usual place;
#'   spread or autofluorescence mismatch widens the SD.
#' * `neg_median_z`, `neg_sd_z`: robust z of `neg_median` and of `log(neg_sd)`
#'   across files, per channel (`(v - median) / (1.4826 * MAD)`).
#' * `neg_share`: share of events in the negative component.
#' * `separation`: distance between the two component means in units of the
#'   root mean component variance. A small value means the channel is not
#'   bimodal in this file, so the "negative" is a split of one population.
#' * `hyper_negative_share`: share of all events below
#'   `neg_median - hyper_k * neg_sd`. Autofluorescence over-correction and
#'   mismatched controls push events there.
#' * `neg_corr`, `neg_corr_z`: Spearman correlation of two channels inside the
#'   events that are negative in both, and its robust z across files, for each
#'   pair in `pairs`. A tilted ("angled") negative shows here while a
#'   whole-file correlation can hide it.
#' * `spill_corr`, `spill_corr_z`: for each pair, in both directions, the
#'   Spearman correlation of the two channels inside the events that are
#'   positive in the first and negative in the second (`"A > B"`), and its
#'   robust z across files. Signal leaking from a bright channel into one that
#'   should be negative, in proportion to the bright signal, shows here: tandem
#'   breakdown into its parent (for example PE-Cy7 into PE) and dye-dye
#'   interactions (for example Spark UV 387 into RB670), as described in O.
#'   Burton's Colibri Cytometry posts. Pass those pairs in `pairs`.
#' * `spread_log_ratio`, `spread_z`: spillover spread, for each pair in both
#'   directions (`"A > B"`), when `spread_quantile` is given. Among the events
#'   negative in the target B, the robust SD of B where the source A is at or
#'   above its `spread_quantile` quantile, over the robust SD of B where A is
#'   at or below its median, as a natural log; and its robust z across files.
#'   Spread widens the target's negative in proportion to the source signal, so
#'   a file (or instrument) with more spread from A into B has a higher ratio.
#'   Fewer than 10 events on either side gives `NA`.
#'
#' @details
#' **Scale is yours.** `data` must be on the analysis scale. The mixture is fit
#' on those values; a transform changes which events are negative.
#'
#' **No defaults.** `n_events`, `seed`, `hyper_k`, `pairs` and
#' `spread_quantile` change the
#' numbers, so each must be passed. Events are drawn per file without
#' replacement with `seed`; a file with fewer events is an error. The caller's
#' RNG state is restored. mclust itself is deterministic.
#'
#' **What can go wrong.** A unimodal channel still gets two components; read
#' `separation` before trusting the other numbers. A fit that fails, a
#' constant channel, or a channel with missing or infinite values gives `NA`
#' with the reason in `detail`.
#'
#' @param data A named list of numeric matrices, one per file, one column per
#'   channel.
#' @param channels Character vector of channels to audit.
#' @param n_events Whole number of events drawn per file.
#' @param seed Whole number for the draws.
#' @param hyper_k Positive number: how many negative robust SDs below the
#'   negative median an event must sit to count as hyper-negative.
#' @param pairs `NULL` (no pair checks), the string `"all"` (every pair of
#'   `channels`), or a two-column character matrix of channel pairs. Each
#'   pair gives `neg_corr` once and `spill_corr` in both directions.
#' @param spread_quantile `NULL` (no spread check), or one number in (0.5, 1):
#'   the source-channel quantile above which events count as high for
#'   `spread_log_ratio`. Needs `pairs`.
#'
#' @return A tibble of class `sw_negative_audit` with columns `file`, `check`,
#'   `channel` (a channel, `"A ~ B"` for a pair, or `"A > B"` for a directed
#'   pair), `value` and `detail`.
#'
#' @seealso [sw_fcs_signal_audit()], [sw_diagnose()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   neg <- stats::runif(1, 0.5, 0.7)
#'   x <- sapply(c("A", "B"), function(ch) {
#'     ifelse(stats::runif(2000) < neg, stats::rnorm(2000, 0, 0.5), stats::rnorm(2000, 3, 0.8))
#'   })
#'   x
#' })
#' names(files) <- paste0("f", 1:6)
#' sw_negative_audit(files, c("A", "B"), n_events = 1000, seed = 1, hyper_k = 4, pairs = "all",
#'   spread_quantile = 0.9)
#' }
#'
#' @export
sw_negative_audit <- function(data, channels, n_events, seed, hyper_k, pairs, spread_quantile) {
  .sw_stop_missing(c(
    data = missing(data),
    channels = missing(channels),
    n_events = missing(n_events),
    seed = missing(seed),
    hyper_k = missing(hyper_k),
    pairs = missing(pairs),
    spread_quantile = missing(spread_quantile)
  ))
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for sw_negative_audit().", call. = FALSE)
  }
  if (
    !is.null(spread_quantile) &&
      (!is.numeric(spread_quantile) ||
        length(spread_quantile) != 1 ||
        !is.finite(spread_quantile) ||
        spread_quantile <= 0.5 ||
        spread_quantile >= 1)
  ) {
    stop("`spread_quantile` must be NULL or one number in (0.5, 1).", call. = FALSE)
  }
  if (!is.null(spread_quantile) && is.null(pairs)) {
    stop("`spread_quantile` needs `pairs`.", call. = FALSE)
  }
  if (!is.list(data) || is.null(names(data)) || length(data) < 2) {
    stop("`data` must be a named list of at least two matrices.", call. = FALSE)
  }
  .sw_check_count(n_events, "n_events")
  .sw_check_count(seed, "seed")
  if (!is.numeric(hyper_k) || length(hyper_k) != 1 || !is.finite(hyper_k) || hyper_k <= 0) {
    stop("`hyper_k` must be one positive number.", call. = FALSE)
  }
  pair_mat <- .sw_negative_pairs(pairs, channels)
  for (f in names(data)) {
    absent <- setdiff(channels, colnames(data[[f]]))
    if (length(absent) > 0) {
      stop("File ", f, " lacks channels: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    if (nrow(data[[f]]) < n_events) {
      stop("File ", f, " has fewer than `n_events` events.", call. = FALSE)
    }
  }

  draws <- .sw_with_seed(seed, {
    lapply(data, \(m) as.matrix(m[sample.int(nrow(m), n_events), channels, drop = FALSE]))
  })

  per_file <- lapply(names(draws), function(f) {
    x <- draws[[f]]
    fits <- lapply(channels, \(ch) .sw_negative_fit(x[, ch]))
    names(fits) <- channels
    rows <- lapply(channels, function(ch) {
      ft <- fits[[ch]]
      v <- x[, ch]
      hyper <- if (is.na(ft$neg_median)) {
        NA_real_
      } else {
        mean(v < ft$neg_median - hyper_k * ft$neg_sd)
      }
      tibble::tibble(
        file = f,
        check = c("neg_median", "neg_sd", "neg_share", "separation", "hyper_negative_share"),
        channel = ch,
        value = c(ft$neg_median, ft$neg_sd, ft$neg_share, ft$separation, hyper),
        detail = ft$detail
      )
    })
    if (!is.null(pair_mat)) {
      pair_rows <- lapply(seq_len(nrow(pair_mat)), function(i) {
        a <- pair_mat[i, 1]
        b <- pair_mat[i, 2]
        both <- fits[[a]]$negative & fits[[b]]$negative
        value <- if (anyNA(both) || sum(both) < 10) {
          NA_real_
        } else {
          suppressWarnings(stats::cor(x[both, a], x[both, b], method = "spearman"))
        }
        # Directed leak: positive in the source, negative in the target.
        spill <- lapply(list(c(a, b), c(b, a)), function(st) {
          src <- fits[[st[1]]]$negative
          tgt <- fits[[st[2]]]$negative
          keep <- !src & tgt
          v <- if (anyNA(keep) || sum(keep) < 10) {
            NA_real_
          } else {
            suppressWarnings(stats::cor(x[keep, st[1]], x[keep, st[2]], method = "spearman"))
          }
          row <- tibble::tibble(
            file = f,
            check = "spill_corr",
            channel = paste(st[1], st[2], sep = " > "),
            value = v,
            detail = paste0("events=", if (anyNA(keep)) NA else sum(keep))
          )
          if (is.null(spread_quantile)) {
            return(row)
          }
          # Spread: target-negative width where the source is high vs low.
          s <- x[, st[1]]
          hi <- tgt & s >= stats::quantile(s, spread_quantile)
          lo <- tgt & s <= stats::median(s)
          ratio <- if (anyNA(tgt) || sum(hi) < 10 || sum(lo) < 10) {
            NA_real_
          } else {
            sd_hi <- stats::mad(x[hi, st[2]])
            sd_lo <- stats::mad(x[lo, st[2]])
            if (sd_hi > 0 && sd_lo > 0) log(sd_hi / sd_lo) else NA_real_
          }
          dplyr::bind_rows(
            row,
            tibble::tibble(
              file = f,
              check = "spread_log_ratio",
              channel = paste(st[1], st[2], sep = " > "),
              value = ratio,
              detail = paste0("high=", if (anyNA(tgt)) NA else sum(hi), ", low=", if (anyNA(tgt)) NA else sum(lo))
            )
          )
        })
        dplyr::bind_rows(
          tibble::tibble(
            file = f,
            check = "neg_corr",
            channel = paste(a, b, sep = " ~ "),
            value = value,
            detail = paste0("joint negatives=", if (anyNA(both)) NA else sum(both))
          ),
          spill
        )
      })
      rows <- c(rows, pair_rows)
    }
    dplyr::bind_rows(rows)
  })
  out <- dplyr::bind_rows(per_file)

  # Cohort robust z per channel (or pair), across files.
  robust_z <- function(v) {
    s <- 1.4826 * stats::mad(v, constant = 1, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(NA_real_, length(v)) else (v - stats::median(v, na.rm = TRUE)) / s
  }
  z_rows <- lapply(
    list(
      c("neg_median", "neg_median_z"),
      c("neg_sd", "neg_sd_z"),
      c("neg_corr", "neg_corr_z"),
      c("spill_corr", "spill_corr_z"),
      c("spread_log_ratio", "spread_z")
    ),
    function(ck) {
      d <- out[out$check == ck[[1]], ]
      if (nrow(d) == 0) {
        return(NULL)
      }
      # A negative SD of 0 (ties) has no log; it becomes NA, not -Inf.
      v <- if (ck[[1]] == "neg_sd") log(ifelse(d$value > 0, d$value, NA_real_)) else d$value
      d$value <- stats::ave(v, d$channel, FUN = robust_z)
      d$check <- ck[[2]]
      d$detail <- NA_character_
      d
    }
  )
  out <- dplyr::bind_rows(c(list(out), z_rows))
  class(out) <- c("sw_negative_audit", class(out))
  out
}

# Two-component one-dimensional Gaussian mixture, unequal variances.
.sw_negative_fit <- function(v) {
  empty <- list(
    neg_median = NA_real_,
    neg_sd = NA_real_,
    neg_share = NA_real_,
    separation = NA_real_,
    negative = rep(NA, length(v)),
    detail = NA_character_
  )
  if (!all(is.finite(v))) {
    empty$detail <- "non-finite values (see nonfinite_share in sw_fcs_signal_audit())"
    return(empty)
  }
  if (stats::sd(v) == 0) {
    empty$detail <- "constant channel"
    return(empty)
  }
  # Mclust() evaluates mclustBIC() in the caller's frame, so bind it here.
  mclustBIC <- mclust::mclustBIC
  fit <- tryCatch(
    mclust::Mclust(v, G = 2, modelNames = "V", verbose = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error") || is.null(fit)) {
    empty$detail <- paste(
      "mixture fit failed:",
      if (inherits(fit, "error")) conditionMessage(fit) else "no model"
    )
    return(empty)
  }
  means <- fit$parameters$mean
  vars <- fit$parameters$variance$sigmasq
  neg_comp <- which.min(means)
  negative <- fit$classification == neg_comp
  if (sum(negative) < 2) {
    empty$detail <- "negative component has fewer than 2 events"
    return(empty)
  }
  list(
    neg_median = stats::median(v[negative]),
    neg_sd = 1.4826 * stats::mad(v[negative], constant = 1),
    neg_share = mean(negative),
    separation = abs(diff(means)) / sqrt(mean(vars)),
    negative = negative,
    detail = sprintf("means %.3g, %.3g", min(means), max(means))
  )
}

# Resolve `pairs` to a two-column character matrix, or NULL.
.sw_negative_pairs <- function(pairs, channels) {
  if (is.null(pairs)) {
    return(NULL)
  }
  if (identical(pairs, "all")) {
    if (length(channels) < 2) {
      stop("`pairs = \"all\"` needs at least two channels.", call. = FALSE)
    }
    return(t(utils::combn(channels, 2)))
  }
  if (!is.matrix(pairs) || ncol(pairs) != 2 || !is.character(pairs)) {
    stop("`pairs` must be NULL, \"all\", or a two-column character matrix.", call. = FALSE)
  }
  absent <- setdiff(as.vector(pairs), channels)
  if (length(absent) > 0) {
    stop("`pairs` names channels not in `channels`: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  pairs
}
