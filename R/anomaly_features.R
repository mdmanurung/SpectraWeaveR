# R/anomaly_features.R
# Per-file features for anomaly detection: quantiles and 1D EMD to a pooled sample.
# Clean-room implementation from a written method description (CytoScan paper:
# quantile features per channel, and earth mover's distance to a pooled sample);
# no upstream code is used. The EMD is the exact 1D Wasserstein-1 distance,
# written as the integral of the gap between two empirical quantile functions.

#' Exact 1D earth mover's distance between two samples
#'
#' The Wasserstein-1 distance of two empirical distributions equals the area
#' between their quantile functions. Both functions are step functions, so the
#' area is a finite sum over the union of their jump points. There are no bins
#' and no tuning values.
#'
#' @param x,y Numeric vectors without `NA`. Sizes may differ.
#'
#' @return One non-negative number, in the units of `x` and `y`.
#'
#' @examples
#' sw_emd_1d(c(0, 1), c(2, 3)) # 2
#' x <- stats::rnorm(100)
#' sw_emd_1d(x, x + 0.5) # exactly 0.5
#'
#' @export
sw_emd_1d <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || anyNA(x) || anyNA(y)) {
    stop("`x` and `y` must be numeric without NA.", call. = FALSE)
  }
  if (length(x) == 0 || length(y) == 0) {
    stop("`x` and `y` must not be empty.", call. = FALSE)
  }
  .sw_emd_sorted(sort(x), sort(y), .sw_emd_plan(length(x), length(y)))
}

# Jump points of both quantile functions on (0, 1], as the interval widths and
# the index each sample uses on each interval. It depends only on the two
# sample sizes, so one plan serves every file of a cohort.
.sw_emd_plan <- function(n, m) {
  breaks <- sort(unique(c(0, seq_len(n) / n, seq_len(m) / m)))
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  list(width = diff(breaks), ix = ceiling(mid * n), iy = ceiling(mid * m))
}

# EMD of two sorted samples given their plan.
.sw_emd_sorted <- function(x_sorted, y_sorted, plan) {
  sum(abs(x_sorted[plan$ix] - y_sorted[plan$iy]) * plan$width)
}

#' Per-file features for anomaly detection
#'
#' @description
#' Turn a set of already-transformed files into one row of features per file.
#' Two feature types are supported, and you name which ones to compute:
#'
#' * `"quantile"`: for each channel, the quantiles at `quantile_step`,
#'   `2 * quantile_step`, and so on up to `1 - quantile_step`. The minimum and
#'   maximum are left out. A step of 0.1 gives 9 features per channel.
#' * `"emd"`: for each channel, the exact 1D earth mover's distance
#'   ([sw_emd_1d()]) between the file and a pooled sample.
#'
#' @details
#' **Scale is yours.** `data` must already be on the analysis scale, for
#' example arcsinh with cofactor 6000 for spectral flow. This function never
#' transforms. `scale` is a text note stored with the result so a reader can
#' see which scale produced the features. It is not checked against the data.
#'
#' **No defaults.** `features`, `quantile_step`, `n_events`, `seed`,
#' `aggregate_from`, `agg_size` and `scale` each change which files look
#' unusual, so each is required. CytoScan uses a quantile step of 0.1 and
#' draws 1000 events per file in its demo; those are conventions, not
#' validated cutoffs. `agg_size` and `reference` are needed only when
#' `"emd"` is requested.
#'
#' **Pooled sample.** `aggregate_from = "reference"` pools the files named in
#' `reference`. `aggregate_from = "all"` pools every file, anomalous ones
#' included, which pulls the pool towards them. Each pooled file contributes
#' `ceiling(agg_size / n_pooled_files)` events, drawn with `seed`.
#'
#' **Events.** Each file contributes `n_events` events drawn without
#' replacement with `seed`. A file with fewer events is an error, so no file
#' is silently down-weighted. The result does not depend on the order of
#' `data`, because files are processed in sorted name order.
#'
#' @param data A named list. Each element is a numeric matrix or data frame
#'   with one column per channel. Names are the file ids.
#' @param channels Character vector of channel (column) names to use. Required
#'   in the sense that there is no automatic choice.
#' @param features Character vector, any of `"quantile"` and `"emd"`.
#' @param quantile_step Number in (0, 0.5]. Needed for `"quantile"`.
#' @param n_events Whole number, events drawn per file.
#' @param seed Whole number for the draws. The caller's RNG state is restored.
#' @param scale One string describing the scale of `data`, stored in the
#'   result.
#' @param aggregate_from `"reference"` or `"all"`. Needed for `"emd"`.
#' @param agg_size Whole number, total events in the pooled sample.
#' @param reference Character vector of file ids to pool when
#'   `aggregate_from = "reference"`.
#' @param reference_group `NULL`, or a named character vector with one group
#'   label per file (names = file ids, covering every file in `data`), for
#'   example the donor. With groups, the EMD of a file is measured against the
#'   pooled events of its own group only (`aggregate_from = "reference"` pools
#'   the `reference` files of that group, and every group needs at least one),
#'   and `agg_size` is the pool size per group. `NULL` pools over all files.
#'   Quantile features do not change.
#'
#' @return A tibble of class `sw_anomaly_features` with `file` and one column
#'   per feature: `<channel>_q<probability>` and `<channel>_emd`. Attribute
#'   `params` echoes every argument except `data`.
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
#' })
#' names(files) <- paste0("f", 1:6)
#' sw_anomaly_features(
#'   files,
#'   channels = c("CD3", "CD4"),
#'   features = c("quantile", "emd"),
#'   quantile_step = 0.25, # illustrative
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   scale = "standard normal (example)",
#'   aggregate_from = "all",
#'   agg_size = 300, # illustrative
#'   reference_group = NULL
#' )
#'
#' @export
sw_anomaly_features <- function(
  data,
  channels,
  features,
  quantile_step,
  n_events,
  seed,
  scale,
  aggregate_from,
  agg_size,
  reference,
  reference_group
) {
  req <- c(
    data = missing(data),
    channels = missing(channels),
    features = missing(features),
    n_events = missing(n_events),
    seed = missing(seed),
    scale = missing(scale),
    reference_group = missing(reference_group)
  )
  .sw_stop_missing(req)

  if (
    !is.list(data) ||
      is.null(names(data)) ||
      anyDuplicated(names(data)) ||
      length(data) < 2
  ) {
    stop("`data` must be a named list of at least two files.", call. = FALSE)
  }
  if (!is.character(features) || !all(features %in% c("quantile", "emd"))) {
    stop("`features` must be \"quantile\", \"emd\" or both.", call. = FALSE)
  }
  .sw_check_count(n_events, "n_events")
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
  if (!is.character(scale) || length(scale) != 1 || !nzchar(scale)) {
    stop("`scale` must be one non-empty string.", call. = FALSE)
  }
  if (!is.character(channels) || length(channels) == 0) {
    stop("`channels` must be a non-empty character vector.", call. = FALSE)
  }
  if ("quantile" %in% features) {
    .sw_stop_missing(c(quantile_step = missing(quantile_step)))
    if (
      !is.numeric(quantile_step) ||
        length(quantile_step) != 1 ||
        is.na(quantile_step) ||
        quantile_step <= 0 ||
        quantile_step > 0.5
    ) {
      stop("`quantile_step` must be one number in (0, 0.5].", call. = FALSE)
    }
  }
  if ("emd" %in% features) {
    .sw_stop_missing(c(
      aggregate_from = missing(aggregate_from),
      agg_size = missing(agg_size)
    ))
    if (
      !is.character(aggregate_from) ||
        length(aggregate_from) != 1 ||
        !aggregate_from %in% c("reference", "all")
    ) {
      stop("`aggregate_from` must be \"reference\" or \"all\".", call. = FALSE)
    }
    .sw_check_count(agg_size, "agg_size")
    if (aggregate_from == "reference") {
      .sw_stop_missing(c(reference = missing(reference)))
      if (!is.character(reference) || !all(reference %in% names(data))) {
        stop("`reference` must name files present in `data`.", call. = FALSE)
      }
    }
  }

  # Sorted file order makes the draws independent of the list order.
  files <- sort(names(data))
  grp <- .sw_check_reference_group(reference_group, files)
  mats <- lapply(files, function(f) {
    m <- as.matrix(data[[f]])
    if (!all(channels %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' lacks channel(s): ",
        paste(setdiff(channels, colnames(m)), collapse = ", "),
        call. = FALSE
      )
    }
    if (nrow(m) < n_events) {
      stop(
        "File '",
        f,
        "' has ",
        nrow(m),
        " events, fewer than `n_events` (",
        n_events,
        ").",
        call. = FALSE
      )
    }
    m <- m[, channels, drop = FALSE]
    if (anyNA(m)) {
      stop(
        "File '",
        f,
        "' has NA values in the chosen channels.",
        call. = FALSE
      )
    }
    m
  })
  names(mats) <- files

  .sw_with_seed(seed, {
    drawn <- lapply(mats, function(m) {
      m[sample.int(nrow(m), n_events), , drop = FALSE]
    })
    pools <- NULL
    if ("emd" %in% features) {
      # One pool for all files, or one pool per reference group: a file is then
      # measured against the pooled events of its own group.
      pool_groups <- if (is.null(grp)) list(files) else split(files, grp)
      pools <- lapply(pool_groups, function(gf) {
        pooled <- if (aggregate_from == "all") gf else sort(intersect(reference, gf))
        if (length(pooled) == 0) {
          stop("A reference group has no file in `reference` to pool.", call. = FALSE)
        }
        per_file <- ceiling(agg_size / length(pooled))
        do.call(
          rbind,
          lapply(pooled, function(f) {
            m <- mats[[f]]
            m[sample.int(nrow(m), min(per_file, nrow(m))), , drop = FALSE]
          })
        )
      })
    }
  })

  cols <- list()
  if ("quantile" %in% features) {
    probs <- seq(quantile_step, 1 - quantile_step, by = quantile_step)
    probs <- probs[probs < 1 - 1e-9]
    for (ch in channels) {
      q <- t(vapply(
        drawn,
        function(m) stats::quantile(m[, ch], probs, names = FALSE, type = 7),
        numeric(length(probs))
      ))
      for (k in seq_along(probs)) {
        cols[[paste0(ch, "_q", format(probs[k], trim = TRUE))]] <- q[, k]
      }
    }
  }
  if ("emd" %in% features) {
    # The plan depends only on the sizes, and each pool is sorted once per channel.
    pool_of <- if (is.null(grp)) rep(1L, length(files)) else match(grp, names(pools))
    plans <- lapply(pools, function(p) .sw_emd_plan(n_events, nrow(p)))
    for (ch in channels) {
      pool_sorted <- lapply(pools, function(p) sort(p[, ch]))
      cols[[paste0(ch, "_emd")]] <- vapply(
        seq_along(drawn),
        function(i) {
          .sw_emd_sorted(sort(drawn[[i]][, ch]), pool_sorted[[pool_of[i]]], plans[[pool_of[i]]])
        },
        numeric(1)
      )
    }
  }

  out <- tibble::tibble(file = files)
  for (nm in names(cols)) {
    out[[nm]] <- unname(cols[[nm]])
  }
  attr(out, "params") <- list(
    channels = channels,
    features = features,
    quantile_step = if (!missing(quantile_step)) quantile_step,
    n_events = n_events,
    seed = seed,
    scale = scale,
    aggregate_from = if (!missing(aggregate_from)) aggregate_from,
    agg_size = if (!missing(agg_size)) agg_size,
    reference = if (!missing(reference)) reference,
    reference_group = reference_group
  )
  class(out) <- c("sw_anomaly_features", class(out))
  out
}

# Whole number >= 1.
.sw_check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != round(x)) {
    stop("`", name, "` must be one whole number >= 1.", call. = FALSE)
  }
}

# Run code with a seed and restore the caller's RNG state afterwards.
.sw_with_seed <- function(seed, code) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  set.seed(seed)
  force(code)
}
