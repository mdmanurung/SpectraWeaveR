# R/cluster_scores.R
# Cluster-level file checks: cluster frequencies and per-cluster marker medians,
# scored against the other files. Targets changes in a rare population or in the
# composition of a file, which file-wide summaries average away. Written from
# scratch with stats::kmeans.

#' Score each file's clusters against the other files
#'
#' @description
#' Cluster the pooled events of all files with k-means, then describe every file
#' by the share of its events in each cluster and by the median of each marker
#' inside each cluster. Each value is turned into a robust z-score against the
#' other files, so a file whose rare population shifted, or whose composition
#' changed, stands out even when its whole-file medians barely move.
#'
#' Two checks (the `check` column):
#'
#' * `freq_z`: robust z of the share of a file's events in a cluster.
#' * `median_z`: robust z of the median of a marker inside a cluster, using
#'   only files with at least `min_cluster_events` events in that cluster.
#'
#' A robust z is `(value - median) / (1.4826 * MAD)` over files.
#'
#' @details
#' **No scaling.** The markers are clustered on the scale of `data`. This
#' function never centres or scales them, because that would be a model-input
#' transform that you approve. Markers on very different scales then weigh very
#' differently in the clustering; put them on a common scale first if that is
#' what you want.
#'
#' **No defaults.** `k`, `nstart`, `n_events`, `seed` and `min_cluster_events`
#' each change the clusters or the values, so each is required. K-means runs
#' with at most 100 iterations. Events are drawn per file without replacement
#' with `seed`, and the caller's RNG state is restored. A file with fewer events
#' than `n_events` is an error.
#'
#' **k-means warnings.** On large pools `stats::kmeans` may warn that its
#' quick-transfer stage reached its step limit. The clusters are still returned.
#' The warning is re-issued once with context, and attribute `kmeans_warnings`
#' keeps its text, so you can see how often it happens.
#'
#' **Reference is the cohort.** The z-scores compare a file with the other
#' files you pass. A fault shared by most files, or a real batch effect, shifts
#' the reference too.
#'
#' **Where the clusters come from.** With `reference = NULL`, k-means runs on
#' the pooled draws of all files, faulty ones included. With `reference` set to
#' file ids (for example control or known-good files), k-means runs on those
#' files' draws only and every file's events are assigned to the nearest
#' centre, so the clusters do not depend on the files being judged. On 36 real
#' files with a +1 SD shift inside the smallest of 6 clusters, this did not
#' raise the AUROC (0.86 with clean reference files, 0.88 without); scoring the
#' injected cluster's own row with all events gave 0.96.
#'
#' @param data A named list. Each element is a numeric matrix or data frame
#'   with one column per marker. Names are the file ids. At least 3 files.
#' @param markers Character vector of marker (column) names to cluster on.
#' @param k Whole number, number of clusters.
#' @param n_events Whole number, events drawn per file.
#' @param seed Whole number for the draws and the k-means start.
#' @param nstart Whole number, random starts for k-means.
#' @param min_cluster_events Whole number, fewest events a file needs in a
#'   cluster for its median to be scored.
#' @param reference `NULL` to fit k-means on all files, or a character vector
#'   of file ids to fit it on; all files are then assigned to the nearest
#'   centre.
#'
#' @return A tibble of class `sw_cluster_scores` with columns `file`, `cluster`
#'   (integer), `marker` (`NA` for `freq_z`), `check`, `value`, `n_events`
#'   (events of the file in the cluster) and `share` (that count over the
#'   file's drawn events). Attribute `centers` holds the cluster centres.
#'
#' @seealso [sw_fcs_signal_audit()], [sw_batch_scores()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   matrix(c(stats::rnorm(300, -2), stats::rnorm(300, 2)), ncol = 2,
#'     dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:6)
#' sw_cluster_scores(
#'   files, c("A", "B"),
#'   k = 2, n_events = 200, seed = 1, nstart = 5, min_cluster_events = 20,
#'   reference = NULL
#' )
#'
#' @export
sw_cluster_scores <- function(
  data,
  markers,
  k,
  n_events,
  seed,
  nstart,
  min_cluster_events,
  reference
) {
  .sw_stop_missing(c(
    data = missing(data),
    markers = missing(markers),
    k = missing(k),
    n_events = missing(n_events),
    seed = missing(seed),
    nstart = missing(nstart),
    min_cluster_events = missing(min_cluster_events),
    reference = missing(reference)
  ))
  if (
    !is.list(data) ||
      is.null(names(data)) ||
      anyDuplicated(names(data)) ||
      length(data) < 3
  ) {
    stop("`data` must be a named list of at least three files.", call. = FALSE)
  }
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  for (nm in c("k", "n_events", "nstart", "min_cluster_events")) {
    .sw_check_count(get(nm), nm)
  }
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
  if (k < 2) {
    stop("`k` must be at least 2.", call. = FALSE)
  }
  if (
    !is.null(reference) &&
      (!is.character(reference) ||
        length(reference) == 0 ||
        anyDuplicated(reference) ||
        !all(reference %in% names(data)))
  ) {
    stop(
      "`reference` must be NULL or distinct file ids from `names(data)`.",
      call. = FALSE
    )
  }

  files <- sort(names(data))
  drawn_in <- lapply(files, function(f) {
    m <- as.matrix(data[[f]])
    if (!all(markers %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' lacks marker(s): ",
        paste(setdiff(markers, colnames(m)), collapse = ", "),
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
    m[, markers, drop = FALSE]
  })
  names(drawn_in) <- files

  res <- .sw_with_seed(seed, {
    drawn <- lapply(drawn_in, function(m) {
      m[sample.int(nrow(m), n_events), , drop = FALSE]
    })
    pool <- do.call(rbind, drawn)
    if (anyNA(pool)) {
      stop("Events must not contain NA in the chosen markers.", call. = FALSE)
    }
    # k-means is fitted on the reference files' draws only, when given.
    fit_pool <- if (is.null(reference)) pool else do.call(rbind, drawn[sort(reference)])
    if (nrow(unique(fit_pool)) < k) {
      stop("`k` is larger than the number of distinct reference events.", call. = FALSE)
    }
    if (nrow(unique(pool)) < k) {
      stop("`k` is larger than the number of distinct events.", call. = FALSE)
    }
    # Keep kmeans' own warnings (for example "Quick-TRANSfer stage steps
    # exceeded maximum" on large pools) visible, but once and with context.
    km_warn <- character(0)
    km <- withCallingHandlers(
      stats::kmeans(fit_pool, centers = k, nstart = nstart, iter.max = 100),
      warning = function(w) {
        km_warn <<- c(km_warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(km = km, warnings = unique(km_warn))
  })
  km <- res$km
  if (length(res$warnings) > 0) {
    warning(
      "k-means reported: ",
      paste(res$warnings, collapse = "; "),
      ". The clusters are still used; a larger `nstart` or a smaller pool ",
      "(`n_events`) can help.",
      call. = FALSE
    )
  }
  # Cluster of every drawn event: k-means' own labels, or the nearest
  # reference centre when the clusters were fitted on reference files.
  labels <- if (is.null(reference)) {
    km$cluster
  } else {
    pool <- do.call(rbind, drawn)
    d2 <- outer(rowSums(pool^2), rowSums(km$centers^2), "+") -
      2 * tcrossprod(pool, km$centers)
    max.col(-d2, ties.method = "first")
  }
  cl <- split(labels, rep(files, each = n_events))[files]

  robust_z <- function(x) {
    # stats::mad() scales by 1.4826 by default; constant = 1 keeps the raw MAD.
    s <- 1.4826 * stats::mad(x, constant = 1, na.rm = TRUE)
    if (sum(!is.na(x)) < 3 || is.na(s) || s == 0) {
      return(rep(NA_real_, length(x)))
    }
    (x - stats::median(x, na.rm = TRUE)) / s
  }

  # Shares and counts: files x clusters.
  counts <- t(vapply(cl, function(v) tabulate(v, nbins = k), numeric(k)))
  shares <- counts / n_events
  rows <- list()
  for (j in seq_len(k)) {
    z <- robust_z(shares[, j])
    rows[[length(rows) + 1L]] <- tibble::tibble(
      file = files,
      cluster = j,
      marker = NA_character_,
      check = "freq_z",
      value = as.numeric(z),
      n_events = as.integer(counts[, j]),
      share = as.numeric(shares[, j])
    )
    for (mk in markers) {
      med <- vapply(
        seq_along(files),
        function(i) {
          v <- drawn[[i]][cl[[i]] == j, mk]
          if (length(v) >= min_cluster_events) stats::median(v) else NA_real_
        },
        numeric(1)
      )
      rows[[length(rows) + 1L]] <- tibble::tibble(
        file = files,
        cluster = j,
        marker = mk,
        check = "median_z",
        value = as.numeric(robust_z(med)),
        n_events = as.integer(counts[, j]),
        share = as.numeric(shares[, j])
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  attr(out, "centers") <- km$centers
  attr(out, "kmeans_warnings") <- res$warnings
  class(out) <- c("sw_cluster_scores", class(out))
  out
}
