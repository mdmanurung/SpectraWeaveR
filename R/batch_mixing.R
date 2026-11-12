# R/batch_mixing.R
# Batch mixing in a k-nearest-neighbour graph of events: LISI and the kBET
# rejection rate. Written from the published definitions only:
# Korsunsky et al. 2019, Nature Methods 16:1289 (LISI) and Buttner et al. 2019,
# Nature Methods 16:43 (kBET). No code from the lisi, harmony or kBET packages
# was opened or used.

#' Batch mixing by LISI and kBET
#'
#' @description
#' Draws `n_events` events per file, pools them, and asks how well the batches
#' mix among nearest neighbours in the space of `markers`. Two numbers per
#' unit (all events, or each `population`):
#'
#' * `lisi_batch`, the local inverse Simpson index of batch labels (Korsunsky
#'   et al. 2019). Around each event, the `k` nearest neighbours get Gaussian
#'   kernel weights whose width is set so the weights have the given
#'   `perplexity`. The index is `1 / sum_b p_b^2`, where `p_b` is the weight
#'   share of batch `b`. It is the effective number of batches in the
#'   neighbourhood: near the number of batches when batches mix evenly, near 1
#'   when they separate.
#' * `kbet_rejection`, the share of events whose `k` nearest neighbours reject
#'   the null of a batch mix equal to the global one (Buttner et al. 2019):
#'   a chi-square test of the batch counts among the neighbours against the
#'   expected counts `k * share_b`, at level `alpha`. Near `alpha` when
#'   batches mix; near 1 when they separate.
#'
#' @details
#' **Neighbours** are the `k` nearest events by Euclidean distance, not
#' counting the event itself. Every event is tested; neighbourhoods overlap, so
#' the tests are not independent and the rejection rate under perfect mixing is
#' near `alpha` only on average. The chi-square approximation is poor when
#' `k * share_b` is below about 5 for a batch.
#'
#' **Rows.** Each unit gives an overall row per check (`batch` is `NA`) and one
#' row per batch (mean LISI, or rejection rate, over that batch's events). The
#' overall mean depends on batch sizes; the per-batch rows do not hide a small
#' batch. Nothing is combined over markers, and there is no pass or fail.
#'
#' **Not assessed.** A unit with fewer than `min_events` drawn events, fewer
#' than `k + 1` events, or fewer than two batches gets `status =
#' "not_assessed"` and `NA` values.
#'
#' **Units.** Files are the replicates. Events are drawn per file, so a batch
#' with more files contributes more events. Cells are not independent
#' replicates; read the rates as descriptions of this draw.
#'
#' **No defaults.** `k`, `perplexity`, `alpha`, `n_events`, `seed` and
#' `min_events` change the numbers and are uncalibrated. `perplexity` must be
#' below `k`. The caller's RNG state is restored.
#'
#' @param data A named list of numeric matrices, one per file, columns named.
#' @param batch A named vector of batch labels, one per file in `data`.
#' @param markers Character vector of columns that define the space.
#' @param k Whole number, neighbours per event (uncalibrated).
#' @param perplexity Number above 1 and below `k`, the effective neighbour
#'   count of the LISI kernel (uncalibrated).
#' @param alpha Number in (0, 1), the kBET test level (uncalibrated).
#' @param n_events Whole number, events drawn per file (uncalibrated).
#' @param seed Whole number for the draws (uncalibrated).
#' @param min_events Whole number, fewest drawn events a unit needs
#'   (uncalibrated).
#' @param population `NULL` (one unit with all events), or a named list with
#'   one character vector per file, giving each event's population (length
#'   `nrow()` of that file). Mixing is then measured inside each population.
#'   Drawn events with an `NA` label are left out with a warning that gives
#'   their number; all-`NA` labels are an error.
#'
#' @return A tibble of class `sw_batch_mixing` with columns `check`
#'   (`lisi_batch` or `kbet_rejection`), `population` (`NA` without
#'   `population`), `batch` (`NA` for the overall row), `value`, `n_events`
#'   (events in the unit, or for a batch row, in that batch), `n_batches` and
#'   `status`.
#'
#' @seealso [sw_batch_scores()], [sw_batch_otd()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:4, function(i) {
#'   matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:4)
#' sw_batch_mixing(
#'   files,
#'   batch = stats::setNames(c("b1", "b1", "b2", "b2"), names(files)),
#'   markers = c("A", "B"),
#'   k = 30, # illustrative
#'   perplexity = 10, # illustrative
#'   alpha = 0.05, # illustrative
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   min_events = 50, # illustrative
#'   population = NULL
#' )
#'
#' @export
sw_batch_mixing <- function(
  data,
  batch,
  markers,
  k,
  perplexity,
  alpha,
  n_events,
  seed,
  min_events,
  population
) {
  .sw_stop_missing(c(
    data = missing(data),
    batch = missing(batch),
    markers = missing(markers),
    k = missing(k),
    perplexity = missing(perplexity),
    alpha = missing(alpha),
    n_events = missing(n_events),
    seed = missing(seed),
    min_events = missing(min_events),
    population = missing(population)
  ))
  .sw_check_count(k, "k")
  .sw_check_count(n_events, "n_events")
  .sw_check_count(seed, "seed")
  .sw_check_count(min_events, "min_events")
  if (
    !is.numeric(perplexity) || length(perplexity) != 1 || !is.finite(perplexity)
  ) {
    stop("`perplexity` must be one number.", call. = FALSE)
  }
  if (perplexity <= 1 || perplexity >= k) {
    stop("`perplexity` must be above 1 and below `k`.", call. = FALSE)
  }
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1 ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop("`alpha` must be one number between 0 and 1.", call. = FALSE)
  }
  files <- .sw_check_batch_files(data, batch, markers)
  batch <- as.character(batch[files])
  names(batch) <- files
  if (length(unique(batch)) < 2) {
    stop("`batch` needs at least two batches.", call. = FALSE)
  }
  if (!is.null(population)) {
    ok <- is.list(population) &&
      all(files %in% names(population)) &&
      all(vapply(
        files,
        \(f) length(population[[f]]) == nrow(data[[f]]),
        logical(1)
      ))
    if (!ok) {
      stop(
        "`population` must be NULL or a named list with one label per event for every file.",
        call. = FALSE
      )
    }
  }
  if (any(vapply(files, \(f) nrow(data[[f]]) < n_events, logical(1)))) {
    stop("A file has fewer events than `n_events`.", call. = FALSE)
  }

  # Draw n_events per file, in sorted file order.
  drawn <- .sw_with_seed(seed, {
    lapply(files, \(f) {
      idx <- sample.int(nrow(data[[f]]), n_events)
      list(
        x = as.matrix(data[[f]][idx, markers, drop = FALSE]),
        pop = if (is.null(population)) {
          rep(NA_character_, n_events)
        } else {
          as.character(population[[f]][idx])
        }
      )
    })
  })
  x_all <- do.call(rbind, lapply(drawn, `[[`, "x"))
  batch_all <- rep(batch, each = n_events)
  pop_all <- unlist(lapply(drawn, `[[`, "pop"), use.names = FALSE)

  if (!is.null(population) && anyNA(pop_all)) {
    if (all(is.na(pop_all))) {
      stop("Every drawn event has an NA `population` label.", call. = FALSE)
    }
    warning(
      sum(is.na(pop_all)),
      " of ",
      length(pop_all),
      " drawn events have an NA `population` label and are left out.",
      call. = FALSE
    )
  }
  units <- if (is.null(population)) {
    NA_character_
  } else {
    sort(unique(pop_all[!is.na(pop_all)]))
  }
  rows <- lapply(units, \(u) {
    keep <- if (is.na(u)) rep(TRUE, length(batch_all)) else which(pop_all == u)
    .sw_mixing_unit(
      x_all[keep, , drop = FALSE],
      batch_all[keep],
      u,
      k = k,
      perplexity = perplexity,
      alpha = alpha,
      min_events = min_events
    )
  })
  out <- dplyr::bind_rows(rows)
  class(out) <- c("sw_batch_mixing", class(out))
  out
}

# Rows for one unit: overall and per-batch LISI and kBET.
.sw_mixing_unit <- function(x, lab, unit, k, perplexity, alpha, min_events) {
  batches <- sort(unique(lab))
  n <- nrow(x)
  checks <- c("lisi_batch", "kbet_rejection")
  thin <- n < max(min_events, k + 1) || length(batches) < 2
  if (thin) {
    return(tibble::tibble(
      check = checks,
      population = unit,
      batch = NA_character_,
      value = NA_real_,
      n_events = n,
      n_batches = length(batches),
      status = "not_assessed"
    ))
  }
  nn <- .sw_knn_index(x, k)
  code <- match(lab, batches)
  nb_code <- matrix(code[nn], nrow = n)
  d2 <- .sw_knn_sqdist(x, nn)

  # LISI: inverse Simpson over the kernel weights of the neighbours.
  p <- .sw_perplexity_weights(d2, perplexity)
  share <- vapply(
    seq_along(batches),
    \(b) rowSums(p * (nb_code == b)),
    numeric(n)
  )
  lisi <- 1 / rowSums(share^2)

  # kBET: chi-square of neighbour batch counts against k * global share.
  global <- tabulate(code, nbins = length(batches)) / n
  counts <- vapply(
    seq_along(batches),
    \(b) rowSums(nb_code == b),
    numeric(n)
  )
  expected <- matrix(k * global, nrow = n, ncol = length(batches), byrow = TRUE)
  stat <- rowSums((counts - expected)^2 / expected)
  reject <- stats::pchisq(stat, df = length(batches) - 1, lower.tail = FALSE) <
    alpha

  per_event <- list(lisi_batch = lisi, kbet_rejection = as.numeric(reject))
  dplyr::bind_rows(lapply(checks, \(chk) {
    v <- per_event[[chk]]
    tibble::tibble(
      check = chk,
      population = unit,
      batch = c(NA_character_, batches),
      value = c(mean(v), vapply(batches, \(b) mean(v[lab == b]), numeric(1))),
      n_events = c(n, as.integer(table(lab)[batches])),
      n_batches = length(batches),
      status = "assessed"
    )
  }))
}

# Indices of the k nearest rows for every row, self excluded. Chunked so the
# distance matrix never exceeds chunk * n.
.sw_knn_index <- function(x, k, chunk = 500) {
  n <- nrow(x)
  sq <- rowSums(x^2)
  out <- matrix(0L, nrow = n, ncol = k)
  for (s in seq(1, n, by = chunk)) {
    r <- s:min(n, s + chunk - 1)
    d <- outer(sq[r], sq, "+") - 2 * tcrossprod(x[r, , drop = FALSE], x)
    for (i in seq_along(r)) {
      o <- order(d[i, ], method = "radix")[seq_len(k + 1)]
      o <- o[o != r[i]]
      out[r[i], ] <- o[seq_len(k)]
    }
  }
  out
}

# Squared Euclidean distance from each row to its neighbours (n x k).
.sw_knn_sqdist <- function(x, nn) {
  d2 <- vapply(
    seq_len(ncol(nn)),
    \(j) rowSums((x - x[nn[, j], , drop = FALSE])^2),
    numeric(nrow(x))
  )
  matrix(d2, nrow = nrow(x))
}

# Gaussian kernel weights over the neighbours of each row, width found by
# bisection so that exp(entropy) = perplexity. All rows are solved together.
# Squared distances are first divided by the row mean, so the search range
# does not depend on the data scale; this does not change the weights.
.sw_perplexity_weights <- function(d2, perplexity) {
  d2 <- d2 - apply(d2, 1, min)
  scale <- rowMeans(d2)
  scale[scale <= 0] <- 1
  s <- d2 / scale
  lo <- rep(-20, nrow(s))
  hi <- rep(20, nrow(s))
  weights <- function(log_beta) {
    e <- exp(-exp(log_beta) * s)
    e / rowSums(e)
  }
  for (it in 1:60) {
    mid <- (lo + hi) / 2
    p <- weights(mid)
    ent <- -rowSums(p * log(pmax(p, 1e-300)))
    too_flat <- ent > log(perplexity) # entropy falls as beta rises
    lo <- ifelse(too_flat, mid, lo)
    hi <- ifelse(too_flat, hi, mid)
  }
  weights((lo + hi) / 2)
}
