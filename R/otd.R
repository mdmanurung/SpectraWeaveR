# R/otd.R
# Debiased Sinkhorn divergence (entropic optimal transport) between event sets,
# and a per-file score against a pooled reference. Written from the standard
# definitions (Cuturi 2013; Feydy et al. 2019) as used by Glehr et al. 2026,
# EBioMedicine 131:106427, for their optimal transport distance (OTD). No code
# from geomloss or otcyto was used.

#' Debiased Sinkhorn divergence between two sets of events
#'
#' @description
#' The entropic optimal transport cost with squared Euclidean cost
#' `C(x, y) = |x - y|^2 / 2`, with the self-transport terms removed:
#'
#' `S(a, b) = OT(a, b) - OT(a, a) / 2 - OT(b, b) / 2`
#'
#' where `OT` is the regularised cost at `epsilon = blur^2`. `S` is 0 when the
#' two sets are equal and grows with how far mass has to move. For a pure
#' translation by a vector `d`, `S = |d|^2 / 2` at any `blur`. Glehr et al.
#' (2026) used this score (via geomloss) with `blur = 0.75` on arcsinh data.
#'
#' @details
#' Solved in the log domain. Epsilon is first annealed from the largest cost
#' down to `blur^2` by a factor 0.7 per step (as in geomloss), then Sinkhorn
#' updates run at `blur^2` until `tol` or `max_iter`. Annealing changes speed,
#' not the converged value. A `blur` that is small against the spacing of the
#' events converges slowly; check the `converged` attribute and raise
#' `max_iter` when it is `FALSE`. Memory grows with `nrow(x) * nrow(y)`, so
#' subsample large files first. Every event has equal weight.
#'
#' **No defaults.** `blur` is in the units of the data and sets how much fine
#' detail counts; `max_iter` and `tol` set when the solver stops. All three
#' change the number.
#'
#' @param x,y Numeric matrices with the same columns, one row per event.
#' @param blur Positive number, the blur scale (`epsilon = blur^2`).
#' @param max_iter Whole number, most Sinkhorn iterations per transport problem.
#' @param tol Positive number. Stop when no dual potential moves by more than
#'   `tol` in one iteration.
#'
#' @return One number, the divergence, with attribute `converged` (`TRUE` when
#'   all three problems met `tol` within `max_iter`).
#'
#' @seealso [sw_otd_scores()], [sw_emd_1d()]
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(400), ncol = 2)
#' sw_sinkhorn_divergence(x, x + 1, blur = 0.5, max_iter = 500, tol = 1e-6) # |d|^2 / 2 = 1
#'
#' @export
sw_sinkhorn_divergence <- function(x, y, blur, max_iter, tol) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (!is.numeric(x) || !is.numeric(y) || ncol(x) != ncol(y)) {
    stop("`x` and `y` must be numeric matrices with the same number of columns.", call. = FALSE)
  }
  if (anyNA(x) || anyNA(y) || nrow(x) < 1 || nrow(y) < 1) {
    stop("`x` and `y` need at least one row and no missing values.", call. = FALSE)
  }
  if (!is.numeric(blur) || length(blur) != 1 || !is.finite(blur) || blur <= 0) {
    stop("`blur` must be one positive number.", call. = FALSE)
  }
  .sw_check_count(max_iter, "max_iter")
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one positive number.", call. = FALSE)
  }
  eps <- blur^2
  xy <- .sw_ot_entropic(x, y, eps, max_iter, tol)
  xx <- .sw_ot_entropic(x, x, eps, max_iter, tol)
  yy <- .sw_ot_entropic(y, y, eps, max_iter, tol)
  out <- xy$value - xx$value / 2 - yy$value / 2
  attr(out, "converged") <- xy$converged && xx$converged && yy$converged
  out
}

#' Optimal transport distance of each file to a pooled reference
#'
#' @description
#' For each file, draw `n_events` events, and compute the debiased Sinkhorn
#' divergence ([sw_sinkhorn_divergence()]) to a reference built by pooling
#' events from all files and drawing `reference_events` of them. This is the
#' per-file score of Glehr et al. (2026), who pooled the training samples and
#' subsampled to 10,000 cells. A file far from the pool has a large value.
#'
#' @details
#' **Scale is yours.** `data` must already be on the analysis scale (for
#' example arcsinh); `blur` is in those units. Scores from two differently
#' transformed versions of the same data are not comparable, because each has
#' its own reference space. Compare a paired difference instead, as the paper
#' did.
#'
#' **No defaults.** `n_events`, `reference_events`, `blur`, `seed`, `max_iter`
#' and `tol` change the numbers. The reference draws the same number of events
#' from each file (`reference_events / length(data)`, rounded down), so every
#' file weighs the same. The caller's RNG state is restored.
#'
#' @param data A named list of numeric matrices, one per file.
#' @param markers Character vector of columns to use jointly.
#' @param n_events Whole number of events drawn per file.
#' @param reference_events Whole number of events in the pooled reference.
#' @param seed Whole number for the draws.
#' @inheritParams sw_sinkhorn_divergence
#'
#' @return A tibble with columns `file`, `otd` and `converged`.
#'
#' @seealso [sw_sinkhorn_divergence()], [sw_batch_scores()]
#'
#' @export
sw_otd_scores <- function(data, markers, n_events, reference_events, blur, seed, max_iter, tol) {
  if (!is.list(data) || is.null(names(data)) || length(data) < 2) {
    stop("`data` must be a named list of at least two matrices.", call. = FALSE)
  }
  .sw_check_count(n_events, "n_events")
  .sw_check_count(reference_events, "reference_events")
  .sw_check_count(seed, "seed")
  per_file_ref <- floor(reference_events / length(data))
  if (per_file_ref < 1) {
    stop("`reference_events` must be at least the number of files.", call. = FALSE)
  }
  for (f in names(data)) {
    absent <- setdiff(markers, colnames(data[[f]]))
    if (length(absent) > 0) {
      stop("File ", f, " lacks markers: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    if (nrow(data[[f]]) < max(n_events, per_file_ref)) {
      stop(
        "File ",
        f,
        " has ",
        nrow(data[[f]]),
        " events; fewer than `n_events` or the reference share.",
        call. = FALSE
      )
    }
  }

  .sw_with_seed(seed, {
    reference <- do.call(
      rbind,
      lapply(data, \(m) as.matrix(m[sample.int(nrow(m), per_file_ref), markers, drop = FALSE]))
    )
    draws <- lapply(data, \(m) as.matrix(m[sample.int(nrow(m), n_events), markers, drop = FALSE]))
  })
  # The reference self-transport term is the same for every file: solve it once.
  eps <- blur^2
  ref_self <- .sw_ot_entropic(reference, reference, eps, max_iter, tol)
  scores <- lapply(names(draws), function(f) {
    xy <- .sw_ot_entropic(draws[[f]], reference, eps, max_iter, tol)
    xx <- .sw_ot_entropic(draws[[f]], draws[[f]], eps, max_iter, tol)
    tibble::tibble(
      file = f,
      otd = xy$value - xx$value / 2 - ref_self$value / 2,
      converged = xy$converged && xx$converged && ref_self$converged
    )
  })
  dplyr::bind_rows(scores)
}

# Entropic OT between uniform weights on the rows of x and y, log domain.
# Returns the dual value <a, f> + <b, g> at convergence. When x and y are the
# same object the symmetric (averaged) update is used, which converges fast.
.sw_ot_entropic <- function(x, y, eps, max_iter, tol) {
  n <- nrow(x)
  m <- nrow(y)
  cost <- (outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)) / 2
  cost[cost < 0] <- 0 # rounding
  log_a <- rep(-log(n), n)
  log_b <- rep(-log(m), m)
  # Soft-min over columns (for f) or rows (for g) of (potential - cost) / eps.
  softmin_rows <- function(g) {
    z <- sweep(-cost, 2, g, "+") / eps
    z <- sweep(z, 2, log_b, "+")
    mx <- apply(z, 1, max)
    -eps * (mx + log(rowSums(exp(z - mx))))
  }
  softmin_cols <- function(f) {
    z <- (-cost + f) / eps + log_a
    mx <- apply(z, 2, max)
    -eps * (mx + log(colSums(exp(sweep(z, 2, mx)))))
  }
  symmetric <- identical(x, y)
  f <- rep(0, n)
  g <- rep(0, m)
  # Epsilon annealing (as in geomloss, factor 0.7): start at the largest cost
  # and shrink to the target, one update per level, warm-starting the
  # potentials. It speeds convergence; the fixed point at the target epsilon
  # is unchanged.
  target_eps <- eps
  schedule <- max(cost) * 0.7^(0:200)
  schedule <- schedule[schedule > target_eps]
  for (e in schedule) {
    eps <- e
    if (symmetric) {
      f <- (f + softmin_rows(f)) / 2
    } else {
      f <- softmin_rows(g)
      g <- softmin_cols(f)
    }
  }
  eps <- target_eps
  if (symmetric) {
    g <- f
  }
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    if (symmetric) {
      f_new <- (f + softmin_rows(f)) / 2
      g_new <- f_new
    } else {
      f_new <- softmin_rows(g)
      g_new <- softmin_cols(f_new)
    }
    change <- max(abs(f_new - f), abs(g_new - g))
    f <- f_new
    g <- g_new
    if (change < tol) {
      converged <- TRUE
      break
    }
  }
  if (symmetric) {
    # Final full update for the self-transport value.
    f <- softmin_rows(f)
    g <- f
  }
  list(value = mean(f) + mean(g), converged = converged)
}
