# Shared-batch-effect variant of the joint REML engine (engine "heb_shared").
#
# Model per marker, on aliquot summaries of reference samples, all populations together:
#   y = donor-by-population (fixed) + study (random, shared over populations)
#       + run/batch (random, shared over populations)
#       + batch-by-population deviation (random, own variance) + error.
# The effect of (batch b, population k) is the BLUP difference to the anchor of
#   u_study(b) + v_b + w_{b,k}.
# A population whose summaries in b were dropped by min_cells still gets
# u + v of b with a zero deviation (w has no data there, so its BLUP is 0).
# Variance components: study, run, deviation, error; one set per marker.
# Written for SpectraWeaveR; see R/corrector_heb.R for the single-population engine.

#' Joint REML engine with a batch effect shared across populations
#'
#' \code{engine = "heb_shared"} in \code{\link{sw_corrector_fit}}. For each
#' marker it fits one mixed model over all populations: donor-by-population
#' effects fixed; study and run effects random and shared by every
#' population; a batch-by-population deviation random with its own variance.
#' The effect of (batch, population) is the BLUP difference to the anchor of
#' the study, run and deviation parts.
#'
#' @details
#' Needs \code{hierarchy} (columns \code{batch}, \code{study}) and
#' \code{variance_floor} > 0, both without default.
#'
#' \strong{Borrowing.} A population with no summary in a batch (every
#' aliquot below \code{min_cells}) still gets the shared study and run effect
#' of that batch, with zero deviation. Its rows carry
#' \code{support = "shared_borrowed"} and \code{n_donors = 0}.
#'
#' \strong{Abstention.} A population with no summary in the anchor batch has
#' no donor-by-population target. The engine writes no effect for it, so the
#' corrector contract abstains on it, as \code{paired_shift} and \code{heb} do.
#' Its data still inform the shared effects. A batch with no summary in any
#' population has no effect either.
#'
#' A batch unseen at fit time is handled as in \code{paired_shift}: its own
#' reference cells against the frozen anchor summaries, without shrinkage.
#' Effects of known batches never change.
#'
#' @name sw_corrector_heb_shared
#' @seealso \code{\link{sw_corrector_heb}}, \code{\link{sw_corrector_fit}}
NULL

# Negative REML log-likelihood (up to a constant); Zq is a list of tcrossprod(Z) matrices.
.sw_hebs_reml <- function(log_theta, y, X, Zq) {
  th <- exp(log_theta)
  V <- th[[length(th)]] * diag(length(y))
  for (i in seq_along(Zq)) V <- V + th[[i]] * Zq[[i]]
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) {
    return(1e10)
  }
  Vinv_X <- backsolve(R, forwardsolve(t(R), X))
  Vinv_y <- backsolve(R, forwardsolve(t(R), y))
  Rx <- tryCatch(chol(crossprod(X, Vinv_X)), error = function(e) NULL)
  if (is.null(Rx)) {
    return(1e10)
  }
  beta <- backsolve(Rx, forwardsolve(t(Rx), crossprod(X, Vinv_y)))
  res <- y - X %*% beta
  Vinv_res <- backsolve(R, forwardsolve(t(R), res))
  0.5 * (2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + sum(res * Vinv_res))
}

# Design for one marker: all populations together. NULL when no contrast to the anchor exists.
.sw_hebs_design <- function(d, hier, anchor) {
  batches <- sort(unique(d$batch))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  study_of <- hier$study[match(batches, hier$batch)]
  studies <- sort(unique(study_of))
  pops <- sort(unique(d$population))
  # Deviation columns for every (batch, population) so a population with no data in b still has one.
  bp <- expand.grid(batch = batches, population = pops, stringsAsFactors = FALSE)
  bp_key <- paste(bp$batch, bp$population, sep = "\r")
  Zs <- outer(hier$study[match(d$batch, hier$batch)], studies, `==`) * 1
  Zb <- outer(d$batch, batches, `==`) * 1
  Zw <- outer(paste(d$batch, d$population, sep = "\r"), bp_key, `==`) * 1
  list(
    d = d, batches = batches, studies = studies, study_of = study_of, pops = pops, bp = bp,
    y = d$value,
    X = stats::model.matrix(~ 0 + factor(paste(d$reference_id, d$population, sep = "\r"))),
    Z = list(Zs, Zb, Zw)
  )
}

.sw_hebs_theta <- function(z, variance_floor) {
  Zq <- lapply(z$Z, tcrossprod)
  scale <- max(stats::var(z$y), 1e-6)
  lo <- rep(log(variance_floor), 4)
  fit <- stats::optim(
    log(scale * c(0.1, 0.1, 0.1, 0.5)),
    function(lt) .sw_hebs_reml(lt, z$y, z$X, Zq),
    method = "L-BFGS-B", lower = lo, upper = lo + 40
  )
  exp(fit$par)
}

# BLUPs of the effect for every (batch, emitted population) at given variance components.
.sw_hebs_blup <- function(z, th, anchor, emit_pops) {
  y <- z$y
  X <- z$X
  Z <- do.call(cbind, z$Z)
  g_diag <- rep(th[1:3], times = vapply(z$Z, ncol, integer(1)))
  V <- Z %*% (g_diag * t(Z)) + th[[4]] * diag(length(y))
  Vinv <- chol2inv(chol(V))
  P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
  GZt <- g_diag * t(Z)
  u <- GZt %*% P %*% y
  pev <- diag(g_diag) - GZt %*% P %*% t(GZt)
  ns <- ncol(z$Z[[1]])
  nb <- ncol(z$Z[[2]])
  col_w <- function(b, k) ns + nb + match(paste(b, k, sep = "\r"), paste(z$bp$batch, z$bp$population, sep = "\r"))
  grid <- expand.grid(batch = z$batches, population = emit_pops, stringsAsFactors = FALSE)
  L <- matrix(0, nrow(grid), ncol(Z))
  for (i in seq_len(nrow(grid))) {
    b <- grid$batch[i]
    k <- grid$population[i]
    sb <- match(z$study_of[match(b, z$batches)], z$studies)
    sa <- match(z$study_of[match(anchor, z$batches)], z$studies)
    L[i, sb] <- L[i, sb] + 1
    L[i, sa] <- L[i, sa] - 1
    L[i, ns + match(b, z$batches)] <- L[i, ns + match(b, z$batches)] + 1
    L[i, ns + match(anchor, z$batches)] <- L[i, ns + match(anchor, z$batches)] - 1
    L[i, col_w(b, k)] <- L[i, col_w(b, k)] + 1
    L[i, col_w(anchor, k)] <- L[i, col_w(anchor, k)] - 1
  }
  grid$estimate <- as.vector(L %*% u)
  grid$se <- sqrt(pmax(rowSums((L %*% pev) * L), 0))
  grid$theta <- list(th)
  grid
}

.sw_hebs_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, variance_floor, ...) {
  if (!is.data.frame(hierarchy) || !all(c("batch", "study") %in% names(hierarchy))) {
    stop("heb_shared needs hierarchy = data frame with columns 'batch' and 'study'.", call. = FALSE)
  }
  if (missing(variance_floor) || !is.numeric(variance_floor) || variance_floor <= 0) {
    stop("heb_shared needs variance_floor > 0 (no default).", call. = FALSE)
  }
  s <- .sw_ps_summaries(xm, cells, min_cells)
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
  donor <- .sw_ps_donor_level(s$summ)
  anc <- donor[donor$batch == anchor, , drop = FALSE]
  if (nrow(anc) == 0) {
    stop(sprintf("Anchor batch '%s' has no reference summary.", anchor), call. = FALSE)
  }
  targets <- tibble::tibble(
    reference_id = anc$reference_id, population = anc$population,
    marker = anc$marker, target = anc$value, n_aliquots_anchor = anc$n_aliquots
  )
  # Populations that may be corrected: a donor-by-population target exists in the anchor.
  emit_pops <- sort(unique(anc$population))
  rows <- lapply(split(seq_len(nrow(s$summ)), s$summ$marker), function(ii) {
    d <- s$summ[ii, , drop = FALSE]
    z <- .sw_hebs_design(d, hierarchy, anchor)
    if (is.null(z)) {
      return(NULL)
    }
    r <- .sw_hebs_blup(z, .sw_hebs_theta(z, variance_floor), anchor, emit_pops)
    n_don <- mapply(function(b, k) length(unique(d$reference_id[d$batch == b & d$population == k])), r$batch, r$population)
    n_aliq <- mapply(function(b, k) length(unique(d$aliquot[d$batch == b & d$population == k])), r$batch, r$population)
    n_anc <- mapply(function(k) length(unique(d$aliquot[d$batch == anchor & d$population == k])), r$population)
    tibble::tibble(
      batch = r$batch, population = r$population, marker = d$marker[[1]],
      estimate = r$estimate, se_donor = r$se,
      n_donors = as.integer(n_don), n_aliquots = as.integer(n_aliq),
      n_aliquots_anchor = as.integer(n_anc),
      support = ifelse(r$batch == anchor, "anchor", ifelse(n_aliq == 0, "shared_borrowed", "joint_reml_shared")),
      source = "fit"
    )
  })
  list(effects = dplyr::bind_rows(rows), targets = targets, excluded = s$excluded)
}

sw_corrector_register_engine(
  "heb_shared",
  fit = .sw_hebs_fit,
  apply = function(model, xm, cells, new_batches) .sw_ps_apply(model, xm, cells, new_batches),
  can_freeze = TRUE,
  description = "Joint REML mixed model per marker over all populations: batch effect shared across populations plus a population-specific deviation; effect = BLUP difference to the anchor."
)
