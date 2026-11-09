# Joint REML hierarchical engine for the corrector contract (HEB-v0).
#
# Model per (population, marker), on aliquot summaries of reference samples:
#   y = donor (fixed) + study (random) + run/batch inside study (random) + error.
# The correction of batch b is the BLUP difference beta_b - beta_anchor, where
# beta_b = study effect + run effect. See the derivation in the tusoskill run
# (derivations/heb_v0_derivation.md). Written for SpectraWeaveR; it does not
# use code from the cyTether bundle's Python prototype.

#' Joint REML hierarchical engine
#'
#' \code{engine = "heb"} in \code{\link{sw_corrector_fit}}. For each
#' (population, marker) it fits a linear mixed model to aliquot summaries of
#' the reference donors: donor effects fixed, study and run effects random,
#' variance components by REML. The effect of batch \eqn{b} is the difference
#' of the best linear unbiased predictions of the batch effects of \eqn{b} and
#' of the anchor. \code{corrected = x - effect}.
#'
#' @details
#' \code{hierarchy} is a data frame with columns \code{batch} and
#' \code{study}, one row per batch. It has no default. Variance components
#' are bounded below by \code{variance_floor} (also required, no default).
#' \code{se_donor} holds the posterior standard error of the effect, from the
#' prediction error covariance. A batch at integrate time that the model has
#' not seen is estimated as in \code{paired_shift}: its own reference cells
#' against the frozen anchor summaries, without shrinkage. Effects of known
#' batches are never changed.
#'
#' @name sw_corrector_heb
#' @seealso \code{\link{sw_corrector_fit}}, \code{\link{sw_corrector_paired_shift}}
NULL

# Negative REML log-likelihood (up to a constant) for log variance components.
.sw_heb_reml <- function(log_theta, y, X, Zs, Zb) {
  th <- exp(log_theta)
  V <- th[[1]] * tcrossprod(Zs) + th[[2]] * tcrossprod(Zb) + th[[3]] * diag(length(y))
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) {
    return(1e10)
  }
  Vinv_X <- backsolve(R, forwardsolve(t(R), X))
  Vinv_y <- backsolve(R, forwardsolve(t(R), y))
  XtVX <- crossprod(X, Vinv_X)
  Rx <- tryCatch(chol(XtVX), error = function(e) NULL)
  if (is.null(Rx)) {
    return(1e10)
  }
  beta <- backsolve(Rx, forwardsolve(t(Rx), crossprod(X, Vinv_y)))
  res <- y - X %*% beta
  Vinv_res <- backsolve(R, forwardsolve(t(R), res))
  0.5 * (2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + sum(res * Vinv_res))
}

# One (population, marker): effect of every batch present, relative to the anchor.
.sw_heb_one <- function(d, hier, anchor, variance_floor) {
  batches <- sort(unique(d$batch))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  studies <- sort(unique(hier$study[match(batches, hier$batch)]))
  y <- d$value
  donor <- factor(d$reference_id)
  X <- stats::model.matrix(~ 0 + donor)
  Zs <- outer(hier$study[match(d$batch, hier$batch)], studies, `==`) * 1
  Zb <- outer(d$batch, batches, `==`) * 1
  start <- log(rep(max(stats::var(y), 1e-6), 3) * c(0.1, 0.1, 0.5))
  lo <- rep(log(variance_floor), 3)
  fit <- stats::optim(
    start, .sw_heb_reml, y = y, X = X, Zs = Zs, Zb = Zb,
    method = "L-BFGS-B", lower = lo, upper = lo + 40
  )
  th <- exp(fit$par)
  V <- th[[1]] * tcrossprod(Zs) + th[[2]] * tcrossprod(Zb) + th[[3]] * diag(length(y))
  Vinv <- chol2inv(chol(V))
  XtVX <- crossprod(X, Vinv %*% X)
  P <- Vinv - Vinv %*% X %*% solve(XtVX, crossprod(X, Vinv))
  Z <- cbind(Zs, Zb)
  G <- diag(c(rep(th[[1]], ncol(Zs)), rep(th[[2]], ncol(Zb))))
  u <- G %*% crossprod(Z, P %*% y)
  pev <- G - G %*% crossprod(Z, P %*% Z) %*% G
  # beta_b = study effect + run effect; contrast against the anchor.
  idx_s <- match(hier$study[match(batches, hier$batch)], studies)
  L <- matrix(0, length(batches), ncol(Z))
  L[cbind(seq_along(batches), idx_s)] <- 1
  L[cbind(seq_along(batches), length(studies) + seq_along(batches))] <- 1
  a <- match(anchor, batches)
  L <- sweep(L, 2, L[a, ], `-`)
  est <- as.vector(L %*% u)
  se <- sqrt(pmax(rowSums((L %*% pev) * L), 0))
  list(batch = batches, estimate = est, se = se, theta = th)
}

.sw_heb_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, variance_floor, ...) {
  if (!is.data.frame(hierarchy) || !all(c("batch", "study") %in% names(hierarchy))) {
    stop("heb needs hierarchy = data frame with columns 'batch' and 'study'.", call. = FALSE)
  }
  if (missing(variance_floor) || !is.numeric(variance_floor) || variance_floor <= 0) {
    stop("heb needs variance_floor > 0 (no default).", call. = FALSE)
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
  key <- paste(s$summ$population, s$summ$marker, sep = "\r")
  rows <- lapply(split(seq_len(nrow(s$summ)), key), function(ii) {
    d <- s$summ[ii, , drop = FALSE]
    r <- .sw_heb_one(d, hierarchy, anchor, variance_floor)
    if (is.null(r)) {
      return(NULL)
    }
    n_don <- vapply(r$batch, function(b) length(unique(d$reference_id[d$batch == b])), integer(1))
    n_aliq <- vapply(r$batch, function(b) length(unique(d$aliquot[d$batch == b])), integer(1))
    tibble::tibble(
      batch = r$batch, population = d$population[[1]], marker = d$marker[[1]],
      estimate = r$estimate, se_donor = r$se, n_donors = n_don, n_aliquots = n_aliq,
      n_aliquots_anchor = n_aliq[match(anchor, r$batch)],
      support = if_else_chr(r$batch == anchor, "anchor", "joint_reml"), source = "fit"
    )
  })
  effects <- dplyr::bind_rows(rows)
  list(effects = effects, targets = targets, excluded = s$excluded)
}

if_else_chr <- function(cond, yes, no) ifelse(cond, yes, no)

sw_corrector_register_engine(
  "heb",
  fit = .sw_heb_fit,
  apply = function(model, xm, cells, new_batches) .sw_ps_apply(model, xm, cells, new_batches),
  can_freeze = TRUE,
  description = "Joint REML mixed model: donor fixed, study and run random; effect = BLUP difference to the anchor."
)
