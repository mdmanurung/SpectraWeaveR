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
# V = th1 * M1 + th2 * M2 + th3 * I, with M1 = Z K1 Z' and M2 = Z K2 Z'.
.sw_heb_reml <- function(log_theta, y, X, M1, M2) {
  th <- exp(log_theta)
  V <- th[[1]] * M1 + th[[2]] * M2 + th[[3]] * diag(length(y))
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

# Design pieces for one (population, marker): NULL when it cannot give a contrast to the anchor.
# anchor_effect = "random": study and run effects of every batch are random and the anchor's own
#   effect is a BLUP; the correction is the BLUP difference to the anchor (derivation v0).
# anchor_effect = "fixed": the anchor-relative effects b_b = beta_b - beta_anchor are the random
#   effects, b_anchor = 0 exactly, with covariance
#   tau_s^2 (S_ij - S_ia - S_ja + 1) + tau_b^2 (delta_ij + 1), S = same-study indicator (derivation v1).
.sw_heb_design <- function(d, hier, anchor, anchor_effect) {
  batches <- sort(unique(d$batch))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  study_of <- hier$study[match(batches, hier$batch)]
  studies <- sort(unique(study_of))
  nb <- length(batches)
  if (anchor_effect == "random") {
    Zs <- outer(hier$study[match(d$batch, hier$batch)], studies, `==`) * 1
    Zb <- outer(d$batch, batches, `==`) * 1
    Z <- cbind(Zs, Zb)
    q1 <- length(studies)
    K1 <- diag(c(rep(1, q1), rep(0, nb)))
    K2 <- diag(c(rep(0, q1), rep(1, nb)))
    L <- matrix(0, nb, ncol(Z))
    L[cbind(seq_len(nb), match(study_of, studies))] <- 1
    L[cbind(seq_len(nb), q1 + seq_len(nb))] <- 1
    L <- sweep(L, 2, L[match(anchor, batches), ], `-`)
  } else {
    na <- batches[batches != anchor]
    s_n <- study_of[match(na, batches)]
    s_a <- study_of[match(anchor, batches)]
    Z <- outer(d$batch, na, `==`) * 1
    ia <- as.numeric(s_n == s_a)
    K1 <- outer(s_n, s_n, `==`) * 1 - outer(ia, rep(1, length(na))) - outer(rep(1, length(na)), ia) + 1
    K2 <- diag(length(na)) + 1
    L <- matrix(0, nb, length(na))
    L[cbind(match(na, batches), seq_along(na))] <- 1
  }
  list(
    d = d, batches = batches, y = d$value,
    X = stats::model.matrix(~ 0 + factor(d$reference_id)),
    Z = Z, K1 = K1, K2 = K2, L = L,
    M1 = Z %*% K1 %*% t(Z), M2 = Z %*% K2 %*% t(Z)
  )
}

# Variance components by REML. Several designs share one set of components when pooled:
# the objective is the sum of their REML criteria.
.sw_heb_theta <- function(designs, variance_floor) {
  scale <- max(mean(vapply(designs, function(z) stats::var(z$y), numeric(1))), 1e-6)
  obj <- function(lt) {
    sum(vapply(designs, function(z) .sw_heb_reml(lt, z$y, z$X, z$M1, z$M2), numeric(1)))
  }
  lo <- rep(log(variance_floor), 3)
  fit <- stats::optim(
    log(rep(scale, 3) * c(0.1, 0.1, 0.5)), obj,
    method = "L-BFGS-B", lower = lo, upper = lo + 40
  )
  exp(fit$par)
}

# BLUPs and prediction error for one design at given variance components.
.sw_heb_blup <- function(z, th, anchor) {
  y <- z$y
  X <- z$X
  V <- th[[1]] * z$M1 + th[[2]] * z$M2 + th[[3]] * diag(length(y))
  Vinv <- chol2inv(chol(V))
  XtVX <- crossprod(X, Vinv %*% X)
  P <- Vinv - Vinv %*% X %*% solve(XtVX, crossprod(X, Vinv))
  G <- th[[1]] * z$K1 + th[[2]] * z$K2
  u <- G %*% crossprod(z$Z, P %*% y)
  pev <- G - G %*% crossprod(z$Z, P %*% z$Z) %*% G
  list(
    batch = z$batches, estimate = as.vector(z$L %*% u),
    se = sqrt(pmax(rowSums((z$L %*% pev) * z$L), 0))
  )
}

.sw_heb_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, variance_floor, pool, anchor_effect, ...) {
  if (!is.data.frame(hierarchy) || !all(c("batch", "study") %in% names(hierarchy))) {
    stop("heb needs hierarchy = data frame with columns 'batch' and 'study'.", call. = FALSE)
  }
  if (missing(variance_floor) || !is.numeric(variance_floor) || variance_floor <= 0) {
    stop("heb needs variance_floor > 0 (no default).", call. = FALSE)
  }
  if (missing(pool) || !pool %in% c("none", "population", "all")) {
    stop("heb needs pool = 'none', 'population' or 'all' (no default).", call. = FALSE)
  }
  if (missing(anchor_effect) || !anchor_effect %in% c("random", "fixed")) {
    stop("heb needs anchor_effect = 'random' or 'fixed' (no default).", call. = FALSE)
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
  designs <- lapply(split(seq_len(nrow(s$summ)), key), function(ii) {
    .sw_heb_design(s$summ[ii, , drop = FALSE], hierarchy, anchor, anchor_effect)
  })
  designs <- designs[!vapply(designs, is.null, logical(1))]
  # Which designs share variance components.
  pool_key <- switch(pool,
    none = names(designs),
    population = vapply(designs, function(z) z$d$marker[[1]], character(1)),
    all = rep("all", length(designs))
  )
  names(pool_key) <- names(designs)
  theta <- lapply(split(designs, pool_key), .sw_heb_theta, variance_floor = variance_floor)
  rows <- lapply(names(designs), function(nm) {
    z <- designs[[nm]]
    r <- .sw_heb_blup(z, theta[[pool_key[[nm]]]], anchor)
    d <- z$d
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
