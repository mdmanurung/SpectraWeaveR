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
#' \code{variance_prior} (required) is \code{"none"} for plain REML or
#' \code{"boundary_avoiding"} for the gamma(2) prior of Chung et al. (2013) on
#' the study and run standard deviations, which keeps them off zero when one
#' donor per batch cannot separate them from the residual.
#' \code{residual_model} (required) is \code{"estimated"} for one residual
#' variance, or \code{"known_plus_extra"}: each summary's residual is its known
#' sampling variance, var(cells) / n, plus one estimated extra variance.
#' \code{anchors} (required) is \code{NULL} or the result of
#' \code{sw_corrector_anchors()}; internal anchors enter as extra rows.
#' \code{anchor_discrepancy} (required) is \code{"none"} or \code{"batch"}:
#' with \code{"batch"}, the anchor rows of each batch share a random
#' discrepancy \eqn{d_b \sim N(0, \tau_d^2)} that is not part of the batch
#' effect, so anchors that disagree with the reference donors lose weight.
#' \code{anchor_batches} (required when \code{anchors} is given) is
#' \code{"all"} or \code{"unreferenced"}: with \code{"unreferenced"}, a
#' (population, marker) design keeps anchor rows only in batches without a
#' reference row and in the anchor batch, so batches with reference donors are
#' fitted from their donors alone.
#' \code{scale_model} (required) is \code{"none"} for a pure shift, or
#' \code{"log_sd"}: the same mixed model is fitted to the log standard
#' deviation of each reference aliquot summary (known variance
#' \eqn{1/(2(n-1))}), and the cells of batch \eqn{b} become
#' \eqn{c + e^{-\delta_b}(x - \mathrm{shift} - c)}, where \eqn{\delta_b} is the
#' log-SD contrast to the anchor and \eqn{c} the batch's shifted reference
#' median. The effects gain \code{scale}, \code{scale_se} (of \eqn{\delta_b})
#' and \code{scale_centre}; a batch without reference cells in a population
#' keeps the shift only. New batches at integrate time get the shift only.
#' \code{se_donor} holds the posterior standard error of the effect, from the
#' prediction error covariance; \code{df} holds its Satterthwaite degrees of
#' freedom from the inverse REML information, with components at the floor
#' treated as known (\code{Inf} when all are). A batch at integrate time that the model has
#' not seen is estimated as in \code{paired_shift}: its own reference cells
#' against the frozen anchor summaries, without shrinkage. Effects of known
#' batches are never changed.
#'
#' @name sw_corrector_heb
#' @seealso \code{\link{sw_corrector_fit}}, \code{\link{sw_corrector_paired_shift}}
NULL

# Negative REML log-likelihood (up to a constant) for log variance components.
# V = th1 * M1 + th2 * M2 + th3 * I + diag(D), with M1 = Z K1 Z' and M2 = Z K2 Z'.
# With anchor_discrepancy = "batch" (plan A11), V gains th5 * M3, M3 = Z K3 Z': the anchor rows of one
# batch share d_b ~ N(0, th5), which never enters the batch contrast.
# D holds known sampling variances of the summaries (residual_model = "known_plus_extra"), else 0.
# ra marks internal-anchor rows (plan A6): their extra residual variance is th4, donor rows keep th3.
.sw_heb_reml <- function(log_theta, y, X, M1, M2, D = 0, ra = 0, M3 = NULL) {
  th <- exp(log_theta)
  V <- th[[1]] * M1 + th[[2]] * M2 + diag(.sw_heb_resid(th, D, ra), length(y))
  if (!is.null(M3)) {
    V <- V + th[[5]] * M3 # plan A11: per-batch anchor discrepancy
  }
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

# Residual variance per row: th3 for donor rows, th4 for anchor rows, plus the known part D.
.sw_heb_resid <- function(th, D, ra) {
  th4 <- if (length(th) >= 4) th[[4]] else 0
  th[[3]] * (1 - ra) + th4 * ra + D
}

# Design pieces for one (population, marker): NULL when it cannot give a contrast to the anchor.
# anchor_effect = "random": study and run effects of every batch are random and the anchor's own
#   effect is a BLUP; the correction is the BLUP difference to the anchor (derivation v0).
# anchor_effect = "fixed": the anchor-relative effects b_b = beta_b - beta_anchor are the random
#   effects, b_anchor = 0 exactly, with covariance
#   tau_s^2 (S_ij - S_ia - S_ja + 1) + tau_b^2 (delta_ij + 1), S = same-study indicator (derivation v1).
.sw_heb_design <- function(d, hier, anchor, anchor_effect, residual_model, anchor_rows = NULL,
                           anchor_discrepancy = "none") {
  # Internal anchors (plan A0): extra rows with their own intercept, same batch random effects.
  na <- if (is.null(anchor_rows)) 0L else nrow(anchor_rows)
  b_all <- c(d$batch, if (na > 0) anchor_rows$batch)
  batches <- sort(unique(b_all))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  study_of <- hier$study[match(batches, hier$batch)]
  studies <- sort(unique(study_of))
  nb <- length(batches)
  if (anchor_effect == "random") {
    Zs <- outer(hier$study[match(b_all, hier$batch)], studies, `==`) * 1
    Zb <- outer(b_all, batches, `==`) * 1
    Z <- cbind(Zs, Zb)
    q1 <- length(studies)
    K1 <- diag(c(rep(1, q1), rep(0, nb)))
    K2 <- diag(c(rep(0, q1), rep(1, nb)))
    L <- matrix(0, nb, ncol(Z))
    L[cbind(seq_len(nb), match(study_of, studies))] <- 1
    L[cbind(seq_len(nb), q1 + seq_len(nb))] <- 1
    L <- sweep(L, 2, L[match(anchor, batches), ], `-`)
  } else {
    nab <- batches[batches != anchor]
    s_n <- study_of[match(nab, batches)]
    s_a <- study_of[match(anchor, batches)]
    Z <- outer(b_all, nab, `==`) * 1
    ia <- as.numeric(s_n == s_a)
    K1 <- outer(s_n, s_n, `==`) * 1 - outer(ia, rep(1, length(nab))) - outer(rep(1, length(nab)), ia) + 1
    K2 <- diag(length(nab)) + 1
    L <- matrix(0, nb, length(nab))
    L[cbind(match(nab, batches), seq_along(nab))] <- 1
  }
  # Donor indicators; outer() also covers a design with one donor, which model.matrix refuses.
  X <- outer(d$reference_id, sort(unique(d$reference_id)), `==`) * 1
  # Known part of each summary's residual variance: var(cells) / n (Fay-Herriot style).
  D <- if (residual_model == "known_plus_extra") d$var_mean else rep(0, nrow(d))
  if (na > 0) {
    X <- rbind(cbind(X, 0), cbind(matrix(0, na, ncol(X)), 1))
    D <- c(D, anchor_rows$var_mean)
  }
  ra <- c(rep(0, nrow(d)), rep(1, na))
  # Plan A11: one discrepancy d_b per batch for its anchor rows, kept out of the contrast L.
  K3 <- NULL
  if (na > 0 && anchor_discrepancy == "batch") {
    q <- ncol(Z)
    pad <- function(K) {
      out <- matrix(0, q + nb, q + nb)
      out[seq_len(q), seq_len(q)] <- K
      out
    }
    Z <- cbind(Z, outer(b_all, batches, `==`) * ra)
    K1 <- pad(K1)
    K2 <- pad(K2)
    K3 <- diag(c(rep(0, q), rep(1, nb)))
    L <- cbind(L, matrix(0, nrow(L), nb))
  }
  list(
    d = d, batches = batches, y = c(d$value, if (na > 0) anchor_rows$value), X = X,
    Z = Z, K1 = K1, K2 = K2, K3 = K3, L = L,
    M1 = Z %*% K1 %*% t(Z), M2 = Z %*% K2 %*% t(Z),
    M3 = if (is.null(K3)) NULL else Z %*% K3 %*% t(Z),
    D = D, ra = ra, anchor_rows = anchor_rows
  )
}

# Variance components by REML. Several designs share one set of components when pooled:
# the objective is the sum of their REML criteria.
#
# variance_prior = "boundary_avoiding" adds Chung et al. (2013, PMID 24092484) gamma(shape 2, rate -> 0)
# priors on the study and run SDs: -log p(sigma_k) = -log sigma_k = -0.5 log theta_k, once per pooled set.
# It keeps those components off zero when the data cannot separate them from the residual (one donor per
# batch). The residual variance is not penalised. "none" is plain REML.
.sw_heb_theta <- function(designs, variance_floor, variance_prior) {
  scale <- max(mean(vapply(designs, function(z) stats::var(z$y), numeric(1))), 1e-6)
  penalty <- if (variance_prior == "boundary_avoiding") function(lt) -0.5 * (lt[[1]] + lt[[2]]) else function(lt) 0
  obj <- function(lt) {
    sum(vapply(designs, function(z) .sw_heb_reml(lt, z$y, z$X, z$M1, z$M2, z$D, z$ra, z$M3), numeric(1))) +
      penalty(lt)
  }
  # A 4th component (anchor-row extra variance) only when a design in this set has anchor rows.
  # A 5th (anchor discrepancy, plan A11) only when a design carries M3; no prior on it.
  k <- if (any(vapply(designs, function(z) !is.null(z$M3), logical(1)))) {
    5
  } else if (any(vapply(designs, function(z) any(z$ra > 0), logical(1)))) {
    4
  } else {
    3
  }
  lo <- rep(log(variance_floor), k)
  fit <- stats::optim(
    log(rep(scale, k) * c(0.1, 0.1, 0.5, 0.5, 0.1)[seq_len(k)]), obj,
    method = "L-BFGS-B", lower = lo, upper = lo + 40
  )
  # Covariance of the log components from the inverse observed REML information.
  # A component pinned at the floor is treated as known (user decision 2026-11-09):
  # its rows and columns are zero. A non-positive curvature direction gets no variance.
  free <- fit$par - lo > 1e-4
  cov_log <- matrix(0, k, k)
  if (any(free)) {
    H <- stats::optimHess(fit$par, obj)[free, free, drop = FALSE]
    e <- eigen((H + t(H)) / 2, symmetric = TRUE)
    pos <- e$values > 1e-10 * max(abs(e$values), 1)
    cov_log[free, free] <- e$vectors[, pos, drop = FALSE] %*%
      (t(e$vectors[, pos, drop = FALSE]) / e$values[pos])
  }
  list(theta = exp(fit$par), cov_log = cov_log)
}

# BLUPs and prediction error variance of each anchor contrast, from Henderson's
# mixed-model equations scaled by sigma^2:
#   C = [X'X, X'Z; Z'X, Z'Z + sigma^2 G^-1],  C [b; u] = [X'y; Z'y],  PEV(u) = sigma^2 (C^-1)_uu.
# Same model as V-based GLS; the form G - G Z'PZ G cancels catastrophically once G / sigma^2
# exceeds about 1e6, which happens when sigma^2 is near the floor.
# With known sampling variances D the residual is heteroscedastic, R = diag(th3 + D), and the
# unscaled equations are used: C = W' R^-1 W + [0, 0; 0, G^-1], PEV(u) = (C^-1)_uu.
.sw_heb_pev <- function(z, th) {
  G <- th[[1]] * z$K1 + th[[2]] * z$K2
  if (!is.null(z$K3)) {
    G <- G + th[[5]] * z$K3
  }
  p <- ncol(z$X)
  W <- cbind(z$X, z$Z)
  if (any(z$D > 0) || any(z$ra > 0)) {
    rinv <- 1 / .sw_heb_resid(th, z$D, z$ra)
    C <- crossprod(W * sqrt(rinv))
    iu <- p + seq_len(ncol(z$Z))
    C[iu, iu] <- C[iu, iu] + solve(G)
    R <- chol(C)
    sol <- backsolve(R, forwardsolve(t(R), crossprod(W, rinv * z$y)))
    Cinv_uu <- chol2inv(R)[iu, iu, drop = FALSE]
    return(list(u = sol[iu], var = pmax(rowSums((z$L %*% Cinv_uu) * z$L), 0)))
  }
  C <- crossprod(W)
  iu <- p + seq_len(ncol(z$Z))
  C[iu, iu] <- C[iu, iu] + th[[3]] * solve(G)
  R <- chol(C)
  sol <- backsolve(R, forwardsolve(t(R), crossprod(W, z$y)))
  Cinv_uu <- chol2inv(R)[iu, iu, drop = FALSE]
  list(u = sol[iu], var = pmax(th[[3]] * rowSums((z$L %*% Cinv_uu) * z$L), 0))
}

# BLUPs, prediction error and Satterthwaite df for one design at given variance components.
# df = 2 v^2 / (g' C g), v the contrast's prediction error variance, g its gradient in the
# log components (central differences) and C their covariance. cov_log = NULL gives df = NA.
.sw_heb_blup <- function(z, th, anchor, cov_log = NULL) {
  p <- .sw_heb_pev(z, th)
  u <- p$u
  df <- rep(NA_real_, length(z$batches))
  if (!is.null(cov_log)) {
    h <- 1e-4
    nk <- length(th)
    g <- vapply(seq_len(nk), function(k) {
      if (all(cov_log[k, ] == 0)) {
        return(rep(0, length(p$var)))
      }
      step <- replace(rep(0, nk), k, h)
      (.sw_heb_pev(z, th * exp(step))$var - .sw_heb_pev(z, th * exp(-step))$var) / (2 * h)
    }, numeric(length(p$var)))
    g <- matrix(g, nrow = length(p$var))
    denom <- rowSums((g %*% cov_log) * g)
    df <- ifelse(denom > 0, 2 * p$var^2 / denom, Inf)
    df[z$batches == anchor] <- NA_real_
  }
  list(
    batch = z$batches, estimate = as.vector(z$L %*% u),
    se = sqrt(p$var), df = df
  )
}

.sw_heb_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, variance_floor, pool, anchor_effect,
                        variance_prior, residual_model, anchors, anchor_discrepancy, anchor_batches = NULL,
                        scale_model, ...) {
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
  if (missing(variance_prior) || !variance_prior %in% c("none", "boundary_avoiding")) {
    stop("heb needs variance_prior = 'none' or 'boundary_avoiding' (no default).", call. = FALSE)
  }
  if (missing(residual_model) || !residual_model %in% c("estimated", "known_plus_extra")) {
    stop("heb needs residual_model = 'estimated' or 'known_plus_extra' (no default).", call. = FALSE)
  }
  if (missing(anchors) || !(is.null(anchors) || (is.list(anchors) && is.data.frame(anchors$modes)))) {
    stop("heb needs anchors = NULL or the result of sw_corrector_anchors() (no default).", call. = FALSE)
  }
  if (missing(anchor_discrepancy) || !(identical(anchor_discrepancy, "none") || identical(anchor_discrepancy, "batch"))) {
    stop("heb needs anchor_discrepancy = 'none' or 'batch' (no default; plan A11).", call. = FALSE)
  }
  # Plan A13: needed only when anchors are supplied; NULL stands for "not given".
  if (!is.null(anchors) &&
    !(identical(anchor_batches, "all") || identical(anchor_batches, "unreferenced"))) {
    stop("heb with anchors needs anchor_batches = 'all' or 'unreferenced' (no default; plan A13).", call. = FALSE)
  }
  if (missing(scale_model) || !(identical(scale_model, "none") || identical(scale_model, "log_sd"))) {
    stop("heb needs scale_model = 'none' or 'log_sd' (no default; plan F12).", call. = FALSE)
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
  anchor_for <- .sw_heb_anchor_rows(anchors, anc, anchor, hierarchy)
  designs <- lapply(split(seq_len(nrow(s$summ)), key), function(ii) {
    d <- s$summ[ii, , drop = FALSE]
    ar <- anchor_for(d$population[[1]], d$marker[[1]])
    # Plan A13 "unreferenced": anchor rows only in batches without a reference row in this design, plus the
    # anchor batch that ties the anchor intercept; none at all when no such batch has anchor rows.
    if (!is.null(ar) && identical(anchor_batches, "unreferenced")) {
      ar <- ar[ar$batch == anchor | !ar$batch %in% d$batch, , drop = FALSE]
      if (!any(ar$batch != anchor)) ar <- NULL
    }
    .sw_heb_design(d, hierarchy, anchor, anchor_effect, residual_model,
                   anchor_rows = ar, anchor_discrepancy = anchor_discrepancy)
  })
  designs <- designs[!vapply(designs, is.null, logical(1))]
  blups <- .sw_heb_fit_designs(designs, pool, anchor, variance_floor, variance_prior)
  rows <- lapply(names(designs), function(nm) {
    z <- designs[[nm]]
    r <- blups[[nm]]
    d <- z$d
    n_don <- vapply(r$batch, function(b) length(unique(d$reference_id[d$batch == b])), integer(1))
    n_aliq <- vapply(r$batch, function(b) length(unique(d$aliquot[d$batch == b])), integer(1))
    tibble::tibble(
      batch = r$batch, population = d$population[[1]], marker = d$marker[[1]],
      estimate = r$estimate, se_donor = r$se, df = r$df, n_donors = n_don, n_aliquots = n_aliq,
      n_aliquots_anchor = n_aliq[match(anchor, r$batch)],
      support = dplyr::case_when(
        r$batch == anchor ~ "anchor",
        n_don == 0 ~ "joint_reml_internal_anchor",
        .default = "joint_reml"
      ), source = "fit"
    )
  })
  effects <- dplyr::bind_rows(rows)
  if (scale_model == "log_sd") {
    effects <- .sw_heb_scale(effects, xm, cells, s$summ, hierarchy, anchor, anchor_effect, residual_model,
                             pool, variance_floor, variance_prior)
  }
  list(effects = effects, targets = targets, excluded = s$excluded)
}

# Fit pooled variance components, then the BLUP contrasts of every design (a named list by design).
.sw_heb_fit_designs <- function(designs, pool, anchor, variance_floor, variance_prior) {
  # Which designs share variance components.
  pool_key <- switch(pool,
    none = names(designs),
    population = vapply(designs, function(z) z$d$marker[[1]], character(1)),
    all = rep("all", length(designs))
  )
  names(pool_key) <- names(designs)
  theta <- lapply(split(designs, pool_key), .sw_heb_theta,
                  variance_floor = variance_floor, variance_prior = variance_prior)
  lapply(stats::setNames(names(designs), names(designs)), function(nm) {
    tk <- theta[[pool_key[[nm]]]]
    .sw_heb_blup(designs[[nm]], tk$theta, anchor, tk$cov_log)
  })
}

# Plan F12, scale_model = "log_sd": the same mixed model on the log SD of each reference aliquot summary
# (known sampling variance 1 / (2 (n - 1))), no internal anchors. delta_b is the BLUP log-SD contrast to
# the anchor; the correction multiplies deviations from the batch's shifted reference median by
# scale = exp(-delta_b). Batches without reference cells in a population get no scale (NA).
.sw_heb_scale <- function(effects, xm, cells, summ, hierarchy, anchor, anchor_effect, residual_model,
                          pool, variance_floor, variance_prior) {
  ref <- cells[cells$is_reference, , drop = FALSE]
  ref_key <- paste(ref$reference_id, ref$aliquot, ref$batch, ref$population, sep = "\r")
  rows_of <- split(ref$row, ref_key)
  summ_key <- paste(summ$reference_id, summ$aliquot, summ$batch, summ$population, sep = "\r")
  ls <- summ
  ls$value <- vapply(seq_len(nrow(summ)), function(i) {
    log(stats::sd(xm[rows_of[[summ_key[[i]]]], summ$marker[[i]]]))
  }, numeric(1))
  ls$var_mean <- 1 / (2 * (summ$n_cells - 1))
  ls <- ls[is.finite(ls$value), , drop = FALSE]
  key <- paste(ls$population, ls$marker, sep = "\r")
  designs <- lapply(split(seq_len(nrow(ls)), key), function(ii) {
    .sw_heb_design(ls[ii, , drop = FALSE], hierarchy, anchor, anchor_effect, residual_model)
  })
  designs <- designs[!vapply(designs, is.null, logical(1))]
  blups <- .sw_heb_fit_designs(designs, pool, anchor, variance_floor, variance_prior)
  sc <- dplyr::bind_rows(lapply(names(designs), function(nm) {
    tibble::tibble(
      batch = blups[[nm]]$batch, population = designs[[nm]]$d$population[[1]],
      marker = designs[[nm]]$d$marker[[1]], log_sd_effect = blups[[nm]]$estimate, scale_se = blups[[nm]]$se
    )
  }))
  # Centre: median of the batch's reference cells for that population, minus the location effect.
  grp <- split(ref$row, paste(ref$batch, ref$population, sep = "\r"))
  centre <- dplyr::bind_rows(lapply(names(grp), function(k) {
    tibble::tibble(
      batch = ref$batch[match(k, paste(ref$batch, ref$population, sep = "\r"))],
      population = ref$population[match(k, paste(ref$batch, ref$population, sep = "\r"))],
      marker = colnames(xm), ref_median = apply(xm[grp[[k]], , drop = FALSE], 2, stats::median)
    )
  }))
  out <- dplyr::left_join(effects, sc, by = c("batch", "population", "marker"), relationship = "one-to-one")
  out <- dplyr::left_join(out, centre, by = c("batch", "population", "marker"), relationship = "one-to-one")
  ok <- !is.na(out$log_sd_effect) & !is.na(out$ref_median)
  out$scale <- unname(ifelse(ok, exp(-out$log_sd_effect), NA_real_))
  out$scale_centre <- unname(ifelse(ok, out$ref_median - out$estimate, NA_real_))
  out$scale_se <- unname(out$scale_se)
  out$scale_se[!ok] <- NA_real_
  out$log_sd_effect <- NULL
  out$ref_median <- NULL
  out
}

if_else_chr <- function(cond, yes, no) ifelse(cond, yes, no)

# Plan A9: an anchor joins the design of a population whose anchor-batch level for that marker lies nearer
# that anchor's mode centre (low or high) in the anchor batch. Returns a function (population, marker) -> rows.
.sw_heb_anchor_rows <- function(anchors, anc, anchor, hierarchy) {
  if (is.null(anchors)) {
    return(function(population, marker) NULL)
  }
  modes <- anchors$modes[anchors$modes$batch %in% hierarchy$batch, , drop = FALSE]
  centre <- dplyr::summarise(modes[modes$batch == anchor, , drop = FALSE], centre = mean(value), .by = c(marker, mode))
  level <- dplyr::summarise(anc, level = mean(value), .by = c(population, marker))
  used <- modes[modes$use, c("batch", "marker", "mode", "value", "var_mean"), drop = FALSE]
  function(population, marker) {
    lv <- level$level[level$population == population & level$marker == marker]
    cl <- centre$centre[centre$marker == marker & centre$mode == "low"]
    ch <- centre$centre[centre$marker == marker & centre$mode == "high"]
    if (length(lv) != 1 || length(cl) != 1 || length(ch) != 1) {
      return(NULL)
    }
    m <- if (abs(lv - cl) <= abs(lv - ch)) "low" else "high"
    rows <- used[used$marker == marker & used$mode == m, , drop = FALSE]
    if (nrow(rows) == 0) NULL else rows
  }
}

sw_corrector_register_engine(
  "heb",
  fit = .sw_heb_fit,
  apply = function(model, xm, cells, new_batches) .sw_ps_apply(model, xm, cells, new_batches),
  can_freeze = TRUE,
  description = "Joint REML mixed model: donor fixed, study and run random; effect = BLUP difference to the anchor."
)
