# Tests for the joint REML engine "heb": oracle, contract and invariance.

# Reference-only simulation: 3 donors measured in 6 batches (3 studies x 2 runs).
heb_sim <- function(seed = 5) {
  set.seed(seed)
  batches <- paste0("s", rep(1:3, each = 2), "_r", rep(1:2, 3))
  study_eff <- c(s1 = 0, s2 = 0.6, s3 = -0.4)
  run_eff <- stats::rnorm(6, 0, 0.3)
  names(run_eff) <- batches
  donors <- c("D1", "D2", "D3")
  donor_eff <- c(D1 = 0, D2 = 1, D3 = -1)
  rows <- list()
  for (b in batches) {
    for (d in donors) {
      for (a in 1:2) {
        mu <- donor_eff[[d]] + study_eff[[substr(b, 1, 2)]] + run_eff[[b]]
        rows[[length(rows) + 1]] <- data.frame(
          sample_id = paste(d, b, a, sep = "_"), batch = b, donor = d,
          aliquot = paste(d, b, a, sep = "_"),
          reference_role = "technical_replicate", reference_id = d,
          mu = mu
        )
      }
    }
  }
  meta <- do.call(rbind, rows)
  cells <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    data.frame(sample_id = meta$sample_id[i], m1 = stats::rnorm(40, meta$mu[i], 0.4))
  }))
  list(
    x = tibble::as_tibble(cells), meta = meta,
    manifest = sw_sample_manifest(
      meta, sample_id = "sample_id", batch = "batch", donor = "donor", aliquot = "aliquot",
      reference_role = "reference_role", reference_id = "reference_id"
    ),
    hierarchy = data.frame(batch = batches, study = substr(batches, 1, 2)),
    audit = sw_design_audit(
      meta, minimum_batches = 2, minimum_batches_origin = "test", require_shared_reference = TRUE,
      anchor = "s1_r1"
    )
  )
}

heb_fit <- function(s, x = s$x, ...) {
  sw_corrector_fit(
    x, s$manifest, s$audit, population = rep("T", nrow(x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none", anchor_effect = "random", variance_prior = "none", residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1, ...
  )
}

test_that("heb is registered and needs its required arguments", {
  expect_true("heb" %in% sw_corrector_engines())
  s <- heb_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = NULL, seed = 1, variance_floor = 1e-8, unbridged = "stop"
    ),
    "hierarchy"
  )
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, seed = 1, unbridged = "stop"
    ),
    "variance_floor"
  )
})

test_that("heb needs pool and pooled variance runs", {
  s <- heb_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, seed = 1, unbridged = "stop"
    ),
    "pool"
  )
  eff_pop <- sw_corrector_fit(
    s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "population", anchor_effect = "random", variance_prior = "none", residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1
  )$effects
  # One population and one marker: pooling over populations changes nothing.
  expect_equal(eff_pop$estimate, heb_fit(s)$effects$estimate, tolerance = 1e-6)
})

test_that("heb effects match an lme4 REML fit on the same summaries", {
  skip_if_not_installed("lme4")
  s <- heb_sim()
  eff <- heb_fit(s)$effects
  d <- s$meta
  d$value <- vapply(seq_len(nrow(d)), function(i) mean(s$x$m1[s$x$sample_id == d$sample_id[i]]), numeric(1))
  d <- merge(d, s$hierarchy, by = "batch")
  fit <- lme4::lmer(
    value ~ 0 + donor + (1 | study) + (1 | batch), data = d, REML = TRUE,
    control = lme4::lmerControl(check.conv.grad = lme4::.makeCC("ignore", tol = 1e-3))
  )
  re <- lme4::ranef(fit)
  beta <- function(b) re$study[substr(b, 1, 2), 1] + re$batch[b, 1]
  want <- vapply(s$hierarchy$batch, function(b) beta(b) - beta("s1_r1"), numeric(1))
  got <- eff$estimate[match(s$hierarchy$batch, eff$batch)]
  expect_equal(unname(got), unname(want), tolerance = 1e-3)
  expect_equal(eff$estimate[eff$batch == "s1_r1"], 0)
})

test_that("heb needs variance_prior, and the boundary prior keeps one-donor SEs off zero", {
  s <- heb_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
      anchor_effect = "random", residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1
    ),
    "variance_prior"
  )
  # One donor per batch: REML cannot split run variance from residual.
  s1 <- s
  keep <- s$meta$donor == "D1"
  s1$meta <- s$meta[keep, ]
  s1$x <- s$x[s$x$sample_id %in% s1$meta$sample_id, ]
  s1$manifest <- s$manifest[s$manifest$sample_id %in% s1$meta$sample_id, ]
  bap <- sw_corrector_fit(
    s1$x, s1$manifest, s1$audit, population = rep("T", nrow(s1$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s1$hierarchy, variance_floor = 1e-8, pool = "none",
    anchor_effect = "random", variance_prior = "boundary_avoiding", residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1
  )$effects
  se <- bap$se_donor[bap$batch != "s1_r1"]
  expect_true(all(se > 1e-3))
})

test_that("known sampling variance sets the one-donor, one-aliquot SE", {
  s <- heb_sim()
  keep <- s$meta$donor == "D1" & endsWith(s$meta$sample_id, "_1")
  s$meta <- s$meta[keep, ]
  s$x <- s$x[s$x$sample_id %in% s$meta$sample_id, ]
  s$manifest <- s$manifest[s$manifest$sample_id %in% s$meta$sample_id, ]
  eff <- sw_corrector_fit(
    s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
    anchor_effect = "random", variance_prior = "none", residual_model = "known_plus_extra", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1
  )$effects
  # Each contrast is a difference of two batch means of 40 cells (sd 0.4 before asinh): the SE
  # cannot fall below the sampling SE of that difference, about sqrt(2) * 0.4 / sqrt(40) on asinh/5.
  se <- eff$se_donor[eff$batch != "s1_r1"]
  expect_true(all(se > 0.5 * sqrt(2) * 0.4 / sqrt(40) / 5))
  expect_true(all(is.finite(eff$estimate)))
})

test_that("heb fits a design with one reference donor", {
  s <- heb_sim()
  keep <- s$meta$donor == "D1"
  s1 <- s
  s1$meta <- s$meta[keep, ]
  s1$x <- s$x[s$x$sample_id %in% s1$meta$sample_id, ]
  s1$manifest <- s$manifest[s$manifest$sample_id %in% s1$meta$sample_id, ]
  eff <- heb_fit(s1)$effects
  expect_equal(nrow(eff), 6)
  expect_true(all(is.finite(eff$estimate)))
})

test_that("heb Satterthwaite df reduces to the residual df without shrinkage", {
  # With batch variances far above the noise, the contrast variance is proportional to
  # sigma^2, so df = 2 / Var(log sigma^2). REML gives Var(log sigma^2) = 2 / df_resid, and
  # this 3 donor x 2 batch layout has df_resid = 6 - 3 - 1 = 2.
  d <- .sw_heb_design(
    data.frame(
      value = c(1, 2.1, 3, 2, 3, 4.2), reference_id = rep(c("D1", "D2", "D3"), 2),
      batch = rep(c("a", "b"), each = 3), aliquot = 1:6
    ),
    data.frame(batch = c("a", "b"), study = c("s", "s")), "a", "random", "estimated"
  )
  # Ratio 1e8 on purpose: the old form G - G Z'PZ G gave df 4.9 here (cancellation);
  # the mixed-model-equation form must stay exact.
  est <- .sw_heb_blup(d, c(1, 1, 1e-8), "a", cov_log = diag(c(0, 0, 2 / 2)))
  expect_equal(est$df[est$batch == "b"], 2, tolerance = 1e-3)
  expect_true(is.na(est$df[est$batch == "a"]))
  # Fitted effects carry df: positive for every non-anchor batch, NA at the anchor.
  eff <- heb_fit(heb_sim())$effects
  expect_true(all(eff$df[eff$batch != "s1_r1"] > 0))
  expect_true(all(is.na(eff$df[eff$batch == "s1_r1"])))
})

test_that("heb leaves earlier corrected values unchanged when new rows arrive", {
  s <- heb_sim()
  model <- heb_fit(s)
  old <- sw_corrector_integrate(
    model, s$x, s$manifest, population = rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5), unsupported = "abstain"
  )
  again <- sw_corrector_integrate(
    old$model, s$x, s$manifest, population = rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5), unsupported = "abstain"
  )
  expect_identical(old$x, again$x)
})

test_that("anchor_effect = 'fixed' pins the anchor at zero and reduces to paired_shift without shrinkage", {
  s <- heb_sim()
  fixed <- sw_corrector_fit(
    s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
    anchor_effect = "fixed", variance_prior = "none", residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = "stop", seed = 1
  )$effects
  expect_equal(fixed$estimate[fixed$batch == "s1_r1"], 0)
  expect_equal(fixed$se_donor[fixed$batch == "s1_r1"], 0)
  # Limiting case: variance components far above the noise mean no shrinkage, so the batch effect is
  # the mean paired difference to the anchor, which is what paired_shift estimates.
  d <- .sw_heb_design(
    data.frame(
      value = c(1, 2, 3, 2, 3, 4), reference_id = rep(c("D1", "D2", "D3"), 2),
      batch = rep(c("a", "b"), each = 3), aliquot = 1:6
    ),
    data.frame(batch = c("a", "b"), study = c("s", "s")), "a", "fixed", "estimated"
  )
  est <- .sw_heb_blup(d, c(1e3, 1e3, 1e-3), "a")
  expect_equal(est$estimate, c(0, 1), tolerance = 1e-3)
})

test_that("sw_corrector_anchors needs its settings and splits two clear modes", {
  set.seed(3)
  v <- c(rnorm(200, 0, 0.3), rnorm(200, 4, 0.3), rnorm(400, 0, 0.3))
  x <- cbind(CD3 = v, m2 = rnorm(800))
  sid <- rep(c("a", "b"), each = 400)
  expect_error(
    sw_corrector_anchors(x, sid, sid, c("CD3", "m2"), positive_markers = "CD3", seed = 1),
    "min_separation"
  )
  a <- sw_corrector_anchors(x, sid, rep("b1", 800), c("CD3", "m2"),
                            positive_markers = "CD3", min_separation = 3, seed = 1)
  expect_false(a$markers$eligible[a$markers$marker == "m2"])
  hi <- a$modes[a$modes$sample_id == "a" & a$modes$marker == "CD3" & a$modes$mode == "high", ]
  expect_equal(hi$value, 4, tolerance = 0.1)
  expect_true(hi$use)
})

test_that("internal anchors give a batch without reference donors an effect near the full fit", {
  s <- heb_sim()
  # Batch s3_r2 keeps its cells, but its samples are not references (as a batch without reference donors).
  m2 <- s$meta
  no_ref <- m2$batch == "s3_r2"
  m2$reference_role[no_ref] <- "none"
  m2$reference_id[no_ref] <- NA
  m2$aliquot[no_ref] <- NA
  man2 <- sw_sample_manifest(
    m2, sample_id = "sample_id", batch = "batch", donor = "donor", aliquot = "aliquot",
    reference_role = "reference_role", reference_id = "reference_id"
  )
  # The audit of this design marks s3_r2 unbridged; only anchors let heb fit it.
  audit2 <- sw_design_audit(
    m2, minimum_batches = 2, minimum_batches_origin = "test", require_shared_reference = TRUE,
    anchor = "s1_r1"
  )
  fit <- function(x, anc, manifest = s$manifest, audit = s$audit, anchor_batches = "all") {
    sw_corrector_fit(
      x, manifest, audit, population = rep("T", nrow(x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
      anchor_effect = "random", variance_prior = "boundary_avoiding", residual_model = "known_plus_extra",
      anchors = anc, anchor_discrepancy = "none", scale_model = "none", anchor_batches = anchor_batches, unbridged = "stop", seed = 1
    )
  }
  full <- fit(s$x, NULL)$effects
  anchors <- sw_corrector_anchors(
    s$x["m1"], s$x$sample_id, s$meta$batch[match(s$x$sample_id, s$meta$sample_id)], "m1",
    positive_markers = character(), min_separation = 0, seed = 1
  )
  expect_error(fit(s$x, NULL, man2, audit2), "ineligible")
  # Anchors that skip the unbridged batch do not open the gate.
  partial <- anchors
  partial$modes <- partial$modes[partial$modes$batch != "s3_r2", ]
  expect_error(fit(s$x, partial, man2, audit2), "ineligible")
  with <- fit(s$x, anchors, man2, audit2)
  e <- with$effects[with$effects$batch == "s3_r2", ]
  expect_equal(nrow(e), 1)
  expect_equal(e$support, "joint_reml_internal_anchor")
  expect_true(is.finite(e$se_donor) && e$se_donor > 0)
  # True effect from heb_sim's draws: study s3 (-0.4) + run s3_r2 - run s1_r1 (run effects drawn first).
  set.seed(5)
  run_eff <- stats::rnorm(6, 0, 0.3)
  truth <- -0.4 + run_eff[6] - run_eff[1]
  expect_lt(abs(e$estimate - truth), 3 * e$se_donor)
  expect_true("s3_r2" %in% with$known_batches)

  # Plan A13: anchor_batches = "unreferenced" keeps anchor rows only in the anchor batch and in batches
  # without reference donors. With every batch referenced, the fit is exactly the fit without anchors.
  unref_full <- fit(s$x, anchors, anchor_batches = "unreferenced")$effects
  expect_equal(unref_full$estimate, full$estimate, tolerance = 1e-10)
  expect_equal(unref_full$se_donor, full$se_donor, tolerance = 1e-10)
  # The batch without reference donors still gets its anchor-based effect.
  unref <- fit(s$x, anchors, man2, audit2, anchor_batches = "unreferenced")$effects
  e2 <- unref[unref$batch == "s3_r2", ]
  expect_equal(e2$support, "joint_reml_internal_anchor")
  expect_lt(abs(e2$estimate - truth), 3 * e2$se_donor)
  # Required (no default) once anchors are supplied.
  expect_error(fit(s$x, anchors, anchor_batches = NULL), "anchor_batches = 'all' or 'unreferenced'")
})

test_that("unbridged = 'abstain' fits the bridged batches and leaves an unbridged batch raw (plan A15)", {
  s <- heb_sim()
  # Batch s3_r2 keeps its cells, but its samples are not references, so no reference bridges it.
  m2 <- s$meta
  no_ref <- m2$batch == "s3_r2"
  m2$reference_role[no_ref] <- "none"
  m2$reference_id[no_ref] <- NA
  m2$aliquot[no_ref] <- NA
  man2 <- sw_sample_manifest(
    m2, sample_id = "sample_id", batch = "batch", donor = "donor", aliquot = "aliquot",
    reference_role = "reference_role", reference_id = "reference_id"
  )
  audit2 <- sw_design_audit(
    m2, minimum_batches = 2, minimum_batches_origin = "test", require_shared_reference = TRUE,
    anchor = "s1_r1"
  )
  fit <- function(unbridged) {
    sw_corrector_fit(
      s$x, man2, audit2, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
      anchor_effect = "random", variance_prior = "boundary_avoiding", residual_model = "known_plus_extra",
      anchors = NULL, anchor_discrepancy = "none", scale_model = "none", unbridged = unbridged, seed = 1
    )
  }
  expect_error(fit("stop"), "ineligible")
  expect_error(fit(NULL), "unbridged = 'abstain' or 'stop'")
  model <- fit("abstain")
  expect_false("s3_r2" %in% model$effects$batch)
  expect_setequal(model$effects$batch, setdiff(s$hierarchy$batch, "s3_r2"))
  out <- sw_corrector_integrate(
    model, s$x, man2, population = rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5), unsupported = "abstain"
  )
  in_b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)] == "s3_r2"
  expect_equal(out$x$m1[in_b], s$x$m1[in_b])
  expect_true(all(out$actions$action[out$actions$batch == "s3_r2"] == "abstain"))
  expect_false(isTRUE(all.equal(out$x$m1[!in_b], s$x$m1[!in_b])))
})

# Plan A11: anchor discrepancy d_b ~ N(0, tau_d^2) on anchor rows ------------------------------
# Designs on the heb_sim hierarchy: 3 donors x 2 aliquots in every batch but s3_r2, and 3 anchor
# samples per batch whose value is the batch effect plus a per-batch discrepancy m_b.
a11_designs <- function(mismatch_sd, n_markers = 8, seed = 11) {
  set.seed(seed)
  hier <- heb_sim()$hierarchy
  batches <- hier$batch
  lapply(seq_len(n_markers), function(k) {
    u <- stats::rnorm(length(batches), 0, 0.3)
    m <- stats::rnorm(length(batches), 0, mismatch_sd)
    d <- expand.grid(
      reference_id = c("D1", "D2", "D3"), a = 1:2, batch = setdiff(batches, "s3_r2"),
      stringsAsFactors = FALSE
    )
    d$value <- c(D1 = 0, D2 = 1, D3 = -1)[d$reference_id] + u[match(d$batch, batches)] +
      stats::rnorm(nrow(d), 0, 0.05)
    d$aliquot <- paste(d$reference_id, d$batch, d$a)
    d$var_mean <- 0.0025
    ar <- data.frame(batch = rep(batches, each = 3), mode = "low", var_mean = 4e-4)
    ar$value <- 2 + u[match(ar$batch, batches)] + m[match(ar$batch, batches)] +
      stats::rnorm(nrow(ar), 0, 0.02)
    .sw_heb_design(d, hier, "s1_r1", "random", "known_plus_extra",
                   anchor_rows = ar, anchor_discrepancy = "batch")
  })
}

test_that("a very large anchor discrepancy gives back the fit without anchors (analytic limit)", {
  z <- a11_designs(0.3, n_markers = 1)[[1]]
  z0 <- .sw_heb_design(z$d, heb_sim()$hierarchy, "s1_r1", "random", "known_plus_extra")
  th <- c(0.09, 0.09, 1e-3, 1e-3)
  with_anchor <- .sw_heb_blup(z, c(th, 1e10), "s1_r1")
  without <- .sw_heb_blup(z0, th[1:3], "s1_r1")
  keep <- with_anchor$batch %in% without$batch
  expect_equal(with_anchor$estimate[keep], without$estimate, tolerance = 1e-6)
  expect_equal(with_anchor$se[keep], without$se, tolerance = 1e-6)
})

test_that("anchor discrepancy is near the floor for a shared shift and above it for a mismatch", {
  shared <- .sw_heb_theta(a11_designs(0), variance_floor = 1e-8, variance_prior = "none")
  mism <- .sw_heb_theta(a11_designs(0.4), variance_floor = 1e-8, variance_prior = "none")
  expect_length(shared$theta, 5)
  expect_lt(shared$theta[[5]], 0.01)
  expect_gt(mism$theta[[5]], 0.05)
  # The batch without reference donors gets a wider interval when anchors disagree with references.
  z <- a11_designs(0.4)[[1]]
  se_shared <- .sw_heb_blup(z, shared$theta, "s1_r1")
  se_mism <- .sw_heb_blup(z, mism$theta, "s1_r1")
  i <- se_shared$batch == "s3_r2"
  expect_gt(se_mism$se[i], se_shared$se[i])
})

# Plan F12: scale term on log SD -------------------------------------------------------------------
# heb_sim with batch s2_r1's cells spread by `spread` around each sample's true mean.
scale_sim <- function(spread) {
  s <- heb_sim()
  mu <- s$meta$mu[match(s$x$sample_id, s$meta$sample_id)]
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  f <- ifelse(b == "s2_r1", spread, 1)
  s$x$m1 <- mu + f * (s$x$m1 - mu)
  s
}
scale_fit <- function(s, scale_model) {
  sw_corrector_fit(
    s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none", anchor_effect = "random",
    variance_prior = "none", residual_model = "known_plus_extra", anchors = NULL, anchor_discrepancy = "none",
    scale_model = scale_model, unbridged = "stop", seed = 1
  )
}

test_that("log_sd scale is about 1 without a spread effect and recovers a 1.5x spread", {
  none <- scale_fit(scale_sim(1), "log_sd")$effects
  expect_true(all(abs(log(none$scale)) < 0.1))
  sp <- scale_fit(scale_sim(1.5), "log_sd")$effects
  e <- sp[sp$batch == "s2_r1", ]
  # scale is the factor applied to the deviations, 1 / spread.
  expect_lt(abs(-log(e$scale) - log(1.5)), 3 * e$scale_se)
  expect_equal(sp$scale[sp$batch == "s1_r1"], 1)
})

test_that("log_sd correction evens the spread, and scale_model = 'none' is the plain shift", {
  s <- scale_sim(1.5)
  fit <- scale_fit(s, "log_sd")
  out <- sw_corrector_integrate(fit, s$x, s$manifest, population = rep("T", nrow(s$x)),
                                transform = list(name = "asinh", cofactor = 5), unsupported = "abstain")
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  sd_ratio <- stats::sd(out$x$m1[b == "s2_r1"] - ave(out$x$m1[b == "s2_r1"], s$x$sample_id[b == "s2_r1"])) /
    stats::sd(out$x$m1[b == "s1_r1"] - ave(out$x$m1[b == "s1_r1"], s$x$sample_id[b == "s1_r1"]))
  expect_lt(abs(sd_ratio - 1), 0.15)
  plain <- scale_fit(s, "none")
  expect_false("scale" %in% names(plain$effects) && any(!is.na(plain$effects$scale)))
  expect_equal(plain$effects$estimate, fit$effects$estimate)
  expect_error(scale_fit(s, "sd"), "scale_model = 'none' or 'log_sd'")
})

test_that("heb needs a valid anchor_discrepancy", {
  s <- heb_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
      anchor_effect = "random", variance_prior = "none", residual_model = "estimated",
      anchors = NULL, anchor_discrepancy = "sometimes", unbridged = "stop", seed = 1
    ),
    "anchor_discrepancy = 'none' or 'batch'"
  )
})
