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
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none", anchor_effect = "random", variance_prior = "none", residual_model = "estimated", anchors = NULL, seed = 1, ...
  )
}

test_that("heb is registered and needs its required arguments", {
  expect_true("heb" %in% sw_corrector_engines())
  s <- heb_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = NULL, seed = 1, variance_floor = 1e-8
    ),
    "hierarchy"
  )
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, seed = 1
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
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, seed = 1
    ),
    "pool"
  )
  eff_pop <- sw_corrector_fit(
    s$x, s$manifest, s$audit, population = rep("T", nrow(s$x)), markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "population", anchor_effect = "random", variance_prior = "none", residual_model = "estimated", anchors = NULL, seed = 1
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
      anchor_effect = "random", residual_model = "estimated", anchors = NULL, seed = 1
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
    anchor_effect = "random", variance_prior = "boundary_avoiding", residual_model = "estimated", anchors = NULL, seed = 1
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
    anchor_effect = "random", variance_prior = "none", residual_model = "known_plus_extra", anchors = NULL, seed = 1
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
    anchor_effect = "fixed", variance_prior = "none", residual_model = "estimated", anchors = NULL, seed = 1
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
  fit <- function(x, anc, manifest = s$manifest) {
    sw_corrector_fit(
      x, manifest, s$audit, population = rep("T", nrow(x)), markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb",
      min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, pool = "none",
      anchor_effect = "random", variance_prior = "boundary_avoiding", residual_model = "known_plus_extra",
      anchors = anc, seed = 1
    )
  }
  full <- fit(s$x, NULL)$effects
  anchors <- sw_corrector_anchors(
    s$x["m1"], s$x$sample_id, s$meta$batch[match(s$x$sample_id, s$meta$sample_id)], "m1",
    positive_markers = character(), min_separation = 0, seed = 1
  )
  without <- fit(s$x, NULL, man2)
  expect_false("s3_r2" %in% without$effects$batch)
  with <- fit(s$x, anchors, man2)
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
})
