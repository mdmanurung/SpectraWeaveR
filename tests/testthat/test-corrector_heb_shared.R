# Tests for engine "heb_shared": lme4 oracle, borrowing, abstention and invariance.

# 3 donors, 6 batches (3 studies x 2 runs), populations P and Q, 2 aliquots each.
# cells_pq[[b]] can cut the cells of one (batch, population) to test borrowing.
hebs_sim <- function(seed = 11, n_cells = 40, few = NULL, drop_anchor_pop = NULL) {
  set.seed(seed)
  batches <- paste0("s", rep(1:3, each = 2), "_r", rep(1:2, 3))
  study_eff <- c(s1 = 0, s2 = 0.6, s3 = -0.4)
  run_eff <- stats::setNames(stats::rnorm(6, 0, 0.3), batches)
  dev_eff <- matrix(stats::rnorm(12, 0, 0.2), 6, 2, dimnames = list(batches, c("P", "Q")))
  donor_eff <- outer(c(D1 = 0, D2 = 1, D3 = -1), c(P = 0, Q = 2), `+`)
  meta <- list()
  cells <- list()
  for (b in batches) for (d in c("D1", "D2", "D3")) for (a in 1:2) for (k in c("P", "Q")) {
    sid <- paste(d, b, a, sep = "_")
    n <- if (!is.null(few) && b == few[[1]] && k == few[[2]]) 3 else n_cells
    if (!is.null(drop_anchor_pop) && b == "s1_r1" && k == drop_anchor_pop) n <- 3
    mu <- donor_eff[d, k] + study_eff[[substr(b, 1, 2)]] + run_eff[[b]] + dev_eff[b, k]
    cells[[length(cells) + 1]] <- data.frame(sample_id = sid, pop = k, m1 = stats::rnorm(n, mu, 0.4))
    if (k == "P" && !any(vapply(meta, function(m) m$sample_id == sid, logical(1)))) {
      meta[[length(meta) + 1]] <- data.frame(
        sample_id = sid, batch = b, donor = d, aliquot = sid,
        reference_role = "technical_replicate", reference_id = d
      )
    }
  }
  meta <- do.call(rbind, meta)
  cells <- do.call(rbind, cells)
  list(
    x = tibble::as_tibble(cells["m1"] |> cbind(sample_id = cells$sample_id)),
    pop = cells$pop, meta = meta,
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

hebs_fit <- function(s, engine = "heb_shared", ...) {
  sw_corrector_fit(
    s$x[c("sample_id", "m1")], s$manifest, s$audit, population = s$pop, markers = "m1",
    transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = engine,
    min_cells = 10, hierarchy = s$hierarchy, variance_floor = 1e-8, seed = 1, ...
  )
}

test_that("heb_shared is registered and needs hierarchy and variance_floor", {
  expect_true("heb_shared" %in% sw_corrector_engines())
  s <- hebs_sim()
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = s$pop, markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb_shared",
      min_cells = 10, hierarchy = NULL, seed = 1, variance_floor = 1e-8
    ),
    "hierarchy"
  )
  expect_error(
    sw_corrector_fit(
      s$x, s$manifest, s$audit, population = s$pop, markers = "m1",
      transform = list(name = "asinh", cofactor = 5), anchor = "s1_r1", engine = "heb_shared",
      min_cells = 10, hierarchy = s$hierarchy, seed = 1
    ),
    "variance_floor"
  )
})

test_that("heb_shared effects match an lme4 REML fit with a batch-by-population term", {
  skip_if_not_installed("lme4")
  s <- hebs_sim()
  eff <- hebs_fit(s)$effects
  d <- stats::aggregate(m1 ~ sample_id + pop, data = cbind(s$x, pop = s$pop), FUN = mean)
  d <- merge(d, s$meta, by = "sample_id")
  d <- merge(d, s$hierarchy, by = "batch")
  d$dp <- paste(d$donor, d$pop)
  d$bp <- paste(d$batch, d$pop)
  fit <- lme4::lmer(
    m1 ~ 0 + dp + (1 | study) + (1 | batch) + (1 | bp), data = d, REML = TRUE,
    control = lme4::lmerControl(check.conv.grad = lme4::.makeCC("ignore", tol = 1e-3))
  )
  re <- lme4::ranef(fit)
  beta <- function(b, k) re$study[substr(b, 1, 2), 1] + re$batch[b, 1] + re$bp[paste(b, k), 1]
  for (k in c("P", "Q")) {
    want <- vapply(s$hierarchy$batch, function(b) beta(b, k) - beta("s1_r1", k), numeric(1))
    e <- eff[eff$population == k, ]
    expect_equal(unname(e$estimate[match(s$hierarchy$batch, e$batch)]), unname(want), tolerance = 2e-3)
    expect_equal(e$estimate[e$batch == "s1_r1"], 0)
  }
})

test_that("a population dropped by min_cells in a batch borrows the shared effect", {
  s <- hebs_sim(few = c("s2_r2", "Q"))
  eff <- hebs_fit(s)$effects
  q <- eff[eff$population == "Q" & eff$batch == "s2_r2", ]
  p <- eff[eff$population == "P" & eff$batch == "s2_r2", ]
  expect_equal(nrow(q), 1)
  expect_equal(q$support, "shared_borrowed")
  expect_equal(q$n_donors, 0L)
  expect_true(is.finite(q$estimate) && is.finite(q$se_donor))
  # Same shared part as P; the deviation is zero, so the intervals are wider than for P.
  expect_gt(q$se_donor, p$se_donor)
  # The single-population engine has no effect there.
  v0 <- hebs_fit(s, engine = "heb", pool = "none", anchor_effect = "random", variance_prior = "none",
                 residual_model = "estimated", anchors = NULL, anchor_discrepancy = "none", scale_model = "none")$effects
  expect_equal(nrow(v0[v0$population == "Q" & v0$batch == "s2_r2", ]), 0)
})

test_that("a population with no anchor summary abstains", {
  s <- hebs_sim(drop_anchor_pop = "Q")
  model <- hebs_fit(s)
  expect_false("Q" %in% model$effects$population)
  expect_true("P" %in% model$effects$population)
  res <- sw_corrector_integrate(
    model, s$x[c("sample_id", "m1")], s$manifest, population = s$pop,
    transform = list(name = "asinh", cofactor = 5), unsupported = "abstain"
  )
  q_act <- res$actions[res$actions$population == "Q", ]
  expect_true(all(q_act$action == "abstain"))
  expect_true(all(res$actions$action[res$actions$population == "P"] == "apply"))
})

test_that("heb_shared leaves earlier corrected values unchanged when rows are integrated again", {
  s <- hebs_sim()
  model <- hebs_fit(s)
  tf <- list(name = "asinh", cofactor = 5)
  old <- sw_corrector_integrate(model, s$x[c("sample_id", "m1")], s$manifest, population = s$pop, transform = tf, unsupported = "abstain")
  again <- sw_corrector_integrate(old$model, s$x[c("sample_id", "m1")], s$manifest, population = s$pop, transform = tf, unsupported = "abstain")
  expect_identical(old$x, again$x)
})
