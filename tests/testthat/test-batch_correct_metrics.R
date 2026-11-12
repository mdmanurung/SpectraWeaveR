# Tests for R/batch_correct_metrics.R (A2c: metric records for candidates).

make_cm <- function(shift = 1, n = 300, seed = 11) {
  set.seed(seed)
  ids <- rep(paste0("s", 1:4), each = n)
  batch <- c(s1 = "b1", s2 = "b1", s3 = "b2", s4 = "b2")
  in_b2 <- ids %in% c("s3", "s4")
  x <- cbind(CD3 = stats::rnorm(4 * n), CD4 = stats::rnorm(4 * n))
  x[in_b2, "CD3"] <- x[in_b2, "CD3"] + shift
  fixed <- x
  fixed[in_b2, "CD3"] <- fixed[in_b2, "CD3"] - shift
  list(x = x, fixed = fixed, ids = ids, batch = batch)
}

cm_args <- function(d, ...) {
  utils::modifyList(
    list(
      candidates = list(fixed = d$fixed),
      uncorrected = d$x,
      sample_id = d$ids,
      batch = d$batch,
      population = NULL,
      markers = c("CD3", "CD4"),
      identity_id = "uncorrected",
      n_quantiles = 50,
      min_cells = 50,
      weights = "samples",
      k = 30,
      perplexity = 10,
      alpha = 0.05,
      n_events = 200,
      seed = 3,
      min_events = 100
    ),
    list(...),
    keep.null = TRUE
  )
}

test_that("a planted batch shift is found by EMD and removed by the fix", {
  d <- make_cm()
  out <- do.call(sw_correct_metrics, cm_args(d))
  emd <- out[out$metric == "emd_within_population" & out$marker == "CD3", ]
  expect_gt(emd$value[emd$candidate_id == "uncorrected"], 0.8)
  expect_lt(emd$value[emd$candidate_id == "fixed"], 0.3)
  cd4 <- out[out$metric == "emd_within_population" & out$marker == "CD4", ]
  expect_lt(max(cd4$value), 0.3)
  expect_true(all(out$calibration == "uncalibrated"))
})

test_that("identity has zero MAD compression and structure distance", {
  d <- make_cm()
  out <- do.call(sw_correct_metrics, cm_args(d))
  idn <- out[out$candidate_id == "uncorrected", ]
  expect_equal(idn$value[idn$metric == "mad_compression_fraction"], c(0, 0))
  expect_equal(idn$value[idn$metric == "joint_structure_distance"], 0)
})

test_that("mixing rows exist and mixing improves after correction", {
  d <- make_cm()
  out <- do.call(sw_correct_metrics, cm_args(d))
  lisi <- out[
    out$metric == "lisi_batch" & out$comparison_id == "all_batches",
  ]
  expect_gt(
    lisi$value[lisi$candidate_id == "fixed"],
    lisi$value[lisi$candidate_id == "uncorrected"]
  )
  expect_true(all(c("lisi_batch", "kbet_rejection") %in% out$metric))
  expect_true(any(out$comparison_id == "batch:b1"))
})

test_that("records round-trip through sw_batch_candidates_from_metrics", {
  d <- make_cm()
  out <- do.call(sw_correct_metrics, cm_args(d))
  cand <- sw_batch_candidates_from_metrics(out, identity_id = "uncorrected")
  expect_setequal(cand$candidate_id, c("uncorrected", "fixed"))
  expect_gt(cand$emd_relative_reduction[cand$candidate_id == "fixed"], 0.5)
  expect_equal(cand$coverage[cand$candidate_id == "fixed"], 1)
})

test_that("the seed reproduces the draw and the caller's RNG is restored", {
  d <- make_cm()
  set.seed(99)
  before <- .Random.seed
  a <- do.call(sw_correct_metrics, cm_args(d))
  expect_identical(.Random.seed, before)
  b <- do.call(sw_correct_metrics, cm_args(d))
  expect_identical(a, b)
  other <- do.call(sw_correct_metrics, cm_args(d, seed = 4))
  m <- a$metric %in% c("lisi_batch", "kbet_rejection")
  expect_false(isTRUE(all.equal(a$value[m], other$value[m])))
})

test_that("populations get their own rows and thin ones are not assessed", {
  d <- make_cm()
  pop <- rep("T", length(d$ids))
  # 295 of 300 cells per sample are "rare"; "T" keeps 5 cells per sample.
  for (s in unique(d$ids)) {
    pop[d$ids == s][1:295] <- "rare"
  }
  out <- do.call(
    sw_correct_metrics,
    cm_args(d, population = pop, min_cells = 100, n_events = 200)
  )
  expect_setequal(unique(out$population), c("T", "rare"))
  thin <- out[
    out$population == "T" & out$metric == "emd_within_population",
  ]
  expect_true(all(thin$status != "ok"))
  expect_true(all(is.na(thin$value)))
  mix <- out[out$population == "T" & out$metric == "lisi_batch", ]
  expect_true(all(mix$status == "not_assessed"))
  expect_true(all(is.na(mix$value)))
  ok <- out[
    out$population == "rare" & out$metric == "emd_within_population",
  ]
  expect_true(all(ok$status == "ok"))
})

test_that("samples under n_events are left out of mixing with a warning", {
  d <- make_cm()
  expect_warning(
    out <- do.call(sw_correct_metrics, cm_args(d, n_events = 301)),
    "s1"
  )
  expect_false(any(out$metric == "lisi_batch"))
  expect_true(any(out$metric == "emd_within_population"))
})

test_that("weights = equal sets every weight to 1", {
  d <- make_cm()
  out <- do.call(sw_correct_metrics, cm_args(d, weights = "equal"))
  expect_true(all(out$weight == 1))
  out2 <- do.call(sw_correct_metrics, cm_args(d))
  expect_true(any(out2$weight != 1))
})

test_that("missing required arguments and bad inputs stop", {
  d <- make_cm()
  full <- cm_args(d)
  for (nm in names(full)) {
    args <- full[setdiff(names(full), nm)]
    expect_error(
      do.call(sw_correct_metrics, args),
      paste0("'", nm, "'"),
      fixed = TRUE
    )
  }
  expect_error(
    do.call(sw_correct_metrics, cm_args(d, weights = "median")),
    "weights"
  )
  expect_error(
    do.call(
      sw_correct_metrics,
      cm_args(d, candidates = list(uncorrected = d$fixed))
    ),
    "reserved"
  )
  expect_error(
    do.call(
      sw_correct_metrics,
      cm_args(d, candidates = list(a = d$fixed[-1, ]))
    ),
    "same rows"
  )
  expect_error(
    do.call(sw_correct_metrics, cm_args(d, markers = c("CD3", "XX"))),
    "XX"
  )
})

test_that("thin samples with a large pool give no ok EMD record of weight 0", {
  d <- make_cm(n = 40)
  out <- suppressWarnings(do.call(
    sw_correct_metrics,
    cm_args(d, n_events = 40)
  ))
  emd <- out[out$metric == "emd_within_population", ]
  expect_false(any(emd$status == "ok" & emd$weight == 0))
  expect_true(all(emd$status == "too_few_cells"))
  expect_true(all(is.na(emd$value)))
})
