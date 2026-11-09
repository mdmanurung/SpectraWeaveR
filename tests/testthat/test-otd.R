# Tests for R/otd.R (debiased Sinkhorn divergence and per-file OTD).

test_that("a pure translation gives |d|^2 / 2 at any blur", {
  set.seed(1)
  x <- matrix(stats::rnorm(300), ncol = 3)
  shifted <- sweep(x, 2, c(1, 0, 0.5), "+")
  for (b in c(0.2, 1)) {
    s <- sw_sinkhorn_divergence(x, shifted, blur = b, max_iter = 2000, tol = 1e-6)
    expect_equal(as.numeric(s), 0.625, tolerance = 1e-4)
    # A blur small against the point spacing converges slowly: the value is
    # already right, but the strict potential tolerance is not met, and the
    # flag says so.
    expect_equal(attr(s, "converged"), b >= 1)
  }
})

test_that("the divergence is zero on equal sets and symmetric", {
  set.seed(2)
  x <- matrix(stats::rnorm(200), ncol = 2)
  y <- matrix(stats::rnorm(240, sd = 1.5), ncol = 2)
  expect_equal(as.numeric(sw_sinkhorn_divergence(x, x, 0.5, 1000, 1e-8)), 0, tolerance = 1e-8)
  expect_equal(
    as.numeric(sw_sinkhorn_divergence(x, y, 0.5, 1000, 1e-8)),
    as.numeric(sw_sinkhorn_divergence(y, x, 0.5, 1000, 1e-8)),
    tolerance = 1e-6
  )
  expect_gt(as.numeric(sw_sinkhorn_divergence(x, y, 0.5, 1000, 1e-8)), 0)
})

test_that("an unconverged solve is reported, not hidden", {
  set.seed(3)
  x <- matrix(stats::rnorm(200), ncol = 2)
  y <- matrix(stats::rnorm(200, 3), ncol = 2)
  s <- sw_sinkhorn_divergence(x, y, blur = 0.05, max_iter = 2, tol = 1e-12)
  expect_false(attr(s, "converged"))
})

test_that("sw_otd_scores ranks a shifted file first and restores the RNG state", {
  set.seed(4)
  files <- lapply(1:6, \(i) matrix(stats::rnorm(600), ncol = 2, dimnames = list(NULL, c("A", "B"))))
  names(files) <- paste0("f", 1:6)
  files$f6[, "A"] <- files$f6[, "A"] + 1.5
  before <- .Random.seed
  sc <- sw_otd_scores(
    files,
    markers = c("A", "B"),
    n_events = 150,
    reference_events = 600,
    blur = 0.5,
    seed = 1,
    max_iter = 1000,
    tol = 1e-6
  )
  expect_identical(.Random.seed, before)
  expect_equal(sc$file[which.max(sc$otd)], "f6")
  expect_true(all(sc$converged))
})

test_that("sw_otd_scores and the divergence stop on bad input", {
  x <- matrix(1:4, 2)
  expect_error(sw_sinkhorn_divergence(x, x, 0, 10, 1e-6), "blur")
  expect_error(sw_sinkhorn_divergence(x, x, 1, 0, 1e-6), "max_iter")
  expect_error(sw_sinkhorn_divergence(x, matrix(1:3, 1), 1, 10, 1e-6), "same number")
  files <- list(a = cbind(A = 1:5), b = cbind(A = 1:5))
  expect_error(sw_otd_scores(files, "A", 10, 4, 1, 1, 10, 1e-6), "fewer")
  expect_error(sw_otd_scores(files, "Z", 2, 4, 1, 1, 10, 1e-6), "lacks")
  expect_error(sw_otd_scores(files["a"], "A", 2, 4, 1, 1, 10, 1e-6), "at least two")
})

test_that("sw_diagnose runs the otd module", {
  set.seed(5)
  files <- lapply(1:4, \(i) matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("A", "B"))))
  names(files) <- paste0("f", 1:4)
  dx <- sw_diagnose(
    files,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = list(
      markers = c("A", "B"),
      n_events = 100,
      reference_events = 400,
      blur = 0.5,
      seed = 1,
      max_iter = 500,
      tol = 1e-6
    )
  )
  expect_equal(attr(dx, "modules"), "otd")
  expect_equal(nrow(dx), 4)
  expect_true(all(dx$text == "converged=TRUE"))
})
