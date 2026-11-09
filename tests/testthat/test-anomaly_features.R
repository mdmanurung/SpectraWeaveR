# tests/testthat/test-anomaly_features.R
# Unit tests for R/anomaly_features.R. Pure matrices; no flow packages needed.

make_files <- function(n_files = 6, n = 200, seed = 1) {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    matrix(
      stats::rnorm(n * 2),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    )
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

feat <- function(data, ...) {
  args <- list(
    data = data,
    channels = c("A", "B"),
    features = c("quantile", "emd"),
    quantile_step = 0.25,
    n_events = 100,
    seed = 1,
    scale = "test",
    aggregate_from = "all",
    agg_size = 300
  )
  args[names(list(...))] <- list(...)
  do.call(sw_anomaly_features, args)
}

test_that("sw_emd_1d is exact on known cases", {
  expect_equal(sw_emd_1d(c(0, 1), c(2, 3)), 2)
  x <- stats::rnorm(50)
  expect_equal(sw_emd_1d(x, x + 0.5), 0.5)
  expect_equal(sw_emd_1d(x, x), 0)
  # Unequal sizes: point mass at 0 vs {0, 2} moves half the mass by 2.
  expect_equal(sw_emd_1d(0, c(0, 2)), 1)
  # Symmetric.
  y <- stats::rnorm(30, 1)
  expect_equal(sw_emd_1d(x, y), sw_emd_1d(y, x))
  expect_error(sw_emd_1d(c(1, NA), 1), "without NA")
})

test_that("required arguments have no default", {
  d <- make_files()
  expect_error(sw_anomaly_features(), "no default")
  expect_error(
    sw_anomaly_features(d, channels = "A", features = "quantile"),
    "n_events"
  )
  expect_error(
    sw_anomaly_features(
      d,
      "A",
      "quantile",
      n_events = 10,
      seed = 1,
      scale = "x"
    ),
    "quantile_step"
  )
  expect_error(
    sw_anomaly_features(
      d,
      "A",
      "emd",
      n_events = 10,
      seed = 1,
      scale = "x"
    ),
    "aggregate_from|agg_size"
  )
  expect_error(
    feat(d, aggregate_from = "reference"),
    "reference"
  )
})

test_that("quantile columns exclude min and max", {
  d <- make_files()
  res <- feat(d, features = "quantile", quantile_step = 0.25)
  expect_s3_class(res, "sw_anomaly_features")
  expect_equal(
    names(res),
    c("file", "A_q0.25", "A_q0.5", "A_q0.75", "B_q0.25", "B_q0.5", "B_q0.75")
  )
  res10 <- feat(d, features = "quantile", quantile_step = 0.1)
  expect_equal(ncol(res10) - 1, 2 * 9)
})

test_that("a shifted file has a larger EMD than clean files", {
  d <- make_files(n_files = 8, n = 400)
  d$f1[, "A"] <- d$f1[, "A"] + 2
  res <- feat(d, features = "emd", n_events = 300, agg_size = 800)
  expect_equal(res$file[which.max(res$A_emd)], "f1")
  expect_true(res$A_emd[res$file == "f1"] > 1)
  expect_true(all(res$A_emd[res$file != "f1"] < 0.5))
})

test_that("results do not depend on list order and restore the RNG", {
  d <- make_files()
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  a <- feat(d)
  after <- stats::runif(1)
  expect_equal(before, after)
  b <- feat(d[rev(names(d))])
  expect_equal(unclass(a)[-1], unclass(b)[-1])
  expect_equal(a$file, b$file)
  c2 <- feat(d, seed = 2)
  expect_false(isTRUE(all.equal(a$A_emd, c2$A_emd)))
})

test_that("reference pooling uses only the named files", {
  d <- make_files(n_files = 6, n = 300)
  d$f6[, "A"] <- d$f6[, "A"] + 5
  ref <- feat(
    d,
    features = "emd",
    aggregate_from = "reference",
    reference = c("f1", "f2", "f3"),
    agg_size = 600
  )
  all_ <- feat(d, features = "emd", aggregate_from = "all", agg_size = 600)
  # An anomalous file in the pool pulls the pool towards it, so f1 looks
  # further from a pool that holds f6 than from a clean pool.
  expect_true(all_$A_emd[all_$file == "f1"] > ref$A_emd[ref$file == "f1"])
})

test_that("a file with too few events is an error", {
  d <- make_files(n = 50)
  expect_error(feat(d, n_events = 100), "fewer than `n_events`")
  d$f1 <- d$f1[, "A", drop = FALSE]
  expect_error(
    feat(
      make_files() |>
        (\(x) {
          x$f1 <- x$f1[, "A", drop = FALSE]
          x
        })()
    ),
    "lacks channel"
  )
})
