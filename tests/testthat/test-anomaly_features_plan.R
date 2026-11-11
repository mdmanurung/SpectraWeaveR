# tests/testthat/test-anomaly_features_plan.R
# The shared EMD plan and the NA check of sw_anomaly_features().

make_two <- function(n_files = 4, n = 300) {
  set.seed(3)
  f <- lapply(seq_len(n_files), function(i) {
    matrix(stats::rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("A", "B")))
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

test_that("NA values in the chosen channels are an error, not dropped", {
  d <- make_two()
  d$f2[3, "A"] <- NA
  expect_error(
    sw_anomaly_features(
      d,
      channels = c("A", "B"),
      features = "emd",
      n_events = 100,
      seed = 1,
      scale = "x",
      aggregate_from = "all",
      agg_size = 200, reference_group = NULL
    ),
    "NA values"
  )
})

test_that("the shared plan gives the same EMD as the direct call", {
  x <- c(0, 1, 5)
  y <- c(0, 2)
  expect_equal(sw_emd_1d(x, y), .sw_emd_sorted(x, y, .sw_emd_plan(3, 2)))
  set.seed(5)
  a <- stats::rnorm(500)
  b <- stats::rnorm(5000)
  expect_equal(
    sw_emd_1d(a, b),
    .sw_emd_sorted(sort(a), sort(b), .sw_emd_plan(500, 5000))
  )
})

test_that("the feature EMD equals sw_emd_1d on the same draws", {
  d <- make_two(n_files = 3, n = 200)
  # One pooled file and n_events = all rows: the draw is a permutation, and EMD
  # does not depend on order, so the feature must equal the direct EMD.
  res <- sw_anomaly_features(
    d,
    channels = "A",
    features = "emd",
    n_events = 200,
    seed = 1,
    scale = "x",
    aggregate_from = "reference",
    agg_size = 200,
    reference = "f1", reference_group = NULL
  )
  expect_equal(res$A_emd[res$file == "f2"], sw_emd_1d(d$f2[, "A"], d$f1[, "A"]))
  expect_equal(res$A_emd[res$file == "f1"], 0)
})
