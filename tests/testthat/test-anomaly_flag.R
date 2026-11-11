# tests/testthat/test-anomaly_flag.R
# Unit tests for R/anomaly_flag.R. Feature tables are simulated directly.

make_feats <- function(n = 40, p = 4, bad = character(0), seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  colnames(x) <- paste0("v", seq_len(p))
  df <- tibble::as_tibble(x)
  df$file <- paste0("f", seq_len(n))
  df <- df[c("file", colnames(x))]
  for (b in bad) {
    df[df$file == b, -1] <- df[df$file == b, -1] + 6
  }
  df
}

test_that("outlier flag: required arguments have no default", {
  skip_if_not_installed("isotree")
  f <- make_feats()
  expect_error(sw_anomaly_flag_outlier(), "no default")
  expect_error(
    sw_anomaly_flag_outlier(
      f,
      ntrees = 100,
      sample_size = 20,
      ndim = 1,
      seed = 1
    ),
    "cutoff"
  )
})

test_that("outlier flag ranks a shifted file first", {
  skip_if_not_installed("isotree")
  f <- make_feats(bad = "f7")
  res <- sw_anomaly_flag_outlier(
    f,
    ntrees = 200,
    sample_size = 32,
    ndim = 1,
    cutoff = 0.5,
    seed = 1
  )
  expect_s3_class(res, "sw_anomaly_flag")
  expect_equal(res$file[which.max(res$score)], "f7")
  expect_true(res$flag[res$file == "f7"])
  # Same seed, same scores.
  res2 <- sw_anomaly_flag_outlier(
    f,
    ntrees = 200,
    sample_size = 32,
    ndim = 1,
    cutoff = 0.5,
    seed = 1
  )
  expect_equal(res$score, res2$score)
})

test_that("outlier flag checks its arguments", {
  skip_if_not_installed("isotree")
  f <- make_feats()
  expect_error(
    sw_anomaly_flag_outlier(
      f,
      100,
      sample_size = 500,
      ndim = 1,
      cutoff = 0.5,
      seed = 1
    ),
    "larger than the number of files"
  )
  expect_error(
    sw_anomaly_flag_outlier(f, 100, 20, ndim = 9, cutoff = 0.5, seed = 1),
    "larger than the number of features"
  )
  expect_error(
    sw_anomaly_flag_outlier(f, 100, 20, 1, cutoff = 1.5, seed = 1),
    "cutoff"
  )
})

test_that("novelty flag finds a shifted test file and leaves reference unflagged", {
  skip_if_not_installed("mclust")
  f <- make_feats(n = 60, p = 3, bad = "f60")
  ref <- paste0("f", 1:50)
  res <- sw_anomaly_flag_novelty(f, reference = ref, quantile = 0.05)
  expect_s3_class(res, "sw_anomaly_flag")
  expect_true(all(is.na(res$flag[res$role == "reference"])))
  expect_true(res$flag[res$file == "f60"])
  # The shifted file has the lowest density; a clean test file may still fall
  # under the reference quantile now and then (about 5% of them).
  expect_equal(res$file[which.min(res$log_density)], "f60")
  expect_lte(sum(res$flag[res$file %in% paste0("f", 51:59)]), 2)
  expect_true(!is.null(attr(res, "model")$model))
})

test_that("novelty flag stops when reference files do not outnumber features", {
  skip_if_not_installed("mclust")
  f <- make_feats(n = 30, p = 8)
  expect_error(
    sw_anomaly_flag_novelty(f, reference = paste0("f", 1:8), quantile = 0.05),
    "8 reference files, 8 features"
  )
})

test_that("novelty flag checks inputs", {
  skip_if_not_installed("mclust")
  f <- make_feats(n = 30, p = 3)
  expect_error(sw_anomaly_flag_novelty(), "no default")
  expect_error(
    sw_anomaly_flag_novelty(f, reference = "nope", quantile = 0.05),
    "must name files"
  )
  expect_error(
    sw_anomaly_flag_novelty(f, reference = paste0("f", 1:20), quantile = 2),
    "quantile"
  )
})
