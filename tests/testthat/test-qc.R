# tests/testthat/test-qc.R
# Unit tests for R/qc.R — Signal quality control (PeacoQC)

test_that("sw_signal_qc rejects non-flowFrame input", {
  skip_if_not_installed("flowCore")

  expect_error(sw_signal_qc(data.frame(x = 1)),
               "flowFrame")
})

test_that("sw_signal_qc validates IT_limit range", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 2,
                dimnames = list(NULL, c("BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)

  expect_error(sw_signal_qc(ff, IT_limit = 0), "between 0 and 1")
  expect_error(sw_signal_qc(ff, IT_limit = 1), "between 0 and 1")
  expect_error(sw_signal_qc(ff, IT_limit = -0.5), "between 0 and 1")
})

test_that("sw_signal_qc validates MAD", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 2,
                dimnames = list(NULL, c("BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)

  expect_error(sw_signal_qc(ff, MAD = 0), "positive number")
  expect_error(sw_signal_qc(ff, MAD = -1), "positive number")
})

test_that("sw_signal_qc_batch rejects empty list", {
  expect_error(sw_signal_qc_batch(list()),
               "non-empty list")
})

test_that("sw_signal_qc_batch rejects non-list input", {
  expect_error(sw_signal_qc_batch("not_a_list"),
               "non-empty list")
})

test_that("sw_qc_summary flags high-removal samples", {
  # Simulate qc_results structure
  qc_results <- list(
    summary = tibble::tibble(
      sample = c("S1", "S2", "S3"),
      n_before = c(1000, 1000, 1000),
      n_after = c(800, 600, 950),
      n_removed = c(200, 400, 50),
      pct_removed = c(20.0, 40.0, 5.0)
    )
  )

  result <- expect_warning(
    sw_qc_summary(qc_results, threshold = 30),
    "exceeded"
  )

  expect_s3_class(result, "tbl_df")
  expect_true("flagged" %in% names(result))
  expect_equal(result$flagged, c(FALSE, TRUE, FALSE))
})

test_that("sw_qc_summary returns no warnings when all below threshold", {
  qc_results <- list(
    summary = tibble::tibble(
      sample = c("S1", "S2"),
      n_before = c(1000, 1000),
      n_after = c(950, 900),
      n_removed = c(50, 100),
      pct_removed = c(5.0, 10.0)
    )
  )

  expect_no_warning(sw_qc_summary(qc_results, threshold = 30))
})

test_that("sw_qc_summary rejects invalid input", {
  expect_error(sw_qc_summary(list()), "output from sw_signal_qc_batch")
  expect_error(sw_qc_summary("not_a_list"), "output from sw_signal_qc_batch")
})

test_that("sw_qc_summary validates threshold", {
  qc_results <- list(
    summary = tibble::tibble(
      sample = "S1",
      n_before = 100,
      n_after = 90,
      n_removed = 10,
      pct_removed = 10.0
    )
  )

  expect_error(sw_qc_summary(qc_results, threshold = -1),
               "between 0 and 100")
  expect_error(sw_qc_summary(qc_results, threshold = 101),
               "between 0 and 100")
})
