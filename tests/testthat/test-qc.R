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

# =========================================================================
# Integration tests — actual PeacoQC execution
# =========================================================================

test_that("sw_signal_qc runs on synthetic flowFrame", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")

  set.seed(42)
  # Create 500 events with 4 channels; include some extreme values
  # that PeacoQC should flag
  n <- 500
  bv421 <- c(abs(rnorm(480, mean = 50000, sd = 15000)),
             rep(262143, 20))
  pe <- c(abs(rnorm(480, mean = 40000, sd = 12000)),
          rep(0, 20))
  mat <- matrix(
    c(abs(rnorm(n, 100000, 30000)),
      abs(rnorm(n, 60000, 20000)),
      bv421, pe),
    ncol = 4,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)

  result <- sw_signal_qc(ff, channels = c("BV421-A", "PE-A"))

  expect_true(methods::is(result$FinalFF, "flowFrame"))
  expect_true(is.logical(result$GoodCells))
  expect_equal(length(result$GoodCells), n)
  expect_true(result$n_removed >= 0)
  expect_true(result$pct_removed >= 0 && result$pct_removed <= 100)
  expect_equal(nrow(result$FinalFF), sum(result$GoodCells))
})

test_that("sw_signal_qc_batch processes multiple samples", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")

  set.seed(42)
  make_ff <- function(n = 200) {
    mat <- matrix(
      abs(rnorm(n * 3, mean = 50000, sd = 15000)),
      ncol = 3,
      dimnames = list(NULL, c("FSC-A", "BV421-A", "PE-A"))
    )
    flowCore::flowFrame(mat)
  }

  ff_list <- list(S1 = make_ff(), S2 = make_ff())
  result <- sw_signal_qc_batch(ff_list, channels = c("BV421-A", "PE-A"))

  expect_true(is.list(result))
  expect_true("cleaned" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_equal(length(result$cleaned), 2)
  expect_s3_class(result$summary, "tbl_df")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("sample", "n_before", "n_after",
                     "n_removed", "pct_removed") %in% names(result$summary)))
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
