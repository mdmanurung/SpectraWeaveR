# tests/testthat/test-batch_correct.R
# Unit tests for R/batch_correct.R — Batch correction (cyCombine)

test_that("sw_prepare_for_correction creates correct tibble", {
  skip_if_not_installed("flowCore")

  # Create mock flowFrames
  mat1 <- matrix(c(1000, 2000, 3000, 4000, 5000, 6000), ncol = 2,
                 dimnames = list(NULL, c("CD3", "CD4")))
  mat2 <- matrix(c(1500, 2500, 3500, 4500, 5500, 6500), ncol = 2,
                 dimnames = list(NULL, c("CD3", "CD4")))

  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)

  ff_list <- list(S1 = ff1, S2 = ff2)

  sample_meta <- tibble::tibble(
    sample = c("S1", "S2"),
    batch = c("B1", "B2"),
    condition = c("Ctrl", "Stim")
  )

  result <- sw_prepare_for_correction(ff_list, sample_meta,
                                      markers = c("CD3", "CD4"),
                                      cofactor = 6000)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("CD3", "CD4", "sample", "batch", "condition", "id") %in%
                    names(result)))
  expect_equal(nrow(result), 6) # 3 cells per sample * 2 samples

  # Check arcsinh transformation was applied
  expect_equal(result$CD3[1], asinh(1000 / 6000), tolerance = 1e-10)
  expect_equal(result$CD4[1], asinh(4000 / 6000), tolerance = 1e-10)

  # Check metadata
  expect_true(all(result$sample[1:3] == "S1"))
  expect_true(all(result$batch[1:3] == "B1"))
})

test_that("sw_prepare_for_correction rejects unnamed list", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_prepare_for_correction(list(ff), tibble::tibble(sample = "S1",
                                                       batch = "B1"),
                              markers = c("CD3")),
    "named list"
  )
})

test_that("sw_prepare_for_correction rejects missing sample_meta columns", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_prepare_for_correction(list(S1 = ff),
                              tibble::tibble(sample = "S1"),
                              markers = c("CD3")),
    "missing required column"
  )
})

test_that("sw_prepare_for_correction rejects missing sample in metadata", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_prepare_for_correction(list(S1 = ff),
                              tibble::tibble(sample = "S2", batch = "B1"),
                              markers = c("CD3")),
    "not found in sample_meta"
  )
})

test_that("sw_prepare_for_correction validates cofactor", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_prepare_for_correction(list(S1 = ff),
                              tibble::tibble(sample = "S1", batch = "B1"),
                              markers = c("CD3"),
                              cofactor = -1),
    "positive number"
  )
})

test_that("sw_batch_correct rejects data without batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_batch_correct(df, markers = "CD3"),
               "batch")
})

test_that("sw_batch_correct rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_batch_correct(df, markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_batch_correct rejects non-data.frame", {
  expect_error(sw_batch_correct("not_a_df", markers = "CD3"),
               "data.frame")
})

test_that("sw_evaluate_correction computes MAD metrics", {
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(50, mean = 1), rnorm(50, mean = 3)),
    CD4 = c(rnorm(50, mean = 2), rnorm(50, mean = 4)),
    batch = rep(c("B1", "B2"), each = 50)
  )

  # Simulate correction: reduce batch difference
  corrected <- tibble::tibble(
    CD3 = c(rnorm(50, mean = 2), rnorm(50, mean = 2.1)),
    CD4 = c(rnorm(50, mean = 3), rnorm(50, mean = 3.1)),
    batch = rep(c("B1", "B2"), each = 50)
  )

  result <- sw_evaluate_correction(uncorrected, corrected,
                                   markers = c("CD3", "CD4"))

  expect_type(result, "list")
  expect_true(all(c("emd", "mad", "improved", "emd_reduction_pct") %in%
                    names(result)))
  expect_s3_class(result$mad, "tbl_df")
  expect_true("marker" %in% names(result$mad))
  expect_true("mad_before" %in% names(result$mad))
  expect_true("mad_after" %in% names(result$mad))

  # Correction should have improved things
  expect_true(result$improved)
})

test_that("sw_evaluate_correction handles single batch", {
  df <- tibble::tibble(
    CD3 = rnorm(20),
    batch = "B1"
  )

  result <- sw_evaluate_correction(df, df, markers = "CD3")
  expect_true(is.na(result$improved))
})

test_that("sw_evaluate_correction rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_evaluate_correction(df, df, markers = "CD3"),
               "batch")
})

test_that("sw_evaluate_correction rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_evaluate_correction(df, df, markers = "MISSING"),
               "not found")
})
