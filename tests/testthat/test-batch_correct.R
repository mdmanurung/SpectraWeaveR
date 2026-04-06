# tests/testthat/test-batch_correct.R
# Unit tests for R/batch_correct.R — Batch correction (cyCombine)

# ===========================================================================
# sw_prepare_for_correction
# ===========================================================================

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

# ===========================================================================
# sw_batch_correct
# ===========================================================================

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

# ===========================================================================
# sw_evaluate_correction
# ===========================================================================

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

# ===========================================================================
# sw_normalize — modular workflow step 1
# ===========================================================================

test_that("sw_normalize rejects non-data.frame", {
  expect_error(sw_normalize("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_normalize rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_normalize(df, markers = "CD3"),
               "batch")
})

test_that("sw_normalize rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_normalize(df, markers = character(0)),
               "non-empty")
})

test_that("sw_normalize rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_normalize(df, markers = "MISSING"),
               "not found")
})

test_that("sw_normalize validates norm_method", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_normalize(df, markers = "CD3", norm_method = "invalid"),
               "arg")
})

# ===========================================================================
# sw_create_som — modular workflow step 2
# ===========================================================================

test_that("sw_create_som rejects non-data.frame", {
  expect_error(sw_create_som("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_create_som rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_create_som(df, markers = character(0)),
               "non-empty")
})

test_that("sw_create_som rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_create_som(df, markers = "MISSING"),
               "not found")
})

# ===========================================================================
# sw_correct_data — modular workflow step 3
# ===========================================================================

test_that("sw_correct_data rejects non-data.frame", {
  expect_error(sw_correct_data("bad", label = 1L, markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_data rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_data(df, label = rep(1L, 10), markers = "CD3"),
               "batch")
})

test_that("sw_correct_data rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_data(df, label = rep(1L, 10), markers = "MISSING"),
               "not found")
})

test_that("sw_correct_data rejects non-numeric label", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_data(df, label = letters[1:10], markers = "CD3"),
               "numeric")
})

test_that("sw_correct_data rejects label length mismatch", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_data(df, label = c(1L, 2L), markers = "CD3"),
               "must equal nrow")
})

# ===========================================================================
# sw_detect_batch_effect
# ===========================================================================

test_that("sw_detect_batch_effect rejects non-data.frame", {
  expect_error(sw_detect_batch_effect("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_detect_batch_effect rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_detect_batch_effect(df, markers = "CD3"),
               "batch")
})

test_that("sw_detect_batch_effect rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_detect_batch_effect(df, markers = "MISSING"),
               "not found")
})

# ===========================================================================
# sw_compute_emd
# ===========================================================================

test_that("sw_compute_emd rejects non-data.frame", {
  expect_error(sw_compute_emd("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_compute_emd rejects missing label column", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_compute_emd(df, markers = "CD3"),
               "label")
})

test_that("sw_compute_emd rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  expect_error(sw_compute_emd(df, markers = "CD3"),
               "batch")
})

test_that("sw_compute_emd rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_compute_emd(df, markers = "MISSING"),
               "not found")
})

test_that("sw_compute_emd rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_compute_emd(df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_evaluate_emd
# ===========================================================================

test_that("sw_evaluate_emd rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_emd("bad", df, markers = "CD3"),
               "data.frame")
  expect_error(sw_evaluate_emd(df, "bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_evaluate_emd rejects missing columns", {
  df_ok <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  df_no_batch <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  df_no_label <- tibble::tibble(CD3 = rnorm(10), batch = "B1")

  expect_error(sw_evaluate_emd(df_no_batch, df_ok, markers = "CD3"),
               "batch")
  expect_error(sw_evaluate_emd(df_ok, df_no_label, markers = "CD3"),
               "label")
})

test_that("sw_evaluate_emd rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_emd(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_evaluate_emd rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_emd(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_evaluate_mad (cyCombine-backed)
# ===========================================================================

test_that("sw_evaluate_mad rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_mad("bad", df, markers = "CD3"),
               "data.frame")
  expect_error(sw_evaluate_mad(df, "bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_evaluate_mad rejects missing columns", {
  df_ok <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  df_no_batch <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  df_no_label <- tibble::tibble(CD3 = rnorm(10), batch = "B1")

  expect_error(sw_evaluate_mad(df_no_batch, df_ok, markers = "CD3"),
               "batch")
  expect_error(sw_evaluate_mad(df_ok, df_no_label, markers = "CD3"),
               "label")
})

test_that("sw_evaluate_mad rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_mad(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_evaluate_mad rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_evaluate_mad(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_plot_batch_densities
# ===========================================================================

test_that("sw_plot_batch_densities rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_plot_batch_densities("bad", df, markers = "CD3"),
               "data.frame")
})

test_that("sw_plot_batch_densities rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_densities(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_plot_batch_densities rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_densities(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_plot_batch_dimred
# ===========================================================================

test_that("sw_plot_batch_dimred rejects non-data.frame", {
  expect_error(sw_plot_batch_dimred("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_plot_batch_dimred rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = "MISSING"),
               "not found")
})

test_that("sw_plot_batch_dimred rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_plot_batch_dimred(df, markers = "CD3"),
               "batch")
})

test_that("sw_plot_batch_dimred rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = character(0)),
               "non-empty")
})

test_that("sw_plot_batch_dimred validates type argument", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = "CD3", type = "invalid"),
               "arg")
})
