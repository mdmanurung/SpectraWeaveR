# tests/testthat/test-unmix_diagnostics.R
# Unit and integration tests for R/unmix_diagnostics.R

# =========================================================================
# sw_unmix_spillover_matrix — input validation
# =========================================================================

test_that("sw_unmix_spillover_matrix rejects non-matrix input", {
  expect_error(
    sw_unmix_spillover_matrix(data.frame(a = 1:3, b = 4:6)),
    "numeric matrix"
  )
})

test_that("sw_unmix_spillover_matrix rejects matrix without row names", {
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  expect_error(
    sw_unmix_spillover_matrix(mat),
    "row names"
  )
})

test_that("sw_unmix_spillover_matrix rejects single-row matrix", {
  mat <- matrix(rnorm(4), nrow = 1,
                dimnames = list("BV421", paste0("D", 1:4)))
  expect_error(
    sw_unmix_spillover_matrix(mat),
    "at least 2 rows"
  )
})

# =========================================================================
# sw_unmix_spillover_matrix — integration
# =========================================================================

test_that("SSM computation produces correct structure", {
  # Create a synthetic 4-fluorophore x 6-detector spectra matrix
  set.seed(42)
  spectra <- matrix(
    c(1.0, 0.1, 0.0, 0.0, 0.0, 0.0,  # BV421: peak in D1
      0.1, 1.0, 0.3, 0.0, 0.0, 0.0,  # FITC: peak in D2
      0.0, 0.0, 0.1, 1.0, 0.2, 0.0,  # PE: peak in D4
      0.0, 0.0, 0.0, 0.1, 0.3, 1.0), # APC: peak in D6
    nrow = 4, ncol = 6, byrow = TRUE,
    dimnames = list(c("BV421", "FITC", "PE", "APC"),
                    paste0("D", 1:6))
  )

  result <- sw_unmix_spillover_matrix(spectra)

  expect_s3_class(result, "sw_ssm")
  expect_true(is.matrix(result$matrix))
  expect_equal(nrow(result$matrix), 4)
  expect_equal(ncol(result$matrix), 4)
  expect_equal(rownames(result$matrix), c("BV421", "FITC", "PE", "APC"))

  # Diagonal should be zero (no self-spreading)
  expect_equal(diag(result$matrix), c(0, 0, 0, 0))

  # All values should be non-negative
  expect_true(all(result$matrix >= 0))

  # All values should be <= 1 (absolute cosine similarity)
  expect_true(all(result$matrix <= 1 + 1e-10))

  # Summary tibble
  expect_s3_class(result$summary, "tbl_df")
  expect_equal(nrow(result$summary), 4)
  expect_true(all(c("fluorophore", "max_spreading", "mean_spreading",
                     "worst_partner") %in% names(result$summary)))
})

test_that("Orthogonal spectra produce zero spreading", {
  # Perfectly orthogonal spectra (identity-like)
  spectra <- diag(4)
  rownames(spectra) <- c("F1", "F2", "F3", "F4")
  colnames(spectra) <- paste0("D", 1:4)

  result <- sw_unmix_spillover_matrix(spectra)

  # All off-diagonal should be 0
  expect_true(all(result$matrix == 0))
  expect_true(all(result$summary$max_spreading == 0))
})

test_that("Identical spectra produce maximum spreading", {
  # Two fluorophores with identical spectra
  spectra <- matrix(c(1, 0.5, 0.2,
                       1, 0.5, 0.2),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(c("F1", "F2"), paste0("D", 1:3)))

  result <- sw_unmix_spillover_matrix(spectra)

  # SSM[1,2] and SSM[2,1] should be ~1
  expect_equal(result$matrix[1, 2], 1.0, tolerance = 1e-10)
  expect_equal(result$matrix[2, 1], 1.0, tolerance = 1e-10)
})

test_that("SSM with unmixed flowFrame incorporates empirical data", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  spectra <- matrix(
    c(1.0, 0.1, 0.0,
      0.1, 1.0, 0.0,
      0.0, 0.1, 1.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("BV421", "FITC", "PE"), paste0("D", 1:3))
  )

  # Create unmixed flowFrame with correlated channels
  n <- 200
  mat <- matrix(
    c(rnorm(n, 1000, 200),
      rnorm(n, 500, 150),
      rnorm(n, 800, 180)),
    ncol = 3,
    dimnames = list(NULL, c("BV421", "FITC", "PE"))
  )
  ff <- flowCore::flowFrame(mat)

  result <- sw_unmix_spillover_matrix(spectra, unmixed_ff = ff)
  expect_s3_class(result, "sw_ssm")
  expect_equal(nrow(result$matrix), 3)
})

# =========================================================================
# sw_plot_spillover_matrix — input validation
# =========================================================================

test_that("sw_plot_spillover_matrix rejects non-sw_ssm input", {
  expect_error(sw_plot_spillover_matrix(list(matrix = diag(3))),
               "sw_ssm")
})

test_that("sw_plot_spillover_matrix runs without error", {
  spectra <- matrix(
    c(1.0, 0.1, 0.0,
      0.1, 1.0, 0.3,
      0.0, 0.3, 1.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("F1", "F2", "F3"), paste0("D", 1:3))
  )
  ssm <- sw_unmix_spillover_matrix(spectra)

  # Plot to temp PDF
  tmp_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp_pdf))

  expect_silent(sw_plot_spillover_matrix(ssm, plot_file = tmp_pdf))
  expect_true(file.exists(tmp_pdf))
})

# =========================================================================
# sw_unmix_quality — input validation
# =========================================================================

test_that("sw_unmix_quality rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_unmix_quality(data.frame(x = 1)), "flowFrame")
})

test_that("sw_unmix_quality rejects invalid cv_threshold", {
  skip_if_not_installed("flowCore")
  mat <- matrix(rnorm(200), ncol = 2,
                dimnames = list(NULL, c("BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_unmix_quality(ff, cv_threshold = -1), "positive number")
})

test_that("sw_unmix_quality rejects missing channels", {
  skip_if_not_installed("flowCore")
  mat <- matrix(rnorm(200), ncol = 2,
                dimnames = list(NULL, c("BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_unmix_quality(ff, channels = "MISSING"), "not found")
})

# =========================================================================
# sw_unmix_quality — integration
# =========================================================================

test_that("sw_unmix_quality computes per-channel metrics", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  # Channel 1: low CV (tight distribution around high mean)
  # Channel 2: high CV (wide distribution around low mean)
  mat <- matrix(
    c(rnorm(300, mean = 10000, sd = 500),     # BV421-A: CV ~ 0.05
      rnorm(300, mean = 100, sd = 200)),       # PE-A: CV ~ 2.0
    ncol = 2,
    dimnames = list(NULL, c("BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)

  result <- sw_unmix_quality(ff, channels = c("BV421-A", "PE-A"),
                                cv_threshold = 0.5)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(c("channel", "mean", "sd", "cv",
                     "pct_negative", "flagged") %in% names(result)))

  # BV421-A should NOT be flagged, PE-A SHOULD be flagged
  bv421_row <- result[result$channel == "BV421-A", ]
  pe_row <- result[result$channel == "PE-A", ]
  expect_false(bv421_row$flagged)
  expect_true(pe_row$flagged)

  # PE-A should have higher pct_negative (negative values in that channel)
  expect_true(pe_row$pct_negative > bv421_row$pct_negative)
})
