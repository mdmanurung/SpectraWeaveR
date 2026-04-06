# tests/testthat/test-unmix.R
# Unit tests for R/unmix.R — Spectral unmixing wrappers

test_that("sw_load_unmixed rejects non-existent directory", {
  expect_error(sw_load_unmixed("/no/such/dir"),
               "does not exist")
})

test_that("sw_load_unmixed rejects non-character input", {
  expect_error(sw_load_unmixed(42),
               "single directory path")
})

test_that("sw_load_unmixed detects empty directory", {
  tmp_dir <- tempdir()
  empty_dir <- file.path(tmp_dir, "empty_fcs_test")
  dir.create(empty_dir, showWarnings = FALSE)
  on.exit(unlink(empty_dir, recursive = TRUE))

  expect_error(sw_load_unmixed(empty_dir),
               "No FCS files")
})

test_that("sw_load_unmixed respects pattern argument", {
  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "pattern_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create non-matching files only

  writeLines("dummy", file.path(test_dir, "file.csv"))
  expect_error(sw_load_unmixed(test_dir, pattern = "\\.fcs$"),
               "No FCS files")
})

test_that("sw_unmix_autospectral rejects non-existent directory", {
  expect_error(sw_unmix_autospectral("/no/such/dir"),
               "does not exist")
})

test_that("sw_unmix_autospectral rejects non-character input", {
  expect_error(sw_unmix_autospectral(42),
               "single directory path")
})

test_that("sw_remove_margins rejects non-flowFrame input", {
  skip_if_not_installed("flowCore")

  expect_error(sw_remove_margins(data.frame(x = 1)),
               "flowFrame")
})

test_that("sw_remove_margins requires PeacoQC", {
  skip_if_not_installed("flowCore")
  skip_if(requireNamespace("PeacoQC", quietly = TRUE),
          message = "PeacoQC is installed; cannot test missing package")

  mat <- matrix(rnorm(100), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)

  expect_error(sw_remove_margins(ff),
               "PeacoQC")
})
