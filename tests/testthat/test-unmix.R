# tests/testthat/test-unmix.R
# Unit tests for R/unmix.R — Spectral unmixing wrappers
# Input validation tests run without AutoSpectral installed.
# Functional tests are gated by skip_if_not_installed("AutoSpectral").

# ============================================================================
# sw_load_unmixed
# ============================================================================

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

# ============================================================================
# sw_remove_margins
# ============================================================================

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

# ============================================================================
# sw_autospectral_setup
# ============================================================================

test_that("sw_autospectral_setup rejects non-character control_dir", {
  skip_if_not_installed("AutoSpectral")
  expect_error(sw_autospectral_setup(42),
               "single directory path")
})

test_that("sw_autospectral_setup rejects non-existent directory", {
  skip_if_not_installed("AutoSpectral")
  expect_error(sw_autospectral_setup("/no/such/dir"),
               "does not exist")
})

test_that("sw_autospectral_setup rejects invalid cytometer", {
  skip_if_not_installed("AutoSpectral")
  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "cytometer_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  expect_error(sw_autospectral_setup(test_dir, cytometer = "invalid"),
               "must be one of")
})

test_that("sw_autospectral_setup rejects non-character cytometer", {
  skip_if_not_installed("AutoSpectral")
  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "cytometer_type_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  expect_error(sw_autospectral_setup(test_dir, cytometer = 123),
               "single character string")
})

test_that("sw_autospectral_setup rejects non-existent control_file", {
  skip_if_not_installed("AutoSpectral")
  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "cf_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  expect_error(
    sw_autospectral_setup(test_dir, control_file = "/no/file.csv"),
    "does not exist"
  )
})

test_that("sw_autospectral_setup rejects non-character output_dir", {
  skip_if_not_installed("AutoSpectral")
  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "outdir_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  expect_error(
    sw_autospectral_setup(test_dir, output_dir = 42),
    "single directory path"
  )
})

test_that("sw_autospectral_setup requires AutoSpectral", {
  skip_if(requireNamespace("AutoSpectral", quietly = TRUE),
          "AutoSpectral is installed; cannot test missing package")

  tmp_dir <- tempdir()
  test_dir <- file.path(tmp_dir, "asp_req_test")
  dir.create(test_dir, showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE))

  expect_error(sw_autospectral_setup(test_dir),
               "AutoSpectral")
})

# ============================================================================
# sw_prepare_controls
# ============================================================================

test_that("sw_prepare_controls rejects invalid setup object", {
  skip_if_not_installed("AutoSpectral")
  expect_error(sw_prepare_controls(list(asp = NULL)),
               "sw_setup object")
})

test_that("sw_prepare_controls rejects invalid clean argument", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(sw_prepare_controls(mock_setup, clean = "yes"),
               "TRUE or FALSE")
})

test_that("sw_prepare_controls rejects invalid af_remove argument", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(sw_prepare_controls(mock_setup, af_remove = "yes"),
               "TRUE or FALSE")
})

test_that("sw_prepare_controls rejects invalid parallel argument", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(sw_prepare_controls(mock_setup, parallel = "yes"),
               "TRUE or FALSE")
})

# ============================================================================
# sw_extract_af_spectra
# ============================================================================

test_that("sw_extract_af_spectra rejects invalid setup", {
  skip_if_not_installed("AutoSpectral")
  expect_error(
    sw_extract_af_spectra("file.fcs", setup = list(), spectra = matrix()),
    "sw_setup object"
  )
})

test_that("sw_extract_af_spectra rejects non-matrix spectra", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra("file.fcs", setup = mock_setup,
                          spectra = data.frame(x = 1)),
    "matrix"
  )
})

test_that("sw_extract_af_spectra rejects non-existent file", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(figures = FALSE)),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra("/no/file.fcs", setup = mock_setup,
                          spectra = matrix(1:4, ncol = 2)),
    "does not exist"
  )
})

test_that("sw_extract_af_spectra rejects invalid refine", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list()),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra("file.fcs", setup = mock_setup,
                          spectra = matrix(1:4, ncol = 2),
                          refine = "yes"),
    "TRUE or FALSE"
  )
})

test_that("sw_extract_af_spectra rejects invalid som_dim", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list()),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra("file.fcs", setup = mock_setup,
                          spectra = matrix(1:4, ncol = 2),
                          som_dim = 1),
    "single integer >= 2"
  )
})

test_that("sw_extract_af_spectra rejects unnamed list", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(figures = FALSE)),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra(list("a.fcs", "b.fcs"),
                          setup = mock_setup,
                          spectra = matrix(1:4, ncol = 2)),
    "non-empty names"
  )
})

test_that("sw_extract_af_spectra rejects invalid input type", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list()),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_af_spectra(42, setup = mock_setup,
                          spectra = matrix(1:4, ncol = 2)),
    "single file path or a named list"
  )
})

test_that("sw_extract_af_spectra requires AutoSpectral", {
  skip_if(requireNamespace("AutoSpectral", quietly = TRUE),
          "AutoSpectral is installed; cannot test missing package")

  expect_error(
    sw_extract_af_spectra("f.fcs", list(), matrix()),
    "AutoSpectral"
  )
})

# ============================================================================
# sw_extract_spectral_variants
# ============================================================================

test_that("sw_extract_spectral_variants rejects invalid setup", {
  skip_if_not_installed("AutoSpectral")
  expect_error(
    sw_extract_spectral_variants(
      setup = list(),
      spectra = matrix(1:4, ncol = 2),
      af_spectra = matrix(1:4, ncol = 2)
    ),
    "sw_setup object"
  )
})

test_that("sw_extract_spectral_variants rejects non-matrix spectra", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_spectral_variants(
      setup = mock_setup,
      spectra = data.frame(x = 1),
      af_spectra = matrix(1:4, ncol = 2)
    ),
    "matrix"
  )
})

test_that("sw_extract_spectral_variants rejects single-row af_spectra", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_spectral_variants(
      setup = mock_setup,
      spectra = matrix(1:4, ncol = 2),
      af_spectra = matrix(1:2, ncol = 2)
    ),
    "at least 2 rows"
  )
})

test_that("sw_extract_spectral_variants rejects invalid refine", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_spectral_variants(
      setup = mock_setup,
      spectra = matrix(1:4, ncol = 2),
      af_spectra = matrix(1:4, nrow = 2),
      refine = "yes"
    ),
    "TRUE or FALSE"
  )
})

test_that("sw_extract_spectral_variants rejects invalid som_dim", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_spectral_variants(
      setup = mock_setup,
      spectra = matrix(1:4, ncol = 2),
      af_spectra = matrix(1:4, nrow = 2),
      som_dim = 0
    ),
    "single integer >= 2"
  )
})

test_that("sw_extract_spectral_variants rejects invalid n_cells", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_extract_spectral_variants(
      setup = mock_setup,
      spectra = matrix(1:4, ncol = 2),
      af_spectra = matrix(1:4, nrow = 2),
      n_cells = 0
    ),
    "positive integer"
  )
})

# ============================================================================
# sw_unmix
# ============================================================================

test_that("sw_unmix rejects empty input", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix(character(0), matrix(1:4, ncol = 2), mock_setup, list()),
    "non-empty character vector"
  )
})

test_that("sw_unmix rejects non-character input", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix(42, matrix(1:4, ncol = 2), mock_setup, list()),
    "non-empty character vector"
  )
})

test_that("sw_unmix rejects non-matrix spectra", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix("file.fcs", data.frame(x = 1), mock_setup, list()),
    "matrix"
  )
})

test_that("sw_unmix rejects invalid setup", {
  skip_if_not_installed("AutoSpectral")
  expect_error(
    sw_unmix("file.fcs", matrix(1:4, ncol = 2), list(), list()),
    "sw_setup object"
  )
})

test_that("sw_unmix rejects invalid flow_control", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix("file.fcs", matrix(1:4, ncol = 2), mock_setup, "not_a_list"),
    "flow.control list"
  )
})

test_that("sw_unmix requires af_spectra for AutoSpectral method", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(), control_file = "f.csv", control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix("file.fcs", matrix(1:4, ncol = 2), mock_setup, list(),
             method = "AutoSpectral", af_spectra = NULL),
    "af_spectra.*required"
  )
})

test_that("sw_unmix rejects non-existent files", {
  skip_if_not_installed("AutoSpectral")
  mock_setup <- structure(
    list(asp = list(unmixed.fcs.dir = tempdir()), control_file = "f.csv",
         control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix("/no/such/file.fcs", matrix(1:4, ncol = 2),
             mock_setup, list(), method = "OLS"),
    "not found"
  )
})

test_that("sw_unmix rejects empty directory", {
  skip_if_not_installed("AutoSpectral")
  tmp_dir <- tempdir()
  empty_dir <- file.path(tmp_dir, "empty_unmix_test")
  dir.create(empty_dir, showWarnings = FALSE)
  on.exit(unlink(empty_dir, recursive = TRUE))

  mock_setup <- structure(
    list(asp = list(unmixed.fcs.dir = tempdir()), control_file = "f.csv",
         control_dir = "."),
    class = "sw_setup"
  )
  expect_error(
    sw_unmix(empty_dir, matrix(1:4, ncol = 2),
             mock_setup, list(), method = "OLS"),
    "No FCS files"
  )
})

# ============================================================================
# sw_unmix_pipeline
# ============================================================================

test_that("sw_unmix_pipeline rejects empty sample_input", {
  expect_error(
    sw_unmix_pipeline(
      control_dir = ".",
      sample_input = character(0)
    ),
    "non-empty character vector"
  )
})

test_that("sw_unmix_pipeline requires unstained for AutoSpectral", {
  expect_error(
    sw_unmix_pipeline(
      control_dir = ".",
      sample_input = "sample.fcs",
      unstained_fcs = NULL,
      method = "AutoSpectral"
    ),
    "unstained_fcs.*required"
  )
})

test_that("sw_unmix_pipeline requires unstained for Automatic method", {
  expect_error(
    sw_unmix_pipeline(
      control_dir = ".",
      sample_input = "sample.fcs",
      unstained_fcs = NULL,
      method = "Automatic"
    ),
    "unstained_fcs.*required"
  )
})

# ============================================================================
# Helper function tests
# ============================================================================

test_that(".valid_cytometers returns expected set", {
  valid <- SpectraWeaveR:::.valid_cytometers()
  expect_true("aurora" %in% valid)
  expect_true("id7000" %in% valid)
  expect_true("s8" %in% valid)
  expect_true("xenith" %in% valid)
  expect_equal(length(valid), 9)
})
