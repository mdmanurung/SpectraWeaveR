# tests/testthat/test-gate.R
# Unit tests for R/gate.R — Automated gating (openCyto)

test_that("sw_build_gating_template creates valid CSV", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  result <- sw_build_gating_template(tmp_file)

  expect_true(file.exists(tmp_file))
  expect_equal(result, tmp_file)

  # Read and validate structure
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 3)

  # Required columns for openCyto gating template
  required_cols <- c("alias", "pop", "parent", "dims", "gating_method")
  expect_true(all(required_cols %in% names(gt)))

  # Check aliases
  expect_true("nonDebris" %in% gt$alias)
  expect_true("singlets" %in% gt$alias)
  expect_true("lymphocytes" %in% gt$alias)
})

test_that("sw_build_gating_template creates parent directory", {
  tmp_dir <- file.path(tempdir(), "nested", "gating_test")
  tmp_file <- file.path(tmp_dir, "template.csv")
  on.exit(unlink(file.path(tempdir(), "nested"), recursive = TRUE))

  sw_build_gating_template(tmp_file)
  expect_true(file.exists(tmp_file))
})

test_that("sw_build_gating_template rejects invalid template_type", {
  expect_error(sw_build_gating_template(tempfile(), template_type = "invalid"),
               "Unknown template_type")
})

test_that("sw_build_gating_template rejects non-character output_file", {
  expect_error(sw_build_gating_template(42),
               "single file path")
})

test_that("sw_gate rejects empty fcs_files", {
  expect_error(sw_gate(character(0), "template.csv"),
               "non-empty character vector")
})

test_that("sw_gate rejects missing template file", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  # Create a dummy FCS file path (won't actually be read due to template check)
  expect_error(sw_gate("dummy.fcs", "/nonexistent/template.csv"),
               "not found")
})

test_that("sw_extract_gated rejects non-GatingSet input", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_extract_gated("not_a_gs", "/singlets"),
               "GatingSet")
})

test_that("sw_extract_gated rejects non-character node", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_extract_gated("not_a_gs", 42),
               "single character string")
})
