# tests/testthat/test-utils.R
# Unit tests for R/utils.R — Format conversion utilities

# ---- Helper: create a mock flowFrame for testing ----
# Since flowCore may not be available, tests are guarded with skip_if_not_installed

test_that("sw_read_fcs rejects empty file vector", {
  expect_error(sw_read_fcs(character(0)),
               "non-empty character vector")
})

test_that("sw_read_fcs rejects non-character input", {
  expect_error(sw_read_fcs(42),
               "non-empty character vector")
})

test_that("sw_read_fcs rejects missing files", {
  expect_error(sw_read_fcs("/nonexistent/file.fcs"),
               "not found")
})

test_that("sw_flowframe_to_tibble converts flowFrame to tibble", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(300), ncol = 3,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A")))
  ff <- flowCore::flowFrame(mat)
  result <- sw_flowframe_to_tibble(ff, sample = "S1", batch = "B1",
                                   condition = "Stim")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 100)
  expect_true(all(c("FSC-A", "SSC-A", "BV421-A", "sample", "batch",
                     "condition") %in% names(result)))
  expect_equal(result$sample[1], "S1")
  expect_equal(result$batch[1], "B1")
  expect_equal(result$condition[1], "Stim")
})

test_that("sw_flowframe_to_tibble works without metadata", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(200), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "BV421-A")))
  ff <- flowCore::flowFrame(mat)
  result <- sw_flowframe_to_tibble(ff)

  expect_s3_class(result, "tbl_df")
  expect_false("sample" %in% names(result))
  expect_false("batch" %in% names(result))
})

test_that("sw_flowframe_to_tibble rejects non-flowFrame", {
  skip_if_not_installed("flowCore")

  expect_error(sw_flowframe_to_tibble(data.frame(x = 1)),
               "flowFrame")
})

test_that("sw_tibble_to_flowframe converts tibble to flowFrame", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(50), CD4 = rnorm(50), sample = "S1")
  ff <- sw_tibble_to_flowframe(df, markers = c("CD3", "CD4"))

  expect_true(methods::is(ff, "flowFrame"))
  expect_equal(nrow(ff), 50)
  expect_equal(ncol(ff), 2)
  expect_true(all(c("CD3", "CD4") %in% flowCore::colnames(ff)))
})

test_that("sw_tibble_to_flowframe uses all numeric cols when markers=NULL", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(30), CD4 = rnorm(30))
  ff <- sw_tibble_to_flowframe(df)

  expect_true(methods::is(ff, "flowFrame"))
  expect_equal(ncol(ff), 2)
})

test_that("sw_tibble_to_flowframe rejects missing markers", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_tibble_to_flowframe(df, markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_tibble_to_flowframe rejects non-data.frame", {
  skip_if_not_installed("flowCore")

  expect_error(sw_tibble_to_flowframe("not_a_df"),
               "data.frame")
})

test_that("sw_exprs_to_tibble converts matrix to tibble", {
  mat <- matrix(1:12, ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  result <- sw_exprs_to_tibble(mat, sample = "S1", batch = 1)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 4)
  expect_true(all(c("A", "B", "C", "sample", "batch") %in% names(result)))
  expect_equal(result$sample[1], "S1")
})

test_that("sw_exprs_to_tibble renames columns", {
  mat <- matrix(1:6, ncol = 2)
  result <- sw_exprs_to_tibble(mat, colnames = c("X", "Y"))

  expect_true(all(c("X", "Y") %in% names(result)))
})

test_that("sw_exprs_to_tibble rejects mismatched colnames", {
  mat <- matrix(1:6, ncol = 2)
  expect_error(sw_exprs_to_tibble(mat, colnames = c("X", "Y", "Z")),
               "ncol")
})

test_that("sw_exprs_to_tibble rejects non-matrix input", {
  expect_error(sw_exprs_to_tibble("not_a_matrix"),
               "matrix or data.frame")
})

test_that("sw_get_fluor_channels excludes scatter and time", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(500), ncol = 5,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "Time",
                                        "BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)
  fluor <- sw_get_fluor_channels(ff)

  expect_equal(sort(fluor), sort(c("BV421-A", "PE-A")))
  expect_false("FSC-A" %in% fluor)
  expect_false("SSC-A" %in% fluor)
  expect_false("Time" %in% fluor)
})

test_that("sw_get_fluor_channels rejects invalid input", {
  skip_if_not_installed("flowCore")

  expect_error(sw_get_fluor_channels("not_a_ff"),
               "flowFrame or flowSet")
})

test_that("sw_set_marker_names renames channels", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(200), ncol = 2,
                dimnames = list(NULL, c("BV421-A", "PE-A")))
  ff <- flowCore::flowFrame(mat)
  marker_map <- c("BV421-A" = "CD3", "PE-A" = "CD4")
  ff2 <- sw_set_marker_names(ff, marker_map)

  pdata <- flowCore::pData(flowCore::parameters(ff2))
  expect_equal(pdata$desc[pdata$name == "BV421-A"], "CD3")
  expect_equal(pdata$desc[pdata$name == "PE-A"], "CD4")
})

test_that("sw_set_marker_names warns on missing channel", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 1,
                dimnames = list(NULL, c("BV421-A")))
  ff <- flowCore::flowFrame(mat)
  marker_map <- c("NONEXISTENT" = "CD3")

  expect_warning(sw_set_marker_names(ff, marker_map),
                 "not found")
})

test_that("sw_set_marker_names rejects unnamed vector", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 1,
                dimnames = list(NULL, c("BV421-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_set_marker_names(ff, c("CD3")),
               "named character vector")
})

# Round-trip test: flowFrame → tibble → flowFrame
test_that("round-trip flowFrame -> tibble -> flowFrame preserves data", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(300), ncol = 3,
                dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  ff_orig <- flowCore::flowFrame(mat)

  tbl <- sw_flowframe_to_tibble(ff_orig)
  ff_back <- sw_tibble_to_flowframe(tbl, markers = c("CD3", "CD4", "CD8"))

  expect_equal(flowCore::exprs(ff_back), flowCore::exprs(ff_orig),
               tolerance = 1e-10)
})
