# tests/testthat/test-dimred.R
# Unit and integration tests for R/dimred.R — Dimensionality reduction

# =========================================================================
# sw_run_dimred — input validation
# =========================================================================

test_that("sw_run_dimred rejects non-data.frame input", {
  expect_error(sw_run_dimred("not_df", markers = "CD3"),
               "data.frame")
})

test_that("sw_run_dimred rejects empty markers", {
  df <- tibble::tibble(CD3 = 1:10)
  expect_error(sw_run_dimred(df, markers = character(0)),
               "non-empty character vector")
})

test_that("sw_run_dimred rejects missing markers", {
  df <- tibble::tibble(CD3 = 1:10)
  expect_error(sw_run_dimred(df, markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_run_dimred rejects invalid method", {
  df <- tibble::tibble(CD3 = 1:10, CD4 = 1:10)
  expect_error(sw_run_dimred(df, markers = c("CD3", "CD4"),
                             method = "invalid"),
               "arg")
})

test_that("sw_run_dimred rejects invalid n_dims", {
  df <- tibble::tibble(CD3 = 1:10, CD4 = 1:10)
  expect_error(sw_run_dimred(df, markers = c("CD3", "CD4"), n_dims = 0),
               "positive integer")
})

# =========================================================================
# sw_run_dimred — PCA integration
# =========================================================================

test_that("sw_run_dimred PCA runs on synthetic data", {
  set.seed(42)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, 1), rnorm(100, 5)),
    CD4 = c(rnorm(100, 2), rnorm(100, 6)),
    CD8 = c(rnorm(100, 3), rnorm(100, 1)),
    batch = rep(c("B1", "B2"), each = 100)
  )

  result <- sw_run_dimred(df, markers = c("CD3", "CD4", "CD8"),
                          method = "pca", seed = 42)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 200)
  expect_true(all(c("dim1", "dim2") %in% names(result)))
  expect_true("batch" %in% names(result))
  expect_equal(attr(result, "method"), "pca")
  expect_true(!is.null(attr(result, "variance_explained")))
  expect_equal(length(attr(result, "variance_explained")), 2)
})

test_that("sw_run_dimred PCA subsamples large data", {
  set.seed(42)
  df <- tibble::tibble(
    CD3 = rnorm(500),
    CD4 = rnorm(500)
  )

  result <- sw_run_dimred(df, markers = c("CD3", "CD4"),
                          method = "pca", max_cells = 100, seed = 42)
  expect_equal(nrow(result), 100)
})

# =========================================================================
# sw_run_dimred — UMAP integration
# =========================================================================

test_that("sw_run_dimred UMAP runs on synthetic data", {
  skip_if_not_installed("uwot")

  set.seed(42)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, 1), rnorm(100, 5)),
    CD4 = c(rnorm(100, 2), rnorm(100, 6)),
    CD8 = c(rnorm(100, 3), rnorm(100, 1)),
    cluster = rep(c(1L, 2L), each = 100)
  )

  result <- sw_run_dimred(df, markers = c("CD3", "CD4", "CD8"),
                          method = "umap", seed = 42)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 200)
  expect_true(all(c("dim1", "dim2") %in% names(result)))
  expect_true("cluster" %in% names(result))
  expect_equal(attr(result, "method"), "umap")
})

# =========================================================================
# sw_plot_dimred — input validation
# =========================================================================

test_that("sw_plot_dimred rejects non-data.frame", {
  expect_error(sw_plot_dimred("not_df"), "data.frame")
})

test_that("sw_plot_dimred rejects missing dim columns", {
  df <- tibble::tibble(x = 1:10, y = 1:10)
  expect_error(sw_plot_dimred(df), "dim1.*dim2")
})

test_that("sw_plot_dimred rejects missing color_by column", {
  df <- tibble::tibble(dim1 = 1:10, dim2 = 1:10)
  expect_error(sw_plot_dimred(df, color_by = "MISSING"), "not found")
})

# =========================================================================
# sw_plot_dimred — integration
# =========================================================================

test_that("sw_plot_dimred returns ggplot object", {
  skip_if_not_installed("ggplot2")

  df <- tibble::tibble(
    dim1 = rnorm(50),
    dim2 = rnorm(50),
    batch = rep(c("A", "B"), 25),
    CD3 = rnorm(50)
  )
  attr(df, "method") <- "umap"

  # Without colour
  p1 <- sw_plot_dimred(df)
  expect_s3_class(p1, "ggplot")

  # With discrete colour
  p2 <- sw_plot_dimred(df, color_by = "batch")
  expect_s3_class(p2, "ggplot")

  # With continuous colour
  p3 <- sw_plot_dimred(df, color_by = "CD3")
  expect_s3_class(p3, "ggplot")
})
