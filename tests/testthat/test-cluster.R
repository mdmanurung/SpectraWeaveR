# tests/testthat/test-cluster.R
# Unit tests for R/cluster.R — Clustering (FlowSOM)

test_that("sw_cluster rejects empty lineage_markers", {
  expect_error(sw_cluster(matrix(1:10, ncol = 2), lineage_markers = character(0)),
               "non-empty character vector")
})

test_that("sw_cluster rejects n_metaclusters < 2", {
  df <- tibble::tibble(CD3 = rnorm(100), CD4 = rnorm(100))
  expect_error(sw_cluster(df, lineage_markers = c("CD3", "CD4"),
                          n_metaclusters = 1),
               "integer >= 2")
})

test_that("sw_cluster rejects invalid input type", {
  expect_error(sw_cluster("not_a_df", lineage_markers = c("CD3")),
               "matrix, data.frame, or tibble")
})

test_that("sw_cluster rejects missing markers in data.frame", {
  df <- tibble::tibble(CD3 = rnorm(50))
  expect_error(sw_cluster(df, lineage_markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_cluster rejects missing markers in matrix", {
  mat <- matrix(rnorm(100), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  expect_error(sw_cluster(mat, lineage_markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_get_cluster_assignments rejects non-FlowSOM", {
  skip_if_not_installed("FlowSOM")
  expect_error(sw_get_cluster_assignments("not_fsom"),
               "FlowSOM object")
})

test_that("sw_cluster_mfis rejects non-FlowSOM", {
  skip_if_not_installed("FlowSOM")
  expect_error(sw_cluster_mfis("not_fsom"),
               "FlowSOM object")
})

test_that("sw_plot_clusters rejects non-FlowSOM", {
  skip_if_not_installed("FlowSOM")
  expect_error(sw_plot_clusters("not_fsom"),
               "FlowSOM object")
})

test_that("sw_plot_clusters rejects invalid plot_file", {
  skip_if_not_installed("FlowSOM")
  expect_error(sw_plot_clusters("not_fsom", plot_file = 42),
               "single file path")
})

test_that("sw_map_new_data rejects non-FlowSOM", {
  skip_if_not_installed("FlowSOM")
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(sw_map_new_data("not_fsom", ff),
               "FlowSOM object")
})

test_that("sw_map_new_data rejects non-flowFrame", {
  skip_if_not_installed("FlowSOM")
  skip_if_not_installed("flowCore")

  expect_error(sw_map_new_data("not_fsom", data.frame(x = 1)),
               "FlowSOM object")
})

# Integration test: full clustering pipeline (requires FlowSOM + flowCore)
test_that("sw_cluster runs end-to-end with FlowSOM", {
  skip_if_not_installed("FlowSOM")
  skip_if_not_installed("flowCore")

  set.seed(123)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, mean = 1), rnorm(100, mean = 5)),
    CD4 = c(rnorm(100, mean = 2), rnorm(100, mean = 6)),
    CD8 = c(rnorm(100, mean = 3), rnorm(100, mean = 1))
  )

  fsom <- sw_cluster(df, lineage_markers = c("CD3", "CD4", "CD8"),
                     xdim = 5, ydim = 5, n_metaclusters = 5, seed = 42)

  expect_true(methods::is(fsom, "FlowSOM"))

  # Test cluster assignments
  assignments <- sw_get_cluster_assignments(fsom)
  expect_equal(length(assignments), 200)
  expect_true(all(assignments >= 1))
  expect_true(all(assignments <= 5))

  # Test MFI table
  mfi <- sw_cluster_mfis(fsom)
  expect_s3_class(mfi, "tbl_df")
  expect_true("metacluster" %in% names(mfi))
  expect_true(all(c("CD3", "CD4", "CD8") %in% names(mfi)))
})
