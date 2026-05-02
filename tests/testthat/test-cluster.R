# tests/testthat/test-cluster.R
# Unit tests for R/cluster.R — Clustering (kohonen SOM / FastPG)

test_that("sw_cluster_run rejects empty lineage_markers", {
  expect_error(sw_cluster_run(matrix(1:10, ncol = 2), lineage_markers = character(0)),
               "non-empty character vector")
})

test_that("sw_cluster_run rejects n_metaclusters < 2 for SOM", {
  df <- tibble::tibble(CD3 = rnorm(100), CD4 = rnorm(100))
  expect_error(sw_cluster_run(df, lineage_markers = c("CD3", "CD4"),
                          method = "som", n_metaclusters = 1),
               "integer >= 2")
})

test_that("sw_cluster_run rejects invalid input type", {
  expect_error(sw_cluster_run("not_a_df", lineage_markers = c("CD3")),
               "matrix, data.frame, or tibble")
})

test_that("sw_cluster_run rejects missing markers in data.frame", {
  df <- tibble::tibble(CD3 = rnorm(50))
  expect_error(sw_cluster_run(df, lineage_markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_cluster_run rejects missing markers in matrix", {
  mat <- matrix(rnorm(100), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  expect_error(sw_cluster_run(mat, lineage_markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_cluster_run rejects invalid method", {
  df <- tibble::tibble(CD3 = rnorm(50))
  expect_error(sw_cluster_run(df, lineage_markers = c("CD3"),
                          method = "invalid"),
               "should be one of")
})

test_that("sw_cluster_assignments rejects non-sw_cluster_result", {
  expect_error(sw_cluster_assignments("not_a_result"),
               "sw_cluster_result")
})

test_that("sw_cluster_mfi rejects non-sw_cluster_result", {
  expect_error(sw_cluster_mfi("not_a_result"),
               "sw_cluster_result")
})

test_that("sw_plot_clusters rejects non-sw_cluster_result", {
  expect_error(sw_plot_clusters("not_a_result"),
               "sw_cluster_result")
})

test_that("sw_plot_clusters rejects invalid plot_file", {
  expect_error(sw_plot_clusters("not_a_result", plot_file = 42),
               "sw_cluster_result")
})

test_that("sw_cluster_predict rejects non-sw_cluster_result", {
  expect_error(sw_cluster_predict("not_a_result", data.frame(x = 1)),
               "sw_cluster_result")
})

# Integration test: full clustering pipeline with kohonen SOM
test_that("sw_cluster_run runs end-to-end with kohonen SOM", {
  skip_if_not_installed("kohonen")

  set.seed(123)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, mean = 1), rnorm(100, mean = 5)),
    CD4 = c(rnorm(100, mean = 2), rnorm(100, mean = 6)),
    CD8 = c(rnorm(100, mean = 3), rnorm(100, mean = 1))
  )

  result <- sw_cluster_run(df, lineage_markers = c("CD3", "CD4", "CD8"),
                       method = "som",
                       xdim = 5, ydim = 5, n_metaclusters = 5, seed = 42)

  expect_s3_class(result, "sw_cluster_result")
  expect_equal(result$method, "som")

  # Test cluster assignments
  assignments <- sw_cluster_assignments(result)
  expect_equal(length(assignments), 200)
  expect_true(all(assignments >= 1))
  expect_true(all(assignments <= 5))

  # Test MFI table
  mfi <- sw_cluster_mfi(result)
  expect_s3_class(mfi, "tbl_df")
  expect_true("cluster" %in% names(mfi))
  expect_true(all(c("CD3", "CD4", "CD8") %in% names(mfi)))
})

# Integration test: predict new data onto trained SOM
test_that("sw_cluster_predict maps new data onto trained SOM", {
  skip_if_not_installed("kohonen")

  set.seed(123)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, mean = 1), rnorm(100, mean = 5)),
    CD4 = c(rnorm(100, mean = 2), rnorm(100, mean = 6)),
    CD8 = c(rnorm(100, mean = 3), rnorm(100, mean = 1))
  )

  result <- sw_cluster_run(df, lineage_markers = c("CD3", "CD4", "CD8"),
                       method = "som",
                       xdim = 5, ydim = 5, n_metaclusters = 5, seed = 42)

  new_data <- tibble::tibble(
    CD3 = rnorm(50, mean = 3),
    CD4 = rnorm(50, mean = 4),
    CD8 = rnorm(50, mean = 2)
  )

  mapped <- sw_cluster_predict(result, new_data)
  expect_s3_class(mapped, "sw_cluster_result")
  expect_equal(length(mapped$assignments), 50)
  expect_true(all(mapped$assignments >= 1))
  expect_true(all(mapped$assignments <= 5))
})

# Integration test: FastPG clustering
test_that("sw_cluster_run runs end-to-end with FastPG", {
  skip_if_not_installed("FastPG")

  set.seed(123)
  df <- tibble::tibble(
    CD3 = c(rnorm(100, mean = 1), rnorm(100, mean = 5)),
    CD4 = c(rnorm(100, mean = 2), rnorm(100, mean = 6)),
    CD8 = c(rnorm(100, mean = 3), rnorm(100, mean = 1))
  )

  result <- sw_cluster_run(df, lineage_markers = c("CD3", "CD4", "CD8"),
                       method = "fastpg", k = 10, seed = 42)

  expect_s3_class(result, "sw_cluster_result")
  expect_equal(result$method, "fastpg")

  assignments <- sw_cluster_assignments(result)
  expect_equal(length(assignments), 200)
  expect_true(all(assignments >= 1))

  mfi <- sw_cluster_mfi(result)
  expect_s3_class(mfi, "tbl_df")
  expect_true("cluster" %in% names(mfi))
})

test_that("sw_cluster_predict rejects fastpg method", {
  skip_if_not_installed("FastPG")

  df <- tibble::tibble(CD3 = rnorm(50), CD4 = rnorm(50))
  result <- sw_cluster_run(df, lineage_markers = c("CD3", "CD4"),
                       method = "fastpg", k = 10)

  expect_error(sw_cluster_predict(result, df),
               "only supported for method = 'som'")
})
