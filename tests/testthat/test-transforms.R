# tests/testthat/test-transforms.R
# Unit tests for R/transforms.R — Scale transformation utilities

# =========================================================================
# sw_estimate_scale_transforms
# =========================================================================

test_that("sw_estimate_scale_transforms rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(
    sw_estimate_scale_transforms("not_a_ff"),
    "flowFrame"
  )
})

test_that("sw_estimate_scale_transforms rejects invalid fluo_method", {
  skip_if_not_installed("flowCore")
  mat <- matrix(rnorm(300), ncol = 3,
                dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_estimate_scale_transforms(ff, fluo_method = "invalid"),
    "arg"
  )
})

test_that("sw_estimate_scale_transforms requires scatter_ref_marker for linearQuantile", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(300)) + 1, ncol = 3,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_estimate_scale_transforms(ff,
                                  fluo_method = "none",
                                  scatter_method = "linearQuantile"),
    "scatter_ref_marker"
  )
})

test_that("sw_estimate_scale_transforms errors if no transforms estimated", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(300)) + 1, ncol = 3,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "Time")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_estimate_scale_transforms(ff, fluo_method = "none",
                                  scatter_method = "none"),
    "No transformations"
  )
})

# =========================================================================
# sw_apply_scale_transforms
# =========================================================================

test_that("sw_apply_scale_transforms rejects non-flowFrame/flowSet", {
  skip_if_not_installed("flowCore")
  expect_error(
    sw_apply_scale_transforms("not_ff", list()),
    "flowFrame or flowSet"
  )
})

test_that("sw_apply_scale_transforms rejects non-transformList", {
  skip_if_not_installed("flowCore")
  mat <- matrix(rnorm(300), ncol = 3,
                dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_apply_scale_transforms(ff, "not_a_translist"),
    "transformList"
  )
})
