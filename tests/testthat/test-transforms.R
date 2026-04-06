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

# =========================================================================
# Integration tests — actual transformation execution
# =========================================================================

test_that("sw_estimate_scale_transforms estimates logicle for fluor channels", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  # Fluorochrome channels need data with range spanning neg/pos for logicle
  mat <- matrix(
    c(abs(rnorm(500, 100000, 30000)),   # FSC-A (scatter)
      abs(rnorm(500, 60000, 20000)),    # SSC-A (scatter)
      rnorm(500, 500, 300),             # BV421-A (fluor, some negative)
      rnorm(500, 1000, 400)),           # PE-A (fluor, some negative)
    ncol = 4,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)

  trans <- sw_estimate_scale_transforms(ff,
                                         fluo_method = "estimateLogicle",
                                         scatter_method = "none")
  expect_true(methods::is(trans, "transformList"))
})

test_that("sw_apply_scale_transforms transforms values correctly", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  mat <- matrix(
    c(abs(rnorm(500, 100000, 30000)),
      abs(rnorm(500, 60000, 20000)),
      rnorm(500, 500, 300),
      rnorm(500, 1000, 400)),
    ncol = 4,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)

  trans <- sw_estimate_scale_transforms(ff,
                                         fluo_method = "estimateLogicle",
                                         scatter_method = "none")
  ff_t <- sw_apply_scale_transforms(ff, trans)

  expect_true(methods::is(ff_t, "flowFrame"))
  expect_equal(nrow(ff_t), nrow(ff))
  expect_equal(ncol(ff_t), ncol(ff))
  # Transformed values should differ from original
  expect_false(identical(flowCore::exprs(ff), flowCore::exprs(ff_t)))
})
