# tests/testthat/test-synref_shift.R
# Unit tests for R/synref_shift.R on simulated asinh-scale cells.

make_cells <- function(n = 400, seed = 1) {
  set.seed(seed)
  data.frame(
    batch = rep(c("b1", "b2"), each = n),
    population = rep(rep(c("CD4T", "CD8T"), each = n / 2), 2),
    CD3 = rnorm(2 * n, 3, 1),
    CD4 = rnorm(2 * n, 1, 0.5)
  )
}

one_target <- function(
  size,
  batch = "b2",
  marker = "CD3",
  population = "CD4T"
) {
  data.frame(
    batch = batch,
    marker = marker,
    population = population,
    size = size
  )
}

shift <- function(x, targets, family) {
  sw_synref_shift(x, targets, family = family, cofactor = 6000, seed = 7)
}

test_that("the inverse recovers the input within 1e-8 for offset, gain and stretch", {
  x <- make_cells()
  for (fam in c("offset_asinh", "gain_linear", "stretch_median")) {
    size <- if (fam == "offset_asinh") 0.7 else 1.4
    out <- shift(x, one_target(size), fam)
    back <- SpectraWeaveR:::.sw_synref_invert(
      out[["data"]],
      out[["truth"]],
      6000
    )
    expect_lt(max(abs(back[["CD3"]] - x[["CD3"]])), 1e-8)
  }
})

test_that("two_mode is undone by the cell signs and splits the cells in half", {
  x <- make_cells()
  out <- shift(x, one_target(0.9), "two_mode")
  back <- SpectraWeaveR:::.sw_synref_invert(out[["data"]], out[["truth"]], 6000)
  expect_lt(max(abs(back[["CD3"]] - x[["CD3"]])), 1e-8)
  signs <- attr(out[["truth"]], "cell_sign")[["sign"]]
  expect_equal(sum(signs == 1), 100)
  expect_equal(sum(signs == -1), 100)
})

test_that("truth equals the mean and median change of shifted vs unshifted cells", {
  x <- make_cells()
  idx <- x[["batch"]] == "b2" & x[["population"]] == "CD4T"
  for (fam in c("offset_asinh", "gain_linear", "stretch_median", "two_mode")) {
    out <- shift(x, one_target(1.3), fam)
    expect_equal(
      out[["truth"]][["mean_change"]],
      mean(out[["data"]][["CD3"]][idx]) - mean(x[["CD3"]][idx])
    )
    expect_equal(
      out[["truth"]][["median_change"]],
      median(out[["data"]][["CD3"]][idx]) - median(x[["CD3"]][idx])
    )
  }
  out <- shift(x, one_target(0.7), "offset_asinh")
  expect_equal(out[["truth"]][["mean_change"]], 0.7)
  expect_equal(out[["truth"]][["median_change"]], 0.7)
})

test_that("only the targeted cells and marker change, and truth stays out of data", {
  x <- make_cells()
  out <- shift(x, one_target(0.7), "offset_asinh")
  idx <- x[["batch"]] == "b2" & x[["population"]] == "CD4T"
  expect_identical(out[["data"]][["CD4"]], x[["CD4"]])
  expect_identical(out[["data"]][["CD3"]][!idx], x[["CD3"]][!idx])
  expect_identical(names(out[["data"]]), names(x))
})

test_that("stretch keeps the median and scales the spread", {
  x <- make_cells()
  out <- shift(x, one_target(2), "stretch_median")
  idx <- x[["batch"]] == "b2" & x[["population"]] == "CD4T"
  expect_equal(median(out[["data"]][["CD3"]][idx]), median(x[["CD3"]][idx]))
  expect_equal(sd(out[["data"]][["CD3"]][idx]), 2 * sd(x[["CD3"]][idx]))
})

test_that("a seed reproduces two_mode and the caller's RNG state is restored", {
  x <- make_cells()
  set.seed(99)
  before <- .Random.seed
  a <- shift(x, one_target(0.9), "two_mode")
  expect_identical(.Random.seed, before)
  b <- shift(x, one_target(0.9), "two_mode")
  expect_identical(a[["data"]], b[["data"]])
  c2 <- sw_synref_shift(
    x,
    one_target(0.9),
    "two_mode",
    cofactor = 6000,
    seed = 8
  )
  expect_false(identical(a[["data"]][["CD3"]], c2[["data"]][["CD3"]]))
})

test_that("missing required arguments stop and name them", {
  x <- make_cells()
  expect_error(sw_synref_shift(x, one_target(1)), "family.*cofactor.*seed")
  expect_error(
    sw_synref_shift(),
    "'x', 'targets', 'family', 'cofactor', 'seed'"
  )
})

test_that("bad input stops instead of dropping a shift", {
  x <- make_cells()
  expect_error(
    shift(x, one_target(1, batch = "b9"), "offset_asinh"),
    "matches no cell"
  )
  expect_error(
    shift(x, one_target(1, marker = "CD99"), "offset_asinh"),
    "not in `x`"
  )
  expect_error(
    shift(x, rbind(one_target(1), one_target(2)), "offset_asinh"),
    "repeats"
  )
  expect_error(shift(x, one_target(-1), "gain_linear"), "above 0")
  expect_error(shift(x, one_target(1), "nonsense"), "`family` must be one of")
  expect_error(
    sw_synref_shift(x, one_target(1), "offset_asinh", cofactor = 0, seed = 1),
    "cofactor"
  )
})
