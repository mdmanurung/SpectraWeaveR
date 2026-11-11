# Tests for R/unmix_screen.R (negative shape and marker screen).

# ---- Mock data ---------------------------------------------------------------

# Paired setups: the same Time key, TNF negatives drawn per setup.
mock_shape_input <- function(sds, n = 400, donors = c("d1", "d2", "d3")) {
  set.seed(7)
  out <- lapply(sds, \(sd) {
    stats::setNames(
      lapply(donors, \(d) {
        cbind(Time = seq_len(n), TNF = stats::rnorm(n, sd = sd))
      }),
      donors
    )
  })
  out
}

shape_args <- function(x, min_negatives = 50) {
  donors <- names(x[[1]])
  list(
    x = x,
    reference = names(x)[1],
    populations = stats::setNames(
      lapply(donors, \(d) list(T = seq_len(nrow(x[[1]][[d]])))),
      donors
    ),
    channels = "TNF",
    cofactor = 6000,
    cut = 1,
    squeeze_frac = 0.25,
    min_negatives = min_negatives,
    key_cols = "Time"
  )
}

# ---- Negative shape ----------------------------------------------------------

test_that("every impactful argument is required", {
  expect_error(
    sw_unmix_negative_shape(x = list()),
    "reference.*populations.*channels.*cofactor.*cut.*squeeze_frac.*min_negatives.*key_cols"
  )
})

test_that("rSD is the Nguyen 84th minus 50th percentile of below-cut events", {
  x <- mock_shape_input(c(ref = 300, alt = 300))
  shape <- do.call(sw_unmix_negative_shape, shape_args(x))
  row <- shape$metrics[
    shape$metrics$setup == "ref" & shape$metrics$sample == "d1",
  ]
  v <- x$ref$d1[, "TNF"]
  v <- v[asinh(v / 6000) < 1]
  expect_equal(row$rsd_neg, unname(diff(stats::quantile(v, c(0.5, 0.84)))))
  expect_equal(row$n_neg, length(v))
})

test_that("the reference against itself has ratio 1", {
  x <- mock_shape_input(c(ref = 300, alt = 150))
  shape <- do.call(sw_unmix_negative_shape, shape_args(x))
  ref_rows <- shape$metrics[shape$metrics$setup == "ref", ]
  expect_equal(ref_rows$rsd_ratio, rep(1, nrow(ref_rows)))
})

test_that("a point mass at zero is fully squeezed", {
  x <- mock_shape_input(c(ref = 300, alt = 300))
  for (d in names(x$alt)) {
    x$alt[[d]][, "TNF"] <- 0
  }
  shape <- do.call(sw_unmix_negative_shape, shape_args(x, min_negatives = 1))
  alt <- shape$metrics[shape$metrics$setup == "alt", ]
  expect_equal(alt$squeezed_share, rep(1, nrow(alt)))
})

test_that("too few negatives in the setup or the reference gives NA, never 0", {
  x <- mock_shape_input(c(ref = 300, alt = 300), n = 40)
  shape <- do.call(sw_unmix_negative_shape, shape_args(x, min_negatives = 100))
  expect_true(all(is.na(shape$metrics$rsd_neg)))
  expect_true(all(is.na(shape$metrics$squeezed_share)))
  # Reference too thin while the setup is fine: ratio and squeeze still NA.
  x2 <- mock_shape_input(c(ref = 300, alt = 300))
  x2$ref$d1[, "TNF"] <- 1e7 # every reference event above the cut
  shape2 <- do.call(sw_unmix_negative_shape, shape_args(x2))
  alt_d1 <- shape2$metrics[
    shape2$metrics$setup == "alt" & shape2$metrics$sample == "d1",
  ]
  expect_false(is.na(alt_d1$rsd_neg))
  expect_true(is.na(alt_d1$rsd_ratio))
  expect_true(is.na(alt_d1$squeezed_share))
})

test_that("unpaired rows stop with the setup named", {
  x <- mock_shape_input(c(ref = 300, alt = 300))
  x$alt$d2[, "Time"] <- rev(x$alt$d2[, "Time"])
  expect_error(
    do.call(sw_unmix_negative_shape, shape_args(x)),
    "Sample 'd2': setup 'alt' rows do not pair"
  )
})

# ---- Marker screen -----------------------------------------------------------

# Screen input with a chosen per-donor squeezed share for the setup and the reference.
mock_screen_shape <- function(d_squeezed, ref_squeezed = 0.1) {
  n <- length(d_squeezed)
  metrics <- tibble::tibble(
    setup = rep(c("ref", "alt"), each = n),
    sample = rep(paste0("d", seq_len(n)), 2),
    pop = "CD8T",
    channel = "TNF",
    n_pop = 100L,
    n_neg = 100L,
    rsd_neg = 1,
    ref_rsd = 1,
    rsd_ratio = c(rep(1, n), rep(0.8, n)),
    squeezed_share = c(rep(ref_squeezed, n), ref_squeezed + d_squeezed),
    pos_share = 0.3
  )
  structure(
    list(metrics = metrics, reference = "ref"),
    class = "sw_negative_shape"
  )
}

screen_flag <- function(d, min_delta = 0.1, min_donors = 5) {
  sw_unmix_marker_screen(
    mock_screen_shape(d),
    min_delta = min_delta,
    min_donors = min_donors
  )$flags$flag
}

test_that("a rise in every donor above min_delta flags", {
  expect_true(screen_flag(c(0.2, 0.3, 0.15, 0.4, 0.25)))
})

test_that("a rise in only 4 of 5 donors does not flag", {
  expect_false(screen_flag(c(0.2, 0.3, -0.01, 0.4, 0.25)))
})

test_that("a fall in every donor does not flag", {
  expect_false(screen_flag(-c(0.2, 0.3, 0.15, 0.4, 0.25)))
})

test_that("a rise in every donor below min_delta does not flag", {
  expect_false(screen_flag(c(0.02, 0.03, 0.05, 0.04, 0.01)))
})

test_that("too few donors gives NA, not FALSE", {
  expect_true(is.na(screen_flag(c(0.2, 0.3, 0.15), min_donors = 5)))
  d <- c(0.2, NA, 0.3, NA, 0.4)
  expect_true(is.na(screen_flag(d, min_donors = 5)))
})

test_that("screen arguments are required and channel_populations restricts rows", {
  shape <- mock_screen_shape(c(0.2, 0.3, 0.15, 0.4, 0.25))
  expect_error(sw_unmix_marker_screen(shape), "min_delta.*min_donors")
  out <- sw_unmix_marker_screen(
    shape,
    min_delta = 0.1,
    min_donors = 5,
    channel_populations = list(TNF = "CD4T")
  )
  expect_equal(nrow(out$flags), 0L)
})
