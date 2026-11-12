# tests/testthat/test-batch_matrix.R
# Unit tests for R/batch_matrix.R on mock score tables.

# Marker A is planted high in b2 and low in b3; marker B is flat.
make_scores <- function() {
  tibble::tibble(
    check = "rest_emd",
    marker = rep(c("A", "B"), each = 3),
    batch = rep(c("b1", "b2", "b3"), 2),
    value = c(0.1, 0.9, 0.2, 0.3, 0.31, 0.29)
  )
}

test_that("ranks batches inside each check x marker column and finds the planted batch", {
  m <- sw_batch_matrix(
    make_scores(),
    check = "rest_emd",
    direction = c(rest_emd = "high")
  )
  expect_s3_class(m, "sw_batch_matrix")
  expect_named(m, c("check", "batch", "marker", "value", "rank"))
  expect_equal(m$batch[m$marker == "A" & m$rank == 1], "b2")
  expect_equal(sort(m$rank[m$marker == "B"]), 1:3)
  expect_equal(sort(m$rank[m$marker == "A"]), 1:3)
})

test_that("direction low and abs change which batch ranks first", {
  s <- make_scores()
  s$value[1:3] <- c(-0.5, 0.2, 0.4)
  low <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "low"))
  abs_ <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "abs"))
  expect_equal(low$batch[low$marker == "A" & low$rank == 1], "b1")
  expect_equal(abs_$batch[abs_$marker == "A" & abs_$rank == 1], "b1")
  high <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "high"))
  expect_equal(high$batch[high$marker == "A" & high$rank == 1], "b3")
})

test_that("diag_key uses the package convention for gain_batch_ratio", {
  s <- tibble::tibble(
    check = "gain_batch_ratio",
    marker = "A",
    batch = c("b1", "b2", "b3"),
    value = c(1, 2, 0.25)
  )
  m <- sw_batch_matrix(s, "gain_batch_ratio", c(gain_batch_ratio = "diag_key"))
  # abs(log(0.25)) > abs(log(2)): the narrowed batch ranks first.
  expect_equal(m$batch[m$rank == 1], "b3")
})

test_that("wide view holds one column per batch and ranks never mix markers", {
  m <- sw_batch_matrix(make_scores(), "rest_emd", c(rest_emd = "high"))
  w <- attr(m, "wide")
  expect_equal(names(w), c("check", "marker", "b1", "b2", "b3"))
  expect_equal(w$b2[w$marker == "A"], 0.9)
  expect_equal(nrow(w), 2)
})

test_that("rows without a batch are listed in dropped, not ranked", {
  s <- rbind(
    make_scores(),
    tibble::tibble(
      check = "rest_emd",
      marker = "A",
      batch = NA_character_,
      value = 5
    )
  )
  m <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "high"))
  expect_equal(nrow(m), 6)
  expect_equal(nrow(attr(m, "dropped")), 1)
})

test_that("thin columns and NA values get rank NA", {
  s <- tibble::tibble(
    check = "rest_emd",
    marker = c("A", "B", "B"),
    batch = c("b1", "b1", "b2"),
    value = c(1, 1, NA)
  )
  m <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "high"))
  expect_true(all(is.na(m$rank)))
})

test_that("population ranks inside population and flag_iqr form works", {
  s <- make_scores()
  p <- rbind(transform(s, population = "T"), transform(s, population = "B"))
  p$value[p$population == "B"] <- rev(p$value[p$population == "B"])
  m <- sw_batch_matrix(p, "rest_emd", c(rest_emd = "high"))
  expect_true("population" %in% names(m))
  expect_equal(sort(m$rank[m$marker == "A" & m$population == "T"]), 1:3)
  f <- transform(s, metric = check, control = "c1")
  f$check <- NULL
  mf <- sw_batch_matrix(f, "rest_emd", c(rest_emd = "high"))
  expect_equal(nrow(mf), 6)
})

test_that("errors on missing arguments, bad direction and repeated rows", {
  s <- make_scores()
  expect_error(sw_batch_matrix(s), "Required argument")
  expect_error(sw_batch_matrix(s, check = "rest_emd"), "'direction'")
  expect_error(
    sw_batch_matrix(s, "rest_emd", c(rest_emd = "up")),
    "one of"
  )
  expect_error(
    sw_batch_matrix(s, "rest_emd", c(other = "high")),
    "named character vector"
  )
  expect_error(
    sw_batch_matrix(s, "nope", c(nope = "high")),
    "not found"
  )
  expect_error(
    sw_batch_matrix(rbind(s, s), "rest_emd", c(rest_emd = "high")),
    "once"
  )
})

test_that("matrix is deterministic and leaves the RNG state alone", {
  set.seed(7)
  before <- .Random.seed
  a <- sw_batch_matrix(make_scores(), "rest_emd", c(rest_emd = "high"))
  expect_identical(.Random.seed, before)
  b <- sw_batch_matrix(make_scores(), "rest_emd", c(rest_emd = "high"))
  expect_identical(a, b)
})
