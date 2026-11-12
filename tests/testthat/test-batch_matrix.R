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

# Real output shapes from the package's own score functions on mock files.
make_files <- function() {
  set.seed(3)
  files <- lapply(1:9, function(i) {
    matrix(
      stats::rnorm(600) + (i > 6) * c(rep(2, 300), rep(0, 300)),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    )
  })
  names(files) <- paste0("f", 1:9)
  files
}

test_that("NA-population rows from sw_batch_scores(population=) are ranked", {
  files <- make_files()
  lab <- lapply(files, \(x) rep(c("T", "M"), length.out = nrow(x)))
  sc <- sw_batch_scores(
    files,
    batch = stats::setNames(rep(c("b1", "b2", "b3"), each = 3), names(files)),
    markers = c("A", "B"),
    n_events = 100,
    seed = 1,
    order = NULL,
    reference_group = NULL,
    population = list(labels = lab, label_source = "gate"),
    min_events = 20,
    group_summary = "median"
  )
  expect_true(any(is.na(sc$population) & sc$check == "rest_emd"))
  m <- sw_batch_matrix(sc, "rest_emd", c(rest_emd = "high"))
  expect_false(any(m$rank == 0, na.rm = TRUE))
  all_ev <- m[is.na(m$population), ]
  expect_gt(nrow(all_ev), 0)
  expect_false(anyNA(all_ev$rank))
  expect_equal(sort(all_ev$rank[all_ev$marker == "A"]), 1:3)
  # Planted effect: b3 is shifted in marker A.
  expect_equal(all_ev$batch[all_ev$marker == "A" & all_ev$rank == 1], "b3")
  # Labelled populations still rank on their own.
  tm <- m[!is.na(m$population) & m$population == "T" & m$marker == "A", ]
  expect_equal(sort(tm$rank), 1:3)
})

test_that("an sw_batch_mixing table works without a marker column", {
  files <- make_files()
  mx <- sw_batch_mixing(
    files,
    batch = stats::setNames(rep(c("b1", "b2", "b3"), each = 3), names(files)),
    markers = c("A", "B"),
    n_events = 60,
    seed = 1,
    k = 10,
    perplexity = 5,
    alpha = 0.05,
    min_events = 20,
    population = NULL
  )
  expect_false("marker" %in% names(mx))
  m <- sw_batch_matrix(mx, "lisi_batch", c(lisi_batch = "low"))
  expect_true(all(m$marker == "all"))
  expect_equal(sort(m$rank), 1:3)
  expect_equal(m$batch[m$rank == 1], "b3")
})

test_that("not_assessed rows are listed, not ranked, and reference_id is kept", {
  s <- rbind(
    transform(make_scores(), reference_id = "r1", status = "assessed"),
    transform(make_scores(), reference_id = "r2", status = "assessed")
  )
  s$value[s$reference_id == "r2" & s$marker == "A"] <- c(0.5, 0.1, 0.2)
  thin <- tibble::tibble(
    check = "rest_emd",
    marker = "A",
    batch = "b1",
    value = 7,
    reference_id = "r1",
    status = "not_assessed",
    file = c("f1", "f2")
  )
  s <- dplyr::bind_rows(s, thin)
  m <- sw_batch_matrix(s, "rest_emd", c(rest_emd = "high"))
  expect_true("reference_id" %in% names(m))
  expect_equal(nrow(m), 12)
  expect_equal(nrow(attr(m, "not_assessed")), 2)
  expect_equal(
    m$batch[m$reference_id == "r1" & m$marker == "A" & m$rank == 1],
    "b2"
  )
  expect_equal(
    m$batch[m$reference_id == "r2" & m$marker == "A" & m$rank == 1],
    "b1"
  )
})

test_that("a pair column with two batches gives rank NA", {
  s <- tibble::tibble(
    check = "pair_emd",
    marker = c("A", "B"),
    batch = "b1 ~ b2",
    value = c(0.4, 0.1)
  )
  m <- sw_batch_matrix(s, "pair_emd", c(pair_emd = "high"))
  expect_true(all(is.na(m$rank)))
})
