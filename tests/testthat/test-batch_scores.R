# tests/testthat/test-batch_scores.R
# Unit tests for R/batch_scores.R on simulated batches.

# Three batches of 5 files; marker A is shifted in batch b2, marker B is not.
make_batches <- function(shift = 1, n = 300, seed = 1) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  b <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  d <- lapply(ids, function(i) {
    cbind(
      A = stats::rnorm(n) + if (b[[i]] == "b2") shift else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(data = d, batch = b)
}

scores <- function(x, ...) {
  args <- list(
    data = x$data,
    batch = x$batch,
    markers = c("A", "B"),
    n_events = 200,
    seed = 1,
    order = NULL,
    reference_group = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_batch_scores, args)
}

get <- function(res, check, marker, batch = NA) {
  r <- res[res$check == check & res$marker == marker, ]
  if (!is.na(batch)) {
    r <- r[r$batch == batch, ]
  }
  r$value
}

test_that("required arguments have no default", {
  x <- make_batches()
  expect_error(sw_batch_scores(), "no default")
  expect_error(
    sw_batch_scores(x$data, x$batch, c("A", "B"), n_events = 100, seed = 1, reference_group = NULL),
    "order"
  )
})

test_that("a shifted batch shows in EMD, median range and eta2 for that marker only", {
  res <- scores(make_batches(shift = 1))
  expect_s3_class(res, "sw_batch_scores")
  # EMD of a pure +1 shift between equal-shape samples is close to 1.
  expect_equal(get(res, "pair_emd", "A", "b1 ~ b2"), 1, tolerance = 0.15)
  expect_lt(get(res, "pair_emd", "A", "b1 ~ b3"), 0.3)
  expect_lt(get(res, "pair_emd", "B", "b1 ~ b2"), 0.3)
  expect_gt(get(res, "median_range", "A"), 0.8)
  expect_lt(get(res, "median_range", "B"), 0.3)
  expect_gt(get(res, "eta2", "A"), 0.8)
  expect_lt(get(res, "eta2", "B"), 0.5)
})

test_that("no batch effect gives small scores", {
  res <- scores(make_batches(shift = 0))
  expect_lt(max(res$value[res$check == "pair_emd"]), 0.3)
  expect_lt(max(res$value[res$check == "median_range"]), 0.3)
  expect_lt(max(abs(log(res$value[res$check == "spread_ratio"]))), 0.1)
})

test_that("a widened batch shows in spread_ratio for that batch and marker only", {
  x <- make_batches(shift = 0)
  for (i in names(x$batch)[x$batch == "b2"]) {
    v <- x$data[[i]][, "A"]
    x$data[[i]][, "A"] <- mean(v) + (v - mean(v)) * 1.5
  }
  res <- scores(x)
  sr <- res[res$check == "spread_ratio", ]
  a_b2 <- sr$value[sr$marker == "A" & sr$batch == "b2"]
  expect_gt(a_b2, 1.3)
  expect_lt(sr$value[sr$marker == "A" & sr$batch == "b1"], 1.1)
  expect_lt(sr$value[sr$marker == "B" & sr$batch == "b2"], 1.1)
  # The median does not move, so median_range stays small for A.
  expect_lt(get(res, "median_range", "A"), 0.3)
})

test_that("order drift inside a batch is found, and small batches are skipped", {
  x <- make_batches(shift = 0)
  # Add a steady drift to marker A across the 5 files of b1.
  ids <- names(x$batch)[x$batch == "b1"]
  for (k in seq_along(ids)) {
    x$data[[ids[k]]][, "A"] <- x$data[[ids[k]]][, "A"] + 0.5 * k
  }
  ord <- stats::setNames(as.numeric(seq_along(x$batch)), names(x$batch))
  res <- scores(x, order = ord)
  expect_equal(get(res, "order_rho", "A", "b1"), 1)
  expect_lt(abs(get(res, "order_rho", "A", "b3")), 1)
  # A batch with fewer than 4 files gets no order_rho.
  y <- make_batches()
  y$batch[names(y$batch)[y$batch == "b3"][1:2]] <- "b4"
  res2 <- scores(
    y,
    order = stats::setNames(as.numeric(seq_along(y$batch)), names(y$batch))
  )
  expect_false("b4" %in% res2$batch[res2$check == "order_rho"])
})

test_that("input problems stop with a message", {
  x <- make_batches()
  expect_error(scores(x, n_events = 5000), "fewer than `n_events`")
  expect_error(scores(x, markers = "Z"), "lacks marker")
  expect_error(scores(x, batch = x$batch[1:3]), "covering every file")
  expect_error(
    scores(x, batch = stats::setNames(rep("b1", 15), names(x$batch))),
    "two batches"
  )
})

test_that("results do not depend on file order and keep the RNG state", {
  x <- make_batches()
  set.seed(5)
  a <- stats::runif(1)
  set.seed(5)
  r1 <- scores(x)
  expect_equal(stats::runif(1), a)
  x2 <- list(
    data = x$data[rev(names(x$data))],
    batch = x$batch[rev(names(x$batch))]
  )
  r2 <- scores(x2)
  expect_equal(r1$value, r2$value)
})

# reference_group ---------------------------------------------------------------

# Four batches of three files, three donors with offsets 0, 4 and 8 on A.
# Donors are spread unevenly over batches (b1: D1 D1 D2, b2: D2 D2 D3,
# b3: D3 D3 D1, b4: D1 D2 D3), so pooled batches differ by donor mix. Batch b4
# carries a real shift of 0.5 in A for every donor.
make_donor_batches <- function(seed = 5, n = 300) {
  set.seed(seed)
  donor_of <- c("D1", "D1", "D2", "D2", "D2", "D3", "D3", "D3", "D1", "D1", "D2", "D3")
  batch_of <- rep(c("b1", "b2", "b3", "b4"), each = 3)
  ids <- paste0("f", seq_along(donor_of))
  d <- lapply(seq_along(ids), function(i) {
    off <- c(D1 = 0, D2 = 4, D3 = 8)[[donor_of[i]]]
    cbind(
      A = stats::rnorm(n) + off + if (batch_of[i] == "b4") 0.5 else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(
    data = d,
    batch = stats::setNames(batch_of, ids),
    donor = stats::setNames(donor_of, ids)
  )
}

test_that("reference_group finds a batch effect hidden by donor differences", {
  x <- make_donor_batches()
  pair_a <- function(res) {
    r <- res[res$check == "pair_emd" & res$marker == "A", ]
    stats::setNames(r$value, r$batch)
  }
  all_files <- pair_a(scores(x))
  by_donor <- pair_a(scores(x, reference_group = x$donor))
  with_b4 <- grepl("b4", names(by_donor))
  # Within donors, every pair that has b4 is farther than every pair without.
  expect_gt(min(by_donor[with_b4]), max(by_donor[!with_b4]))
  # Over all files the donor mix swamps it.
  b4_all <- grepl("b4", names(all_files))
  expect_lt(min(all_files[b4_all]), max(all_files[!b4_all]))
  # All score types are still reported.
  res <- scores(x, reference_group = x$donor)
  expect_setequal(unique(res$check), unique(scores(x)$check))
})

test_that("batch reference_group is required, covers every file and needs 2 batches in a group", {
  x <- make_donor_batches()
  expect_error(
    sw_batch_scores(x$data, x$batch, c("A", "B"), n_events = 100, seed = 1, order = NULL),
    "no default"
  )
  expect_error(scores(x, reference_group = x$donor[-1]), "label for every file")
  # One file per group and batch: no group has two batches.
  one_each <- stats::setNames(names(x$data), names(x$data))
  expect_error(
    scores(x, reference_group = one_each),
    "at least two batches"
  )
})

# reference_group, two-way (data frame) --------------------------------------------------

test_that("a two-way reference_group is an error for batch scores", {
  x <- make_donor_batches()
  design <- data.frame(
    file = names(x$data),
    donor = unname(x$donor),
    batch = unname(x$batch),
    stringsAsFactors = FALSE
  )
  expect_error(scores(x, reference_group = design), "remove the batch effect")
})
