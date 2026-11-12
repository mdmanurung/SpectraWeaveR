# Tests for R/batch_otd.R (joint batch distance by Sinkhorn divergence).

make_otd_files <- function(shift_b3 = c(3, 0), seed = 11) {
  set.seed(seed)
  shifts <- list(c(0, 0), c(0, 0), c(0, 0), c(0, 0), shift_b3, shift_b3)
  files <- lapply(shifts, \(s) {
    m <- matrix(stats::rnorm(600), ncol = 2, dimnames = list(NULL, c("A", "B")))
    sweep(m, 2, s, "+")
  })
  names(files) <- paste0("f", 1:6)
  list(
    files = files,
    batch = stats::setNames(c("b1", "b1", "b2", "b2", "b3", "b3"), names(files))
  )
}

otd_args <- function(d, mode) {
  list(
    d$files,
    batch = d$batch,
    markers = c("A", "B"),
    n_events = 120,
    blur = 1,
    seed = 3,
    max_iter = 300,
    tol = 1e-5,
    mode = mode
  )
}

test_that("rest mode ranks the shifted batch first with a plausible size", {
  d <- make_otd_files()
  res <- do.call(sw_batch_otd, otd_args(d, "rest"))
  expect_s3_class(res, "sw_batch_otd")
  expect_equal(res[["batch"]][which.max(res[["otd"]])], "b3")
  # b3 is moved by (3, 0); against a pool that is two thirds unshifted, the
  # divergence is well above the two unshifted batches.
  expect_gt(
    min(res[["otd"]][res[["batch"]] == "b3"]),
    2 * max(res[["otd"]][res[["batch"]] != "b3"])
  )
})

test_that("pairwise mode gives one row per pair and the planted pairs are largest", {
  d <- make_otd_files()
  res <- do.call(sw_batch_otd, otd_args(d, "pairwise"))
  expect_equal(nrow(res), 3)
  expect_setequal(res[["batch"]], c("b1 ~ b2", "b1 ~ b3", "b2 ~ b3"))
  expect_equal(res[["batch"]][which.min(res[["otd"]])], "b1 ~ b2")
  # Pure translation by (3, 0): |d|^2 / 2 = 4.5, up to sampling noise.
  big <- res[["otd"]][res[["batch"]] != "b1 ~ b2"]
  expect_equal(big, c(4.5, 4.5), tolerance = 0.25)
})

test_that("the same seed reproduces and the RNG state is restored", {
  d <- make_otd_files()
  set.seed(99)
  before <- .Random.seed
  a <- do.call(sw_batch_otd, otd_args(d, "rest"))
  expect_identical(.Random.seed, before)
  b <- do.call(sw_batch_otd, otd_args(d, "rest"))
  expect_equal(a[["otd"]], b[["otd"]])
  other <- otd_args(d, "rest")
  other[["seed"]] <- 4
  expect_false(isTRUE(all.equal(
    a[["otd"]],
    do.call(sw_batch_otd, other)[["otd"]]
  )))
})

test_that("missing required arguments and bad values stop", {
  d <- make_otd_files()
  expect_error(
    sw_batch_otd(d$files, batch = d$batch, markers = c("A", "B")),
    "Required argument.*n_events.*blur.*seed.*max_iter.*tol.*mode"
  )
  a <- otd_args(d, "rest")
  a[["mode"]] <- "all"
  expect_error(do.call(sw_batch_otd, a), "mode")
  a <- otd_args(d, "rest")
  a[["blur"]] <- 0
  expect_error(do.call(sw_batch_otd, a), "blur")
  a <- otd_args(d, "rest")
  a[["n_events"]] <- 5000
  expect_error(do.call(sw_batch_otd, a), "fewer events")
  a <- otd_args(d, "rest")
  a[["batch"]] <- rep("b1", 6) |> stats::setNames(names(d$files))
  expect_error(do.call(sw_batch_otd, a), "two batches")
})

test_that("rest scores do not depend on which file of a batch is shifted", {
  # b1 holds one shifted and one unshifted file; the order of the two must not
  # change what the other batches see in their rest set.
  mk <- function(shifted_first) {
    d <- make_otd_files(shift_b3 = c(0, 0))
    s <- c(6, 0)
    target <- if (shifted_first) "f1" else "f2"
    d$files[[target]] <- sweep(d$files[[target]], 2, s, "+")
    d
  }
  r1 <- do.call(sw_batch_otd, otd_args(mk(TRUE), "rest"))
  r2 <- do.call(sw_batch_otd, otd_args(mk(FALSE), "rest"))
  other1 <- r1[["otd"]][r1[["batch"]] %in% c("b2", "b3")]
  other2 <- r2[["otd"]][r2[["batch"]] %in% c("b2", "b3")]
  # The shifted file enters every other batch's rest set either way.
  expect_gt(min(other1), 1)
  expect_gt(min(other2), 1)
  expect_equal(other1, other2, tolerance = 0.35)
})

test_that("n_files and n_events report uneven batches", {
  d <- make_otd_files()
  d$batch[["f6"]] <- "b2"
  res <- do.call(sw_batch_otd, otd_args(d, "rest"))
  expect_equal(res[["n_files"]][res[["batch"]] == "b2"], 3L)
  expect_equal(res[["n_files"]][res[["batch"]] == "b3"], 1L)
  # b2: 3 files x floor(120 / 3) = 120 events; rest of b2: 60 from each of b1
  # (2 files x 30) and b3 (1 file x 60).
  expect_equal(res[["n_events"]][res[["batch"]] == "b2"], "120 ~ 120")
  expect_equal(res[["n_events"]][res[["batch"]] == "b1"], "120 ~ 120")
})
