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
