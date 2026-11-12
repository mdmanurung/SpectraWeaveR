# Tests for R/batch_mixing.R (LISI and kBET rejection rate).

make_mix_files <- function(
  n_batches = 3,
  shift = 0,
  per_batch = 2,
  n = 300,
  seed = 21
) {
  set.seed(seed)
  b <- rep(seq_len(n_batches), each = per_batch)
  files <- lapply(b, \(i) {
    m <- matrix(
      stats::rnorm(n * 2),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    )
    m[, 1] <- m[, 1] + shift * i
    m
  })
  names(files) <- paste0("f", seq_along(files))
  list(files = files, batch = stats::setNames(paste0("b", b), names(files)))
}

mix_args <- function(d, ...) {
  utils::modifyList(
    list(
      d$files,
      batch = d$batch,
      markers = c("A", "B"),
      k = 45,
      perplexity = 15,
      alpha = 0.05,
      n_events = 200,
      seed = 5,
      min_events = 100,
      population = NULL
    ),
    list(...)
  )
}

value_of <- function(res, chk, pop = NA_character_) {
  res[["value"]][
    res[["check"]] == chk &
      is.na(res[["batch"]]) &
      (res[["population"]] %in% pop)
  ]
}

test_that("LISI is near the number of batches when mixed and near 1 when separated", {
  mixed <- do.call(sw_batch_mixing, mix_args(make_mix_files()))
  expect_s3_class(mixed, "sw_batch_mixing")
  expect_gt(value_of(mixed, "lisi_batch"), 2.4)
  expect_lte(value_of(mixed, "lisi_batch"), 3)
  apart <- do.call(sw_batch_mixing, mix_args(make_mix_files(shift = 50)))
  expect_lt(value_of(apart, "lisi_batch"), 1.05)
})

test_that("kBET rejection is near alpha when mixed and near 1 when separated", {
  mixed <- do.call(sw_batch_mixing, mix_args(make_mix_files()))
  expect_lt(value_of(mixed, "kbet_rejection"), 0.15)
  apart <- do.call(sw_batch_mixing, mix_args(make_mix_files(shift = 50)))
  expect_gt(value_of(apart, "kbet_rejection"), 0.95)
})

test_that("a partly separated design lies between the two extremes, per batch rows included", {
  d <- make_mix_files(shift = 1.5)
  res <- do.call(sw_batch_mixing, mix_args(d))
  lisi <- value_of(res, "lisi_batch")
  expect_gt(lisi, 1.05)
  expect_lt(lisi, 2.4)
  per_batch <- res[["batch"]][
    !is.na(res[["batch"]]) & res[["check"]] == "lisi_batch"
  ]
  expect_setequal(per_batch, c("b1", "b2", "b3"))
})

test_that("population restricts mixing, and a thin population is not assessed", {
  d <- make_mix_files(n_batches = 2, per_batch = 2)
  set.seed(8)
  pop <- lapply(d$files, \(m) {
    p <- rep("big", nrow(m))
    p[1:12] <- "tiny"
    p
  })
  # Separate the batches only inside "big" by shifting those events.
  for (f in names(d$files)) {
    if (d$batch[[f]] == "b2") {
      big <- pop[[f]] == "big"
      d$files[[f]][big, 1] <- d$files[[f]][big, 1] + 50
    }
  }
  res <- do.call(
    sw_batch_mixing,
    mix_args(d, population = pop, min_events = 100, k = 30, perplexity = 10)
  )
  big <- res[res[["population"]] == "big", ]
  tiny <- res[res[["population"]] == "tiny", ]
  expect_lt(
    big[["value"]][big[["check"]] == "lisi_batch" & is.na(big[["batch"]])],
    1.05
  )
  expect_true(all(tiny[["status"]] == "not_assessed"))
  expect_true(all(is.na(tiny[["value"]])))
  expect_true(all(big[["status"]] == "assessed"))
})

test_that("the same seed reproduces and the RNG state is restored", {
  d <- make_mix_files()
  set.seed(77)
  before <- .Random.seed
  a <- do.call(sw_batch_mixing, mix_args(d))
  expect_identical(.Random.seed, before)
  b <- do.call(sw_batch_mixing, mix_args(d))
  expect_equal(a[["value"]], b[["value"]])
  c2 <- do.call(sw_batch_mixing, mix_args(d, seed = 6))
  expect_false(isTRUE(all.equal(a[["value"]], c2[["value"]])))
})

test_that("missing required arguments and bad values stop", {
  d <- make_mix_files()
  expect_error(
    sw_batch_mixing(d$files, batch = d$batch, markers = c("A", "B")),
    "Required argument.*k.*perplexity.*alpha.*n_events.*seed.*min_events.*population"
  )
  expect_error(
    do.call(sw_batch_mixing, mix_args(d, perplexity = 45)),
    "perplexity"
  )
  expect_error(do.call(sw_batch_mixing, mix_args(d, alpha = 1.5)), "alpha")
  expect_error(
    do.call(sw_batch_mixing, mix_args(d, n_events = 5000)),
    "fewer events"
  )
  expect_error(
    do.call(sw_batch_mixing, mix_args(d, population = list(a = "x"))),
    "population"
  )
})
