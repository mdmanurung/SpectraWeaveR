# tests/testthat/test-ref_noise_floor.R
# sw_ref_noise_floor() on mock data.

make_nf <- function(shift_a2 = 0) {
  set.seed(1)
  d <- list(
    r1 = cbind(A = rnorm(800), B = rnorm(800)),
    r2a = cbind(A = rnorm(500), B = rnorm(500)),
    r2b = cbind(A = rnorm(500) + shift_a2, B = rnorm(500)),
    r3 = cbind(A = rnorm(100), B = rnorm(100))
  )
  man <- data.frame(
    sample_id = names(d),
    batch = c("b1", "b2", "b2", "b3"),
    reference_role = "technical_replicate",
    reference_id = "ref1",
    aliquot = c("a1", "a1", "a2", "a1")
  )
  list(d = d, man = man)
}

nf <- function(x, ...) {
  args <- list(
    data = x$d,
    manifest = x$man,
    markers = c("A", "B"),
    n_events = 200,
    n_splits = 4,
    n_quantiles = NULL,
    seed = 1,
    roles = "technical_replicate"
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_noise_floor, args)
}

test_that("required arguments have no default", {
  x <- make_nf()
  expect_error(sw_ref_noise_floor(), "no default")
  expect_error(
    sw_ref_noise_floor(
      x$d,
      x$man,
      c("A", "B"),
      n_events = 100,
      n_splits = 2,
      seed = 1,
      roles = "technical_replicate"
    ),
    "n_quantiles"
  )
  expect_error(nf(x, n_splits = 0), "n_splits")
})

test_that("the cell split is labelled sampling only and sits near zero", {
  res <- nf(make_nf())
  expect_s3_class(res, "sw_ref_noise_floor")
  s <- res[res$floor_type == "sampling only" & res$status == "assessed", ]
  expect_setequal(s$unit, c("r1", "r2a", "r2b"))
  expect_equal(nrow(s[s$unit == "r1" & s$marker == "A", ]), 4)
  # Two halves of one N(0,1) sample of 200 events differ by about 0.12 on
  # average; single splits of one file have a tail.
  expect_lt(mean(s$value), 0.2)
  expect_lt(max(s$value), 0.6)
  expect_gt(min(s$value), 0)
  expect_false(any(c("pass", "verdict", "flag") %in% names(res)))
})

test_that("an aliquot floor needs two aliquots in a batch and finds a planted shift", {
  res <- nf(make_nf())
  al <- res[res$floor_type == "aliquot", ]
  expect_equal(unique(al$batch), "b2")
  expect_equal(unique(al$unit), "a1 ~ a2")
  expect_lt(max(al$value), 0.35)
  shifted <- nf(make_nf(shift_a2 = 1))
  is_al <- shifted$floor_type == "aliquot"
  expect_equal(
    mean(shifted$value[is_al & shifted$marker == "A"]),
    1,
    tolerance = 0.2
  )
  expect_lt(max(shifted$value[is_al & shifted$marker == "B"]), 0.35)
  # Without an aliquot column there is no aliquot floor.
  x <- make_nf()
  x$man$aliquot <- NULL
  expect_false("aliquot" %in% nf(x)$floor_type)
})

test_that("thin units are not_assessed and the quantile grid form agrees roughly", {
  res <- nf(make_nf())
  thin <- res[res$unit == "r3", ]
  expect_true(all(thin$status == "not_assessed" & is.na(thin$value)))
  expect_match(thin$reason[1], "2 \\* n_events")
  big <- nf(make_nf(), n_events = 600)
  expect_true(all(big$status[big$floor_type == "aliquot"] == "not_assessed"))
  grid <- nf(make_nf(), n_quantiles = 100)
  expect_equal(
    mean(grid$value, na.rm = TRUE),
    mean(res$value, na.rm = TRUE),
    tolerance = 0.2
  )
})

test_that("the seed reproduces the table and the caller's RNG state is restored", {
  x <- make_nf()
  set.seed(5)
  before <- .Random.seed
  a <- nf(x)
  expect_identical(.Random.seed, before)
  expect_equal(a, nf(x))
  expect_false(isTRUE(all.equal(a$value, nf(x, seed = 2)$value)))
})
