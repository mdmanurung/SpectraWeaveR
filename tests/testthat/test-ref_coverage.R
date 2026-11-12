# tests/testthat/test-ref_coverage.R

make_cov <- function(seed = 1, n = 3000) {
  set.seed(seed)
  list(
    ref = cbind(A = stats::rnorm(n), B = stats::rnorm(n)),
    study = list(
      s1 = cbind(A = stats::rnorm(n), B = stats::rnorm(n)),
      s2 = cbind(A = stats::rnorm(n, mean = 2), B = stats::rnorm(n))
    ),
    batch = c(s1 = "b1", s2 = "b2")
  )
}

cov_run <- function(x, ...) {
  args <- list(
    reference = x$ref,
    study = x$study,
    batch = x$batch,
    markers = c("A", "B"),
    probs = c(0.01, 0.99),
    reference_labels = NULL,
    study_labels = NULL,
    min_events = 100,
    transform = "asinh",
    cofactor = 6000
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_coverage, args)
}

test_that("a planted range shift raises share_above only for the shifted batch and marker", {
  res <- cov_run(make_cov())
  g <- \(b, m, col) res[[col]][res$batch == b & res$marker == m]
  expect_lt(g("b1", "A", "share_above"), 0.03)
  expect_gt(g("b2", "A", "share_above"), 0.25)
  expect_lt(g("b2", "A", "share_below"), 0.005)
  expect_lt(g("b2", "B", "share_above"), 0.03)
  expect_equal(unique(res$expected_above), 0.01)
  expect_true(all(res$n_reference_cells == 3000))
})

test_that("populations are scored separately and thin ones are not_assessed", {
  x <- make_cov()
  ref_lab <- rep(c("T", "B"), c(2900, 100))
  st_lab <- list(
    s1 = rep(c("T", "B"), c(2900, 100)),
    s2 = rep(c("T", "B"), c(2990, 10))
  )
  res <- cov_run(x, reference_labels = ref_lab, study_labels = st_lab)
  expect_setequal(res$population, c("T", "B"))
  thin <- res[res$batch == "b2" & res$population == "B", ]
  expect_true(all(thin$status == "not_assessed"))
  expect_true(all(is.na(thin$share_above)))
  ok <- res[res$batch == "b1" & res$population == "B", ]
  expect_true(all(ok$status == "ok"))
})

test_that("it is deterministic and leaves the RNG alone", {
  x <- make_cov()
  set.seed(4)
  before <- .Random.seed
  r1 <- cov_run(x)
  expect_identical(.Random.seed, before)
  expect_identical(r1, cov_run(x))
})

test_that("required arguments and bad inputs stop", {
  expect_error(sw_ref_coverage(), "Required argument")
  x <- make_cov()
  expect_error(cov_run(x, probs = c(0.9, 0.1)), "probs")
  expect_error(cov_run(x, probs = 0.99), "probs")
  expect_error(cov_run(x, reference_labels = rep("a", 3000)), "both")
  expect_error(cov_run(x, markers = "Z"), "all `markers`")
  expect_error(cov_run(x, batch = c(s1 = "b1")), "covering")
})

test_that("study files with reordered marker columns are bound by name", {
  x <- make_cov()
  x$study$s3 <- x$study$s2[, c("B", "A")]
  x$batch <- c(x$batch, s3 = "b2")
  res <- cov_run(x)
  g <- \(m) res$share_above[res$batch == "b2" & res$marker == m]
  expect_gt(g("A"), 0.25)
  expect_lt(g("B"), 0.03)
})
