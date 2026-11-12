# tests/testthat/test-synref_checks.R
# sw_synref_check_overlap() on mock data.

make_ref <- function() {
  set.seed(1)
  d <- list(
    f1 = cbind(M1 = stats::rnorm(100)),
    f2 = cbind(M1 = stats::rnorm(100)),
    f3 = cbind(M1 = stats::rnorm(100))
  )
  man <- data.frame(
    sample_id = c("f1", "f2", "f3"),
    batch = "b1",
    donor = c("d1", "d1", "d2")
  )
  ref <- sw_synref_downsample(
    d,
    man,
    groups = list(r = "f1"),
    strata = "file",
    labels = NULL,
    n_events = 30,
    n_draws = 1,
    seed = 1
  )
  list(ref = ref, man = man)
}

test_that("a disjoint file and donor pass", {
  x <- make_ref()
  out <- sw_synref_check_overlap(
    x$ref,
    against = "f3",
    manifest = x$man,
    max_overlap_share = 0
  )
  expect_equal(out$unit, c("event", "file", "donor"))
  expect_equal(out$status, c("not_assessed", "ok", "ok"))
  expect_equal(out$share[2:3], c(0, 0))
})

test_that("a shared file stops", {
  x <- make_ref()
  expect_error(
    sw_synref_check_overlap(
      x$ref,
      against = "f1",
      manifest = x$man,
      max_overlap_share = 0
    ),
    "r_d01 \\(file"
  )
})

test_that("a sibling file of the same donor stops on the donor unit", {
  x <- make_ref()
  expect_error(
    sw_synref_check_overlap(
      x$ref,
      against = "f2",
      manifest = x$man,
      max_overlap_share = 0
    ),
    "\\(donor, share 1.000\\)"
  )
  # Without the manifest the donor unit is not assessed, so the sibling passes.
  out <- sw_synref_check_overlap(
    x$ref,
    against = "f2",
    manifest = NULL,
    max_overlap_share = 0
  )
  expect_equal(out$status[out$unit == "donor"], "not_assessed")
})

test_that("event overlap is found per row and respects max_overlap_share", {
  x <- make_ref()
  p <- x$ref$provenance
  # Plant 6 of the 30 drawn events in the scoring set.
  planted <- data.frame(sample_id = "f3", row = 1:3)
  none <- sw_synref_check_overlap(x$ref, planted, x$man, max_overlap_share = 0)
  expect_equal(none$n_shared[none$unit == "event"], 0L)
  hit <- data.frame(sample_id = p$source_sample[1:6], row = p$source_row[1:6])
  expect_error(
    sw_synref_check_overlap(x$ref, hit, x$man, 0),
    "event, share 0.200"
  )
  ok <- sw_synref_check_overlap(
    x$ref,
    hit[0, ],
    x$man,
    max_overlap_share = 0.25
  )
  expect_equal(ok$n_shared[ok$unit == "event"], 0L)
  # A file-level hit stays above 0.25, so it stops even when events are within it.
  expect_error(sw_synref_check_overlap(x$ref, hit, x$man, 0.25), "file")
})

test_that("required arguments and bad limits stop", {
  x <- make_ref()
  expect_error(
    sw_synref_check_overlap(x$ref, "f3", x$man),
    "'max_overlap_share'"
  )
  expect_error(
    sw_synref_check_overlap(x$ref, "f3", max_overlap_share = 0),
    "'manifest'"
  )
  expect_error(sw_synref_check_overlap(x$ref, "f3", x$man, 2), "0 to 1")
  expect_error(sw_synref_check_overlap(list(), "f3", x$man, 0), "sw_synref")
})
