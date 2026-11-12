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

test_that("an unlisted against file stops instead of passing the donor guard", {
  x <- make_ref()
  expect_error(
    sw_synref_check_overlap(
      x$ref,
      against = "zz",
      manifest = x$man,
      max_overlap_share = 0
    ),
    "missing from 'manifest' or without a donor: zz"
  )
  man_na <- x$man
  man_na$donor[man_na$sample_id == "f2"] <- NA
  expect_error(
    sw_synref_check_overlap(
      x$ref,
      against = "f2",
      manifest = man_na,
      max_overlap_share = 0
    ),
    "without a donor: f2"
  )
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

# sw_synref_spread() and sw_synref_fidelity() ---------------------------------------

make_spread_data <- function(n = 6000, seed = 8) {
  set.seed(seed)
  lab <- sample(c("T", "B"), n, replace = TRUE, prob = c(0.8, 0.2))
  list(
    data = list(
      a = cbind(M1 = stats::rnorm(n) + 4 * (lab == "B"), M2 = stats::rnorm(n))
    ),
    labels = list(a = lab),
    man = data.frame(sample_id = "a", batch = "b1")
  )
}

draw_n <- function(x, n_events, n_draws = 8, seed = 1) {
  sw_synref_downsample(
    x$data,
    x$man,
    groups = list(r = "a"),
    strata = "file",
    labels = NULL,
    n_events = n_events,
    n_draws = n_draws,
    seed = seed
  )
}

test_that("spread shrinks as events grow", {
  x <- make_spread_data()
  small <- sw_synref_spread(draw_n(x, 100), probs = c(0.5, 0.9))
  large <- sw_synref_spread(draw_n(x, 3000), probs = c(0.5, 0.9))
  expect_equal(nrow(small), 2 * 2)
  expect_equal(
    names(small),
    c(
      "group",
      "marker",
      "prob",
      "sd_across_draws",
      "range_across_draws",
      "n_draws",
      "mean_events",
      "status"
    )
  )
  expect_true(all(small$status == "ok"))
  expect_true(all(small$n_draws == 8))
  expect_true(all(small$mean_events == 100))
  expect_true(all(large$sd_across_draws < small$sd_across_draws))
  expect_true(all(small$range_across_draws >= small$sd_across_draws))
})

test_that("spread matches a hand computation", {
  x <- make_spread_data()
  ref <- draw_n(x, 200, n_draws = 4)
  s <- sw_synref_spread(ref, probs = 0.5)
  med <- vapply(ref$events, \(m) stats::median(m[, "M1"]), numeric(1))
  row <- s[s$marker == "M1", ]
  expect_equal(row$sd_across_draws, stats::sd(med))
  expect_equal(row$range_across_draws, diff(range(med)))
})

test_that("a group with one draw is not_assessed", {
  x <- make_spread_data()
  s <- sw_synref_spread(draw_n(x, 100, n_draws = 1), probs = 0.5)
  expect_true(all(s$status == "not_assessed"))
  expect_true(all(is.na(s$sd_across_draws)))
})

test_that("spread errors on missing or invalid arguments", {
  x <- make_spread_data()
  ref <- draw_n(x, 100)
  expect_error(sw_synref_spread(probs = 0.5), "Required argument")
  expect_error(sw_synref_spread(ref), "Required argument")
  expect_error(sw_synref_spread(ref, probs = 1.5), "probs")
  expect_error(sw_synref_spread(list(), probs = 0.5), "sw_synref")
})

test_that("fidelity is near zero for a random draw and flags a mix", {
  x <- make_spread_data()
  ref <- draw_n(x, 3000, n_draws = 1)
  f <- sw_synref_fidelity(
    ref,
    x$data,
    labels = x$labels,
    n_events = NULL,
    seed = 1
  )
  emd <- f[f$check == "emd_to_source", ]
  expect_equal(emd$item, c("M1", "M2"))
  expect_true(all(emd$emd_to_source < 0.1))
  pop <- f[f$check == "share_ratio", ]
  expect_equal(pop$item, c("B", "T"))
  expect_true(all(abs(pop$share_ratio - 1) < 0.1))
  # Source shares equal the label shares.
  expect_equal(pop$share_source[pop$item == "B"], mean(x$labels$a == "B"))

  # A 50:50 mix doubles the B share and moves M1, the planted B marker.
  mix <- sw_synref_mix(
    x$data,
    x$man,
    groups = list(r = "a"),
    labels = x$labels,
    proportions = c(T = 0.5, B = 0.5),
    n_events = 1000,
    n_draws = 1,
    seed = 2,
    short_stratum = "stop"
  )
  fm <- sw_synref_fidelity(
    mix,
    x$data,
    labels = x$labels,
    n_events = NULL,
    seed = 1
  )
  m1 <- fm$emd_to_source[fm$check == "emd_to_source" & fm$item == "M1"]
  m2 <- fm$emd_to_source[fm$check == "emd_to_source" & fm$item == "M2"]
  expect_gt(m1, 5 * m2)
  ratio_b <- fm$share_ratio[fm$check == "share_ratio" & fm$item == "B"]
  expect_equal(ratio_b, 0.5 / mean(x$labels$a == "B"), tolerance = 1e-8)
})

test_that("fidelity without labels has marker rows only and restores the RNG", {
  x <- make_spread_data()
  ref <- draw_n(x, 500, n_draws = 1)
  set.seed(3)
  before <- .Random.seed
  f <- sw_synref_fidelity(
    ref,
    x$data,
    labels = NULL,
    n_events = 1000,
    seed = 5
  )
  expect_identical(.Random.seed, before)
  expect_true(all(f$check == "emd_to_source"))
  again <- sw_synref_fidelity(
    ref,
    x$data,
    labels = NULL,
    n_events = 1000,
    seed = 5
  )
  expect_identical(f, again)
  other <- sw_synref_fidelity(
    ref,
    x$data,
    labels = NULL,
    n_events = 1000,
    seed = 6
  )
  expect_false(identical(f, other))
})

test_that("fidelity errors on missing or invalid arguments", {
  x <- make_spread_data()
  ref <- draw_n(x, 100, n_draws = 1)
  expect_error(
    sw_synref_fidelity(ref, x$data, labels = NULL, n_events = NULL),
    "Required argument"
  )
  expect_error(
    sw_synref_fidelity(ref, x$data, n_events = NULL, seed = 1),
    "Required argument"
  )
  expect_error(
    sw_synref_fidelity(ref, list(z = x$data$a), NULL, NULL, 1),
    "not in 'data'"
  )
  expect_error(
    sw_synref_fidelity(ref, x$data, list(a = "T"), NULL, 1),
    "labels"
  )
})
