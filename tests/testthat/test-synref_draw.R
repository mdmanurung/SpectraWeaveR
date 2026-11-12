# tests/testthat/test-synref_draw.R
# sw_synref_downsample() and sw_synref_split() on mock data.

make_files <- function(n_files = 4, n = 200, seed = 1) {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    cbind(M1 = stats::rnorm(n, mean = i), M2 = stats::rnorm(n))
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

make_man <- function() {
  data.frame(
    sample_id = paste0("f", 1:4),
    batch = c("b1", "b1", "b1", "b2"),
    donor = c("d1", "d1", "d2", "d3")
  )
}

down <- function(d, ..., man = make_man()) {
  args <- list(
    data = d,
    manifest = man,
    groups = list(r1 = c("f1", "f2", "f3")),
    strata = "file",
    labels = NULL,
    n_events = 50,
    n_draws = 2,
    seed = 7
  )
  args[names(list(...))] <- list(...)
  do.call(sw_synref_downsample, args)
}

test_that("realised counts equal the request and every row is a source row", {
  d <- make_files()
  ref <- down(d)
  expect_s3_class(ref, "sw_synref")
  expect_equal(ref$certainty, "diagnostic_only")
  expect_equal(ref$scope, "within_batch")
  expect_true(all(ref$realised$realised == 50))
  expect_equal(nrow(ref$realised), 2 * 3)
  expect_equal(names(ref$events), c("r1_d01", "r1_d02"))
  expect_equal(nrow(ref$events[["r1_d01"]]), 150)
  # Each drawn row equals its source row.
  p <- ref$provenance[ref$provenance$synref_id == "r1_d01", ]
  rebuilt <- t(mapply(
    function(s, r) d[[s]][r, ],
    p$source_sample,
    p$source_row
  ))
  expect_equal(unname(rebuilt), unname(ref$events[["r1_d01"]]))
  # Without replacement: no repeated source row inside a draw.
  expect_false(anyDuplicated(paste(p$source_sample, p$source_row)) > 0)
  expect_equal(
    names(ref$provenance),
    c(
      "synref_id",
      "draw",
      "source_sample",
      "source_row",
      "batch",
      "donor",
      "population",
      "stratum"
    )
  )
})

test_that("the same seed gives the same draw and the caller's RNG is restored", {
  d <- make_files()
  set.seed(99)
  before <- .Random.seed
  a <- down(d)
  expect_identical(.Random.seed, before)
  b <- down(d)
  expect_identical(a$events, b$events)
  expect_identical(a$digest, b$digest)
  c2 <- down(d, seed = 8)
  expect_false(identical(a$events, c2$events))
})

test_that("a short stratum stops and names the stratum", {
  d <- make_files(n = 40)
  expect_error(down(d, n_events = 50), "Stratum 'f1'.*40 events; 50 requested")
})

test_that("a group that spans batches stops", {
  d <- make_files()
  expect_error(
    down(d, groups = list(r1 = c("f1", "f4"))),
    "within-batch only"
  )
})

test_that("donor and population strata draw per unit", {
  d <- make_files()
  ref <- down(d, strata = "donor", n_events = 100)
  expect_setequal(ref$realised$stratum, c("d1", "d2"))
  expect_true(all(ref$realised$realised == 100))
  expect_true(all(ref$provenance$donor == ref$provenance$stratum))

  labs <- lapply(d, \(m) rep(c("A", "B"), length.out = nrow(m)))
  refp <- down(d, strata = "population", labels = labs, n_events = 80)
  expect_setequal(refp$realised$stratum, c("A", "B"))
  expect_true(all(refp$provenance$population == refp$provenance$stratum))
  expect_error(down(d, strata = "population", labels = NULL), "needs 'labels'")
})

test_that("strata none pools the group into one stratum", {
  d <- make_files()
  ref <- down(d, strata = "none", n_events = 120, n_draws = 1)
  expect_equal(ref$realised$stratum, "all")
  expect_equal(nrow(ref$events[[1]]), 120)
})

test_that("strata none draws a pooled, unequal share per file", {
  d <- make_files()
  ref <- down(d, strata = "none", n_events = 120, n_draws = 1)
  per_file <- table(ref$provenance$source_sample)
  expect_equal(sum(per_file), 120)
  expect_gt(length(per_file), 1)
})

test_that("a repeated sample_id in a group stops, so no source event repeats", {
  d <- make_files()
  expect_error(
    down(d, groups = list(r1 = c("f1", "f1"))),
    "more than once: r1"
  )
  for (st in c("none", "file", "donor")) {
    ref <- down(d, strata = st, n_events = 40, n_draws = 1)
    key <- paste(ref$provenance$source_sample, ref$provenance$source_row)
    expect_false(anyDuplicated(key) > 0)
  }
})

test_that("missing required arguments stop and name them", {
  d <- make_files()
  expect_error(
    sw_synref_downsample(d, make_man()),
    "Required argument\\(s\\) with no default"
  )
  expect_error(
    sw_synref_downsample(
      d,
      make_man(),
      list(r1 = "f1"),
      "file",
      NULL,
      n_events = 10,
      n_draws = 1
    ),
    "'seed'"
  )
  expect_error(down(d, n_events = 0), "n_events")
  expect_error(down(d, seed = NA_real_), "seed")
})

test_that("split gives disjoint parts of equal size labelled sampling only", {
  d <- make_files(n = 300)
  sp <- sw_synref_split(
    d,
    make_man(),
    source = "f2",
    n_parts = 4,
    events_per_part = 60,
    seed = 3
  )
  expect_equal(sp$floor, "sampling only")
  expect_length(sp$events, 4)
  expect_true(all(vapply(sp$events, nrow, integer(1)) == 60))
  rows <- sp$provenance$source_row
  expect_equal(anyDuplicated(rows), 0)
  expect_true(all(sp$provenance$source_sample == "f2"))
  expect_true(all(sp$realised$requested == sp$realised$realised))
  set.seed(5)
  before <- .Random.seed
  a <- sw_synref_split(d, make_man(), "f2", 4, 60, 3)
  expect_identical(.Random.seed, before)
  expect_identical(a$digest, sp$digest)
})

test_that("split stops when the parts need more events than the source has", {
  d <- make_files(n = 100)
  expect_error(
    sw_synref_split(
      d,
      make_man(),
      "f1",
      n_parts = 3,
      events_per_part = 40,
      seed = 1
    ),
    "100 events; 3 parts of 40 need 120"
  )
  expect_error(sw_synref_split(d, make_man(), "f1", 3, 40), "'seed'")
})

test_that("pseudo-aliquots of one sample sit near the sampling floor", {
  d <- make_files(n = 4000)
  sp <- sw_synref_split(
    d,
    make_man(),
    "f1",
    n_parts = 2,
    events_per_part = 1500,
    seed = 2
  )
  m <- sp$events
  expect_lt(
    abs(stats::median(m[[1]][, "M1"]) - stats::median(m[[2]][, "M1"])),
    0.15
  )
})
