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

# sw_synref_mix() and sw_synref_goal() ----------------------------------------------

# Two files in one batch; T is 70% and B 30% of events, B is shifted in M1.
make_mix <- function(n = 400, seed = 3) {
  set.seed(seed)
  lab <- rep(c("T", "B"), c(0.7 * n, 0.3 * n))
  one <- function(shift) {
    cbind(
      M1 = stats::rnorm(n, mean = ifelse(lab == "B", 5, 0)) + shift,
      M2 = stats::rnorm(n)
    )
  }
  list(
    data = list(a = one(0), b = one(0.1), c = one(2)),
    labels = list(a = lab, b = lab, c = lab),
    man = data.frame(
      sample_id = c("a", "b", "c"),
      batch = c("b1", "b1", "b2"),
      donor = c("d1", "d2", "d3")
    )
  )
}

mix_args <- function(x, ...) {
  args <- list(
    data = x$data,
    manifest = x$man,
    groups = list(r = c("a", "b")),
    labels = x$labels,
    proportions = c(T = 0.5, B = 0.5),
    n_events = 100,
    n_draws = 3,
    seed = 11,
    short_stratum = "stop"
  )
  args[names(list(...))] <- list(...)
  args
}

test_that("mix realises the requested counts and proportions", {
  x <- make_mix()
  ref <- do.call(sw_synref_mix, mix_args(x))
  expect_s3_class(ref, "sw_synref")
  expect_equal(ref$certainty, "diagnostic_only")
  expect_equal(ref$scope, "within_batch")
  expect_equal(nrow(ref$realised), 3 * 2)
  expect_equal(ref$realised$realised, ref$realised$requested)
  expect_equal(unique(ref$realised$stratum), c("T", "B"))
  expect_equal(ref$realised$realised_prop, ref$realised$requested_prop)
  # Planted B shift: a 50:50 mix has a higher M1 mean than the 30% source.
  src <- mean(unlist(lapply(x$data[c("a", "b")], \(m) m[, "M1"])))
  expect_gt(mean(ref$events[[1]][, "M1"]), src + 0.5)
  # Every drawn row equals its source row and carries its label.
  p <- ref$provenance[ref$provenance$synref_id == "r_d01", ]
  rebuilt <- t(mapply(\(s, r) x$data[[s]][r, ], p$source_sample, p$source_row))
  expect_equal(unname(rebuilt), unname(ref$events[["r_d01"]]))
  lab_of <- mapply(\(s, r) x$labels[[s]][r], p$source_sample, p$source_row)
  expect_equal(unname(lab_of), p$population)
})

test_that("mix uses largest-remainder counts that sum to n_events", {
  x <- make_mix()
  ref <- do.call(
    sw_synref_mix,
    mix_args(x, proportions = c(T = 0.333, B = 0.667), n_events = 101)
  )
  expect_true(all(
    tapply(ref$realised$realised, ref$realised$synref_id, sum) == 101
  ))
})

test_that("mix is seed reproducible and restores the RNG state", {
  x <- make_mix()
  set.seed(5)
  before <- .Random.seed
  a <- do.call(sw_synref_mix, mix_args(x))
  expect_identical(.Random.seed, before)
  b <- do.call(sw_synref_mix, mix_args(x))
  d <- do.call(sw_synref_mix, mix_args(x, seed = 12))
  expect_identical(a$digest, b$digest)
  expect_false(identical(a$digest, d$digest))
})

test_that("a short population stops or is capped and reported", {
  x <- make_mix()
  # 100% B asks for 200 B events; the group holds 240.
  expect_error(
    do.call(
      sw_synref_mix,
      mix_args(x, proportions = c(T = 0.1, B = 0.9), n_events = 400)
    ),
    "Population 'B'"
  )
  ref <- do.call(
    sw_synref_mix,
    mix_args(
      x,
      proportions = c(T = 0.1, B = 0.9),
      n_events = 400,
      short_stratum = "cap",
      n_draws = 1
    )
  )
  b <- ref$realised[ref$realised$stratum == "B", ]
  expect_equal(b$requested, 360)
  expect_equal(b$realised, 240)
  expect_lt(b$realised_prop, b$requested_prop)
  expect_equal(nrow(ref$events[[1]]), 40 + 240)
})

test_that("mix errors on missing and invalid arguments", {
  x <- make_mix()
  args <- mix_args(x)
  for (nm in names(args)) {
    expect_error(
      do.call(sw_synref_mix, args[setdiff(names(args), nm)]),
      "Required argument"
    )
  }
  expect_error(
    do.call(sw_synref_mix, mix_args(x, short_stratum = "drop")),
    "short_stratum"
  )
  expect_error(do.call(sw_synref_mix, mix_args(x, labels = NULL)), "labels")
  expect_error(
    do.call(sw_synref_mix, mix_args(x, proportions = c(T = 0.5, B = 0.4))),
    "sum to 1"
  )
  expect_error(
    do.call(sw_synref_mix, mix_args(x, groups = list(r = c("a", "c")))),
    "within-batch"
  )
  expect_error(
    do.call(
      sw_synref_mix,
      mix_args(x, proportions = c(T = 0.5, Z = 0.5))
    ),
    "not in the labels"
  )
})

make_goal_data <- function() {
  set.seed(2)
  list(
    data = list(
      a = cbind(M1 = stats::rnorm(300, 0), M2 = stats::rnorm(300)),
      b = cbind(M1 = stats::rnorm(300, 0), M2 = stats::rnorm(300)),
      c = cbind(M1 = stats::rnorm(300, 1), M2 = stats::rnorm(300)),
      d = cbind(M1 = stats::rnorm(300, 4), M2 = stats::rnorm(300))
    ),
    man = data.frame(
      sample_id = c("a", "b", "c", "d"),
      batch = c("b1", "b1", "b2", "b3")
    )
  )
}

goal_call <- function(x, ...) {
  args <- list(
    data = x$data,
    manifest = x$man,
    probs = c(0.25, 0.5, 0.75),
    goal_statistic = "median",
    n_events = 200,
    seed = 4
  )
  args[names(list(...))] <- list(...)
  do.call(sw_synref_goal, args)
}

test_that("goal is the mean or median of batch quantiles and has role goal", {
  x <- make_goal_data()
  g <- goal_call(x)
  expect_s3_class(g, "sw_synref_goal")
  expect_equal(g$role, "goal")
  expect_equal(g$certainty, "diagnostic_only")
  expect_false(inherits(g, "sw_synref"))
  expect_equal(nrow(g$goal), 2 * 3)
  expect_equal(nrow(g$batch_quantiles), 3 * 2 * 3)
  bq <- g$batch_quantiles
  one <- bq[bq$marker == "M1" & bq$prob == 0.5, ]
  expect_equal(
    g$goal$goal[g$goal$marker == "M1" & g$goal$prob == 0.5],
    stats::median(one$value)
  )
  gm <- goal_call(x, goal_statistic = "mean")
  expect_equal(
    gm$goal$goal[gm$goal$marker == "M1" & gm$goal$prob == 0.5],
    mean(gm$batch_quantiles$value[
      gm$batch_quantiles$marker == "M1" & gm$batch_quantiles$prob == 0.5
    ])
  )
  # The planted M1 batch difference shows as a larger spread than in M2.
  sd_m1 <- g$goal$spread_sd[g$goal$marker == "M1" & g$goal$prob == 0.5]
  sd_m2 <- g$goal$spread_sd[g$goal$marker == "M2" & g$goal$prob == 0.5]
  expect_gt(sd_m1, 5 * sd_m2)
  expect_equal(
    bq$delta,
    bq$value -
      g$goal$goal[match(
        paste(bq$marker, bq$prob),
        paste(g$goal$marker, g$goal$prob)
      )]
  )
  expect_output(print(g), "role: goal")
})

test_that("goal is seed reproducible and restores the RNG state", {
  x <- make_goal_data()
  set.seed(9)
  before <- .Random.seed
  a <- goal_call(x)
  expect_identical(.Random.seed, before)
  expect_identical(a$digest, goal_call(x)$digest)
  expect_false(identical(a$digest, goal_call(x, seed = 5)$digest))
})

test_that("goal errors on missing or invalid arguments", {
  x <- make_goal_data()
  args <- list(
    data = x$data,
    manifest = x$man,
    probs = 0.5,
    goal_statistic = "mean",
    n_events = 100,
    seed = 1
  )
  for (nm in names(args)) {
    expect_error(
      do.call(sw_synref_goal, args[setdiff(names(args), nm)]),
      "Required argument"
    )
  }
  expect_error(goal_call(x, goal_statistic = "mode"), "goal_statistic")
  expect_error(goal_call(x, probs = c(0, 0.5)), "probs")
  # b1 holds 600 events, b2 and b3 hold 300.
  expect_error(goal_call(x, n_events = 400), "Batch 'b2'")
  one_batch <- x$man
  one_batch$batch <- "b1"
  expect_error(goal_call(x, manifest = one_batch), "two batches")
})

test_that("cap warns for an absent population and stops when all are absent", {
  x <- make_mix()
  expect_warning(
    ref <- do.call(
      sw_synref_mix,
      mix_args(
        x,
        proportions = c(T = 0.5, B = 0.3, Z = 0.2),
        short_stratum = "cap",
        n_draws = 1
      )
    ),
    "Z"
  )
  z <- ref$realised[ref$realised$stratum == "Z", ]
  expect_equal(z$realised, 0)
  expect_error(
    suppressWarnings(do.call(
      sw_synref_mix,
      mix_args(x, proportions = c(Z = 1), short_stratum = "cap")
    )),
    "none of the requested populations"
  )
})
