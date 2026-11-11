# tests/testthat/test-anomaly_features.R
# Unit tests for R/anomaly_features.R. Pure matrices; no flow packages needed.

make_files <- function(n_files = 6, n = 200, seed = 1) {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    matrix(
      stats::rnorm(n * 2),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    )
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

feat <- function(data, ...) {
  args <- list(
    data = data,
    channels = c("A", "B"),
    features = c("quantile", "emd"),
    quantile_step = 0.25,
    n_events = 100,
    seed = 1,
    scale = "test",
    aggregate_from = "all",
    agg_size = 300,
    reference_group = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_anomaly_features, args)
}

test_that("sw_emd_1d is exact on known cases", {
  expect_equal(sw_emd_1d(c(0, 1), c(2, 3)), 2)
  x <- stats::rnorm(50)
  expect_equal(sw_emd_1d(x, x + 0.5), 0.5)
  expect_equal(sw_emd_1d(x, x), 0)
  # Unequal sizes: point mass at 0 vs {0, 2} moves half the mass by 2.
  expect_equal(sw_emd_1d(0, c(0, 2)), 1)
  # Symmetric.
  y <- stats::rnorm(30, 1)
  expect_equal(sw_emd_1d(x, y), sw_emd_1d(y, x))
  expect_error(sw_emd_1d(c(1, NA), 1), "without NA")
})

test_that("required arguments have no default", {
  d <- make_files()
  expect_error(sw_anomaly_features(), "no default")
  expect_error(
    sw_anomaly_features(d, channels = "A", features = "quantile", reference_group = NULL),
    "n_events"
  )
  expect_error(
    sw_anomaly_features(
      d,
      "A",
      "quantile",
      n_events = 10,
      seed = 1,
      scale = "x", reference_group = NULL
    ),
    "quantile_step"
  )
  expect_error(
    sw_anomaly_features(
      d,
      "A",
      "emd",
      n_events = 10,
      seed = 1,
      scale = "x", reference_group = NULL
    ),
    "aggregate_from|agg_size"
  )
  expect_error(
    feat(d, aggregate_from = "reference"),
    "reference"
  )
})

test_that("quantile columns exclude min and max", {
  d <- make_files()
  res <- feat(d, features = "quantile", quantile_step = 0.25)
  expect_s3_class(res, "sw_anomaly_features")
  expect_equal(
    names(res),
    c("file", "A_q0.25", "A_q0.5", "A_q0.75", "B_q0.25", "B_q0.5", "B_q0.75")
  )
  res10 <- feat(d, features = "quantile", quantile_step = 0.1)
  expect_equal(ncol(res10) - 1, 2 * 9)
})

test_that("a shifted file has a larger EMD than clean files", {
  d <- make_files(n_files = 8, n = 400)
  d$f1[, "A"] <- d$f1[, "A"] + 2
  res <- feat(d, features = "emd", n_events = 300, agg_size = 800)
  expect_equal(res$file[which.max(res$A_emd)], "f1")
  expect_true(res$A_emd[res$file == "f1"] > 1)
  expect_true(all(res$A_emd[res$file != "f1"] < 0.5))
})

test_that("results do not depend on list order and restore the RNG", {
  d <- make_files()
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  a <- feat(d)
  after <- stats::runif(1)
  expect_equal(before, after)
  b <- feat(d[rev(names(d))])
  expect_equal(unclass(a)[-1], unclass(b)[-1])
  expect_equal(a$file, b$file)
  c2 <- feat(d, seed = 2)
  expect_false(isTRUE(all.equal(a$A_emd, c2$A_emd)))
})

test_that("reference pooling uses only the named files", {
  d <- make_files(n_files = 6, n = 300)
  d$f6[, "A"] <- d$f6[, "A"] + 5
  ref <- feat(
    d,
    features = "emd",
    aggregate_from = "reference",
    reference = c("f1", "f2", "f3"),
    agg_size = 600
  )
  all_ <- feat(d, features = "emd", aggregate_from = "all", agg_size = 600)
  # An anomalous file in the pool pulls the pool towards it, so f1 looks
  # further from a pool that holds f6 than from a clean pool.
  expect_true(all_$A_emd[all_$file == "f1"] > ref$A_emd[ref$file == "f1"])
})

test_that("a file with too few events is an error", {
  d <- make_files(n = 50)
  expect_error(feat(d, n_events = 100), "fewer than `n_events`")
  d$f1 <- d$f1[, "A", drop = FALSE]
  expect_error(
    feat(
      make_files() |>
        (\(x) {
          x$f1 <- x$f1[, "A", drop = FALSE]
          x
        })()
    ),
    "lacks channel"
  )
})

# reference_group ---------------------------------------------------------------

# Three donors of four files; donor offsets 0, 4 and 8 on A. File d2_4 is
# shifted by 0.5 on top of its donor.
make_donor_files <- function(seed = 6, n = 200, per = 4) {
  set.seed(seed)
  ids <- paste0("d", rep(1:3, each = per), "_", rep(seq_len(per), 3))
  d <- lapply(ids, function(i) {
    donor <- as.integer(substr(i, 2, 2))
    cbind(
      A = stats::rnorm(n) + 4 * (donor - 1) + if (i == "d2_4") 0.5 else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(data = d, donor = stats::setNames(substr(ids, 1, 2), ids))
}

test_that("reference_group measures the EMD against the file's own group", {
  x <- make_donor_files()
  all_files <- feat(x$data, features = "emd", agg_size = 600)
  by_donor <- feat(x$data, features = "emd", agg_size = 600, reference_group = x$donor)
  # Over all files, the EMD mostly records the donor offset.
  expect_gt(all_files$A_emd[all_files$file == "d1_1"], 3)
  expect_lt(by_donor$A_emd[by_donor$file == "d1_1"], 0.5)
  # Within the group, the shifted file is the one farthest from its pool.
  expect_equal(by_donor$file[which.max(by_donor$A_emd)], "d2_4")
  # Quantile features are unchanged.
  q_all <- feat(x$data, features = "quantile")
  q_grp <- feat(x$data, features = "quantile", reference_group = x$donor)
  expect_equal(q_grp, q_all, ignore_attr = TRUE)
})

test_that("with aggregate_from = 'reference' every group needs a reference file", {
  x <- make_donor_files()
  ok <- feat(
    x$data,
    features = "emd",
    aggregate_from = "reference",
    reference = c("d1_1", "d1_2", "d2_1", "d2_2", "d3_1", "d3_2"),
    reference_group = x$donor
  )
  expect_equal(nrow(ok), 12)
  expect_error(
    feat(
      x$data,
      features = "emd",
      aggregate_from = "reference",
      reference = c("d1_1", "d2_1"),
      reference_group = x$donor
    ),
    "no file in `reference`"
  )
})

test_that("anomaly reference_group is required and must cover every file", {
  x <- make_donor_files()
  expect_error(
    sw_anomaly_features(x$data, "A", "quantile", 0.25, 100, 1, "test"),
    "no default"
  )
  expect_error(feat(x$data, reference_group = x$donor[-1]), "label for every file")
  expect_error(feat(x$data, reference_group = unname(x$donor)), "label for every file")
})

test_that("sw_diagnose takes emd_z within groups", {
  skip_if_not_installed("isotree")
  x <- make_donor_files(per = 6)
  an <- list(
    channels = c("A", "B"), features = "emd", quantile_step = NULL, n_events = 200,
    seed = 1, scale = "test", aggregate_from = "all", agg_size = 600,
    reference = NULL, ntrees = 50, sample_size = 8, ndim = 1
  )
  ez <- function(rg) {
    dx <- sw_diagnose(
      x$data, NULL, NULL, c(an, list(reference_group = rg)),
      NULL, NULL, NULL, NULL, NULL, NULL
    )
    dx[dx$check == "emd_z" & dx$item == "A", ]
  }
  by_donor <- ez(x$donor)
  expect_equal(by_donor$subject[which.max(abs(by_donor$value))], "d2_4")
  all_files <- ez(NULL)
  expect_false(all_files$subject[which.max(abs(all_files$value))] == "d2_4")
})

# reference_group, two-way (data frame) --------------------------------------------------

test_that("a two-way reference_group pools all files and medpolishes the EMD values in sw_diagnose", {
  skip_if_not_installed("isotree")
  x <- make_donor_files(per = 6)
  design <- data.frame(
    file = names(x$data),
    donor = unname(x$donor),
    batch = rep(paste0("b", 1:6), 3),
    stringsAsFactors = FALSE
  )
  # Features are the same as with NULL: the EMD reference is the pool of all files.
  expect_equal(
    feat(x$data, features = "emd", agg_size = 600, reference_group = design),
    feat(x$data, features = "emd", agg_size = 600, reference_group = NULL),
    ignore_attr = TRUE
  )
  an <- list(
    channels = c("A", "B"), features = "emd", quantile_step = NULL, n_events = 200,
    seed = 1, scale = "test", aggregate_from = "all", agg_size = 600,
    reference = NULL, ntrees = 50, sample_size = 8, ndim = 1
  )
  dx <- sw_diagnose(
    x$data, NULL, NULL, c(an, list(reference_group = design)),
    NULL, NULL, NULL, NULL, NULL, NULL
  )
  ez <- dx[dx$check == "emd_z" & dx$item == "A", ]
  emd <- feat(x$data, features = "emd", agg_size = 600, n_events = 200, reference_group = NULL)
  expected <- .sw_robust_z(
    emd$A_emd,
    structure(
      list(a = design$donor[match(emd$file, design$file)], b = design$batch[match(emd$file, design$file)]),
      class = "sw_two_way"
    )
  )
  expect_equal(ez$value[match(emd$file, ez$subject)], expected)
  expect_error(
    feat(x$data, reference_group = design[c("donor", "batch")]),
    "file. column"
  )
  expect_error(feat(x$data, reference_group = design[-1, ]), "covers every file")
})

test_that("medpolish residual z is NA with fewer than 3 finite values and handles empty cells", {
  a <- c("d1", "d1", "d2", "d2", "d3")
  b <- c("b1", "b2", "b1", "b3", "b2")
  tw <- structure(list(a = a, b = b), class = "sw_two_way")
  z <- .sw_robust_z(c(1, 2, 3.1, 4, 5.2), tw)
  expect_length(z, 5)
  expect_true(any(is.finite(z)))
  expect_true(all(is.na(.sw_robust_z(c(1, 2, NA, NA, NA), tw))))
})
