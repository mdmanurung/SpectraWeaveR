# tests/testthat/test-fcs_signal_audit.R
# Unit tests for R/fcs_signal_audit.R on small simulated files.

make_sig <- function(n_files = 5, n = 400, seed = 1) {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    cbind(
      Time = sort(stats::runif(n, 0, 100)),
      A = stats::rnorm(n),
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

audit <- function(d, ...) {
  args <- list(
    data = d,
    channels = c("A", "B", "C"),
    range_max = NULL,
    range_min = NULL,
    time_channel = NULL,
    time_bin = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_fcs_signal_audit, args)
}

val <- function(res, file, check, channel = NA) {
  r <- res[res$file == file & res$check == check, ]
  if (!is.na(channel)) {
    r <- r[r$channel == channel, ]
  }
  r$value
}

test_that("required arguments have no default", {
  expect_error(sw_fcs_signal_audit(), "no default")
  d <- make_sig()
  expect_error(
    sw_fcs_signal_audit(d, channels = c("A", "B")),
    "range_max"
  )
})

test_that("clean files give the base checks only", {
  res <- audit(make_sig())
  expect_s3_class(res, "sw_fcs_signal_audit")
  expect_setequal(
    unique(res$check),
    c(
      "n_events",
      "nonfinite_share",
      "zero_variance",
      "negative_share",
      "corr_shift"
    )
  )
  expect_equal(val(res, "f1", "n_events"), 400)
  expect_equal(val(res, "f1", "zero_variance", "A"), 0)
})

test_that("a dead channel, non-finite values and saturation are found", {
  d <- make_sig()
  d$f2[, "B"] <- 0
  d$f3[1:40, "C"] <- NA
  d$f4[1:100, "A"] <- 1000
  res <- audit(d, range_max = 1000, range_min = -1000)
  expect_equal(val(res, "f2", "zero_variance", "B"), 1)
  expect_equal(val(res, "f1", "zero_variance", "B"), 0)
  expect_equal(val(res, "f3", "nonfinite_share", "C"), 0.1)
  expect_equal(val(res, "f4", "saturation_share", "A"), 0.25)
  expect_equal(val(res, "f1", "saturation_share", "A"), 0)
  expect_equal(val(res, "f1", "floor_share", "A"), 0)
})

test_that("rate_cv counts events per bin, so thinning a window raises it", {
  # Regression: rate_cv once counted occupied bins, not events, so a window
  # with 80% fewer events (no empty bin) looked like a clean file.
  d <- make_sig(n = 2000)
  win <- d$f5[, "Time"] > 40 & d$f5[, "Time"] < 70
  set.seed(9)
  drop <- win & stats::runif(nrow(d$f5)) < 0.8
  d$f5 <- d$f5[!drop, ]
  res <- audit(d, time_channel = "Time", time_bin = 5)
  clean <- c(val(res, "f1", "rate_cv"), val(res, "f2", "rate_cv"))
  expect_true(val(res, "f5", "rate_cv") > 2 * max(clean))
  # Exact value on a hand-made file: bins of 1 unit hold 3, 1, 0, 2 events.
  m <- cbind(Time = c(0.1, 0.2, 0.3, 1.5, 3.2, 3.9), A = 1:6, B = 6:1)
  d2 <- list(x = m, y = m)
  r2 <- sw_fcs_signal_audit(d2, c("A", "B"), NULL, NULL, "Time", 1)
  counts <- c(3, 1, 0, 2)
  expect_equal(
    r2$value[r2$check == "rate_cv"][1],
    stats::sd(counts) / mean(counts)
  )
})

test_that("a marker drifting with time shows in time_trend for that channel only", {
  d <- make_sig(n = 600)
  d$f2[, "A"] <- d$f2[, "A"] + 0.03 * d$f2[, "Time"] # about 3 SD over the run
  res <- audit(d, time_channel = "Time", time_bin = 5)
  # A ramp 0 to 3 on unit noise correlates at about sqrt(0.75 / 1.75) = 0.65.
  expect_equal(val(res, "f2", "time_trend", "A"), 0.65, tolerance = 0.1)
  expect_lt(abs(val(res, "f2", "time_trend", "B")), 0.2)
  expect_lt(abs(val(res, "f1", "time_trend", "A")), 0.2)
  # Without a time channel the check is absent.
  expect_false("time_trend" %in% audit(d)$check)
})

test_that("a gap in event times raises rate_cv and max_gap", {
  d <- make_sig()
  keep <- !(d$f5[, "Time"] > 40 & d$f5[, "Time"] < 70)
  d$f5 <- d$f5[keep, ]
  res <- audit(d, time_channel = "Time", time_bin = 5)
  expect_true(
    val(res, "f5", "rate_cv") >
      2 * max(val(res, "f1", "rate_cv"), val(res, "f2", "rate_cv"))
  )
  expect_true(val(res, "f5", "max_gap") > 10 * val(res, "f1", "max_gap"))
})

test_that("a changed channel correlation is named by pair", {
  d <- make_sig(n = 600)
  d$f3[, "B"] <- d$f3[, "A"] * 0.9 + stats::rnorm(600, sd = 0.3)
  res <- audit(d)
  r <- res[res$check == "corr_shift", ]
  expect_equal(r$file[which.max(r$value)], "f3")
  expect_equal(r$detail[r$file == "f3"], "A ~ B")
  expect_true(r$value[r$file == "f3"] > 0.5)
})

test_that("limits can be per channel and must cover every channel", {
  d <- make_sig()
  expect_error(audit(d, range_max = c(A = 1, B = 1)), "name every channel")
  res <- audit(d, range_max = c(A = 1000, B = 1000, C = 1000))
  expect_true("saturation_share" %in% res$check)
  expect_error(audit(d, time_channel = "Time"), "time_bin")
  expect_error(audit(d, time_channel = "T2", time_bin = 1), "lacks column")
})
