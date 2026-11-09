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
