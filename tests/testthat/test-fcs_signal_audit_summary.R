# tests/testthat/test-fcs_signal_audit_summary.R
# The `summary` argument of sw_fcs_signal_audit(): per file x channel counts,
# location and quantiles, with thin files reported as not_assessed.

make_sum <- function() {
  set.seed(3)
  d <- lapply(1:4, function(i) {
    cbind(
      A = stats::rnorm(200, i),
      B = stats::rnorm(200),
      C = stats::rnorm(200)
    )
  })
  names(d) <- paste0("f", 1:4)
  d$f2[1:10, "A"] <- NA # 10 non-finite values
  d$f3[1:5, "B"] <- Inf
  d
}

run_sum <- function(d, summary) {
  sw_fcs_signal_audit(
    d,
    channels = c("A", "B", "C"),
    range_max = NULL,
    range_min = NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL,
    summary = summary
  )
}

sum_set <- function(...) {
  args <- list(probs = c(0.25, 0.75), min_events = 50, stat = "median")
  args[names(list(...))] <- list(...)
  args
}

test_that("a missing summary stops with no default", {
  expect_error(
    sw_fcs_signal_audit(
      make_sum(),
      c("A", "B", "C"),
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "no default.*summary"
  )
})

test_that("summary = NULL leaves the output unchanged", {
  d <- make_sum()
  off <- run_sum(d, NULL)
  on <- run_sum(d, sum_set())
  expect_false(any(grepl("^summary_", off$check)))
  expect_null(attr(off, "provenance"))
  base <- on[!grepl("^summary_", on$check), ]
  expect_equal(
    as.data.frame(base),
    as.data.frame(off),
    ignore_attr = TRUE
  )
  expect_identical(
    names(off)[1:6],
    c("file", "check", "channel", "value", "detail", "status")
  )
  expect_true("reason" %in% names(off))
})

test_that("summary rows hold counts, location and quantiles", {
  d <- make_sum()
  res <- run_sum(d, sum_set())
  get <- function(f, ch, chk) {
    res$value[res$file == f & res$channel == ch & res$check == chk]
  }
  expect_equal(get("f1", "A", "summary_n"), 200)
  expect_equal(get("f2", "A", "summary_nonfinite"), 10)
  expect_equal(get("f3", "B", "summary_nonfinite"), 5)
  expect_equal(get("f1", "A", "summary_nonfinite"), 0)
  expect_equal(
    get("f2", "A", "summary_median"),
    stats::median(d$f2[11:200, "A"])
  )
  expect_equal(
    get("f1", "C", "summary_q0.75"),
    unname(stats::quantile(d$f1[, "C"], 0.75))
  )
  # one row per file and channel for each summary check
  expect_equal(sum(res$check == "summary_median"), 12)
  expect_equal(sum(res$check == "summary_q0.25"), 12)
})

test_that("stat = mean gives the mean and no median row", {
  d <- make_sum()
  res <- run_sum(d, sum_set(stat = "mean"))
  expect_equal(
    res$value[
      res$file == "f1" & res$channel == "B" & res$check == "summary_mean"
    ],
    mean(d$f1[, "B"])
  )
  expect_false("summary_median" %in% res$check)
})

test_that("files below min_events become not_assessed rows with a reason", {
  d <- make_sum()
  d$f4 <- d$f4[1:30, ]
  res <- run_sum(d, sum_set(min_events = 50))
  thin <- res[res$file == "f4" & grepl("^summary_(median|q)", res$check), ]
  expect_equal(nrow(thin), 9) # median + 2 quantiles x 3 channels
  expect_true(all(is.na(thin$value)))
  expect_true(all(grepl("^not_assessed: ", thin$detail)))
  expect_match(thin$detail[1], "30")
  # `reason` carries the text without the `not_assessed: ` prefix
  expect_false(anyNA(thin$reason))
  expect_match(thin$reason[1], "^30 finite events")
  # the counts are still reported, and no file or channel is dropped
  expect_equal(
    res$value[res$file == "f4" & res$channel == "A" & res$check == "summary_n"],
    30
  )
  expect_setequal(unique(res$file[grepl("^summary_", res$check)]), names(d))
  # a file that passes is not marked
  ok <- res[res$file == "f1" & grepl("^summary_(median|q)", res$check), ]
  expect_true(all(is.na(ok$detail)))
  expect_true(all(is.na(ok$reason)))
})

test_that("a channel with too few finite values is not_assessed even in a big file", {
  d <- make_sum()
  d$f1[1:180, "C"] <- NA
  res <- run_sum(d, sum_set(min_events = 50))
  r <- res[
    res$file == "f1" & res$channel == "C" & res$check == "summary_median",
  ]
  expect_true(is.na(r$value))
  expect_match(r$detail, "^not_assessed")
  r2 <- res[
    res$file == "f1" & res$channel == "A" & res$check == "summary_median",
  ]
  expect_false(is.na(r2$value))
})

test_that("summary carries a provenance attribute", {
  res <- run_sum(make_sum(), sum_set())
  p <- attr(res, "provenance")
  expect_true(all(
    c("package_version", "r_version", "seed", "summary") %in% names(p)
  ))
  expect_equal(p$summary$min_events, 50)
  expect_true(is.na(p$seed))
})

test_that("summary is validated before any heavy work", {
  d <- make_sum()
  expect_error(run_sum(d, 1), "summary")
  expect_error(run_sum(d, list(probs = 0.5, min_events = 5)), "stat")
  expect_error(run_sum(d, sum_set(probs = NULL)), "probs")
  expect_error(run_sum(d, sum_set(probs = c(0.5, 1.2))), "probs")
  expect_error(run_sum(d, sum_set(min_events = -1)), "min_events")
  expect_error(run_sum(d, sum_set(min_events = NULL)), "min_events")
  expect_error(run_sum(d, sum_set(stat = "mode")), "stat")
  expect_error(run_sum(d, sum_set(probs = c(0.5, 0.5))), "probs")
})

test_that("quantile check names do not depend on the other probs", {
  d <- make_sum()
  q_checks <- function(probs) {
    unique(grep(
      "^summary_q",
      run_sum(d, sum_set(probs = probs))$check,
      value = TRUE
    ))
  }
  expect_identical(q_checks(0.5), "summary_q0.5")
  expect_true("summary_q0.5" %in% q_checks(c(0.5, 0.25)))
  expect_true("summary_q0.5" %in% q_checks(c(0.5, 0.123)))
  expect_true("summary_q0.123" %in% q_checks(c(0.5, 0.123)))
})

test_that("provenance starts with the shared keys (signal_audit)", {
  pv <- attr(run_sum(make_sum(), sum_set()), "provenance")
  expect_equal(
    head(names(pv), 8),
    c(
      "package_version",
      "flowcore_version",
      "r_version",
      "platform",
      "datetime",
      "seed",
      "session",
      "args"
    )
  )
  expect_true(all(
    c("running", "locale", "attached") %in% names(pv$session)
  ))
})
