# tests/testthat/test-signal_reference.R
# Frozen signal reference: fit, save and load, and use in sw_fcs_signal_audit().

make_cohort <- function(n_files, n = 300, seed = 1, shift = 0, prefix = "f") {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    cbind(
      A = stats::rnorm(n, 0 + shift),
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
  })
  names(f) <- paste0(prefix, seq_len(n_files))
  f
}

audit_ref <- function(d, ref, ...) {
  args <- list(
    data = d,
    channels = c("A", "B", "C"),
    range_max = NULL,
    range_min = NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = ref
  )
  args[names(list(...))] <- list(...)
  do.call(sw_fcs_signal_audit, args)
}

zval <- function(res, check, channel) {
  r <- res[res$check == check & res$channel %in% channel, ]
  stats::setNames(r$value, r$file)
}

test_that("fit stores centre and scale of file-level medians and SDs", {
  d <- make_cohort(6, seed = 2, prefix = "r")
  ref <- sw_signal_reference_fit(d, c("A", "B"))
  expect_s3_class(ref, "sw_signal_reference")
  med <- vapply(d, \(m) stats::median(m[, "A"]), numeric(1))
  sds <- vapply(d, \(m) stats::sd(m[, "A"]), numeric(1))
  a <- ref$channels[ref$channels$channel == "A", ]
  expect_equal(a$median_centre, stats::median(med))
  expect_equal(a$median_scale, 1.4826 * stats::mad(med, constant = 1))
  expect_equal(a$sd_centre, stats::median(sds))
  expect_equal(a$sd_scale, 1.4826 * stats::mad(sds, constant = 1))
  expect_equal(ref$files, sort(names(d)))
  expect_equal(ref$n_files, 6)
  expect_false(is.na(ref$digest))
  expect_equal(ref$schema_version, "sw_signal_reference/1")
})

test_that("frozen z equals a hand computation on a toy cohort", {
  ref_files <- make_cohort(8, seed = 3, prefix = "r")
  ref <- sw_signal_reference_fit(ref_files, c("A", "B", "C"))
  new <- make_cohort(5, seed = 4, prefix = "n")
  res <- audit_ref(new, ref)
  for (ch in c("A", "B", "C")) {
    rc <- ref$channels[ref$channels$channel == ch, ]
    med <- vapply(new, \(m) stats::median(m[, ch]), numeric(1))
    sds <- vapply(new, \(m) stats::sd(m[, ch]), numeric(1))
    expect_equal(
      zval(res, "median_z", ch)[names(new)],
      (med - rc$median_centre) / rc$median_scale
    )
    expect_equal(
      zval(res, "sd_z", ch)[names(new)],
      (sds - rc$sd_centre) / rc$sd_scale
    )
  }
  # Other checks are the cohort ones, unchanged.
  base <- audit_ref(new, NULL)
  keep <- !base$check %in% c("median_z", "sd_z")
  expect_equal(res[keep, ], base[keep, ])
})

test_that("a fault shared by every scored file shows frozen, not cohort", {
  ref <- sw_signal_reference_fit(
    make_cohort(12, seed = 5, prefix = "r"),
    c("A", "B", "C")
  )
  new <- make_cohort(10, seed = 6, shift = 1, prefix = "n") # all files shifted on A
  frozen <- audit_ref(new, ref)
  cohort <- audit_ref(new, NULL)
  expect_true(all(abs(zval(frozen, "median_z", "A")) > 5))
  # The cohort z only sees spread among the shifted files: noise-sized.
  expect_gt(
    min(abs(zval(frozen, "median_z", "A"))),
    2 * max(abs(zval(cohort, "median_z", "A")))
  )
  # An unshifted channel stays unremarkable under the frozen reference.
  expect_true(all(abs(zval(frozen, "median_z", "B")) < 5))
})

test_that("save and load round-trip identically", {
  ref <- sw_signal_reference_fit(make_cohort(5, seed = 7), c("A", "B"))
  path <- tempfile(fileext = ".rds")
  expect_identical(sw_signal_reference_save(ref, path), path)
  back <- sw_signal_reference_load(path)
  expect_identical(back, ref)
  expect_s3_class(back, "sw_signal_reference")
})

test_that("a changed file fails the digest, a wrong class fails too", {
  ref <- sw_signal_reference_fit(make_cohort(5, seed = 8), c("A", "B"))
  tampered <- ref
  tampered$channels$median_centre[1] <- tampered$channels$median_centre[1] + 1
  path <- tempfile(fileext = ".rds")
  saveRDS(tampered, path)
  expect_error(sw_signal_reference_load(path), "Digest mismatch")
  expect_error(sw_signal_reference_save(tampered, path), "Digest mismatch")
  saveRDS(list(a = 1), path)
  expect_error(sw_signal_reference_load(path), "sw_signal_reference")
  old <- ref
  old$schema_version <- "sw_signal_reference/0"
  saveRDS(old, path)
  expect_error(sw_signal_reference_load(path), "schema_version")
})

test_that("fit needs 3 files, named data and its arguments", {
  expect_error(sw_signal_reference_fit(), "no default")
  expect_error(
    sw_signal_reference_fit(make_cohort(2), c("A", "B")),
    "at least 3 files"
  )
  expect_error(
    sw_signal_reference_fit(unname(make_cohort(4)), c("A", "B")),
    "named list"
  )
  expect_error(
    sw_signal_reference_fit(make_cohort(4), c("A", "Z")),
    "lacks column"
  )
})

test_that("reference with reference_group stops; missing channel stops", {
  ref <- sw_signal_reference_fit(make_cohort(6, seed = 9), c("A", "B"))
  new <- make_cohort(6, seed = 10, prefix = "n")
  grp <- stats::setNames(rep(c("x", "y"), 3), names(new))
  expect_error(
    audit_ref(new, ref, reference_group = grp, channels = c("A", "B")),
    "cannot be used together"
  )
  expect_error(audit_ref(new, ref), "missing from `reference`")
  expect_error(
    sw_fcs_signal_audit(
      new,
      c("A", "B"),
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "'reference'"
  )
  expect_error(
    audit_ref(new, list(a = 1), channels = c("A", "B")),
    "sw_signal_reference"
  )
})

test_that("a zero reference scale gives NA z", {
  const <- lapply(1:4, \(i) cbind(A = rep(i * 0 + 1, 20), B = stats::rnorm(20)))
  names(const) <- paste0("r", 1:4)
  ref <- sw_signal_reference_fit(const, c("A", "B"))
  new <- make_cohort(4, seed = 11, prefix = "n")
  res <- audit_ref(new, ref, channels = c("A", "B"))
  expect_true(all(is.na(zval(res, "median_z", "A"))))
  expect_true(all(is.finite(zval(res, "median_z", "B"))))
})
