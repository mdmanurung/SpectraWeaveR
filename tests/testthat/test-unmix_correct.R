# Tests for R/unmix_correct.R. AutoSpectral is not required: argument checks
# run first, and the AutoSpectral call is replaced by a fake in the mock tests.

mk_spectra <- function() {
  matrix(
    c(
      1,
      0,
      0,
      0,
      0,
      1,
      0.5,
      0,
      0.2,
      0.3,
      1,
      0.1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("FITC", "PE", "AF"), c("d1", "d2", "d3", "d4"))
  )
}
tmp_fcs <- function() {
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  f
}

# Replace the three AutoSpectral seams for the calling test only. Works when
# the file is sourced into globalenv() and under devtools::test().
mock_autospectral <- function(fake_call, env = parent.frame()) {
  fakes <- list(
    .sw_correct_require_autospectral = function() invisible(TRUE),
    .sw_correct_call = fake_call,
    .sw_correct_version = function() "1.8.8-mock"
  )
  home <- environment(sw_unmix_correct)
  if (identical(home, globalenv())) {
    old <- mget(names(fakes), envir = home)
    for (n in names(fakes)) {
      assign(n, fakes[[n]], envir = home)
    }
    withr::defer(
      for (n in names(old)) {
        assign(n, old[[n]], envir = home)
      },
      envir = env
    )
  } else {
    do.call(local_mocked_bindings, c(fakes, list(.env = env)))
  }
}

# ---- sw_unmix_compare_spectra ------------------------------------------------

test_that("compare_spectra gives 0 for identical and 90 for orthogonal rows", {
  orig <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("A", "B"), c("d1", "d2"))
  )
  corr <- orig
  corr["B", ] <- c(1, 0)
  out <- sw_unmix_compare_spectra(orig, corr)
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("fluorophore", "angle_deg", "max_abs_delta", "changed"))
  expect_equal(out[["fluorophore"]], c("A", "B"))
  expect_identical(out[["angle_deg"]][1], 0)
  expect_equal(out[["angle_deg"]][2], 90)
  expect_equal(out[["max_abs_delta"]], c(0, 1))
  expect_equal(out[["changed"]], c(FALSE, TRUE))
})

test_that("compare_spectra peak-normalises and treats rescaling as unchanged", {
  orig <- matrix(c(1, 0.5), nrow = 1, dimnames = list("A", c("d1", "d2")))
  corr <- orig * 3
  out <- sw_unmix_compare_spectra(orig, corr)
  expect_equal(out[["max_abs_delta"]], 0)
  expect_equal(out[["angle_deg"]], 0, tolerance = 1e-6)

  corr2 <- matrix(c(2, 2), nrow = 1, dimnames = list("A", c("d1", "d2")))
  out2 <- sw_unmix_compare_spectra(orig, corr2)
  expect_equal(out2[["max_abs_delta"]], 0.5)
  expect_equal(
    out2[["angle_deg"]],
    acos(sum(c(1, 0.5) * c(2, 2)) / sqrt(1.25 * 8)) * 180 / pi
  )
})

test_that("compare_spectra aligns rows and columns by name", {
  orig <- mk_spectra()
  corr <- orig
  corr["PE", "d3"] <- 0.6
  shuffled <- corr[c("AF", "PE", "FITC"), c("d4", "d2", "d1", "d3")]
  expect_equal(
    sw_unmix_compare_spectra(orig, shuffled),
    sw_unmix_compare_spectra(orig, corr)
  )
  out <- sw_unmix_compare_spectra(orig, shuffled)
  expect_equal(out[["fluorophore"]], rownames(orig))
  expect_equal(out[["changed"]], c(FALSE, TRUE, FALSE))
  expect_equal(out[["max_abs_delta"]][2], 0.1)
})

test_that("compare_spectra errors on missing names, never aligns by position", {
  orig <- mk_spectra()
  expect_error(
    sw_unmix_compare_spectra(orig, orig[c("FITC", "PE"), ]),
    "fluorophore names differ: AF",
    fixed = TRUE
  )
  renamed <- orig
  rownames(renamed)[3] <- "BV421"
  expect_error(
    sw_unmix_compare_spectra(orig, renamed),
    "AF, BV421",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_compare_spectra(orig, orig[, c("d1", "d2", "d3")]),
    "detector names differ: d4",
    fixed = TRUE
  )
})

test_that("spectra checks reject bad matrices", {
  m <- mk_spectra()
  expect_error(sw_unmix_compare_spectra(as.data.frame(m), m), "numeric matrix")
  no_names <- unname(m)
  expect_error(sw_unmix_compare_spectra(no_names, m), "row names")
  dup <- m
  rownames(dup) <- c("FITC", "FITC", "AF")
  expect_error(sw_unmix_compare_spectra(dup, m), "unique")
  dup_col <- m
  colnames(dup_col) <- c("d1", "d1", "d3", "d4")
  expect_error(sw_unmix_compare_spectra(dup_col, m), "unique")
  bad <- m
  bad[1, 1] <- NA
  expect_error(sw_unmix_compare_spectra(bad, m), "finite")
})

# ---- Wrapper validation ------------------------------------------------------

test_that("sw_unmix_correct_signatures validates input before the dependency", {
  s <- mk_spectra()
  raw <- matrix(1, 5, 4, dimnames = list(NULL, colnames(s)))
  f <- tmp_fcs()
  fc <- list(a = 1)
  expect_error(
    sw_unmix_correct_signatures(raw, unname(s), NULL),
    "row names (fluorophores)",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct_signatures(raw, s, asp = "x"),
    "`asp` must be a list"
  )
  expect_error(
    sw_unmix_correct_signatures(raw[, 4:1], s, NULL),
    "columns must equal `colnames(spectra)`",
    fixed = TRUE
  )
  expect_error(sw_unmix_correct_signatures(1:3, s, NULL), "numeric matrix")
  expect_error(
    sw_unmix_correct_signatures("no_such_file.fcs", s, NULL),
    "not found: no_such_file.fcs",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct_signatures(c(f, f), s, NULL),
    "reads one fully stained file"
  )
  expect_error(
    sw_unmix_correct_signatures(f, s, NULL),
    "needs `unstained.sample`"
  )
  expect_error(
    sw_unmix_correct_signatures(
      f,
      s,
      NULL,
      unstained.sample = "missing.fcs",
      flow.control = fc
    ),
    "unstained.sample` file(s) not found",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct_signatures(f, s, NULL, unstained.sample = f),
    "needs `flow.control`"
  )
  expect_error(
    sw_unmix_correct_signatures(raw, s, NULL, raw.data = raw),
    "Do not pass `raw.data`"
  )
  expect_error(sw_unmix_correct_signatures(raw, s, NULL, 5), "must be named")
})

test_that("sw_unmix_correct_spillover validates input before the dependency", {
  s <- mk_spectra()
  f <- tmp_fcs()
  fc <- list(a = 1)
  expect_error(sw_unmix_correct_spillover(f, s, NULL), "`asp` must be a list")
  expect_error(
    sw_unmix_correct_spillover(matrix(1, 2, 4), s, list()),
    "reads files only"
  )
  expect_error(
    sw_unmix_correct_spillover(c(f, "gone.fcs"), s, list()),
    "not found: gone.fcs",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct_spillover(f, s, list()),
    "needs `unstained.sample`"
  )
  expect_error(
    sw_unmix_correct_spillover(f, s, list(), unstained.sample = f),
    "needs `flow.control`"
  )
  expect_error(
    sw_unmix_correct_spillover(
      f,
      s,
      list(),
      unstained.sample = f,
      flow.control = fc,
      fully.stained.sample = f
    ),
    "Do not pass `fully.stained.sample`"
  )
  expect_error(
    sw_unmix_correct_spillover(f, s[, c(1, 1, 2, 3)], list()),
    "unique"
  )
})

test_that("sw_unmix_correct requires explicit, ordered steps", {
  s <- mk_spectra()
  f <- tmp_fcs()
  expect_error(sw_unmix_correct(f, s, list()), "`steps` is required")
  expect_error(sw_unmix_correct(f, s, list(), steps = "fix"), "unique values")
  expect_error(
    sw_unmix_correct(f, s, list(), steps = c("signatures", "signatures")),
    "unique values"
  )
  expect_error(
    sw_unmix_correct(f, s, list(), steps = c("spillover", "signatures")),
    "Run \"signatures\" before \"spillover\"",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct(f, s, list(), steps = "spillover", spillover_args = "x"),
    "must be lists"
  )
  expect_error(
    sw_unmix_correct(f, s, list(), steps = "spillover"),
    "needs `unstained.sample`"
  )
})

test_that("missing AutoSpectral gives a clear message after validation", {
  skip_if(
    requireNamespace("AutoSpectral", quietly = TRUE),
    "AutoSpectral is installed"
  )
  s <- mk_spectra()
  f <- tmp_fcs()
  raw <- matrix(1, 5, 4, dimnames = list(NULL, colnames(s)))
  fc <- list(a = 1)
  msg <- "Package 'AutoSpectral' is required."
  expect_error(.sw_correct_require_autospectral(), msg, fixed = TRUE)
  expect_error(sw_unmix_correct_signatures(raw, s, NULL), msg, fixed = TRUE)
  expect_error(
    sw_unmix_correct_spillover(
      f,
      s,
      list(),
      unstained.sample = f,
      flow.control = fc
    ),
    msg,
    fixed = TRUE
  )
  expect_error(
    sw_unmix_correct(
      f,
      s,
      list(),
      steps = c("signatures", "spillover"),
      signatures_args = list(unstained.sample = f, flow.control = fc),
      spillover_args = list(unstained.sample = f, flow.control = fc)
    ),
    msg,
    fixed = TRUE
  )
})

# ---- Mocked AutoSpectral -----------------------------------------------------

test_that("signatures result is packaged by name from a fake AutoSpectral", {
  s <- mk_spectra()
  raw <- matrix(1, 5, 4, dimnames = list(NULL, colnames(s)))
  seen <- new.env()
  mock_autospectral(function(fn, args) {
    seen$fn <- fn
    seen$args <- args
    out <- args$spectra[c("AF", "PE", "FITC"), ]
    out["PE", ] <- c(0, 1, 0.4, 0)
    list(
      spectra = out,
      fit.log = data.frame(fluorophore = c("FITC", "PE"), iter = 1L),
      accepted = c(PE = TRUE, FITC = FALSE)
    )
  })
  res <- sw_unmix_correct_signatures(raw, s, NULL, min.gain = 0.01)

  expect_equal(seen$fn, "correct.unmixing.signatures")
  expect_identical(seen$args$raw.data, raw)
  expect_equal(seen$args$min.gain, 0.01)
  expect_null(seen$args$fully.stained.sample)

  expect_s3_class(res, "sw_unmix_correction")
  expect_named(
    res,
    c("spectra", "original", "changes", "fit_logs", "note", "provenance")
  )
  expect_equal(rownames(res$spectra), rownames(s))
  expect_identical(res$original, s)
  expect_equal(res$spectra["PE", "d3"], 0.4)
  ch <- res$changes
  expect_named(
    ch,
    c(
      "fluorophore",
      "angle_deg",
      "max_abs_delta",
      "changed",
      "accepted",
      "steps_run"
    )
  )
  expect_equal(ch$changed, c(FALSE, TRUE, FALSE))
  expect_equal(ch$accepted, c(FALSE, TRUE, NA))
  expect_equal(ch$steps_run, rep("signatures", 3))
  expect_equal(ch$max_abs_delta[2], 0.1)
  expect_equal(res$fit_logs$autospectral$signatures$fit.log$iter, c(1L, 1L))
  expect_equal(res$note, "rerun spectral variants with corrected spectra")
  expect_equal(res$provenance$autospectral_version, "1.8.8-mock")
  expect_equal(res$provenance$steps, "signatures")
  expect_setequal(
    res$provenance$call_args$signatures,
    c("spectra", "asp", "min.gain", "raw.data")
  )
  expect_output(print(res), "changed: 1; accepted: 1; accepted NA: 1")
  expect_output(print(res), "rerun spectral variants")
})

test_that("spillover reads accepted from signature.log by fluorophore name", {
  s <- mk_spectra()
  f <- tmp_fcs()
  fc <- list(a = 1)
  mock_autospectral(function(fn, args) {
    expect_equal(fn, "fix.my.unmix")
    expect_equal(args$fully.stained.sample, c(f, f))
    out <- args$spectra
    out["FITC", "d2"] <- 0.05
    list(
      spectra = out,
      signature.log = data.frame(
        fluorophore = c("PE", "FITC"),
        accepted = c(FALSE, TRUE)
      )
    )
  })
  res <- sw_unmix_correct_spillover(
    c(f, f),
    s,
    list(),
    unstained.sample = f,
    flow.control = fc,
    figures = FALSE
  )
  expect_equal(res$changes$accepted, c(TRUE, FALSE, NA))
  expect_equal(res$changes$changed, c(TRUE, FALSE, FALSE))
  expect_equal(res$changes$steps_run, rep("spillover", 3))
  expect_true(is.list(res$fit_logs$autospectral$spillover))
})

test_that("two-step run chains spectra and combines accepted", {
  s <- mk_spectra()
  f <- tmp_fcs()
  fc <- list(a = 1)
  calls <- character()
  mock_autospectral(function(fn, args) {
    calls <<- c(calls, fn)
    out <- args$spectra
    if (fn == "correct.unmixing.signatures") {
      out["PE", "d3"] <- 0.45
      list(spectra = out, accepted = c(PE = TRUE, FITC = FALSE))
    } else {
      # Step two must see step one's spectra.
      expect_equal(args$spectra["PE", "d3"], 0.45)
      out["PE", "d3"] <- 0.4
      list(
        spectra = out,
        signature.log = data.frame(
          fluorophore = c("FITC", "PE"),
          accepted = c(FALSE, TRUE)
        )
      )
    }
  })
  res <- sw_unmix_correct(
    f,
    s,
    list(),
    steps = c("signatures", "spillover"),
    signatures_args = list(unstained.sample = f, flow.control = fc),
    spillover_args = list(unstained.sample = f, flow.control = fc, save = FALSE)
  )
  expect_equal(calls, c("correct.unmixing.signatures", "fix.my.unmix"))
  expect_equal(res$spectra["PE", "d3"], 0.4)
  expect_equal(res$changes$max_abs_delta[2], 0.1)
  expect_equal(res$changes$accepted, c(FALSE, TRUE, NA))
  expect_equal(res$changes$steps_run, rep("signatures+spillover", 3))
  expect_named(res$fit_logs$autospectral, c("signatures", "spillover"))
  expect_true("save" %in% res$provenance$call_args$spillover)
  expect_output(print(res), "signatures -> spillover")
})
