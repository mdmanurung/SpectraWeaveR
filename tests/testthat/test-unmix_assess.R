# Tests for R/unmix_assess.R. AutoSpectral is not required: argument checks run first.

mk_mat <- function(rows = c("A", "B", "C"), cols = c("d1", "d2", "d3", "d4")) {
  m <- matrix(
    runif(length(rows) * length(cols)),
    nrow = length(rows),
    dimnames = list(rows, cols)
  )
  m
}
mk_variants <- function() {
  list(A = mk_mat(c("ref", "v1")), B = mk_mat(c("ref", "v1")))
}
tmp_fcs <- function() {
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  f
}

test_that("helper error text for missing AutoSpectral is exact", {
  skip_if(
    requireNamespace("AutoSpectral", quietly = TRUE),
    "AutoSpectral is installed"
  )
  expect_error(
    .sw_assess_require_autospectral(),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_condition_number(mk_mat()),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_assess_mismatch(mk_variants(), mk_variants(), metric = "angle"),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_assess_variability(mk_variants(), what = "mad"),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_assess_alignment(mk_mat(), mk_mat()),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
  expect_error(
    sw_unmix_assess_mismatch_clusters(mk_mat(), mk_mat()),
    "Package 'AutoSpectral' is required.",
    fixed = TRUE
  )
})

test_that("sw_unmix_assess_mismatch validates arguments", {
  v <- mk_variants()
  expect_error(sw_unmix_assess_mismatch(1, v, metric = "cosine"), "named list")
  expect_error(
    sw_unmix_assess_mismatch(unname(v), v, metric = "cosine"),
    "named list"
  )
  expect_error(
    sw_unmix_assess_mismatch(list(A = 1), v, metric = "cosine"),
    "detector column names"
  )
  expect_error(sw_unmix_assess_mismatch(v, v), "`metric`")
  expect_error(sw_unmix_assess_mismatch(v, v, metric = "euclid"), "`metric`")
  other <- list(Z = mk_mat(c("ref", "v1")))
  expect_error(
    sw_unmix_assess_mismatch(v, other, metric = "cosine"),
    "share no fluorophore"
  )
})

test_that("sw_unmix_assess_mismatch_clusters validates arguments", {
  m <- mk_mat()
  expect_error(
    sw_unmix_assess_mismatch_clusters(as.data.frame(m), m),
    "numeric matrix"
  )
  expect_error(sw_unmix_assess_mismatch_clusters(unname(m), m), "row names")
  other <- mk_mat(c("X", "Y"), c("e1", "e2"))
  expect_error(
    sw_unmix_assess_mismatch_clusters(m, other),
    "share no row or column"
  )
})

test_that("sw_unmix_assess_variability validates arguments", {
  v <- mk_variants()
  expect_error(sw_unmix_assess_variability(v), "`what`")
  expect_error(sw_unmix_assess_variability(v, what = "x"), "`what`")
  expect_error(
    sw_unmix_assess_variability(v, what = "mad", peak_height = 2),
    "peak_height"
  )
  expect_error(
    sw_unmix_assess_variability(v, what = "mad", peak_height = "a"),
    "peak_height"
  )
  expect_error(sw_unmix_assess_variability(list(), what = "mad"), "named list")
})

test_that("sw_unmix_assess_alignment validates arguments", {
  m <- mk_mat()
  expect_error(sw_unmix_assess_alignment(m, 1:3), "numeric matrix")
  expect_error(
    sw_unmix_assess_alignment(m, mk_mat(c("X", "Y"))),
    "share no fluorophore"
  )
})

test_that("sw_unmix_condition_number validates arguments", {
  expect_error(sw_unmix_condition_number(1:3), "numeric matrix")
  expect_error(sw_unmix_condition_number(unname(mk_mat())), "row names")
  expect_error(sw_unmix_condition_number(mk_mat("A")), "at least 2")
  bad <- mk_mat()
  bad[1, 1] <- NA
  expect_error(sw_unmix_condition_number(bad), "NA, NaN or Inf")
})

test_that("sw_unmix_check_controls validates arguments", {
  f <- tmp_fcs()
  d <- tempdir()
  expect_error(
    sw_unmix_check_controls(file.path(d, "nope"), f, list()),
    "not a directory"
  )
  expect_error(
    sw_unmix_check_controls(d, file.path(d, "nope.csv"), list()),
    "does not exist"
  )
  expect_error(sw_unmix_check_controls(d, f, "asp"), "parameter list")
  expect_error(
    sw_unmix_check_controls(d, f, list(), spectra = mk_mat("A")),
    "at least 2"
  )
  expect_error(
    sw_unmix_check_controls(d, f, list(), gate_list = list(a = 1)),
    "together"
  )
  expect_error(
    sw_unmix_check_controls(d, f, list(), spectral_channels = 1),
    "character vector"
  )
})

test_that("sw_unmix_compare validates arguments", {
  f <- tmp_fcs()
  s <- mk_mat()
  sp <- c(d1 = 0.1, d2 = 0.2)
  call_with <- function(...) {
    args <- modifyList(
      list(
        single_stained_fcs = f,
        unstained_fcs = f,
        fluorophore = "A",
        spectra = s,
        ref_spectrum = sp,
        test_spectrum = sp,
        cytometer = "aurora",
        plot_dir = tempfile("figs")
      ),
      list(...)
    )
    do.call(sw_unmix_compare, args)
  }
  expect_error(call_with(single_stained_fcs = "missing.fcs"), "does not exist")
  expect_error(call_with(fluorophore = c("A", "B")), "single non-empty string")
  expect_error(call_with(spectra = 1), "numeric matrix")
  expect_error(call_with(ref_spectrum = c(0.1, 0.2)), "detector names")
  expect_error(
    call_with(test_spectrum = c(zz = 0.1)),
    "not columns of `spectra`"
  )
  expect_error(call_with(cytometer = NA_character_), "single non-empty string")
})

test_that("sw_unmix_compare_data validates arguments", {
  a <- data.frame(x = 1:3, y = 1:3, z = 1:3)
  expect_error(
    sw_unmix_compare_data(list(a), c("x", "y"), list()),
    "two or more"
  )
  expect_error(
    sw_unmix_compare_data(list(a, 1), c("x", "y"), list()),
    "column names"
  )
  expect_error(sw_unmix_compare_data(list(a, a), "x", list()), "at least two")
  expect_error(
    sw_unmix_compare_data(list(a, a), c("x", "q"), list()),
    "Fewer than two"
  )
  expect_error(
    sw_unmix_compare_data(list(a, a), c("x", "y"), "asp"),
    "parameter list"
  )
})

test_that("sw_unmix_compare_af validates arguments", {
  f <- tmp_fcs()
  s <- mk_mat()
  af <- mk_mat(c("af1", "af2"))
  expect_error(
    sw_unmix_compare_af("missing.fcs", s, list(af)),
    "does not exist"
  )
  expect_error(sw_unmix_compare_af(f, 1, list(af)), "numeric matrix")
  expect_error(sw_unmix_compare_af(f, s, list()), "non-empty list")
  expect_error(
    sw_unmix_compare_af(f, s, list(mk_mat(cols = c("q1", "q2")))),
    "columns match"
  )
})

test_that("real AutoSpectral calls return tibbles", {
  skip_if_not_installed("AutoSpectral")
  v <- list(
    A = rbind(
      ref = c(d1 = 1, d2 = 0.5, d3 = 0.1),
      v1 = c(d1 = 1, d2 = 0.4, d3 = 0.1)
    ),
    B = rbind(
      ref = c(d1 = 0.1, d2 = 0.5, d3 = 1),
      v1 = c(d1 = 0.1, d2 = 0.6, d3 = 1)
    )
  )
  cos_tbl <- sw_unmix_assess_mismatch(v, v, metric = "cosine")
  expect_s3_class(cos_tbl, "tbl_df")
  expect_named(cos_tbl, c("fluorophore", "cosine_similarity"))
  ang_tbl <- sw_unmix_assess_mismatch(v, v, metric = "angle")
  expect_named(ang_tbl, c("fluorophore", "angle_degrees"))
  # identical inputs give ~0 degrees (AutoSpectral returns 0.0023 from acos rounding)
  expect_true(all(ang_tbl$angle_degrees < 0.01))

  var_tbl <- sw_unmix_assess_variability(v, what = "summary")
  expect_true("fluorophore" %in% names(var_tbl))
  mad_tbl <- sw_unmix_assess_variability(v, what = "mad")
  # two toy variants are too few for the MAD routine: expect a tibble with a
  # fluorophore column, not a fixed width
  expect_s3_class(mad_tbl, "tbl_df")
  expect_true("fluorophore" %in% names(mad_tbl))

  spectra <- rbind(
    A = c(d1 = 1, d2 = 0.5, d3 = 0.1),
    B = c(d1 = 0.1, d2 = 0.5, d3 = 1)
  )
  expect_type(sw_unmix_condition_number(spectra), "double")
})

test_that("sw_unmix_compare requires plot_dir", {
  expect_error(
    sw_unmix_compare(
      "a.fcs", "b.fcs", "FITC", matrix(1, 1, 1, dimnames = list("FITC", "D1")),
      ref_spectrum = "FITC", test_spectrum = "FITC", cytometer = "Aurora"
    ),
    "plot_dir"
  )
})
