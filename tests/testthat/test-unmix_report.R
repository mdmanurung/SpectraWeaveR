# Tests for R/unmix_report.R. Mock data only; flowCore and AutoSpectral are
# not needed.

# ---- Fixtures ----------------------------------------------------------------

# 3 fluorophores x 4 detectors with known cosines: A-B 1/sqrt(2), C orthogonal.
mock_spectra <- function() {
  matrix(
    c(1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), paste0("D", 1:4))
  )
}

mock_af <- function() {
  matrix(
    c(1, 1, 1, 1, 1, 0, 1, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("AF1", "AF2"), paste0("D", 1:4))
  )
}

mock_variants <- function() {
  s <- mock_spectra()
  list(
    A = s[c(1, 1, 1), , drop = FALSE],
    B = s[c(2, 2), , drop = FALSE]
  )
}

mock_data <- function(n = 100) {
  set.seed(11)
  list(
    s1 = matrix(
      c(rnorm(n, 50, 10), rnorm(n, -5, 1), rnorm(n, 10, 2)),
      ncol = 3,
      dimnames = list(NULL, c("A-A", "B-A", "FSC-A"))
    )
  )
}

status_of <- function(dx, metric, item = NULL) {
  ck <- dx$checks
  keep <- ck$metric == metric
  if (!is.null(item)) {
    keep <- keep & ck$item == item
  }
  ck$status[keep]
}


# ---- Constructor -------------------------------------------------------------

test_that("constructor accepts hand-built pieces", {
  res <- sw_unmix_result(
    spectra = mock_spectra(),
    af_spectra = mock_af(),
    spectra_variants = mock_variants(),
    unmixed_files = character(),
    seed = 42,
    args = list(method = "AutoSpectral"),
    run_id = "r1"
  )
  expect_s3_class(res, "sw_unmix_result")
  expect_identical(res$run_id, "r1")
  expect_named(res$af_spectra, "af")
  expect_named(res$variants, c("A", "B"))
  expect_identical(res$provenance$seed, 42)
  expect_identical(res$provenance$args, list(method = "AutoSpectral"))
  expect_output(print(res), "3 fluorophores x 4 detectors")
})

test_that("constructor wraps a pipeline-shaped list", {
  pipe <- list(
    setup = NULL,
    flow_control = list(),
    spectra = mock_spectra(),
    af_spectra = list(blood = mock_af()),
    spectra_variants = structure(
      list(thresholds = c(A = 1), variants = mock_variants()),
      class = "sw_variants"
    ),
    unmixed = list(output_dir = "x", output_files = character(), method = "OLS")
  )
  res <- sw_unmix_result(pipe, run_id = "p1")
  expect_identical(res$method, "OLS")
  expect_named(res$af_spectra, "blood")
  expect_named(res$variants, c("A", "B"))
  expect_s3_class(res$spectra_variants, "sw_variants")
  expect_error(
    sw_unmix_result(pipe, spectra = mock_spectra()),
    "both in `x` and as an argument: spectra"
  )
  expect_error(sw_unmix_result(list(setup = NULL)), "`spectra` element")
})

test_that("constructor rejects bad shapes and names", {
  s <- mock_spectra()
  expect_error(sw_unmix_result(), "`spectra` is required")
  expect_error(sw_unmix_result(spectra = as.data.frame(s)), "numeric matrix")
  no_rows <- s
  rownames(no_rows) <- NULL
  expect_error(sw_unmix_result(spectra = no_rows), "row names")
  no_cols <- s
  colnames(no_cols) <- NULL
  expect_error(sw_unmix_result(spectra = no_cols), "column names")
  dup <- s
  rownames(dup) <- c("A", "A", "C")
  expect_error(sw_unmix_result(spectra = dup), "duplicated")
  bad <- s
  bad[1, 1] <- NA
  expect_error(sw_unmix_result(spectra = bad), "NA, NaN or Inf")

  af_bad <- mock_af()
  colnames(af_bad)[4] <- "D9"
  expect_error(
    sw_unmix_result(spectra = s, af_spectra = af_bad),
    "detectors do not match"
  )
  expect_error(
    sw_unmix_result(spectra = s, spectra_variants = list(Z = s)),
    "not in `rownames\\(spectra\\)`: Z"
  )
  expect_error(
    sw_unmix_result(spectra = s, af_spectra = list(mock_af())),
    "unique, non-empty names"
  )
  expect_error(sw_unmix_result(spectra = s, seed = "a"), "`seed`")
  expect_error(sw_unmix_result(spectra = s, args = list(1)), "named list")
  expect_error(sw_unmix_result(spectra = s, unmixed_files = 1), "character")
})

test_that("AF and variant matrices are reordered to the spectra detectors", {
  af <- mock_af()[, c(4, 3, 2, 1)]
  res <- sw_unmix_result(spectra = mock_spectra(), af_spectra = af)
  expect_identical(colnames(res$af_spectra$af), paste0("D", 1:4))
  expect_identical(res$af_spectra$af, mock_af())
})

test_that("a bare variant vector counts as zero variants", {
  s <- mock_spectra()
  res <- sw_unmix_result(spectra = s, spectra_variants = list(A = s[1, ]))
  dx <- sw_unmix_diagnose(res)
  ck <- dx$checks[dx$checks$metric == "n_variants", ]
  expect_identical(ck$value[ck$item == "A"], 0)
  expect_identical(ck$status[ck$item == "B"], "not_run")
  expect_match(ck$reason[ck$item == "B"], "No variant matrix")
})


# ---- Provenance --------------------------------------------------------------

test_that("provenance records versions, inventory and baseline", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "in.fcs")
  writeLines("abc", f)
  res <- sw_unmix_result(
    spectra = mock_spectra(),
    inputs = c(f, file.path(dir, "missing.fcs"))
  )
  prov <- res$provenance
  expect_setequal(
    prov$packages$package,
    c("SpectraWeaveR", "AutoSpectral", "AutoSpectralRcpp", "flowCore")
  )
  expect_identical(prov$r_version, R.version.string)
  inv <- prov$files
  expect_identical(inv$md5[inv$path == f], unname(tools::md5sum(f)))
  expect_identical(inv$size_bytes[inv$path == f], as.numeric(file.size(f)))
  expect_false(inv$exists[2])
  expect_true(is.na(inv$md5[2]))
  expect_identical(
    prov$baseline$reference,
    "AutoSpectral 1.8.8 @ 3c504c4560cf1f0fd5dce298eaabaafaebfae9be"
  )
  expect_type(prov$baseline$matches_reference, "logical")
})

test_that("baseline matches only on version and commit together", {
  pin <- "3c504c4560cf1f0fd5dce298eaabaafaebfae9be"
  info <- function(v, s) {
    tibble::tibble(package = "AutoSpectral", version = v, remote_sha = s)
  }
  ok <- .sw_baseline_check(info("1.8.8", pin))
  expect_true(ok$matches_reference)
  absent <- .sw_baseline_check(info(NA_character_, NA_character_))
  expect_false(absent$matches_reference)
  expect_match(absent$message, "not installed")
  old <- .sw_baseline_check(info("1.8.7", pin))
  expect_false(old$matches_reference)
  expect_match(old$message, "1.8.7 differs")
  no_sha <- .sw_baseline_check(info("1.8.8", NA_character_))
  expect_false(no_sha$matches_reference)
  expect_match(no_sha$message, "no RemoteSha")
  other <- .sw_baseline_check(info("1.8.8", "abc"))
  expect_false(other$matches_reference)
  expect_match(other$message, "commit abc differs")
})

test_that("matches_reference is FALSE when AutoSpectral is absent", {
  skip_if(requireNamespace("AutoSpectral", quietly = TRUE))
  res <- sw_unmix_result(spectra = mock_spectra())
  expect_false(res$provenance$baseline$matches_reference)
})


# ---- Diagnostics: numbers ----------------------------------------------------

test_that("max cosine per fluorophore on a hand-made 3x4 matrix", {
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = mock_spectra()))
  ck <- dx$checks[dx$checks$metric == "max_cosine", ]
  expect_equal(ck$value, c(1 / sqrt(2), 1 / sqrt(2), 0))
  expect_identical(ck$item, c("A", "B", "C"))
  expect_match(ck$reason[1], "Closest spectrum: B")
  expect_match(ck$reason[2], "Closest spectrum: A")
  expect_true(all(ck$status == "not_judged"))
})

test_that("condition number matches a known matrix and base kappa", {
  s <- matrix(
    c(3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), paste0("D", 1:4))
  )
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = s))
  v <- dx$checks$value[dx$checks$metric == "condition_number"]
  expect_equal(v, 3)
  expect_equal(v, kappa(s, exact = TRUE))
})

test_that("condition number agrees with AutoSpectral up to its rounding", {
  skip_if_not_installed("AutoSpectral")
  s <- mock_spectra() + 0.1
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = s))
  v <- dx$checks$value[dx$checks$metric == "condition_number"]
  expect_equal(round(v, 2), sw_unmix_condition_number(s))
})

test_that("AF library count and minimum pairwise cosine", {
  dx <- sw_unmix_diagnose(sw_unmix_result(
    spectra = mock_spectra(),
    af_spectra = mock_af()
  ))
  ck <- dx$checks[dx$checks$check_id == "af_library", ]
  expect_identical(ck$value[ck$metric == "n_af_spectra"], 2)
  expect_equal(ck$value[ck$metric == "min_af_cosine"], 2 / (2 * sqrt(2)))
})

test_that("variant counts exclude the reference row", {
  dx <- sw_unmix_diagnose(sw_unmix_result(
    spectra = mock_spectra(),
    spectra_variants = mock_variants()
  ))
  ck <- dx$checks[dx$checks$metric == "n_variants", ]
  expect_identical(ck$value, c(2, 1, NA))
  expect_identical(ck$status, c("not_judged", "not_judged", "not_run"))
})

test_that("per-file statistics match on fluorophore or fluorophore-A", {
  dat <- mock_data()
  dx <- sw_unmix_diagnose(
    sw_unmix_result(spectra = mock_spectra()),
    unmixed_data = dat
  )
  ck <- dx$checks
  expect_identical(ck$value[ck$metric == "n_events"], 100)
  a <- dat$s1[, "A-A"]
  expect_equal(ck$value[ck$metric == "mean" & ck$item == "s1 | A-A"], mean(a))
  expect_equal(ck$value[ck$metric == "sd" & ck$item == "s1 | A-A"], sd(a))
  expect_equal(
    ck$value[ck$metric == "cv" & ck$item == "s1 | A-A"],
    sd(a) / abs(mean(a))
  )
  b <- dat$s1[, "B-A"]
  expect_equal(
    ck$value[ck$metric == "pct_negative" & ck$item == "s1 | B-A"],
    100 * mean(b < 0)
  )
  # C has no channel, and FSC-A is not a fluorophore channel.
  miss <- ck[ck$item == "s1 | C", ]
  expect_identical(miss$status, "not_run")
  expect_match(miss$reason, "No channel named C or C-A")
  expect_false(any(grepl("FSC", ck$item)))
})


# ---- Diagnostics: limits and status ------------------------------------------

test_that("no limit gives not_judged with the value kept", {
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = mock_spectra()))
  run <- dx$checks[dx$checks$status != "not_run", ]
  expect_true(all(run$status == "not_judged"))
  expect_true(all(is.na(run$limit)))
  expect_false(anyNA(run$value))
  expect_identical(
    dx$summary$status,
    c("pass", "warn", "fail", "not_judged", "not_run")
  )
  expect_identical(sum(dx$summary$n), nrow(dx$checks))
  expect_true(any(grepl("No cutoff is built in", dx$notes)))
  expect_output(print(dx), "nothing judged")
})

test_that("max limits judge at the boundaries; equal passes", {
  s <- matrix(
    c(3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), paste0("D", 1:4))
  )
  res <- sw_unmix_result(spectra = s)
  judge <- function(lim) {
    status_of(
      sw_unmix_diagnose(res, limits = list(max_condition_number = lim)),
      "condition_number"
    )
  }
  expect_identical(judge(c(warn = 3, fail = 4)), "pass")
  expect_identical(judge(c(warn = 2, fail = 3)), "warn")
  expect_identical(judge(c(warn = 1, fail = 2.99)), "fail")
  expect_identical(judge(c(warn = 2.99)), "warn")
  expect_identical(judge(c(fail = 3)), "pass")
  expect_identical(judge(2.99), "fail")
  dx <- sw_unmix_diagnose(
    res,
    limits = list(max_condition_number = c(warn = 2, fail = 4))
  )
  expect_identical(
    dx$checks$limit[dx$checks$metric == "condition_number"],
    "warn if > 2; fail if > 4"
  )
  # Limits on one metric leave other metrics unjudged.
  expect_true(all(status_of(dx, "max_cosine") == "not_judged"))
})

test_that("max_cosine limits judge each fluorophore", {
  res <- sw_unmix_result(spectra = mock_spectra())
  dx <- sw_unmix_diagnose(
    res,
    limits = list(max_cosine = c(warn = 0.5, fail = 0.9))
  )
  expect_identical(status_of(dx, "max_cosine"), c("warn", "warn", "pass"))
  dx <- sw_unmix_diagnose(res, limits = list(max_cosine = c(fail = 0.7)))
  expect_identical(status_of(dx, "max_cosine"), c("fail", "fail", "pass"))
})

test_that("min limits judge at the boundaries; equal passes", {
  res <- sw_unmix_result(spectra = mock_spectra())
  judge <- function(lim) {
    status_of(
      sw_unmix_diagnose(
        res,
        limits = list(min_events = lim),
        unmixed_data = mock_data()
      ),
      "n_events"
    )
  }
  expect_identical(judge(c(warn = 100, fail = 50)), "pass")
  expect_identical(judge(c(warn = 101, fail = 100)), "warn")
  expect_identical(judge(c(warn = 200, fail = 101)), "fail")
  expect_identical(judge(101), "fail")
})

test_that("an NA value is not judged even with a limit", {
  s <- mock_spectra()
  s["C", ] <- 0
  dx <- sw_unmix_diagnose(
    sw_unmix_result(spectra = s),
    limits = list(max_cosine = 0.9)
  )
  ck <- dx$checks[dx$checks$metric == "max_cosine" & dx$checks$item == "C", ]
  expect_identical(ck$status, "not_judged")
  expect_match(ck$reason, "limit was not applied")
})

test_that("bad limits error", {
  res <- sw_unmix_result(spectra = mock_spectra())
  expect_error(
    sw_unmix_diagnose(res, limits = list(max_cosin = 0.9)),
    "Unknown limit name\\(s\\): max_cosin"
  )
  expect_error(sw_unmix_diagnose(res, limits = list(0.9)), "must be named")
  expect_error(
    sw_unmix_diagnose(res, limits = c(max_cosine = 0.9)),
    "named list"
  )
  expect_error(
    sw_unmix_diagnose(
      res,
      limits = list(max_cosine = c(warn = 0.99, fail = 0.9))
    ),
    "warn must be <= fail"
  )
  expect_error(
    sw_unmix_diagnose(
      res,
      limits = list(min_events = c(warn = 10, fail = 100))
    ),
    "warn must be >= fail"
  )
  expect_error(
    sw_unmix_diagnose(res, limits = list(max_cosine = c(0.9, 0.95))),
    "named `warn` and/or `fail`"
  )
  expect_error(
    sw_unmix_diagnose(res, limits = list(max_cosine = NA_real_)),
    "non-missing"
  )
  expect_error(sw_unmix_diagnose(list()), "sw_unmix_result")
  expect_error(
    sw_unmix_diagnose(res, unmixed_data = list(matrix(1))),
    "uniquely named"
  )
})


# ---- Diagnostics: not_run reasons --------------------------------------------

test_that("missing inputs give not_run with a reason", {
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = mock_spectra()))
  ck <- dx$checks
  reason <- function(id) unique(ck$reason[ck$check_id == id])
  expect_true(all(ck$status[ck$check_id == "controls"] == "not_run"))
  expect_match(reason("controls"), "No setup")
  expect_match(reason("af_library"), "No AF spectra")
  expect_match(reason("variants"), "No spectral variants")
  expect_match(reason("unmixed_file"), "No unmixed files or data")

  one <- mock_spectra()[1, , drop = FALSE]
  dx1 <- sw_unmix_diagnose(sw_unmix_result(spectra = one))
  expect_match(
    dx1$checks$reason[dx1$checks$metric == "condition_number"],
    "at least 2 spectra"
  )
  expect_identical(
    dx1$checks$status[dx1$checks$metric == "max_cosine"],
    "not_run"
  )
})

test_that("controls are not_run without AutoSpectral", {
  skip_if(requireNamespace("AutoSpectral", quietly = TRUE))
  dir <- withr::local_tempdir()
  cf <- file.path(dir, "controls.csv")
  writeLines("x", cf)
  setup <- structure(
    list(asp = list(), control_file = cf, control_dir = dir),
    class = "sw_setup"
  )
  dx <- sw_unmix_diagnose(sw_unmix_result(
    spectra = mock_spectra(),
    setup = setup
  ))
  ck <- dx$checks[dx$checks$check_id == "controls", ]
  expect_true(all(ck$status == "not_run"))
  expect_match(ck$reason[1], "AutoSpectral is not installed")
})

test_that("unmixed files are not_run without flowCore", {
  skip_if(requireNamespace("flowCore", quietly = TRUE))
  dx <- sw_unmix_diagnose(
    sw_unmix_result(spectra = mock_spectra()),
    unmixed_files = c("a.fcs", "b.fcs")
  )
  ck <- dx$checks[dx$checks$check_id == "unmixed_file", ]
  expect_identical(ck$item, c("a.fcs", "b.fcs"))
  expect_match(ck$reason[1], "flowCore is not installed")
})


# ---- Report ------------------------------------------------------------------

test_that("report writes the files and refuses an existing path", {
  root <- withr::local_tempdir()
  path <- file.path(root, "rep")
  res <- sw_unmix_result(spectra = mock_spectra(), af_spectra = mock_af())
  expect_error(sw_unmix_report(res), "`path` has no default")
  out <- sw_unmix_report(
    res,
    path = path,
    notes = "Not yet run on real FCS data."
  )
  expect_true(all(file.exists(out)))
  expect_setequal(
    basename(out),
    c(
      "checks.csv",
      "spectra.csv",
      "provenance.csv",
      "files.csv",
      "report.json",
      "report.html"
    )
  )
  expect_error(sw_unmix_report(res, path = path), "already exists")
  stale <- file.path(path, "other.txt")
  writeLines("keep", stale)
  expect_silent(sw_unmix_report(
    res,
    path = path,
    overwrite = TRUE,
    notes = "Not yet run on real FCS data."
  ))
  # Overwrite replaces only the report's own files.
  expect_true(file.exists(stale))

  html <- paste(readLines(out[["html"]]), collapse = "\n")
  expect_match(html, "no cutoff is built in", fixed = TRUE)
  expect_match(html, res$provenance$baseline$message, fixed = TRUE)
  expect_match(
    html,
    "Real-data run status: Not yet run on real FCS data.",
    fixed = TRUE
  )
  expect_false(grepl("<script", html, fixed = TRUE))
  expect_false(grepl("<img", html, fixed = TRUE))
  # Section order.
  heads <- regmatches(html, gregexpr("<h2>[^<]+</h2>", html))[[1]]
  expect_identical(
    heads,
    paste0(
      "<h2>",
      c(
        "Provenance and reference check",
        "Inputs",
        "Spectra summary",
        "AF library",
        "Variants",
        "Diagnostics",
        "Limitations",
        "Output files"
      ),
      "</h2>"
    )
  )
})

test_that("report HTML escapes a hostile fluorophore name", {
  s <- mock_spectra()
  rownames(s)[1] <- "<script>alert('x')</script>"
  res <- sw_unmix_result(spectra = s)
  out <- sw_unmix_report(res, path = file.path(withr::local_tempdir(), "rep"))
  html <- paste(readLines(out[["html"]]), collapse = "\n")
  expect_false(grepl("<script>", html, fixed = TRUE))
  expect_match(
    html,
    "&lt;script&gt;alert(&#39;x&#39;)&lt;/script&gt;",
    fixed = TRUE
  )
})

test_that("report refuses a path that overlaps an input folder", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "in.fcs")
  writeLines("abc", f)
  res <- sw_unmix_result(spectra = mock_spectra(), inputs = f)
  expect_error(sw_unmix_report(res, path = file.path(dir, "rep")), "overlaps")
  expect_error(sw_unmix_report(res, path = dir, overwrite = TRUE), "overlaps")
  expect_error(
    sw_unmix_report(res, path = dirname(dir), overwrite = TRUE),
    "overlaps"
  )
  expect_true(file.exists(f))
})

test_that("report CSV files round trip", {
  res <- sw_unmix_result(spectra = mock_spectra(), seed = 7, run_id = "rt")
  dx <- sw_unmix_diagnose(
    res,
    limits = list(max_cosine = c(warn = 0.5)),
    unmixed_data = mock_data()
  )
  out <- sw_unmix_report(
    res,
    dx,
    path = file.path(withr::local_tempdir(), "rep")
  )

  ck <- utils::read.csv(out[["checks"]], stringsAsFactors = FALSE)
  expect_identical(names(ck), names(dx$checks))
  expect_equal(ck$value, dx$checks$value, tolerance = 1e-12)
  expect_identical(ck$status, dx$checks$status)
  expect_identical(ck$item, dx$checks$item)

  sp <- utils::read.csv(out[["spectra"]], stringsAsFactors = FALSE)
  back <- matrix(
    sp$value,
    nrow = 3,
    dimnames = list(unique(sp$fluorophore), unique(sp$detector))
  )
  expect_equal(back, mock_spectra())

  pv <- utils::read.csv(out[["provenance"]], stringsAsFactors = FALSE)
  expect_identical(pv$value[pv$field == "run_id"], "rt")
  expect_identical(pv$value[pv$field == "seed"], "7")
  expect_identical(
    pv$value[pv$field == "matches_reference"],
    as.character(res$provenance$baseline$matches_reference)
  )

  fl <- utils::read.csv(out[["files"]], stringsAsFactors = FALSE)
  expect_identical(names(fl), c("role", "path", "exists", "size_bytes", "md5"))
})

test_that("report JSON parses when jsonlite is present", {
  skip_if_not_installed("jsonlite")
  res <- sw_unmix_result(spectra = mock_spectra(), run_id = "js")
  out <- sw_unmix_report(res, path = file.path(withr::local_tempdir(), "rep"))
  js <- jsonlite::read_json(out[["json"]])
  expect_identical(js$run_id, "js")
  expect_length(js$checks, nrow(sw_unmix_diagnose(res)$checks))
})

test_that("report skips JSON and says so without jsonlite", {
  # Mock in the package namespace, or in the env the file was sourced into.
  env <- environment(sw_unmix_report)
  if (isNamespace(env)) {
    local_mocked_bindings(.sw_has_jsonlite = function() FALSE)
  } else {
    old <- get(".sw_has_jsonlite", envir = env)
    assign(".sw_has_jsonlite", function() FALSE, envir = env)
    withr::defer(assign(".sw_has_jsonlite", old, envir = env))
  }
  res <- sw_unmix_result(spectra = mock_spectra())
  out <- sw_unmix_report(res, path = file.path(withr::local_tempdir(), "rep"))
  expect_false("json" %in% names(out))
  expect_false(file.exists(file.path(dirname(out[[1]]), "report.json")))
  html <- paste(readLines(out[["html"]]), collapse = "\n")
  expect_match(html, "report.json was not written", fixed = TRUE)
})

test_that("report rejects diagnostics from another run", {
  a <- sw_unmix_result(spectra = mock_spectra(), run_id = "a")
  b <- sw_unmix_result(spectra = mock_spectra(), run_id = "b")
  expect_error(
    sw_unmix_report(
      a,
      sw_unmix_diagnose(b),
      path = file.path(withr::local_tempdir(), "rep")
    ),
    "different run_id"
  )
})

test_that("HTML cells show 4 significant digits without scientific notation", {
  expect_identical(
    .sw_cell_text(c(2e6, 0.707106781, 123456.7, NA)),
    c("2000000", "0.7071", "123500", "NA")
  )
})

# ---- SIF metric ----------------------------------------------------------

# Two spectra with cosine c have S = [[1, c], [c, 1]]. Its inverse is
# [[1, -c], [-c, 1]] / (1 - c^2), so the diagonal of the hotspot matrix is
# sqrt(1 / (1 - c^2)) and the off-diagonal is sqrt(c / (1 - c^2)).
test_that("SIF matches the closed form for two spectra", {
  cc <- 0.99
  sp <- rbind(
    A = c(1, 0, 0),
    B = c(cc, sqrt(1 - cc^2), 0),
    C = c(0, 0, 1)
  )
  colnames(sp) <- paste0("D", 1:3)
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = sp))
  ck <- dx$checks[dx$checks$metric == "max_sif", ]
  expect_equal(ck$value[ck$item == "A"], sqrt(1 / (1 - cc^2)), tolerance = 1e-6)
  expect_equal(ck$value[ck$item == "C"], 1, tolerance = 1e-6)
  expect_true(all(ck$status == "not_judged"))
  # The partner is the other member of the close pair.
  expect_match(ck$reason[ck$item == "A"], "B")
})

test_that("AutoSpectral SIF limits judge 4 to 6 as warn, above 6 as fail", {
  expect_identical(sw_unmix_limits_autospectral(), list(max_sif = c(warn = 4, fail = 6)))
  mk <- function(cc) {
    sp <- rbind(A = c(1, 0, 0), B = c(cc, sqrt(1 - cc^2), 0), C = c(0, 0, 1))
    colnames(sp) <- paste0("D", 1:3)
    sw_unmix_diagnose(
      sw_unmix_result(spectra = sp),
      limits = sw_unmix_limits_autospectral()
    )
  }
  # cc = 0.99 gives 7.09 (fail); 0.97 gives 4.11 (warn); 0.9 gives 2.29 (pass).
  expect_identical(status_of(mk(0.99), "max_sif"), c("fail", "fail", "pass"))
  expect_identical(status_of(mk(0.97), "max_sif"), c("warn", "warn", "pass"))
  expect_identical(status_of(mk(0.9), "max_sif"), c("pass", "pass", "pass"))
})

# Parity with AutoSpectral: the SIF is the diagonal of its hotspot matrix.
test_that("max_sif equals the diagonal of AutoSpectral's hotspot matrix", {
  skip_if_not_installed("AutoSpectral")
  # Fixed matrix: 6 spectra by 20 detectors, distinct row names, seed 1.
  set.seed(1)
  spectra <- matrix(abs(rnorm(6 * 20)), nrow = 6, ncol = 20)
  rownames(spectra) <- paste0("F", 1:6)
  colnames(spectra) <- paste0("D", 1:20)

  hotspot <- AutoSpectral::calculate.hotspot.matrix(spectra)
  dx <- sw_unmix_diagnose(sw_unmix_result(spectra = spectra))
  ck <- dx$checks[dx$checks$metric == "max_sif", ]

  expect_equal(
    unname(ck$value[match(rownames(spectra), ck$item)]),
    unname(diag(hotspot)[rownames(spectra)]),
    tolerance = 1e-6
  )
})
