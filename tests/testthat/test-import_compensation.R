# tests/testthat/test-import_compensation.R
# Unit tests for R/import_compensation.R. Writes tiny FCS files to a temp folder.
# flowCore::write.FCS() cannot write "$SPILLOVER" (it expects a matrix), so the
# files use "$SPILL", "$COMP" and a made-up apply key; the function takes any key.

# Write one FCS file with two channels and the given header keywords.
make_comp_fcs <- function(dir, name, kw = list()) {
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(40),
    ncol = 2,
    dimnames = list(NULL, c("B1-A", "V1-A"))
  )
  ff <- flowCore::flowFrame(m)
  desc <- flowCore::keyword(ff)
  for (k in names(kw)) {
    desc[[k]] <- kw[[k]]
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

ID2 <- "2,B1-A,V1-A,1,0,0,1"
REAL2 <- "2,B1-A,V1-A,1,0.12,0.03,1"
MK <- c("$SPILLOVER", "$SPILL", "$COMP")
AK <- c("APPLYCOMP")
# Test vocabulary for the apply flag; the function has no built-in one.
TV <- c("true", "t", "yes", "y", "1")
FV <- c("false", "f", "no", "n", "0")

make_comp_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  c(
    id_yes = make_comp_fcs(
      d,
      "id_yes.fcs",
      list("$SPILL" = ID2, "APPLYCOMP" = "TRUE")
    ),
    id_no = make_comp_fcs(
      d,
      "id_no.fcs",
      list("$SPILL" = ID2, "APPLYCOMP" = "FALSE")
    ),
    real_yes = make_comp_fcs(
      d,
      "real_yes.fcs",
      list("$SPILL" = REAL2, "APPLYCOMP" = "true")
    ),
    real_no = make_comp_fcs(
      d,
      "real_no.fcs",
      list("$SPILL" = REAL2, "APPLYCOMP" = "FALSE")
    ),
    none = make_comp_fcs(d, "none.fcs", list()),
    near = make_comp_fcs(
      d,
      "near.fcs",
      list("$SPILL" = "2,B1-A,V1-A,1,0.0001,0,1")
    ),
    comp_old = make_comp_fcs(
      d,
      "comp_old.fcs",
      list("$COMP" = "2,1,0.2,0,1")
    ),
    both = make_comp_fcs(
      d,
      "both.fcs",
      list("$SPILL" = ID2, "$COMP" = "2,1,0.2,0,1")
    ),
    odd_flag = make_comp_fcs(
      d,
      "odd_flag.fcs",
      list("$SPILL" = ID2, "APPLYCOMP" = "maybe")
    ),
    bad_matrix = make_comp_fcs(
      d,
      "bad_matrix.fcs",
      list("$SPILL" = "2,B1-A,V1-A,1,0,0")
    )
  )
}

run_comp <- function(paths, ...) {
  args <- list(
    paths = unname(paths),
    matrix_keywords = MK,
    apply_keywords = AK,
    apply_true = TV,
    apply_false = FV,
    zero_tol = 1e-6,
    expected_state = NULL
  )
  over <- list(...)
  args[names(over)] <- over
  do.call(sw_fcs_compensation_state, args)
}

row_of <- function(res, f) res[res$file == f, , drop = FALSE]

test_that("every required argument has no default and is named", {
  base <- list(
    paths = "a.fcs",
    matrix_keywords = MK,
    apply_keywords = AK,
    apply_true = TV,
    apply_false = FV,
    zero_tol = 1e-6,
    expected_state = NULL
  )
  for (arg in names(base)) {
    args <- base[setdiff(names(base), arg)]
    expect_error(
      do.call(sw_fcs_compensation_state, args),
      paste0("no default.*'", arg, "'")
    )
  }
  # Several missing at once are all named.
  expect_error(sw_fcs_compensation_state(), "'zero_tol'")
})

test_that("arguments are validated before any file or package is touched", {
  ok <- list(
    paths = c("x/a.fcs", "y/b.fcs"),
    matrix_keywords = MK,
    apply_keywords = AK,
    apply_true = TV,
    apply_false = FV,
    zero_tol = 1e-6,
    expected_state = NULL
  )
  bad <- function(...) {
    a <- ok
    over <- list(...)
    a[names(over)] <- over
    do.call(sw_fcs_compensation_state, a)
  }
  expect_error(bad(paths = character()), "`paths`")
  expect_error(bad(paths = c("x/a.fcs", "y/a.fcs")), "unique")
  expect_error(bad(matrix_keywords = character()), "`matrix_keywords`")
  expect_error(bad(apply_keywords = NA_character_), "`apply_keywords`")
  expect_error(bad(paths = c("a.fcs", NA)), "`paths`")
  expect_error(bad(matrix_keywords = c("$SPILL", "$SPILL")), "unique")
  expect_error(bad(apply_keywords = c("A", "A")), "unique")
  expect_error(bad(apply_true = character()), "`apply_true`")
  expect_error(bad(apply_false = NA_character_), "`apply_false`")
  expect_error(bad(apply_true = c("yes", "No"), apply_false = "no"), "both")
  expect_error(bad(zero_tol = -1), "`zero_tol`")
  expect_error(bad(zero_tol = c(0, 1)), "`zero_tol`")
  expect_error(bad(zero_tol = "0"), "`zero_tol`")
  expect_error(bad(expected_state = "x"), "`expected_state`")
  expect_error(
    bad(expected_state = data.frame(file = "a.fcs")),
    "expected_state"
  )
  expect_error(
    bad(
      expected_state = data.frame(file = "a.fcs", matrix_present = "yes")
    ),
    "logical"
  )
  expect_error(
    bad(
      expected_state = data.frame(
        file = c("a.fcs", "a.fcs"),
        matrix_present = TRUE
      )
    ),
    "once"
  )
})

test_that("matrix, identity and apply flag are separate columns", {
  p <- make_comp_set()
  res <- run_comp(p)
  expect_s3_class(res, "sw_fcs_compensation_state")
  expect_equal(res$file, basename(unname(p)))

  a <- row_of(res, "id_yes.fcs")
  expect_true(a$matrix_present)
  expect_true(a$matrix_identity)
  expect_equal(a$matrix_key, "$SPILL")
  expect_equal(a$matrix_n, 2)
  expect_equal(a$apply_raw, "TRUE")
  expect_true(a$apply_parsed)
  expect_equal(a$apply_key, "APPLYCOMP")

  # Same identity matrix, flag says not applied: the two facts stay apart.
  b <- row_of(res, "id_no.fcs")
  expect_true(b$matrix_identity)
  expect_false(b$apply_parsed)

  # Real matrix, flag says applied (lower-case raw text parses).
  c1 <- row_of(res, "real_yes.fcs")
  expect_false(c1$matrix_identity)
  expect_true(c1$apply_parsed)
  expect_equal(c1$apply_raw, "true")
  expect_equal(c1$matrix_max_dev, 0.12)

  d <- row_of(res, "real_no.fcs")
  expect_false(d$matrix_identity)
  expect_false(d$apply_parsed)

  # No combined verdict and no pass/fail column.
  expect_false(any(grepl("pass|fail|verdict|state$", names(res))))
})

test_that("a file with no matrix or flag says so and keeps its row", {
  res <- run_comp(make_comp_set())
  n <- row_of(res, "none.fcs")
  expect_false(n$matrix_present)
  expect_true(is.na(n$matrix_identity))
  expect_equal(n$matrix_status, "not_present")
  expect_false(n$apply_present)
  expect_true(is.na(n$apply_parsed))
  expect_true(is.na(n$apply_raw))
  expect_equal(n$header_status, "read")
})

test_that("zero_tol decides what counts as identity", {
  p <- make_comp_set()
  loose <- row_of(run_comp(p, zero_tol = 1e-3), "near.fcs")
  tight <- row_of(run_comp(p, zero_tol = 1e-6), "near.fcs")
  expect_true(loose$matrix_identity)
  expect_false(tight$matrix_identity)
  expect_equal(loose$matrix_max_dev, 1e-4)
  # The no-name $COMP layout parses too.
  old <- row_of(run_comp(p), "comp_old.fcs")
  expect_equal(old$matrix_key, "$COMP")
  expect_false(old$matrix_identity)
  expect_equal(old$matrix_max_dev, 0.2)
})

test_that("two matching keys are listed, not ranked", {
  res <- run_comp(make_comp_set())
  b <- row_of(res, "both.fcs")
  expect_true(b$matrix_present)
  expect_equal(b$matrix_keys_matched, "$SPILL; $COMP")
  expect_true(is.na(b$matrix_key))
  expect_true(is.na(b$matrix_identity))
  expect_equal(b$matrix_status, "not_assessed")
  expect_match(b$matrix_reason, "several matrix keys")
  # Order of matrix_keywords does not change the outcome.
  rev <- row_of(
    run_comp(make_comp_set(), matrix_keywords = rev(MK)),
    "both.fcs"
  )
  expect_true(is.na(rev$matrix_identity))
  expect_setequal(
    strsplit(rev$matrix_keys_matched, "; ")[[1]],
    c("$SPILL", "$COMP")
  )
})

test_that("an unrecognised flag keeps its raw text and parses to NA", {
  res <- run_comp(make_comp_set())
  o <- row_of(res, "odd_flag.fcs")
  expect_equal(o$apply_raw, "maybe")
  expect_true(is.na(o$apply_parsed))
  expect_equal(o$apply_status, "not_assessed")
  expect_match(o$apply_reason, "maybe")
  # The matrix of that file is still assessed.
  expect_true(o$matrix_identity)
})

test_that("a malformed matrix is not_assessed with a reason, not a crash", {
  res <- run_comp(make_comp_set())
  m <- row_of(res, "bad_matrix.fcs")
  expect_true(m$matrix_present)
  expect_equal(m$matrix_status, "not_assessed")
  expect_true(is.na(m$matrix_identity))
  expect_match(m$matrix_reason, "expected 7 values")
  expect_equal(m$header_status, "read")
})

test_that("the spillover parser rejects each malformed shape with a reason", {
  chk <- function(x) .sw_parse_spillover(x)
  ok <- chk("2,A,B,1,0,0,1")
  expect_true(ok$ok)
  expect_equal(ok$names, c("A", "B"))
  expect_equal(ok$matrix, diag(2), ignore_attr = TRUE)
  # Row-major order: values 1..4 fill row 1 first.
  rm <- chk("2,A,B,1,2,3,4")
  expect_equal(rm$matrix[1, 2], 2)
  expect_equal(rm$matrix[2, 1], 3)
  # No-name layout.
  un <- chk("2,1,0,0,1")
  expect_true(un$ok)
  expect_null(un$names)
  expect_true(chk(" 1 , A , 1 ")$ok)

  bad <- list(
    empty = "",
    blank = "   ",
    na = NA_character_,
    n_text = "x,A,B,1,0,0,1",
    n_zero = "0",
    n_neg = "-2,A,B,1,0,0,1",
    n_frac = "2.5,A,B,1,0,0,1",
    too_few = "2,A,B,1,0,0",
    too_many = "2,A,B,1,0,0,1,1",
    nonnum = "2,A,B,1,zz,0,1",
    nan = "2,A,B,1,NaN,0,1",
    inf = "2,A,B,1,Inf,0,1",
    trailing = "2,A,B,1,0,0,1,",
    empty_name = "2,A,,1,0,0,1",
    dup_name = "2,A,A,1,0,0,1",
    semicolon = "2;A;B;1;0;0;1"
  )
  for (nm in names(bad)) {
    r <- chk(bad[[nm]])
    expect_false(r$ok, info = nm)
    expect_true(is.character(r$reason) && nzchar(r$reason), info = nm)
    expect_null(r$matrix, info = nm)
  }
  expect_match(chk("2,A,B,1,0,0")$reason, "expected 7 values")
})

test_that("an unreadable or missing file is a not_assessed row", {
  p <- make_comp_set()
  d <- withr::local_tempdir()
  junk <- file.path(d, "junk.fcs")
  writeLines("this is not an fcs file", junk)
  gone <- file.path(d, "gone.fcs")
  res <- run_comp(c(p, junk = junk, gone = gone))
  expect_equal(nrow(res), length(p) + 2)
  for (f in c("junk.fcs", "gone.fcs")) {
    r <- row_of(res, f)
    expect_equal(r$header_status, "not_assessed", info = f)
    expect_true(nzchar(r$header_reason), info = f)
    expect_true(is.na(r$matrix_present), info = f)
    expect_true(is.na(r$matrix_identity), info = f)
    expect_true(is.na(r$apply_parsed), info = f)
    expect_equal(r$matrix_status, "not_assessed", info = f)
  }
  expect_equal(row_of(res, "junk.fcs")$header_reason, "not an FCS file")
  expect_equal(row_of(res, "gone.fcs")$header_reason, "file not found")
  # The good files are unaffected by the bad ones.
  expect_true(row_of(res, "id_yes.fcs")$matrix_identity)
})

test_that("files on disk are never changed", {
  p <- make_comp_set()
  before <- unname(tools::md5sum(p))
  invisible(run_comp(p))
  expect_equal(unname(tools::md5sum(p)), before)
})

test_that("expected_state adds comparison columns and no verdict", {
  p <- make_comp_set()
  exp <- data.frame(
    file = c("id_yes.fcs", "real_yes.fcs", "none.fcs", "odd_flag.fcs"),
    matrix_present = c(TRUE, TRUE, TRUE, TRUE),
    matrix_identity = c(TRUE, TRUE, NA, TRUE),
    apply_parsed = c(TRUE, FALSE, TRUE, TRUE)
  )
  res <- run_comp(p, expected_state = exp)
  expect_true(all(
    c(
      "expected_matrix_present",
      "expected_matrix_identity",
      "expected_apply_parsed",
      "matrix_present_expected_vs_observed",
      "matrix_identity_expected_vs_observed",
      "apply_parsed_expected_vs_observed"
    ) %in%
      names(res)
  ))
  g <- function(f, col) row_of(res, f)[[col]]
  expect_equal(g("id_yes.fcs", "matrix_identity_expected_vs_observed"), "same")
  expect_equal(g("id_yes.fcs", "apply_parsed_expected_vs_observed"), "same")
  expect_equal(
    g("real_yes.fcs", "matrix_identity_expected_vs_observed"),
    "differs"
  )
  expect_equal(
    g("real_yes.fcs", "apply_parsed_expected_vs_observed"),
    "differs"
  )
  # Expected present but absent in the file.
  expect_equal(g("none.fcs", "matrix_present_expected_vs_observed"), "differs")
  # Expectation NA: no expectation, not a match.
  expect_equal(
    g("none.fcs", "matrix_identity_expected_vs_observed"),
    "not_expected"
  )
  # Observed NA (flag unparsed): cannot compare.
  expect_equal(
    g("odd_flag.fcs", "apply_parsed_expected_vs_observed"),
    "not_assessed"
  )
  # A file with no row in expected_state stays, flagged not_expected.
  expect_equal(
    g("id_no.fcs", "matrix_present_expected_vs_observed"),
    "not_expected"
  )
  expect_equal(nrow(res), length(p))
  expect_false(any(grepl("pass|fail|verdict", names(res))))
})

test_that("expected_state rows with no file in paths are kept as rows", {
  p <- make_comp_set()
  exp <- data.frame(
    file = c("id_yes.fcs", "ghost.fcs"),
    matrix_present = c(TRUE, TRUE)
  )
  res <- run_comp(p, expected_state = exp)
  expect_equal(nrow(res), length(p) + 1)
  g <- row_of(res, "ghost.fcs")
  expect_equal(g$header_status, "not_assessed")
  expect_match(g$header_reason, "not in `paths`")
  expect_true(g$expected_matrix_present)
  expect_equal(g$matrix_present_expected_vs_observed, "not_assessed")
})

test_that("NULL expected_state adds no expectation columns", {
  res <- run_comp(make_comp_set(), expected_state = NULL)
  expect_false(any(grepl("expected", names(res))))
})

test_that("the result carries provenance", {
  res <- run_comp(make_comp_set(), zero_tol = 1e-4)
  pv <- attr(res, "provenance")
  expect_true(is.list(pv))
  expect_true(all(
    c("r_version", "datetime", "packages", "seed", "args") %in% names(pv)
  ))
  expect_equal(pv$args$zero_tol, 1e-4)
  expect_equal(pv$args$matrix_keywords, MK)
  expect_true(is.null(pv$seed))
  expect_true("flowCore" %in% pv$packages$package)
})

test_that("the apply-flag vocabulary is the caller's, case and padding ignored", {
  f <- function(x, t = TV, fl = FV) .sw_parse_flag(x, t, fl)
  for (x in c("1", "yes", "Y", "t", " TRUE ")) {
    expect_true(f(x), info = x)
  }
  for (x in c("0", "no", "N", "f", " False ")) {
    expect_false(f(x), info = x)
  }
  expect_true(is.na(f("maybe")))
  # A narrower vocabulary changes the reading: "1" is not assumed to mean TRUE.
  expect_true(is.na(f("1", t = "true", fl = "false")))
  expect_true(f("on", t = c("ON", "true"), fl = "off"))
  # Through the public call, with a vendor word the default list would not know.
  p <- make_comp_set()
  r <- row_of(
    run_comp(p, apply_true = "applied", apply_false = "none"),
    "id_yes.fcs"
  )
  expect_equal(r$apply_raw, "TRUE")
  expect_true(is.na(r$apply_parsed))
  expect_equal(r$apply_status, "not_assessed")
})

test_that("identity uses absolute distance and a closed boundary", {
  chk <- function(m, tol) {
    d <- withr::local_tempdir(.local_envir = parent.frame())
    p <- make_comp_fcs(d, "m.fcs", list("$SPILL" = m))
    run_comp(p, zero_tol = tol)
  }
  neg <- chk("2,B1-A,V1-A,1,-0.2,0,1", 1e-6)
  expect_false(neg$matrix_identity)
  expect_equal(neg$matrix_max_dev, 0.2)
  # Exactly at the tolerance counts as identity (<=).
  edge <- chk("2,B1-A,V1-A,1,0.5,0,1", 0.5)
  expect_equal(edge$matrix_max_dev, 0.5)
  expect_true(edge$matrix_identity)
})

test_that("repeated keys are refused, padded header values are trimmed", {
  p <- make_comp_set()
  expect_error(run_comp(p, matrix_keywords = c("$SPILL", "$SPILL")), "unique")
  skip_if_not_installed("testthat", "3.2.0")
  testthat::local_mocked_bindings(
    read.FCSheader = function(...) {
      list(c(
        "$PAR" = "2",
        "$SPILLOVER" = paste0("  ", ID2, " "),
        "APPLYCOMP" = "  TRUE  "
      ))
    },
    .package = "flowCore"
  )
  r <- run_comp(p[["id_yes"]])
  # Real key name of the FCS 3.1 standard, which write.FCS() cannot write.
  expect_equal(r$matrix_key, "$SPILLOVER")
  expect_true(r$matrix_identity)
  expect_equal(r$apply_raw, "TRUE")
  expect_true(r$apply_parsed)
})

test_that("a huge channel count is a not_assessed row, not a crash", {
  for (x in c("100000,A", "1e5,A", "1e9,A,B", "50000,A")) {
    r <- .sw_parse_spillover(x)
    expect_false(r$ok, info = x)
    expect_match(r$reason, "expected", info = x)
  }
  d <- withr::local_tempdir()
  big <- make_comp_fcs(d, "big.fcs", list("$SPILL" = "1000000000,A,B"))
  good <- make_comp_fcs(d, "good.fcs", list("$SPILL" = ID2))
  res <- run_comp(c(big, good))
  expect_equal(row_of(res, "big.fcs")$matrix_status, "not_assessed")
  expect_true(row_of(res, "good.fcs")$matrix_identity)
})

test_that("a truncated file starting with FCS leaks no connection", {
  d <- withr::local_tempdir()
  bad <- file.path(d, "bad.fcs")
  writeLines(paste0("FCS3.1", strrep(" ", 20)), bad)
  short <- file.path(d, "short.fcs")
  writeBin(charToRaw("FCS.."), short)
  before <- length(getAllConnections())
  res <- run_comp(c(bad, short))
  expect_equal(length(getAllConnections()), before)
  # A leaked connection shows up as "closing unused connection" at gc time.
  expect_no_warning(gc())
  expect_equal(res$header_status, c("not_assessed", "not_assessed"))
  # Reason is one line and does not repeat the temp path.
  expect_false(any(grepl("\n", res$header_reason)))
  expect_false(any(grepl(d, res$header_reason, fixed = TRUE)))
})

test_that("a connection left open by a failed header read is closed", {
  skip_if_not_installed("testthat", "3.2.0")
  p <- make_comp_set()[["id_yes"]]
  # Stand-in for flowCore: opens the file, then fails without closing it.
  testthat::local_mocked_bindings(
    read.FCSheader = function(files, ...) {
      con <- file(files, "rb")
      stop("boom at ", files)
    },
    .package = "flowCore"
  )
  before <- length(getAllConnections())
  res <- run_comp(p)
  expect_equal(length(getAllConnections()), before)
  expect_equal(res$header_status, "not_assessed")
  expect_match(res$header_reason, "boom")
  expect_false(grepl(p, res$header_reason, fixed = TRUE))
})

test_that("a combined `reason` column joins the three specific reasons", {
  p <- file.path(withr::local_tempdir(), "gone.fcs")
  res <- run_comp(p)
  expect_true("reason" %in% names(res))
  expect_true(all(
    c("header_reason", "matrix_reason", "apply_reason") %in% names(res)
  ))
  expect_match(res$reason, res$header_reason, fixed = TRUE)
  expect_match(res$reason, "header not read", fixed = TRUE)
})

test_that("one status and reason per row: not_assessed when any part is", {
  d <- make_comp_set()
  res <- run_comp(d)
  expect_true(all(c("status", "reason") %in% names(res)))
  expect_true(all(res$status %in% c("assessed", "not_assessed")))
  expect_equal(res$status[res$file == "id_yes.fcs"], "assessed")
  expect_true(is.na(res$reason[res$file == "id_yes.fcs"]))
  # A file that is not a readable FCS keeps its row, as not_assessed with a reason
  bad <- withr::local_tempfile(fileext = ".fcs")
  writeLines("not an fcs file", bad)
  res2 <- run_comp(c(d, bad = bad))
  b <- res2[res2$file == basename(bad), ]
  expect_equal(b$status, "not_assessed")
  expect_false(is.na(b$reason))
})

test_that("provenance starts with the shared keys (compensation)", {
  pv <- attr(run_comp(make_comp_set(), zero_tol = 1e-4), "provenance")
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
