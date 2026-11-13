# tests/testthat/test-import_integrity.R
# Unit tests for R/import_integrity.R. Writes tiny FCS files to a temp folder.

# Write one FCS file with a Time channel, one marker channel and chosen keywords.
make_int_fcs <- function(dir, name, n = 20, kw = list()) {
  skip_if_not_installed("flowCore")
  m <- cbind(Time = seq_len(n), "B1-A" = stats::rnorm(n))
  ff <- flowCore::flowFrame(m)
  flowCore::pData(flowCore::parameters(ff))$desc <- c("Time", "CD3")
  desc <- flowCore::keyword(ff)
  for (k in names(kw)) {
    desc[[k]] <- kw[[k]]
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

# Copy a 20-row file and change its $TOT text from 20 to 25 (same width).
# The audit must report this; it must never do it.
make_bad_tot <- function(src, dir, name) {
  b <- readBin(src, "raw", file.size(src))
  pat <- charToRaw("|$TOT|20|")
  at <- which(vapply(
    seq_len(length(b) - length(pat) + 1),
    function(i) all(b[i:(i + length(pat) - 1)] == pat),
    logical(1)
  ))
  stopifnot(length(at) == 1)
  b[at + 7] <- charToRaw("5")
  out <- file.path(dir, name)
  writeBin(b, out)
  out
}

make_int_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  kw <- list("$BTIM" = "10:00:00", "$ETIM" = "10:05:30")
  good <- make_int_fcs(d, "good.fcs", kw = kw)
  c(
    good = good,
    bad = make_bad_tot(good, d, "bad_tot.fcs"),
    junk = {
      writeLines("not an fcs file", p <- file.path(d, "junk.fcs"))
      p
    },
    gone = file.path(d, "gone.fcs")
  )
}

audit_ok <- function(paths, ...) {
  args <- list(
    paths = unname(paths),
    tot_tolerance = 0,
    time_keyword = "Time",
    header_only = FALSE
  )
  args[names(list(...))] <- list(...)
  do.call(sw_fcs_integrity_audit, args)
}

test_that("every required argument stops with 'no default' and is named", {
  expect_error(sw_fcs_integrity_audit(), "no default")
  full <- list(
    paths = "a.fcs",
    tot_tolerance = 0,
    time_keyword = "Time",
    header_only = FALSE
  )
  for (arg in names(full)) {
    expect_error(
      do.call(sw_fcs_integrity_audit, full[setdiff(names(full), arg)]),
      paste0("no default.*'", arg, "'")
    )
  }
})

test_that("bad argument values stop before any file or package is touched", {
  expect_error(audit_ok("a.fcs", tot_tolerance = -1), "tot_tolerance")
  expect_error(audit_ok("a.fcs", tot_tolerance = NA_real_), "tot_tolerance")
  expect_error(audit_ok("a.fcs", header_only = NA), "header_only")
  expect_error(audit_ok("a.fcs", time_keyword = character()), "time_keyword")
  expect_error(audit_ok(character()), "paths")
  expect_error(audit_ok(c("a/x.fcs", "b/x.fcs")), "unique")
})

test_that("a clean file gets one row per item with the measured numbers", {
  p <- make_int_set()
  res <- audit_ok(p[["good"]], time_keyword = "Time")
  expect_s3_class(res, "sw_fcs_integrity_audit")
  expect_equal(
    res$item,
    c("readable", "tot_vs_rows", "time_span", "channel_set")
  )
  expect_true(all(res$file == "good.fcs"))
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$value, 0)
  expect_equal(tot$status, "match")
  span <- res[res$item == "time_span", ]
  expect_equal(span$value, 19) # channel Time runs 1..20
  chans <- res[res$item == "channel_set", ]
  expect_equal(chans$value, 2)
  expect_match(chans$detail, "B1-A")
  expect_match(chans$detail, "CD3")
  # No pass or fail column and no per-file score.
  expect_false(any(c("pass", "fail", "score", "ok") %in% names(res)))
})

test_that("a header keyword pair gives the time span in seconds", {
  p <- make_int_set()
  res <- audit_ok(p[["good"]], time_keyword = c("$BTIM", "$ETIM"))
  span <- res[res$item == "time_span", ]
  expect_equal(span$value, 330)
  expect_equal(span$status, "assessed")
})

test_that("a time channel or keyword that is absent is not_assessed with a reason", {
  p <- make_int_set()
  res <- audit_ok(p[["good"]], time_keyword = "NoSuchChannel")
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "not_assessed")
  expect_match(span$reason, "NoSuchChannel")
  expect_true(is.na(span$value))
})

test_that("a $TOT mismatch is reported and the file is never edited", {
  p <- make_int_set()
  before <- unname(tools::md5sum(p[["bad"]]))
  res <- audit_ok(p[["bad"]])
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$value, 5) # header says 25, 20 rows read
  expect_equal(tot$status, "differs")
  expect_match(tot$detail, "25")
  # flowCore returns no events for this file, so the channel time span cannot
  # be computed; that must show as a reason, not a missing row.
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "not_assessed")
  expect_match(span$reason, "TOT")
  expect_equal(res$value[res$item == "readable"], 1)
  expect_equal(unname(tools::md5sum(p[["bad"]])), before)
  expect_equal(flowCore::read.FCSheader(p[["bad"]])[[1]][["$TOT"]], "25")
})

test_that("tot_tolerance decides match or differs, and is recorded", {
  p <- make_int_set()
  loose <- audit_ok(p[["bad"]], tot_tolerance = 5)
  expect_equal(loose$status[loose$item == "tot_vs_rows"], "match")
  expect_equal(loose$value[loose$item == "tot_vs_rows"], 5)
  strict <- audit_ok(p[["bad"]], tot_tolerance = 4)
  expect_equal(strict$status[strict$item == "tot_vs_rows"], "differs")
  expect_equal(attr(strict, "provenance")$args$tot_tolerance, 4)
})

test_that("header_only = TRUE skips the $TOT check and says so", {
  p <- make_int_set()
  res <- audit_ok(p[["bad"]], header_only = TRUE)
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_true(is.na(tot$value))
  expect_match(tot$reason, "header_only")
  expect_true(attr(res, "header_only"))
  # A time channel needs events, so it is not_assessed too, with the same reason.
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "not_assessed")
  expect_match(span$reason, "header_only")
  # Header-only items still work.
  expect_equal(res$status[res$item == "channel_set"], "assessed")
  keyed <- audit_ok(
    p[["good"]],
    header_only = TRUE,
    time_keyword = c("$BTIM", "$ETIM")
  )
  expect_equal(keyed$value[keyed$item == "time_span"], 330)
  expect_false(attr(audit_ok(p[["good"]]), "header_only"))
})

test_that("unreadable and missing files become not_assessed rows with the error", {
  p <- make_int_set()
  res <- audit_ok(p)
  # Nothing dropped: four files times four items.
  expect_equal(nrow(res), 16)
  expect_setequal(res$file, basename(p))
  for (f in c("junk.fcs", "gone.fcs")) {
    r <- res[res$file == f, ]
    expect_equal(r$status, rep("not_assessed", 4))
    expect_true(all(nzchar(r$reason) & !is.na(r$reason)))
    expect_equal(r$value[r$item == "readable"], 0)
  }
  # Same short reason as every other import function; no raw R error, no path.
  expect_equal(res$reason[res$file == "junk.fcs"], rep("not an FCS file", 4))
  expect_equal(res$reason[res$file == "gone.fcs"], rep("file not found", 4))
  # The readable files in the same call are unaffected.
  expect_equal(
    res$status[res$file == "good.fcs" & res$item == "tot_vs_rows"],
    "match"
  )
  expect_equal(
    res$status[res$file == "bad_tot.fcs" & res$item == "tot_vs_rows"],
    "differs"
  )
})

test_that("provenance carries version, seed slot, session and arguments", {
  p <- make_int_set()
  res <- audit_ok(p[["good"]])
  pv <- attr(res, "provenance")
  expect_true(all(
    c(
      "package_version",
      "r_version",
      "platform",
      "datetime",
      "seed",
      "session",
      "args"
    ) %in%
      names(pv)
  ))
  expect_true(is.na(pv$seed)) # the audit draws no random numbers
  expect_equal(pv$args$time_keyword, "Time")
  expect_false(pv$args$header_only)
})

# ---- Rule-bearing branches ---------------------------------------------------

# Replace one keyword value in the header of a copy, keeping the same width.
patch_bytes <- function(src, dir, name, from, to) {
  stopifnot(nchar(from) == nchar(to))
  b <- readBin(src, "raw", file.size(src))
  pat <- charToRaw(from)
  at <- which(vapply(
    seq_len(length(b) - length(pat) + 1),
    function(i) all(b[i:(i + length(pat) - 1)] == pat),
    logical(1)
  ))
  stopifnot(length(at) >= 1)
  b[at[[1]] + seq_along(pat) - 1] <- charToRaw(to)
  out <- file.path(dir, name)
  writeBin(b, out)
  out
}

test_that("header_only marks readable not_assessed, because no event was read", {
  p <- make_int_set()
  res <- audit_ok(p[["good"]], header_only = TRUE)
  r <- res[res$item == "readable", ]
  expect_equal(r$status, "not_assessed")
  expect_true(is.na(r$value))
  expect_match(r$reason, "header_only")
})

test_that("an event read refusal with unrelated wording is not_assessed with the text", {
  p <- make_int_set()
  testthat::local_mocked_bindings(
    sw_io_read_fcs = function(...) stop("disk exploded", call. = FALSE)
  )
  res <- audit_ok(p[["good"]])
  r <- res[res$item == "readable", ]
  expect_equal(r$value, 0)
  expect_equal(r$status, "not_assessed")
  expect_match(r$reason, "disk exploded")
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_match(tot$reason, "events unreadable")
  expect_match(tot$reason, "disk exploded")
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "not_assessed")
  expect_match(span$reason, "events unreadable")
  # Header items are still filled.
  expect_equal(res$status[res$item == "channel_set"], "assessed")
})

test_that("a $TOT refusal with changed flowCore wording keeps the error text", {
  p <- make_int_set()
  testthat::local_mocked_bindings(
    sw_io_read_fcs = function(...) stop("rows differ from TOT, new wording")
  )
  res <- audit_ok(p[["bad"]])
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_match(tot$reason, "new wording")
})

test_that("the $TOT refusal gives readable = 1 with the counted rows in detail", {
  p <- make_int_set()
  r <- audit_ok(p[["bad"]])
  r <- r[r$item == "readable", ]
  expect_equal(r$status, "assessed")
  expect_match(r$detail, "20 rows")
  expect_match(r$detail, "TOT")
})

test_that("a valid header whose events cannot be decoded is not readable", {
  p <- make_int_set()
  d <- dirname(p[["good"]])
  # $DATATYPE Q is not one flowCore can decode, and it is not the $TOT wording.
  cut <- patch_bytes(
    p[["good"]],
    d,
    "bad_type.fcs",
    "|$DATATYPE|F|",
    "|$DATATYPE|Q|"
  )
  res <- audit_ok(cut)
  r <- res[res$item == "readable", ]
  expect_equal(r$value, 0)
  expect_equal(r$status, "not_assessed")
  expect_match(r$reason, "DATATYPE")
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_true(is.na(tot$value))
  expect_match(tot$reason, "events unreadable")
  expect_equal(res$status[res$item == "channel_set"], "assessed")
})

test_that("a missing or non-numeric $TOT is not_assessed with a reason", {
  p <- make_int_set()
  d <- dirname(p[["good"]])
  testthat::local_mocked_bindings(
    sw_io_read_fcs = function(...) {
      flowCore::flowFrame(cbind(Time = 1:20, "B1-A" = stats::rnorm(20)))
    }
  )
  abc <- patch_bytes(p[["good"]], d, "tot_abc.fcs", "|$TOT|20|", "|$TOT|ab|")
  res <- audit_ok(abc)
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_match(tot$reason, "no readable \\$TOT")
  none <- patch_bytes(p[["good"]], d, "tot_none.fcs", "|$TOT|20|", "|$TOZ|20|")
  res <- audit_ok(none)
  tot <- res[res$item == "tot_vs_rows", ]
  expect_equal(tot$status, "not_assessed")
  expect_match(tot$reason, "no readable \\$TOT")
})

test_that("events are read without truncating to $PnR", {
  p <- make_int_set()
  d <- dirname(p[["good"]])
  # Time runs 1..20; make the header range 10, so truncation would give 9.
  hdr <- flowCore::read.FCSheader(p[["good"]])[[1]]
  rng <- hdr[["$P1R"]]
  expect_true(nchar(rng) >= 2)
  small <- patch_bytes(
    p[["good"]],
    d,
    "small_range.fcs",
    paste0("|$P1R|", rng, "|"),
    paste0("|$P1R|", strrep("1", nchar(rng)), "|")
  )
  res <- audit_ok(small, time_keyword = "Time")
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "assessed")
  expect_equal(span$value, 19)
})

test_that("a zero-row file gives a not_assessed time span, not a crash", {
  p <- make_int_set()
  testthat::local_mocked_bindings(
    sw_io_read_fcs = function(...) {
      ff <- flowCore::flowFrame(cbind(Time = 1:20, "B1-A" = stats::rnorm(20)))
      suppressWarnings(ff[integer(0), ])
    }
  )
  res <- audit_ok(p[["good"]])
  expect_equal(nrow(res), 4)
  expect_equal(
    res$item,
    c("readable", "tot_vs_rows", "time_span", "channel_set")
  )
  expect_equal(res$detail[res$item == "readable"], "0 rows read")
  expect_equal(res$value[res$item == "tot_vs_rows"], 20)
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "not_assessed")
  expect_true(is.na(span$value))
  expect_match(span$reason, "no finite value")
})

test_that("non-finite values in the time channel are counted in detail", {
  p <- make_int_set()
  testthat::local_mocked_bindings(
    sw_io_read_fcs = function(...) {
      ff <- flowCore::flowFrame(cbind(Time = 1:20, "B1-A" = stats::rnorm(20)))
      e <- flowCore::exprs(ff)
      e[19, "Time"] <- NA
      e[20, "Time"] <- Inf
      flowCore::exprs(ff) <- e
      ff
    }
  )
  res <- audit_ok(p[["good"]])
  span <- res[res$item == "time_span", ]
  expect_equal(span$status, "assessed")
  expect_equal(span$value, 17)
  expect_match(span$detail, "2 non-finite")
})

test_that("clock text: midnight, tails and junk", {
  p <- make_int_set()
  d <- dirname(p[["good"]])
  # Rewrite $BTIM and $ETIM on a copy through the mock of the header reader.
  mk <- function(b, e) {
    kw <- list("$BTIM" = b, "$ETIM" = e)
    make_int_fcs(d, paste0("t_", sample.int(1e8, 1), ".fcs"), kw = kw)
  }
  get <- function(b, e) {
    r <- audit_ok(mk(b, e), time_keyword = c("$BTIM", "$ETIM"))
    r[r$item == "time_span", ]
  }
  mid <- get("23:59:00", "00:01:00")
  expect_equal(mid$status, "not_assessed")
  expect_true(is.na(mid$value))
  expect_match(mid$reason, "earlier")
  expect_match(mid$reason, "no day is added")
  hund <- get("10:00:00.50", "10:00:01.50")
  expect_equal(hund$value, 1)
  hund2 <- get("10:00:00.25", "10:00:01.75")
  expect_equal(hund2$value, 1.5)
  sixty <- get("10:00:00:30", "10:00:01:60")
  # 30/60 = 0.5 s and 60/60 = 1 s: begin 0.5, end 2
  expect_equal(sixty$value, 1.5)
  tail_only <- get("23:59:59.50", "23:59:59.75")
  expect_equal(tail_only$value, 0.25)
  junk <- get("abc", "10:00:00")
  expect_equal(junk$status, "not_assessed")
  expect_match(junk$reason, "cannot read time")
  expect_match(junk$reason, "abc")
  bad_tail <- get("10:00:00", "10:00:00.5x")
  expect_equal(bad_tail$status, "not_assessed")
  expect_match(bad_tail$reason, "cannot read time")
  expect_equal(.sw_int_clock("23:59:59.50"), 86399.5)
})

test_that("the shared reason helper gives one short line without the path", {
  d <- withr::local_tempdir()
  expect_equal(.sw_fcs_read_reason(file.path(d, "gone.fcs")), "file not found")
  expect_equal(.sw_fcs_read_reason(d), "file not found")
  junk <- file.path(d, "junk.fcs")
  writeLines("plain text", junk)
  expect_equal(.sw_fcs_read_reason(junk), "not an FCS file")
  # Starts with FCS but is not a valid header: first line of the error, no path.
  odd <- file.path(d, "odd.fcs")
  writeLines("FCS3.1 and then nothing useful", odd)
  e <- tryCatch(.sw_read_fcs_header(odd), error = function(e) e)
  skip_if_not(inherits(e, "error"))
  why <- .sw_fcs_read_reason(odd, e)
  expect_match(why, "^header unreadable: [^\n]+$")
  expect_false(grepl(d, why, fixed = TRUE))
  expect_false(grepl("^header unreadable: Error in", why))
})

test_that("provenance starts with the shared keys (integrity)", {
  pv <- attr(audit_ok(make_int_set()[["good"]]), "provenance")
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
