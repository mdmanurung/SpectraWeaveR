# tests/testthat/test-import_batch_date.R
# Unit tests for R/import_batch_date.R. Writes tiny FCS files to a temp folder.

# Write one FCS file with chosen date, begin time and serial keywords.
make_dt_fcs <- function(
  dir,
  name,
  date = "2026-03-01",
  btim = "09:00:00",
  sn = "SN1"
) {
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(20),
    ncol = 2,
    dimnames = list(NULL, c("B1-A", "V1-A"))
  )
  ff <- flowCore::flowFrame(m)
  desc <- flowCore::keyword(ff)
  if (!is.null(date)) {
    desc[["$DATE"]] <- date
  }
  if (!is.null(btim)) {
    desc[["$BTIM"]] <- btim
  }
  if (!is.null(sn)) {
    desc[["$CYTSN"]] <- sn
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

# Six files: b1 spans two dates and a 10 h gap; b2 is one morning, one serial.
make_dt_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  paths <- c(
    s1 = make_dt_fcs(d, "s1.fcs", "2026-03-01", "09:00:00", "SN1"),
    s2 = make_dt_fcs(d, "s2.fcs", "2026-03-01", "10:00:00", "SN1"),
    s3 = make_dt_fcs(d, "s3.fcs", "2026-03-01", "20:00:00", "SN1"),
    s4 = make_dt_fcs(d, "s4.fcs", "2026-03-02", "09:30:00", "SN2"),
    s5 = make_dt_fcs(d, "s5.fcs", "2026-03-05", "08:00:00", "SN1"),
    s6 = make_dt_fcs(d, "s6.fcs", "2026-03-05", "08:20:00", "SN1")
  )
  man <- sw_sample_manifest(
    data.frame(
      id = names(paths),
      b = c("b1", "b1", "b1", "b1", "b2", "b2"),
      r = c("r1", "r1", "r2", "r2", "r3", "r3")
    ),
    sample_id = "id",
    batch = "b",
    run = "r"
  )
  list(paths = paths, man = man, dir = d)
}

args_ok <- function(set, ...) {
  a <- list(
    paths = set$paths,
    manifest = set$man,
    date_keywords = "$DATE",
    time_keywords = "$BTIM",
    serial_keyword = "$CYTSN",
    date_format = "%Y-%m-%d",
    time_zone = "UTC",
    run_gap_hours = 4
  )
  over <- list(...)
  a[names(over)] <- over
  a
}

run_ok <- function(set, ...) do.call(sw_fcs_batch_date_check, args_ok(set, ...))

test_that("every required argument has no default and is named", {
  set <- make_dt_set()
  a <- args_ok(set)
  for (nm in names(a)) {
    expect_error(
      do.call(sw_fcs_batch_date_check, a[setdiff(names(a), nm)]),
      "no default",
      info = nm
    )
    expect_error(
      do.call(sw_fcs_batch_date_check, a[setdiff(names(a), nm)]),
      sprintf("'%s'", nm),
      fixed = TRUE,
      info = nm
    )
  }
  expect_error(sw_fcs_batch_date_check(), "no default")
})

test_that("bad argument values stop before any file is read", {
  set <- make_dt_set()
  expect_error(run_ok(set, run_gap_hours = -1), "run_gap_hours")
  expect_error(run_ok(set, run_gap_hours = c(1, 2)), "run_gap_hours")
  expect_error(run_ok(set, time_zone = "Nowhere/Land"), "time_zone")
  expect_error(
    run_ok(set, manifest = data.frame(sample_id = "a", batch = "b")),
    "sw_sample_manifest"
  )
  expect_error(run_ok(set, serial_keyword = c("$A", "$B")), "serial_keyword")
  expect_error(run_ok(set, date_keywords = character(0)), "date_keywords")
  unnamed <- unname(set$paths)
  expect_error(
    run_ok(set, paths = c(unnamed, file.path(tempdir(), "x", "s1.fcs"))),
    "unique"
  )
})

test_that("files table parses date, time and serial and keeps the raw text", {
  set <- make_dt_set()
  res <- run_ok(set)
  expect_s3_class(res, "sw_fcs_batch_date_check")
  f <- res$files
  expect_equal(nrow(f), 6)
  expect_equal(f$sample_id, names(set$paths))
  expect_equal(f$batch, c("b1", "b1", "b1", "b1", "b2", "b2"))
  expect_equal(f$run, c("r1", "r1", "r2", "r2", "r3", "r3"))
  expect_equal(f$date_raw[1], "2026-03-01")
  expect_equal(f$time_raw[2], "10:00:00")
  expect_equal(f$serial, c("SN1", "SN1", "SN1", "SN2", "SN1", "SN1"))
  expect_equal(f$date_keyword[1], "$DATE")
  expect_equal(
    format(f$datetime[4], "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "2026-03-02 09:30:00"
  )
  expect_equal(attr(f$datetime, "tzone"), "UTC")
  expect_true(all(f$status == "assessed"))
})

test_that("crosstab counts files per batch and date and keeps empty cells", {
  set <- make_dt_set()
  ct <- run_ok(set)$crosstab
  expect_true(all(c("batch", "date", "n_files") %in% names(ct)))
  n <- function(b, d) ct$n_files[ct$batch == b & ct$date %in% as.Date(d)]
  expect_equal(n("b1", "2026-03-01"), 3)
  expect_equal(n("b1", "2026-03-02"), 1)
  expect_equal(n("b2", "2026-03-05"), 2)
  expect_equal(n("b2", "2026-03-01"), 0)
  expect_equal(sum(ct$n_files), 6)
})

test_that("gaps, dates and serials are counted per batch and per run", {
  set <- make_dt_set()
  g <- run_ok(set)$groups
  b1 <- g[g$level == "batch" & g$batch == "b1", ]
  expect_equal(b1$n_files, 4)
  expect_equal(b1$n_dates, 2)
  expect_equal(b1$n_serials, 2)
  # Order: 09:00, 10:00, 20:00, next day 09:30. Gaps 1, 10, 13.5 h; limit 4.
  expect_equal(b1$n_gaps_over, 2)
  expect_equal(b1$max_gap_hours, 13.5)
  expect_equal(b1$n_sessions, 3)
  expect_equal(b1$run_gap_hours, 4)
  b2 <- g[g$level == "batch" & g$batch == "b2", ]
  expect_equal(
    c(b2$n_dates, b2$n_serials, b2$n_gaps_over, b2$n_sessions),
    c(1, 1, 0, 1)
  )
  # Declared run r1 holds two files one hour apart: one session.
  r1 <- g[g$level == "run" & g$run == "r1", ]
  expect_equal(c(r1$n_files, r1$n_sessions, r1$n_gaps_over), c(2, 1, 0))
  # Declared run r2 spans a 13.5 h gap inside one declared run.
  r2 <- g[g$level == "run" & g$run == "r2", ]
  expect_equal(c(r2$n_sessions, r2$n_gaps_over), c(2, 1))
  expect_equal(run_ok(set)$files$session, c(1, 1, 2, 3, 1, 1))
})

test_that("an unparseable date is a not_assessed row with the raw text", {
  d <- withr::local_tempdir()
  paths <- c(
    a = make_dt_fcs(d, "a.fcs", "2026-03-01", "09:00:00"),
    b = make_dt_fcs(d, "b.fcs", "01/03/26 oops", "09:00:00"),
    c = make_dt_fcs(d, "c.fcs", NULL, "09:00:00")
  )
  man <- sw_sample_manifest(
    data.frame(id = names(paths), b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  res <- sw_fcs_batch_date_check(
    paths,
    man,
    "$DATE",
    "$BTIM",
    "$CYTSN",
    "%Y-%m-%d",
    "UTC",
    4
  )
  f <- res$files
  expect_equal(nrow(f), 3)
  expect_equal(f$status, c("assessed", "not_assessed", "not_assessed"))
  expect_equal(f$date_raw[2], "01/03/26 oops")
  expect_match(f$reason[2], "date")
  expect_true(is.na(f$datetime[2]))
  expect_equal(f$date_status[3], "missing_keyword")
  # The two rows still count in the crosstab, under an NA date.
  ct <- res$crosstab
  expect_equal(sum(ct$n_files), 3)
  expect_equal(ct$n_files[is.na(ct$date)], 2)
  # And in the batch group row.
  g <- res$groups
  expect_equal(g$n_files, 3)
  expect_equal(g$n_not_assessed, 2)
})

test_that("a missing time leaves the date but gives no datetime", {
  d <- withr::local_tempdir()
  paths <- c(
    a = make_dt_fcs(d, "a.fcs", "2026-03-01", "09:00:00"),
    b = make_dt_fcs(d, "b.fcs", "2026-03-01", NULL)
  )
  man <- sw_sample_manifest(
    data.frame(id = names(paths), b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  res <- sw_fcs_batch_date_check(
    paths,
    man,
    "$DATE",
    "$BTIM",
    "$CYTSN",
    "%Y-%m-%d",
    "UTC",
    4
  )
  f <- res$files
  expect_equal(f$status, c("assessed", "date_only"))
  expect_false(is.na(f$date[2]))
  expect_true(is.na(f$datetime[2]))
  expect_equal(
    res$crosstab$n_files[res$crosstab$date %in% as.Date("2026-03-01")],
    2
  )
})

test_that("an unreadable file is a not_assessed row, not a crash", {
  set <- make_dt_set()
  junk <- file.path(set$dir, "junk.fcs")
  writeLines("not an fcs file", junk)
  gone <- file.path(set$dir, "gone.fcs")
  paths <- c(set$paths, junk = junk, gone = gone)
  man <- sw_sample_manifest(
    data.frame(
      id = names(paths),
      b = c("b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2")
    ),
    sample_id = "id",
    batch = "b"
  )
  res <- run_ok(set, paths = paths, manifest = man)
  f <- res$files
  expect_equal(nrow(f), 8)
  bad <- f[f$sample_id %in% c("junk", "gone"), ]
  expect_equal(bad$status, c("not_assessed", "not_assessed"))
  expect_equal(bad$reason[bad$sample_id == "junk"], "not an FCS file")
  expect_equal(bad$reason[bad$sample_id == "gone"], "file not found")
  expect_true(all(is.na(bad$datetime)))
  expect_equal(sum(res$crosstab$n_files), 8)
})

test_that("files without a manifest row and rows without a file are listed", {
  set <- make_dt_set()
  extra <- make_dt_fcs(set$dir, "extra.fcs")
  paths <- c(set$paths[1:5], extra = extra)
  man <- sw_sample_manifest(
    data.frame(id = c(names(set$paths)[1:5], "ghost"), b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  res <- run_ok(set, paths = paths, manifest = man)
  u <- res$unmatched
  expect_setequal(
    u$side,
    c("file_without_manifest_row", "manifest_row_without_file")
  )
  expect_equal(u$sample_id[u$side == "file_without_manifest_row"], "extra")
  expect_equal(u$sample_id[u$side == "manifest_row_without_file"], "ghost")
  expect_equal(
    u$path[u$side == "file_without_manifest_row"],
    unname(extra)
  )
  # The extra file stays in the files table with no declared batch.
  f <- res$files
  expect_equal(nrow(f), 6)
  expect_false(f$in_manifest[f$sample_id == "extra"])
  expect_true(is.na(f$batch[f$sample_id == "extra"]))
  expect_equal(sum(res$crosstab$n_files), 6)
})

test_that("unnamed paths take the file name without extension as sample id", {
  set <- make_dt_set()
  res <- run_ok(set, paths = unname(set$paths))
  expect_equal(res$files$sample_id, names(set$paths))
  expect_equal(nrow(res$unmatched), 0)
})

test_that("a manifest without a run column still works", {
  set <- make_dt_set()
  man <- sw_sample_manifest(
    data.frame(
      id = names(set$paths),
      b = c("b1", "b1", "b1", "b1", "b2", "b2")
    ),
    sample_id = "id",
    batch = "b"
  )
  res <- run_ok(set, manifest = man)
  expect_false("run" %in% res$groups$level)
  expect_true(all(is.na(res$files$run)))
})

test_that("files are never edited and the result has no verdict column", {
  set <- make_dt_set()
  before <- unname(tools::md5sum(set$paths))
  res <- run_ok(set)
  expect_equal(unname(tools::md5sum(set$paths)), before)
  bad <- c("pass", "fail", "verdict", "ok", "flag", "score")
  for (tb in list(res$files, res$crosstab, res$groups, res$unmatched)) {
    expect_false(any(tolower(names(tb)) %in% bad))
  }
})

test_that("the result carries provenance", {
  set <- make_dt_set()
  p <- run_ok(set)$provenance
  expect_true(all(
    c("package_version", "r_version", "datetime", "args", "seed") %in% names(p)
  ))
  expect_true(is.na(p$seed))
  expect_equal(p$args$run_gap_hours, 4)
  expect_equal(p$args$time_zone, "UTC")
})

test_that("the first date keyword that holds a value is used and named", {
  d <- withr::local_tempdir()
  paths <- c(a = make_dt_fcs(d, "a.fcs", "2026-03-01", "09:00:00"))
  man <- sw_sample_manifest(
    data.frame(id = "a", b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  res <- sw_fcs_batch_date_check(
    paths,
    man,
    c("$NOPE", "$DATE"),
    "$BTIM",
    "$CYTSN",
    "%Y-%m-%d",
    "UTC",
    4
  )
  expect_equal(res$files$date_keyword, "$DATE")
})

# Run the check on a small set of (date, time) pairs in one batch.
run_pairs <- function(dates, times, ..., names = letters[seq_along(dates)]) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  paths <- vapply(
    seq_along(dates),
    function(i) make_dt_fcs(d, paste0(names[i], ".fcs"), dates[i], times[i]),
    character(1)
  )
  names(paths) <- names
  man <- sw_sample_manifest(
    data.frame(id = names, b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  a <- list(
    paths = paths,
    manifest = man,
    date_keywords = "$DATE",
    time_keywords = "$BTIM",
    serial_keyword = "$CYTSN",
    date_format = "%Y-%m-%d",
    time_zone = "UTC",
    run_gap_hours = 4
  )
  over <- list(...)
  a[names(over)] <- over
  do.call(sw_fcs_batch_date_check, a)
}

test_that("a two-digit year is unparseable, not year 0026", {
  res <- run_pairs(
    c("26-03-01", "2026-03-01"),
    c("09:00:00", "09:00:00")
  )
  f <- res$files
  expect_equal(f$date_status, c("unparseable", "parsed"))
  expect_equal(f$status, c("not_assessed", "assessed"))
  expect_true(is.na(f$date[1]))
  expect_equal(f$date_raw[1], "26-03-01")
  expect_match(f$reason[1], "year")
  # Same with a month-name format.
  res2 <- run_pairs("01-Mar-26", "09:00:00", date_format = "%d-%b-%Y")
  expect_equal(res2$files$status, "not_assessed")
  expect_true(is.na(res2$files$date))
  # No bad year reaches the crosstab.
  expect_false(any(res$crosstab$date < as.Date("1000-01-01"), na.rm = TRUE))
})

test_that("a bad time gives date_only; a fraction is ignored", {
  res <- run_pairs(
    rep("2026-03-01", 5),
    c("25:00:00", "abc", "09:00:00.5", "09:60:00", "09:00:61")
  )
  f <- res$files
  expect_equal(
    f$time_status,
    c("unparseable", "unparseable", "parsed", "unparseable", "unparseable")
  )
  expect_equal(
    f$status,
    c("date_only", "date_only", "assessed", "date_only", "date_only")
  )
  expect_true(all(is.na(f$datetime[c(1, 2, 4, 5)])))
  expect_match(f$reason[1], "time")
  expect_equal(
    format(f$datetime[3], "%H:%M:%S", tz = "UTC"),
    "09:00:00"
  )
})

test_that("time_zone sets the zone and the clock time of the datetime", {
  res <- run_pairs("2026-07-01", "09:00:00", time_zone = "Europe/Amsterdam")
  dt <- res$files$datetime
  expect_equal(attr(dt, "tzone"), "Europe/Amsterdam")
  expect_equal(format(dt, "%H:%M:%S"), "09:00:00")
  # 09:00 CEST is 07:00 UTC.
  expect_equal(format(dt, "%H:%M:%S", tz = "UTC"), "07:00:00")
})

test_that("a gap exactly equal to run_gap_hours does not start a session", {
  res <- run_pairs(
    rep("2026-03-01", 3),
    c("08:00:00", "12:00:00", "17:00:01")
  )
  f <- res$files
  # 4 h (equal to the limit: same session), then 5 h 1 s (over: new session).
  expect_equal(f$session, c(1, 1, 2))
  expect_equal(f$gap_before_hours[2], 4)
  g <- res$groups
  expect_equal(g$n_gaps_over, 1)
  expect_equal(g$n_sessions, 2)
  res2 <- run_pairs(rep("2026-03-01", 2), c("08:00:00", "12:00:00"))
  expect_equal(res2$groups$n_gaps_over, 0)
  expect_equal(res2$files$session, c(1, 1))
})

test_that("keyword names match without regard to case", {
  res <- run_pairs("2026-03-01", "09:00:00", date_keywords = "$date")
  expect_equal(res$files$date_status, "parsed")
  expect_equal(res$files$date_keyword, "$date")
})

test_that("the first listed keyword with a value wins when two hold values", {
  d <- withr::local_tempdir()
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(20),
    ncol = 2,
    dimnames = list(NULL, c("B1-A", "V1-A"))
  )
  ff <- flowCore::flowFrame(m)
  desc <- flowCore::keyword(ff)
  desc[["$DATE"]] <- "2026-03-01"
  desc[["$ETIM"]] <- "2026-03-09"
  desc[["$BTIM"]] <- "09:00:00"
  desc[["$CYTSN"]] <- "SN1"
  flowCore::keyword(ff) <- desc
  p <- file.path(d, "a.fcs")
  flowCore::write.FCS(ff, p)
  man <- sw_sample_manifest(
    data.frame(id = "a", b = "b1"),
    sample_id = "id",
    batch = "b"
  )
  call <- function(kw) {
    sw_fcs_batch_date_check(
      c(a = p),
      man,
      kw,
      "$BTIM",
      "$CYTSN",
      "%Y-%m-%d",
      "UTC",
      4
    )$files
  }
  expect_equal(call(c("$DATE", "$ETIM"))$date, as.Date("2026-03-01"))
  expect_equal(call(c("$ETIM", "$DATE"))$date, as.Date("2026-03-09"))
})

test_that("sessions follow time order, not the order of paths", {
  fwd <- run_pairs(
    rep("2026-03-01", 3),
    c("08:00:00", "09:00:00", "20:00:00")
  )
  rev <- run_pairs(
    rep("2026-03-01", 3),
    c("20:00:00", "09:00:00", "08:00:00")
  )
  expect_equal(rev$files$session, c(2, 1, 1))
  expect_true(all(rev$files$gap_before_hours >= 0, na.rm = TRUE))
  expect_equal(rev$files$gap_before_hours, c(11, 1, NA))
  expect_equal(rev$groups$n_sessions, fwd$groups$n_sessions)
  expect_equal(rev$groups$n_gaps_over, fwd$groups$n_gaps_over)
  expect_equal(rev$groups$max_gap_hours, fwd$groups$max_gap_hours)
})

test_that("argument errors come before any file is read", {
  set <- make_dt_set()
  ghosts <- c(a = "/no/such/a.fcs", b = "/no/such/b.fcs")
  expect_error(
    run_ok(set, paths = ghosts, run_gap_hours = -1),
    "run_gap_hours"
  )
  # Valid arguments with the same ghost paths give not_assessed rows instead.
  res <- run_ok(set, paths = ghosts)
  expect_true(all(res$files$status == "not_assessed"))
})

test_that("files table leads with file (basename) and path, like the other import tables", {
  set <- make_dt_set()
  res <- run_ok(set)
  expect_identical(names(res$files)[1:3], c("file", "path", "sample_id"))
  expect_identical(res$files$file, basename(unname(set$paths)))
  expect_identical(res$files$file, paste0(res$files$sample_id, ".fcs"))
})

test_that("provenance starts with the shared keys (batch_date)", {
  pv <- {
    res <- run_ok(make_dt_set())
    expect_identical(attr(res, "provenance"), res$provenance)
    res$provenance
  }
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
