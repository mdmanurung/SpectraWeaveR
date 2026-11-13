# tests/testthat/test-fcs_header_audit.R
# Unit tests for R/fcs_header_audit.R. Writes tiny FCS files to a temp folder.

# Write one FCS file with chosen channels, markers, voltages and keywords.
make_fcs <- function(dir, name, chan, marker, volt = NULL, kw = list()) {
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(20 * length(chan)),
    ncol = length(chan),
    dimnames = list(NULL, chan)
  )
  ff <- flowCore::flowFrame(m)
  flowCore::pData(flowCore::parameters(ff))$desc <- marker
  desc <- flowCore::keyword(ff)
  for (i in seq_along(chan)) {
    if (!is.null(volt)) desc[[paste0("$P", i, "V")]] <- as.character(volt[i])
  }
  for (k in names(kw)) {
    desc[[k]] <- kw[[k]]
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

make_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  chan <- c("B1-A", "V1-A")
  marker <- c("CD3", "CD4")
  c(
    ref = make_fcs(d, "ref.fcs", chan, marker, c(500, 600), list("$CYT" = "A")),
    same = make_fcs(
      d,
      "same.fcs",
      chan,
      marker,
      c(500, 600),
      list("$CYT" = "A")
    ),
    volt = make_fcs(
      d,
      "volt.fcs",
      chan,
      marker,
      c(500, 650),
      list("$CYT" = "A")
    ),
    inst = make_fcs(
      d,
      "inst.fcs",
      chan,
      marker,
      c(500, 600),
      list("$CYT" = "B")
    ),
    chan = make_fcs(
      d,
      "chan.fcs",
      c("B1-A", "V2-A"),
      marker,
      c(500, 600),
      list("$CYT" = "A")
    )
  )
}

test_that("required arguments have no default", {
  expect_error(sw_fcs_header_audit(), "`paths` is required")
  expect_error(sw_fcs_header_audit(c("a", "b")), "`keywords` is required")
  expect_error(
    sw_fcs_header_audit(c("a", "b"), keywords = "$CYT"),
    "`reference` is required"
  )
})

test_that("reference must be one of paths and names must be unique", {
  expect_error(
    sw_fcs_header_audit(c("a/x.fcs", "b/y.fcs"), "$CYT", "z.fcs"),
    "exactly one of"
  )
  expect_error(
    sw_fcs_header_audit(c("a/x.fcs", "b/x.fcs"), "$CYT", "a/x.fcs"),
    "unique"
  )
})

test_that("a changed keyword, voltage and channel are each reported", {
  p <- make_set()
  res <- sw_fcs_header_audit(
    unname(p),
    c("$CYT", "$PnV"),
    reference = p[["ref"]]
  )
  expect_s3_class(res, "sw_fcs_header_audit")
  expect_equal(attr(res, "reference"), "ref.fcs")

  # Reference and identical file: every row matches.
  expect_true(all(res$status[res$file == "ref.fcs"] == "match"))
  expect_true(all(res$status[res$file == "same.fcs"] == "match"))

  # Instrument keyword differs only in inst.fcs.
  inst <- res[res$file == "inst.fcs" & res$kind == "keyword", ]
  expect_equal(inst$status, "differs")
  expect_equal(c(inst$value, inst$reference_value), c("B", "A"))

  # Voltage differs for one detector only.
  volt <- res[res$file == "volt.fcs" & res$kind == "parameter_keyword", ]
  expect_equal(volt$channel[volt$status == "differs"], "V1-A")
  expect_equal(sum(volt$status == "match"), 1)

  # Renamed channel: one missing from the file, one absent in the reference.
  ch <- res[res$file == "chan.fcs" & res$kind == "channel", ]
  expect_setequal(
    ch$status[ch$key %in% c("V1-A", "V2-A")],
    c("missing_in_file", "absent_in_reference")
  )
})

test_that("a keyword the reference lacks gives not_assessed rows", {
  p <- make_set()
  res <- sw_fcs_header_audit(
    unname(p),
    c("$CYT", "$NOPE"),
    reference = p[["ref"]]
  )
  na <- res[res$key == "$NOPE", ]
  expect_equal(nrow(na), 5)
  expect_true(all(na$status == "not_assessed"))
  expect_true(all(na$reason == "reference lacks keyword"))
  expect_true(all(res$status[res$key == "$CYT"] %in% c("match", "differs")))
})

test_that("an unreadable file gives a not_assessed row and one path works", {
  p <- make_set()
  d <- dirname(p[["ref"]])
  junk <- file.path(d, "junk.fcs")
  writeLines("not an fcs file", junk)
  empty <- file.path(d, "empty.fcs")
  file.create(empty)
  expect_no_warning(
    res <- sw_fcs_header_audit(
      c(p[["ref"]], p[["same"]], junk, empty),
      "$CYT",
      reference = p[["ref"]]
    )
  )
  bad <- res[res$file %in% c("junk.fcs", "empty.fcs"), ]
  expect_equal(nrow(bad), 2)
  expect_true(all(bad$status == "not_assessed" & bad$kind == "file"))
  expect_false(anyNA(bad$reason))
  s <- sw_fcs_header_summary(res)
  expect_equal(s$not_assessed[s$file == "junk.fcs"], 1)
  one <- sw_fcs_header_audit(p[["ref"]], "$CYT", reference = p[["ref"]])
  expect_true(all(one$status == "match"))
  expect_error(
    sw_fcs_header_audit(junk, "$CYT", reference = junk),
    "Reference file could not be read"
  )
  expect_error(sw_fcs_header_audit(character(0), "$CYT", "a"), "one or more")
})

test_that("duplicated and out-of-range $PnN keys are ignored and reported", {
  h <- c(
    "$PAR" = "2",
    "$P1N" = "A-A",
    "$P1S" = "CD3",
    "$P2N" = "B-A",
    "$P2S" = "CD4",
    "$P3N" = "Width",
    "$P1N" = "Z-A",
    "$CYT" = "X"
  )
  it <- .sw_hdr_items(h, "$CYT")
  expect_equal(it$key[it$kind == "channel"], c("A-A", "B-A"))
  expect_match(attr(it, "notes")$reason, "repeated")
  clean <- .sw_hdr_items(h[1:5], "$CYT")
  expect_equal(nrow(attr(clean, "notes")), 0)
})

test_that("summary counts statuses and puts differing files first", {
  p <- make_set()
  res <- sw_fcs_header_audit(
    unname(p),
    c("$CYT", "$PnV"),
    reference = p[["ref"]]
  )
  s <- sw_fcs_header_summary(res)
  expect_equal(nrow(s), 5)
  expect_true(all(
    c("file", "match", "differs", "missing_in_file", "absent_in_reference") %in%
      names(s)
  ))
  expect_equal(s$differs[s$file == "same.fcs"], 0)
  expect_true(s$file[nrow(s)] %in% c("ref.fcs", "same.fcs"))
  expect_error(sw_fcs_header_summary(data.frame()), "must come from")
})
