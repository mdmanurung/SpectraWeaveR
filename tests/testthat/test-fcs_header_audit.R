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

test_that("a keyword the reference lacks is an error", {
  p <- make_set()
  expect_error(
    sw_fcs_header_audit(unname(p), "$NOPE", reference = p[["ref"]]),
    "Reference file lacks keyword"
  )
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
