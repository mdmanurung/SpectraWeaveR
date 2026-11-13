# tests/testthat/test-import_signatures.R
# Unit tests for R/import_signatures.R. Writes tiny FCS files to a temp folder.

# Write one FCS file with chosen channels, markers, voltages and keywords.
sig_fcs <- function(dir, name, chan, marker, volt = NULL, kw = list()) {
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(20 * length(chan)),
    ncol = length(chan),
    dimnames = list(NULL, chan)
  )
  ff <- flowCore::flowFrame(m)
  flowCore::pData(flowCore::parameters(ff))$desc <- marker
  desc <- flowCore::keyword(ff)
  if (!is.null(volt)) {
    for (i in seq_along(chan)) {
      desc[[paste0("$P", i, "V")]] <- as.character(volt[i])
    }
  }
  for (k in names(kw)) {
    desc[[k]] <- kw[[k]]
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

# a1 and a2 share a signature; b1 differs by $CYT; c1 by channel set; v1 by voltage.
sig_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  chan <- c("B1-A", "V1-A")
  marker <- c("CD3", "CD4")
  c(
    a1 = sig_fcs(d, "a1.fcs", chan, marker, c(500, 600), list("$CYT" = "A")),
    a2 = sig_fcs(d, "a2.fcs", chan, marker, c(500, 600), list("$CYT" = "A")),
    b1 = sig_fcs(d, "b1.fcs", chan, marker, c(500, 600), list("$CYT" = "B")),
    c1 = sig_fcs(
      d,
      "c1.fcs",
      c("B1-A", "V2-A"),
      marker,
      c(500, 600),
      list("$CYT" = "A")
    ),
    v1 = sig_fcs(d, "v1.fcs", chan, marker, c(500, 650), list("$CYT" = "A"))
  )
}

test_that("required arguments have no default", {
  expect_error(sw_fcs_header_signatures(), "no default.*'paths'")
  expect_error(
    sw_fcs_header_signatures(paths = "a.fcs", panel = NULL),
    "no default.*'keywords'"
  )
  expect_error(
    sw_fcs_header_signatures(paths = "a.fcs", keywords = "$CYT"),
    "no default.*'panel'"
  )
})

test_that("bad arguments stop before any file is read", {
  expect_error(
    sw_fcs_header_signatures(c("a/x.fcs", "b/x.fcs"), "$CYT", NULL),
    "unique"
  )
  expect_error(
    sw_fcs_header_signatures("x.fcs", character(), NULL),
    "keywords"
  )
  expect_error(
    sw_fcs_header_signatures("x.fcs", "$CYT", list(c("CD3"))),
    "panel"
  )
  expect_error(
    sw_fcs_header_signatures("x.fcs", "$CYT", list(T = 1:2)),
    "panel"
  )
  expect_error(
    sw_fcs_header_signatures("x.fcs", "$CYT", list(T = "CD3", T = "CD4")),
    "panel"
  )
})

test_that("identical keyword values and channel sets share a signature", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p, c("$CYT", "$PnV"), panel = NULL)
  f <- res$files
  expect_s3_class(res, "sw_fcs_header_signatures")
  expect_equal(f$file, basename(unname(p)))
  expect_equal(f$signature[f$file == "a1.fcs"], f$signature[f$file == "a2.fcs"])
  others <- f$signature[f$file != "a1.fcs" & f$file != "a2.fcs"]
  expect_false(any(others %in% f$signature[f$file == "a1.fcs"]))
  expect_equal(length(unique(f$signature)), 4)
  # Counts add up to the files and say who is in each signature.
  expect_equal(sum(res$counts$n_files), length(p))
  top <- res$counts[1, ]
  expect_equal(top$n_files, 2)
  expect_match(top$files, "a1.fcs")
  expect_match(top$files, "a2.fcs")
})

test_that("channel set alone splits signatures when no keyword differs", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p[c("a1", "c1")], "$CYT", panel = NULL)
  expect_equal(length(unique(res$files$signature)), 2)
})

test_that("a keyword absent from a file is its own value, not a crash", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p[c("a1", "b1")], "$NOSUCH", panel = NULL)
  expect_equal(length(unique(res$files$signature)), 1)
  expect_true(all(res$files$status == "assessed"))
})

test_that("an unreadable file is a not_assessed row and the rest still run", {
  p <- sig_set()
  d <- dirname(p[[1]])
  junk <- file.path(d, "junk.fcs")
  writeLines("this is not an fcs file", junk)
  gone <- file.path(d, "gone.fcs")
  res <- sw_fcs_header_signatures(
    c(p[c("a1", "a2")], junk, gone),
    "$CYT",
    panel = NULL
  )
  f <- res$files
  expect_equal(nrow(f), 4)
  bad <- f[f$file %in% c("junk.fcs", "gone.fcs"), ]
  expect_true(all(bad$status == "not_assessed"))
  expect_true(all(nzchar(bad$reason)))
  expect_true(all(is.na(bad$signature)))
  expect_true(all(f$status[f$file %in% c("a1.fcs", "a2.fcs")] == "assessed"))
  # Not-assessed files are counted, not hidden.
  expect_equal(sum(res$counts$n_files), 2)
  expect_equal(res$n_not_assessed, 2)
})

test_that("panel match lists missing and extra markers per file and entry", {
  p <- sig_set()
  panel <- list(T = c("CD3", "CD4"), B = c("CD3", "CD19"))
  res <- sw_fcs_header_signatures(p[c("a1", "c1")], "$CYT", panel = panel)
  pm <- res$panel_match
  expect_equal(nrow(pm), 2 * 2)
  row <- pm[pm$file == "a1.fcs" & pm$panel == "T", ]
  expect_equal(row$n_missing, 0)
  expect_equal(row$n_extra, 0)
  row <- pm[pm$file == "a1.fcs" & pm$panel == "B", ]
  expect_equal(row$missing, "CD19")
  expect_equal(row$extra, "CD4")
  expect_equal(row$n_matched, 1)
  expect_true(all(pm$status == "assessed"))
})

test_that("panel match keeps a row for every file, also unreadable ones", {
  p <- sig_set()
  junk <- file.path(dirname(p[[1]]), "junk.fcs")
  writeLines("nope", junk)
  res <- sw_fcs_header_signatures(
    c(p["a1"], junk),
    "$CYT",
    panel = list(T = c("CD3", "CD4"), B = "CD19")
  )
  pm <- res$panel_match
  expect_equal(nrow(pm), 4)
  expect_true(all(pm$status[pm$file == "junk.fcs"] == "not_assessed"))
  expect_true(all(is.na(pm$n_missing[pm$file == "junk.fcs"])))
})

test_that("channels without a marker are counted, not treated as extra markers", {
  d <- withr::local_tempdir()
  f <- sig_fcs(d, "t.fcs", c("FSC-A", "B1-A"), c(NA, "CD3"))
  res <- sw_fcs_header_signatures(f, "$CYT", panel = list(T = "CD3"))
  expect_equal(res$panel_match$n_no_marker, 1)
  expect_equal(res$panel_match$extra, "")
})

test_that("panel = NULL skips the match and says so", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p[c("a1", "a2")], "$CYT", panel = NULL)
  expect_null(res$panel_match)
  expect_null(res$provenance$args$panel)
})

test_that("result has provenance and no pass/fail or per-file score", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p, "$CYT", panel = list(T = "CD3"))
  prov <- res$provenance
  expect_true(all(
    c("package_version", "r_version", "seed", "session") %in% names(prov)
  ))
  expect_null(prov$seed)
  all_names <- c(names(res$files), names(res$counts), names(res$panel_match))
  expect_false(any(grepl(
    "pass|fail|score|ok_flag",
    all_names,
    ignore.case = TRUE
  )))
})

test_that("files on disk are not changed", {
  p <- sig_set()
  before <- unname(tools::md5sum(p))
  sw_fcs_header_signatures(p, c("$CYT", "$PnV"), panel = list(T = "CD3"))
  expect_equal(unname(tools::md5sum(p)), before)
})

test_that("channel names that contain commas do not merge with split channels", {
  d <- withr::local_tempdir()
  x <- sig_fcs(d, "x.fcs", "A,B", "CD3")
  y <- sig_fcs(d, "y.fcs", c("A", "B"), c("CD3", "CD4"))
  res <- sw_fcs_header_signatures(c(x, y), "$CYT", panel = NULL)
  expect_equal(length(unique(res$files$signature)), 2)
  expect_equal(nrow(res$counts), 2)
  expect_equal(sort(res$counts$n_channels), c(1L, 2L))
  expect_equal(res$counts$n_files, c(1L, 1L))
})

test_that("channel order does not change the signature", {
  d <- withr::local_tempdir()
  a <- sig_fcs(d, "a.fcs", c("B1-A", "V1-A"), c("CD3", "CD4"))
  b <- sig_fcs(d, "b.fcs", c("V1-A", "B1-A"), c("CD4", "CD3"))
  res <- sw_fcs_header_signatures(c(a, b), "$CYT", panel = NULL)
  expect_equal(length(unique(res$files$signature)), 1)
})

test_that("signature ids run from the largest group, ties by path order", {
  p <- sig_set()
  # The singleton b1 is listed first; the pair a1, a2 must still be sig01.
  res <- sw_fcs_header_signatures(p[c("b1", "a1", "a2")], "$CYT", NULL)
  f <- res$files
  expect_equal(
    f$signature[f$file %in% c("a1.fcs", "a2.fcs")],
    c("sig01", "sig01")
  )
  expect_equal(f$signature[f$file == "b1.fcs"], "sig02")
  expect_equal(res$counts$n_files, c(2L, 1L))
  # Three groups of three sizes: ids follow n_files descending.
  res3 <- sw_fcs_header_signatures(p[c("b1", "c1", "a1", "a2")], "$CYT", NULL)
  expect_true(all(diff(res3$counts$n_files) <= 0))
  # Equal sizes: the group seen first in `paths` takes the lower id.
  tie <- sw_fcs_header_signatures(p[c("c1", "b1")], "$CYT", NULL)
  expect_equal(tie$files$signature, c("sig01", "sig02"))
  tie2 <- sw_fcs_header_signatures(p[c("b1", "c1")], "$CYT", NULL)
  expect_equal(tie2$files$signature, c("sig01", "sig02"))
  expect_match(tie2$counts$files[1], "b1.fcs")
})

test_that("all files unreadable gives an empty counts table and a row per file", {
  res <- sw_fcs_header_signatures(
    c("/nonexistent/a.fcs", "/nonexistent/b.fcs"),
    "$CYT",
    list(T = "CD3")
  )
  expect_equal(nrow(res$counts), 0)
  expect_equal(res$n_not_assessed, 2)
  expect_equal(nrow(res$files), 2)
  expect_true(all(res$files$status == "not_assessed"))
  expect_equal(res$files$reason, rep("file not found", 2))
  expect_equal(nrow(res$panel_match), 2)
})

test_that("panel match is case sensitive: a case mismatch is missing plus extra", {
  d <- withr::local_tempdir()
  f <- sig_fcs(d, "t.fcs", c("B1-A", "V1-A"), c("cd4", "CD3"))
  res <- sw_fcs_header_signatures(f, "$CYT", list(T = c("CD3", "CD4")))
  expect_equal(res$panel_match$missing, "CD4")
  expect_equal(res$panel_match$extra, "cd4")
})

test_that("panel entries are trimmed before matching", {
  d <- withr::local_tempdir()
  f <- sig_fcs(d, "t.fcs", c("B1-A", "V1-A"), c("CD3", "CD4"))
  res <- sw_fcs_header_signatures(f, "$CYT", list(T = c(" CD3 ", "CD4")))
  expect_equal(res$panel_match$n_missing, 0)
})

test_that("keyword values are trimmed before grouping", {
  d <- withr::local_tempdir()
  a <- sig_fcs(d, "a.fcs", "B1-A", "CD3", kw = list("$CYT" = "A"))
  b <- sig_fcs(d, "b.fcs", "B1-A", "CD3", kw = list("$CYT" = "  A "))
  res <- sw_fcs_header_signatures(c(a, b), "$CYT", NULL)
  expect_equal(length(unique(res$files$signature)), 1)
})

test_that("a non-zero $NEXTDATA flags the file and the flowCore warning is muffled", {
  p <- sig_set()
  h <- flowCore::read.FCSheader(p[["a1"]])
  h[[1]][["$NEXTDATA"]] <- "13432728"
  testthat::local_mocked_bindings(
    read.FCSheader = function(...) {
      warning(
        "The file contains 1 additional data segment. The default is to read it."
      )
      h
    },
    .package = "flowCore"
  )
  expect_no_warning(
    res <- sw_fcs_header_signatures(p[c("a1", "a2")], "$CYT", NULL)
  )
  expect_equal(res$files$n_datasets, c(2L, 2L))
  expect_match(res$files$reason, "more than one data set")
  expect_true(all(res$files$status == "assessed"))
})

test_that("single data set files report n_datasets = 1", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p[c("a1", "c1")], "$CYT", NULL)
  expect_equal(res$files$n_datasets, c(1L, 1L))
  expect_equal(res$files$reason, c("", ""))
})

test_that("a keyword found in no file is reported in the keywords table", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(p[c("a1", "b1")], c("$CYT", "$NOPE"), NULL)
  kt <- res$keywords
  expect_equal(kt$keyword, c("$CYT", "$NOPE"))
  expect_equal(kt$n_files_present, c(2L, 0L))
  expect_equal(kt$n_files_absent, c(0L, 2L))
  expect_match(kt$reason[kt$keyword == "$NOPE"], "no file")
  expect_true(is.na(kt$reason[kt$keyword == "$CYT"]))
  expect_false("note" %in% names(kt))
})

test_that("keywords table counts unreadable files apart from absent ones", {
  p <- sig_set()
  res <- sw_fcs_header_signatures(
    c(p["a1"], gone = "/nonexistent/g.fcs"),
    c("$CYT", "$PnV"),
    NULL
  )
  kt <- res$keywords
  expect_equal(kt$n_not_assessed, c(1L, 1L))
  expect_equal(kt$n_files_present, c(1L, 1L))
  expect_equal(kt$n_files_absent, c(0L, 0L))
})

test_that("empty or blank keywords are rejected", {
  expect_error(sw_fcs_header_signatures("x.fcs", "", NULL), "keywords")
  expect_error(
    sw_fcs_header_signatures("x.fcs", c("$CYT", " "), NULL),
    "keywords"
  )
})

test_that("provenance starts with the shared keys (signatures)", {
  pv <- {
    res <- sw_fcs_header_signatures(sig_set(), "$CYT", panel = list(T = "CD3"))
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
