# tests/testthat/test-import_dye_map.R
# Unit tests for R/import_dye_map.R. Writes tiny FCS files to a temp folder.

# Write one FCS file; `marker` becomes $PnS.
make_dye_fcs <- function(dir, name, chan, marker) {
  skip_if_not_installed("flowCore")
  m <- matrix(
    stats::rnorm(20 * length(chan)),
    ncol = length(chan),
    dimnames = list(NULL, chan)
  )
  ff <- flowCore::flowFrame(m)
  flowCore::pData(flowCore::parameters(ff))$desc <- marker
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

# a, b in batch B1; c in B2 where CD4 carries another dye. d has an unsplit name.
make_dye_set <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  chan <- c("B1-A", "V1-A")
  p <- c(
    a = make_dye_fcs(d, "a.fcs", chan, c("CD3 BV421", "CD4 BB515")),
    b = make_dye_fcs(d, "b.fcs", chan, c("CD3 BV421", "CD4 BB515")),
    c = make_dye_fcs(d, "c.fcs", chan, c("CD3 BV421", "CD4 PE")),
    d = make_dye_fcs(d, "d.fcs", chan, c("CD3 BV421", "CD8-APC"))
  )
  p
}
dye_batch <- c(
  a.fcs = "B1",
  b.fcs = "B1",
  c.fcs = "B2",
  d.fcs = "B2"
)

test_that("required arguments have no default", {
  expect_error(sw_fcs_marker_dye_map(), "no default.*'paths'")
  p <- make_dye_set()
  expect_error(
    sw_fcs_marker_dye_map(unname(p), separator = " "),
    "no default.*'batch'"
  )
  expect_error(
    sw_fcs_marker_dye_map(unname(p), batch = dye_batch),
    "no default.*'separator'"
  )
})

test_that("batch must cover every file by name, with no extras", {
  p <- make_dye_set()
  expect_error(
    sw_fcs_marker_dye_map(unname(p), unname(dye_batch), " "),
    "named"
  )
  expect_error(
    sw_fcs_marker_dye_map(unname(p), dye_batch[-1], " "),
    "a.fcs"
  )
  expect_error(
    sw_fcs_marker_dye_map(unname(p), c(dye_batch, z.fcs = "B3"), " "),
    "z.fcs"
  )
  expect_error(
    sw_fcs_marker_dye_map(unname(p), replace(dye_batch, 1, NA), " "),
    "NA"
  )
  expect_error(
    sw_fcs_marker_dye_map(unname(p), dye_batch, c(" ", "-")),
    "separator"
  )
})

test_that("names split into marker and dye; the rest are listed, not guessed", {
  p <- make_dye_set()
  res <- sw_fcs_marker_dye_map(unname(p), dye_batch, separator = " ")
  expect_s3_class(res, "sw_fcs_marker_dye_map")
  long <- res$map
  # Every channel of every file has a row.
  expect_equal(nrow(long), 8)
  expect_true(all(
    c("file", "batch", "channel", "name", "marker", "dye", "status") %in%
      names(long)
  ))
  r <- long[long$file == "a.fcs" & long$channel == "B1-A", ]
  expect_equal(c(r$marker, r$dye, r$batch), c("CD3", "BV421", "B1"))
  # CD8-APC does not split on a space: listed, marker and dye stay NA.
  bad <- long[long$file == "d.fcs" & long$channel == "V1-A", ]
  expect_equal(bad$status, "not_split")
  expect_true(is.na(bad$marker) && is.na(bad$dye))
  expect_equal(res$unsplit$name, "CD8-APC")
  expect_equal(res$unsplit$file, "d.fcs")
})

test_that("marker x batch table counts dyes and flags a marker with two dyes", {
  p <- make_dye_set()
  res <- sw_fcs_marker_dye_map(unname(p), dye_batch, " ")
  mb <- res$marker_batch
  expect_true(all(
    c(
      "marker",
      "batch",
      "dyes",
      "n_files",
      "n_dyes_in_batch",
      "n_dyes_marker"
    ) %in%
      names(mb)
  ))
  cd4 <- mb[mb$marker == "CD4", ]
  expect_setequal(cd4$dyes, c("BB515", "PE"))
  expect_equal(unique(cd4$n_dyes_marker), 2L)
  expect_equal(cd4$n_files[cd4$batch == "B1"], 2L)
  cd3 <- mb[mb$marker == "CD3", ]
  expect_equal(unique(cd3$n_dyes_marker), 1L)
  # No verdict column.
  expect_false(any(grepl("pass|fail|ok", names(mb), ignore.case = TRUE)))
})

test_that("an unreadable file becomes a not_assessed row, not a crash", {
  p <- make_dye_set()
  junk <- file.path(dirname(p[["a"]]), "junk.fcs")
  writeLines("not an fcs file", junk)
  paths <- c(unname(p), junk, file.path(dirname(junk), "gone.fcs"))
  b <- c(dye_batch, junk.fcs = "B1", gone.fcs = "B2")
  res <- sw_fcs_marker_dye_map(paths, b, " ")
  na_rows <- res$map[res$map$status == "not_assessed", ]
  expect_setequal(na_rows$file, c("junk.fcs", "gone.fcs"))
  expect_equal(na_rows$reason[na_rows$file == "junk.fcs"], "not an FCS file")
  expect_equal(na_rows$reason[na_rows$file == "gone.fcs"], "file not found")
  expect_true(all(c("junk.fcs", "gone.fcs") %in% res$unsplit$file))
  # Good files still processed.
  expect_equal(sum(res$map$status == "split"), 7)
})

test_that("every channel is accounted for and the file is not edited", {
  p <- make_dye_set()
  before <- unname(tools::md5sum(p))
  res <- sw_fcs_marker_dye_map(unname(p), dye_batch, " ")
  expect_equal(unname(tools::md5sum(p)), before)
  expect_equal(
    nrow(res$map),
    nrow(res$unsplit) + sum(res$map$status == "split")
  )
})

test_that("provenance is attached", {
  p <- make_dye_set()
  res <- sw_fcs_marker_dye_map(unname(p), dye_batch, " ")
  prov <- attr(res, "provenance")
  expect_true(all(
    c("package_version", "separator", "seed", "r_version") %in% names(prov)
  ))
  expect_equal(prov$separator, " ")
  expect_null(prov$seed)
})

test_that("provenance holds a session summary and the call arguments", {
  p <- make_dye_set()
  res <- sw_fcs_marker_dye_map(unname(p), dye_batch, " ")
  prov <- attr(res, "provenance")
  expect_true(all(c("packages", "args", "files") %in% names(prov)))
  expect_true("flowCore" %in% prov$packages$package)
  expect_equal(prov$args$batch, dye_batch)
  expect_equal(nrow(prov$files), 4)
})

test_that("a channel with an empty $PnS is listed, not split on its $PnN", {
  d <- withr::local_tempdir()
  p <- make_dye_fcs(d, "e.fcs", c("FSC-A", "B1-A"), c(NA, "CD3-BV421"))
  res <- sw_fcs_marker_dye_map(p, c(e.fcs = "B1"), "-")
  fsc <- res$map[res$map$channel == "FSC-A", ]
  expect_equal(fsc$source, "$PnN")
  expect_equal(fsc$status, "not_split")
  expect_true(is.na(fsc$marker) && is.na(fsc$dye))
  expect_match(fsc$reason, "PnS")
  expect_false("FSC" %in% res$marker_batch$marker)
  expect_equal(res$unsplit$channel, "FSC-A")
  # The other channel still splits from $PnS.
  expect_equal(res$marker_batch$marker, "CD3")
})

test_that("path checks stop on empty, NA, non-character and repeated names", {
  p <- make_dye_set()
  expect_error(sw_fcs_marker_dye_map(character(0), dye_batch, " "), "paths")
  expect_error(sw_fcs_marker_dye_map(NA_character_, dye_batch, " "), "paths")
  expect_error(sw_fcs_marker_dye_map(1:2, dye_batch, " "), "paths")
  sub <- file.path(dirname(p[["a"]]), "sub")
  dir.create(sub)
  file.copy(p[["a"]], file.path(sub, "a.fcs"))
  expect_error(
    sw_fcs_marker_dye_map(
      c(unname(p), file.path(sub, "a.fcs")),
      dye_batch,
      " "
    ),
    "unique"
  )
})

test_that("repeated names in batch stop", {
  p <- make_dye_set()
  expect_error(
    sw_fcs_marker_dye_map(unname(p), c(dye_batch, a.fcs = "B9"), " "),
    "unique"
  )
})

test_that("an invalid regular expression stops before any file is read", {
  p <- make_dye_set()
  expect_error(
    expect_no_warning(sw_fcs_marker_dye_map(unname(p), dye_batch, "(")),
    "separator.*regular expression"
  )
})

test_that("batch labels follow the file name, not the position", {
  p <- make_dye_set()
  res <- sw_fcs_marker_dye_map(unname(p), rev(dye_batch), " ")
  got <- unique(res$map[, c("file", "batch")])
  expect_equal(setNames(got$batch, got$file)[names(dye_batch)], dye_batch)
})

test_that("labels need exactly two non-empty parts; parts are trimmed", {
  d <- withr::local_tempdir()
  chan <- c("C1-A", "C2-A", "C3-A", "C4-A")
  p <- make_dye_fcs(
    d,
    "l.fcs",
    chan,
    c("CD3 BV421 X", "CD3  BV421", "-BV421", "CD3 - BV421")
  )
  m <- sw_fcs_marker_dye_map(p, c(l.fcs = "B1"), " ")$map
  expect_equal(m$status, rep("not_split", 4))
  expect_match(m$reason[1], "3 part")
  expect_match(m$reason[2], "3 part")
  expect_true(all(is.na(m$marker[1:2]) & is.na(m$dye[1:2])))
  m2 <- sw_fcs_marker_dye_map(p, c(l.fcs = "B1"), "-")$map
  # "-BV421" has an empty marker: two parts, one empty, so not split.
  expect_equal(m2$status[3], "not_split")
  expect_true(is.na(m2$marker[3]))
  # White space around the parts is trimmed.
  expect_equal(m2$status[4], "split")
  expect_equal(c(m2$marker[4], m2$dye[4]), c("CD3", "BV421"))
})

test_that("a multi-dataset file reports the datasets it did not read", {
  p <- make_dye_set()
  h1 <- flowCore::read.FCSheader(p[["a"]])[[1]]
  testthat::local_mocked_bindings(
    read.FCSheader = function(...) list(h1, h1, h1),
    .package = "flowCore"
  )
  res <- sw_fcs_marker_dye_map(p[["a"]], c(a.fcs = "B1"), " ")
  skipped <- res$map[res$map$status == "not_assessed", ]
  expect_equal(nrow(skipped), 2)
  expect_match(skipped$reason, "dataset")
  expect_equal(sum(res$map$status == "split"), 2)
  expect_equal(nrow(res$unsplit), 2)
})

test_that("channel_index follows $Pn, not the order of keys in the header", {
  p <- make_dye_set()[["a"]]
  # $P10N is stored before $P1N, as in alphabetically written CyTOF files.
  fake <- list(c(
    "$P10N" = "Z-A",
    "$P10S" = "CD8 APC",
    "$P1N" = "Time",
    "$P1S" = "",
    "$P2N" = "Y-A",
    "$P2S" = "CD3 BV421"
  ))
  testthat::local_mocked_bindings(.sw_read_fcs_header = \(path, ...) fake)
  res <- sw_fcs_marker_dye_map(p, c(a.fcs = "B1"), " ")
  expect_equal(res$map$channel, c("Time", "Y-A", "Z-A"))
  expect_equal(res$map$channel_index, c(1L, 2L, 10L))
  expect_equal(res$map$marker, c(NA, "CD3", "CD8"))
})

test_that("the batch error names the expected names", {
  p <- make_dye_set()
  expect_error(
    sw_fcs_marker_dye_map(unname(p), unname(dye_batch), " "),
    "basename"
  )
})

test_that("provenance starts with the shared keys (dye_map)", {
  pv <- attr(
    sw_fcs_marker_dye_map(unname(make_dye_set()), dye_batch, " "),
    "provenance"
  )
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
