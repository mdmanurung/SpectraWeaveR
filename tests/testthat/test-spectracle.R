# tests/testthat/test-spectracle.R
# spectracle is not needed for validation or pure-converter tests.

# Mock of spectracle() output: metadata columns, numeric detector columns,
# hash.md5, and an appended AF row. Detector names use flowstate aliases.
mock_spectracle <- function() {
  df <- data.frame(
    deriving.function = "spectracle_0.2.0",
    sample.id = c("CD3 BV421 (Cells)", "CD4 FITC (Cells)", "Unstained (Cells)"),
    tissue.type = factor("Cells"),
    N = c("BV421", "FITC", "AF"),
    S = c("CD3", "CD4", NA),
    detector = factor(c("V1", "B1", "V1")),
    laser = factor(c("V", "B", "V")),
    V1 = c(1, 0.1, 0.4),
    V2 = c(0.5, 0, 0.2),
    B1 = c(0, 1, 1),
    hash.md5 = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  df[["$CYT"]] <- "Aurora"
  df
}

# ============================================================================
# sw_unmix_spectracle validation
# ============================================================================

test_that("sw_unmix_spectracle validates raw_controls", {
  expect_error(sw_unmix_spectracle(1, af = "drop"), "raw_controls")
  expect_error(sw_unmix_spectracle(character(0), af = "drop"), "raw_controls")
  expect_error(
    sw_unmix_spectracle("/no/such/dir", af = "drop"),
    "existing directory"
  )
  expect_error(
    sw_unmix_spectracle(c("/no/a.fcs", "/no/b.fcs"), af = "drop"),
    "not found"
  )
})

test_that("sw_unmix_spectracle refuses a single FCS file", {
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  on.exit(unlink(f))
  expect_error(sw_unmix_spectracle(f, af = "drop"), "existing directory")
})

test_that("sw_unmix_spectracle requires af", {
  expect_error(sw_unmix_spectracle(tempdir()), "'af' is required")
  expect_error(sw_unmix_spectracle(tempdir(), af = "yes"), "'af' is required")
  expect_error(
    sw_unmix_spectracle(tempdir(), af = c("drop", "keep_as_row")),
    "'af' is required"
  )
})

test_that("sw_unmix_spectracle validates optional arguments", {
  expect_error(
    sw_unmix_spectracle(tempdir(), "drop", unstained_ids = "x"),
    "named character"
  )
  expect_error(
    sw_unmix_spectracle(tempdir(), "drop", filter_top_expressing = 1),
    "character vector"
  )
  expect_error(
    sw_unmix_spectracle(
      tempdir(),
      "drop",
      top_expressing_override = c(a = "1")
    ),
    "named numeric"
  )
  expect_error(
    sw_unmix_spectracle(tempdir(), "drop", rename = c("a", "b")),
    "named character"
  )
})

test_that("sw_unmix_spectracle reports a missing spectracle package", {
  skip_if(requireNamespace("spectracle", quietly = TRUE))
  expect_error(
    sw_unmix_spectracle(tempdir(), af = "drop"),
    "'spectracle' is required"
  )
})

test_that("sw_unmix_spectracle runs on real controls", {
  skip_if_not_installed("spectracle")
  skip("needs a directory of raw Aurora single-color controls")
})

# ============================================================================
# sw_io_spectracle_to_spectra validation
# ============================================================================

test_that("converter validates x", {
  expect_error(sw_io_spectracle_to_spectra(1, af = "drop"), "data.frame")
  expect_error(
    sw_io_spectracle_to_spectra(data.frame(V1 = 1), af = "drop"),
    "'N'"
  )
})

test_that("converter requires af", {
  expect_error(
    sw_io_spectracle_to_spectra(mock_spectracle()),
    "'af' is required"
  )
  expect_error(
    sw_io_spectracle_to_spectra(mock_spectracle(), af = "keep"),
    "'af' is required"
  )
})

test_that("converter validates detectors, rename and tissue_type", {
  x <- mock_spectracle()
  expect_error(
    sw_io_spectracle_to_spectra(x, "drop", detectors = "nope"),
    "not found"
  )
  expect_error(
    sw_io_spectracle_to_spectra(x, "drop", detectors = "S"),
    "not numeric"
  )
  expect_error(
    sw_io_spectracle_to_spectra(x, "drop", detectors = c("V1", "V1")),
    "unique"
  )
  expect_error(
    sw_io_spectracle_to_spectra(x, "drop", rename = c("a", "b")),
    "named character"
  )
  expect_error(
    sw_io_spectracle_to_spectra(x, "drop", tissue_type = "Beads"),
    "No rows"
  )
  expect_error(
    sw_io_spectracle_to_spectra(
      x[, names(x) != "tissue.type"],
      "drop",
      tissue_type = "Cells"
    ),
    "tissue.type"
  )
})

test_that("converter warns on values outside [0, 1]", {
  x <- mock_spectracle()
  x$V1[1] <- 4e6
  expect_warning(sw_io_spectracle_to_spectra(x, "drop"), "\\[0, 1\\]")
  x <- mock_spectracle()
  x$V2[2] <- NA
  expect_error(sw_io_spectracle_to_spectra(x, "drop"), "finite")
})

test_that("converter rejects all-zero spectra", {
  x <- mock_spectracle()
  x[1, c("V1", "V2", "B1")] <- 0
  expect_error(sw_io_spectracle_to_spectra(x, "drop"), "all-zero")
})

# ============================================================================
# sw_io_spectracle_to_spectra: the three af modes
# ============================================================================

test_that("af = 'drop' returns fluorophores only, order preserved", {
  m <- sw_io_spectracle_to_spectra(mock_spectracle(), af = "drop")
  expect_true(is.matrix(m))
  expect_equal(rownames(m), c("BV421", "FITC"))
  expect_equal(colnames(m), c("V1", "V2", "B1"))
  expect_equal(unname(m["BV421", ]), c(1, 0.5, 0))
})

test_that("af = 'keep_as_row' keeps AF as the last row", {
  m <- sw_io_spectracle_to_spectra(mock_spectracle(), af = "keep_as_row")
  expect_equal(rownames(m), c("BV421", "FITC", "AF"))
  expect_equal(unname(m["AF", ]), c(0.4, 0.2, 1))
})

test_that("af = 'separate' returns spectra and af", {
  out <- sw_io_spectracle_to_spectra(mock_spectracle(), af = "separate")
  expect_named(out, c("spectra", "af"))
  expect_equal(rownames(out$spectra), c("BV421", "FITC"))
  expect_equal(rownames(out$af), "AF")
  expect_identical(colnames(out$af), colnames(out$spectra))
})

test_that("af modes other than drop need an AF row", {
  x <- mock_spectracle()[1:2, ]
  expect_error(sw_io_spectracle_to_spectra(x, "keep_as_row"), "no row")
  expect_error(sw_io_spectracle_to_spectra(x, "separate"), "no row")
  expect_equal(nrow(sw_io_spectracle_to_spectra(x, "drop")), 2)
})

test_that("only-AF input fails for drop and separate", {
  x <- mock_spectracle()[3, ]
  expect_error(sw_io_spectracle_to_spectra(x, "drop"), "No fluorophore rows")
  expect_error(sw_io_spectracle_to_spectra(x, "separate"), "No fluorophore")
})

test_that("several AF rows are named by tissue type", {
  x <- rbind(mock_spectracle(), mock_spectracle()[3, ])
  x$tissue.type <- factor(c("Cells", "Cells", "Cells", "Beads"))
  m <- sw_io_spectracle_to_spectra(x, "keep_as_row")
  expect_equal(rownames(m), c("BV421", "FITC", "AF_Cells", "AF_Beads"))
})

# ============================================================================
# names, order and scale
# ============================================================================

test_that("rename maps by name and keeps order", {
  m <- sw_io_spectracle_to_spectra(
    mock_spectracle(),
    "drop",
    rename = c(B1 = "B1-A", V2 = "V2-A", V1 = "V1-A")
  )
  expect_equal(colnames(m), c("V1-A", "V2-A", "B1-A"))
  expect_equal(unname(m["FITC", ]), c(0.1, 0, 1))
})

test_that("rename must cover every detector", {
  expect_error(
    sw_io_spectracle_to_spectra(
      mock_spectracle(),
      "drop",
      rename = c(V1 = "V1-A")
    ),
    "no entry"
  )
})

test_that("rename may not create duplicate names", {
  expect_error(
    sw_io_spectracle_to_spectra(
      mock_spectracle(),
      "drop",
      rename = c(V1 = "X", V2 = "X", B1 = "B1-A")
    ),
    "unique"
  )
})

test_that("explicit detectors fix the column order", {
  m <- sw_io_spectracle_to_spectra(
    mock_spectracle(),
    "drop",
    detectors = c("B1", "V1")
  )
  expect_equal(colnames(m), c("B1", "V1"))
  expect_equal(unname(m["FITC", ]), c(1, 0.1))
})

test_that("duplicate fluorophore names fail unless tissue_type filters", {
  x <- mock_spectracle()
  dup <- x[1, ]
  dup$tissue.type <- "Beads"
  x$tissue.type <- as.character(x$tissue.type)
  x <- rbind(x, dup)
  expect_error(sw_io_spectracle_to_spectra(x, "drop"), "Duplicate row names")
  m <- sw_io_spectracle_to_spectra(x, "drop", tissue_type = "Cells")
  expect_equal(rownames(m), c("BV421", "FITC"))
})

test_that("rows with a peak below 1 are rescaled to peak 1 and reported", {
  x <- mock_spectracle()
  x[1, c("V1", "V2", "B1")] <- c(0.5, 0.25, 0)
  expect_message(m <- sw_io_spectracle_to_spectra(x, "drop"), "rescaled")
  expect_equal(unname(m["BV421", ]), c(1, 0.5, 0))
  expect_equal(unname(apply(abs(m), 1, max)), c(1, 1))
})

test_that("rows already at peak 1 are unchanged and silent", {
  expect_no_message(m <- sw_io_spectracle_to_spectra(mock_spectracle(), "drop"))
  expect_equal(unname(m["FITC", ]), c(0.1, 0, 1))
})

test_that("a tibble input works", {
  m <- sw_io_spectracle_to_spectra(tibble::as_tibble(mock_spectracle()), "drop")
  expect_equal(dim(m), c(2L, 3L))
})
