# tests/testthat/test-flowstate_io.R
# flowstate, flowCore are not needed for validation or pure-converter tests.
# Real calls are gated with skip_if_not_installed().

# Mock with the real flowstate shape: S3 list of data.frames (data.tables in
# the real object); columns carry "N"/"alias" attributes; sample.id is a factor.
mock_flowstate <- function(transformed = FALSE) {
  x <- c(1e6, 4e6, 2e5, 3e3)
  attr(x, "N") <- "UV1-A"
  attr(x, "alias") <- c(N = "UV1-A", N.alias = "UV1", S.alias = "UV1")
  if (transformed) {
    attr(x, "transformed") <- list(transform.func = "asinh", cofactor = 6000)
  }
  structure(
    list(
      data = data.frame(
        `FSC_A` = c(1, 2, 3, 4),
        UV1 = x,
        sample.id = factor(c("a", "a", "b", "b"))
      ),
      parameters = data.frame(par = c("$P1", "$P2"), N = c("FSC-A", "UV1-A")),
      keywords = data.frame(sample.id = factor(c("a", "b")))
    ),
    class = "flowstate"
  )
}

# ============================================================================
# sw_io_read_flowstate
# ============================================================================

test_that("sw_io_read_flowstate validates files", {
  expect_error(
    sw_io_read_flowstate(character(0), colnames_type = "N"),
    "non-empty character"
  )
  expect_error(
    sw_io_read_flowstate(1, colnames_type = "N"),
    "non-empty character"
  )
  expect_error(
    sw_io_read_flowstate("/no/such.fcs", colnames_type = "N"),
    "not found"
  )
})

test_that("sw_io_read_flowstate requires colnames_type", {
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  on.exit(unlink(f))
  expect_error(sw_io_read_flowstate(f), "colnames_type")
  expect_error(sw_io_read_flowstate(f, colnames_type = "X"), "colnames_type")
  expect_error(
    sw_io_read_flowstate(f, colnames_type = c("N", "S")),
    "colnames_type"
  )
})

test_that("sw_io_read_flowstate validates concatenate", {
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  on.exit(unlink(f))
  expect_error(sw_io_read_flowstate(f, "N", concatenate = NA), "TRUE or FALSE")
  expect_error(sw_io_read_flowstate(f, "N", concatenate = TRUE), "at least two")
})

test_that("sw_io_read_flowstate reports a missing flowstate package", {
  skip_if(requireNamespace("flowstate", quietly = TRUE))
  f <- tempfile(fileext = ".fcs")
  file.create(f)
  on.exit(unlink(f))
  expect_error(sw_io_read_flowstate(f, "N"), "'flowstate' is required")
})

test_that("sw_io_read_flowstate reads a bundled example file", {
  skip_if_not_installed("flowstate")
  files <- list.files(
    system.file("extdata", package = "flowstate"),
    pattern = "BLOCK.*fcs$",
    full.names = TRUE
  )
  skip_if(length(files) == 0)
  fs <- sw_io_read_flowstate(files[1], colnames_type = "S")
  expect_s3_class(fs, "flowstate")
})

# ============================================================================
# flowstate validation
# ============================================================================

test_that("flowstate functions reject non-flowstate input", {
  expect_error(sw_io_flowstate_to_tibble(data.frame(a = 1)), "flowstate")
  expect_error(sw_io_flowstate_to_ff(list()), "flowstate")
  expect_error(sw_io_flowstate_transform(1, cofactor = 6000), "flowstate")
  expect_error(sw_io_write_flowstate("x", file_dir = tempdir()), "flowstate")
})

test_that("validation names missing flowstate elements", {
  bad <- mock_flowstate()
  bad$keywords <- NULL
  expect_error(sw_io_flowstate_to_tibble(bad), "keywords")
  bad <- mock_flowstate()
  bad$data$sample.id <- NULL
  expect_error(sw_io_flowstate_to_tibble(bad), "sample.id")
  bad <- mock_flowstate()
  bad$data <- bad$data[0, ]
  expect_error(sw_io_flowstate_to_tibble(bad), "no events")
})

# ============================================================================
# sw_io_flowstate_transform
# ============================================================================

test_that("sw_io_flowstate_transform requires a positive cofactor", {
  fs <- mock_flowstate()
  expect_error(sw_io_flowstate_transform(fs), "cofactor")
  expect_error(sw_io_flowstate_transform(fs, cofactor = 0), "positive")
  expect_error(sw_io_flowstate_transform(fs, cofactor = "6000"), "positive")
  expect_error(sw_io_flowstate_transform(fs, cofactor = c(1, 2)), "positive")
})

test_that("sw_io_flowstate_transform validates j before dependency check", {
  fs <- mock_flowstate()
  expect_error(sw_io_flowstate_transform(fs, 6000, j = 1), "character vector")
  expect_error(sw_io_flowstate_transform(fs, 6000, j = "nope"), "nope")
})

test_that("sw_io_flowstate_transform reports a missing flowstate package", {
  skip_if(requireNamespace("flowstate", quietly = TRUE))
  expect_error(
    sw_io_flowstate_transform(mock_flowstate(), 6000),
    "'flowstate' is required"
  )
})

# ============================================================================
# sw_io_write_flowstate
# ============================================================================

test_that("sw_io_write_flowstate validates file_dir and suffix", {
  fs <- mock_flowstate()
  expect_error(sw_io_write_flowstate(fs), "file_dir")
  expect_error(sw_io_write_flowstate(fs, file_dir = 1), "file_dir")
  expect_error(
    sw_io_write_flowstate(fs, file_dir = "/no/such/dir"),
    "does not exist"
  )
  expect_error(
    sw_io_write_flowstate(fs, tempdir(), fil_suffix = 1),
    "fil_suffix"
  )
  expect_error(
    sw_io_write_flowstate(fs, tempdir(), endianness = "middle"),
    "should be one of"
  )
})

test_that("sw_io_write_flowstate reports a missing flowstate package", {
  skip_if(requireNamespace("flowstate", quietly = TRUE))
  expect_error(
    sw_io_write_flowstate(mock_flowstate(), tempdir()),
    "'flowstate' is required"
  )
})

# ============================================================================
# sw_io_flowstate_to_tibble (pure)
# ============================================================================

test_that("sw_io_flowstate_to_tibble converts shape and drops attributes", {
  out <- sw_io_flowstate_to_tibble(mock_flowstate())
  expect_s3_class(out, "tbl_df")
  expect_equal(names(out), c("FSC_A", "UV1", "sample"))
  expect_equal(nrow(out), 4)
  expect_type(out$sample, "character")
  expect_null(attr(out$UV1, "alias"))
  expect_null(attr(out, "sw_transform"))
})

test_that("sw_io_flowstate_to_tibble keeps Aurora-scale values exactly", {
  out <- sw_io_flowstate_to_tibble(mock_flowstate())
  expect_equal(out$UV1, c(1e6, 4e6, 2e5, 3e3))
})

test_that("sw_io_flowstate_to_tibble records applied transforms", {
  out <- sw_io_flowstate_to_tibble(mock_flowstate(transformed = TRUE))
  tr <- attr(out, "sw_transform")
  expect_named(tr, "UV1")
  expect_equal(tr$UV1$cofactor, 6000)
  expect_equal(tr$UV1$transform.func, "asinh")
})

# ============================================================================
# sw_io_flowstate_to_ff
# ============================================================================

test_that("sw_io_flowstate_to_ff reports a missing flowCore package", {
  skip_if(requireNamespace("flowCore", quietly = TRUE))
  expect_error(
    sw_io_flowstate_to_ff(mock_flowstate()),
    "'flowCore' is required"
  )
})

test_that("sw_io_flowstate_to_ff builds frames without truncation", {
  skip_if_not_installed("flowCore")
  fs <- mock_flowstate()
  fset <- sw_io_flowstate_to_ff(fs)
  expect_s4_class(fset, "flowSet")
  expect_equal(sort(flowCore::sampleNames(fset)), c("a", "b"))
  ff <- fset[[1]]
  expect_equal(flowCore::colnames(ff), c("FSC_A", "UV1"))
  expect_equal(max(flowCore::exprs(fset[["b"]])[, "UV1"], na.rm = TRUE), 2e5)

  one <- fs
  one$data <- one$data[one$data$sample.id == "a", ]
  one$data$sample.id <- droplevels(one$data$sample.id)
  ff1 <- sw_io_flowstate_to_ff(one)
  expect_s4_class(ff1, "flowFrame")
  expect_equal(max(flowCore::exprs(ff1)[, "UV1"]), 4e6)
})
