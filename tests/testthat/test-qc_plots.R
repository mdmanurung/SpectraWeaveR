# tests/testthat/test-qc_plots.R
# Tests for R/qc_plots.R. Argument checks run without PeacoQC; two tests run the
# real PeacoQC on a small synthetic file.

make_qc_ff <- function(name = "s1.fcs", n = 3000) {
  set.seed(1)
  m <- cbind(Time = seq_len(n), A = stats::rnorm(n), B = stats::rnorm(n))
  ff <- flowCore::flowFrame(m)
  flowCore::keyword(ff)[["FILENAME"]] <- name
  ff
}

test_that("sw_qc_plot requires output_dir and a PeacoQC result", {
  expect_error(sw_qc_plot(), "output_dir")
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  ff <- make_qc_ff()
  expect_error(
    sw_qc_plot(
      ff,
      c("A", "B"),
      qc = list(x = 1),
      output_dir = withr::local_tempdir()
    ),
    "must come from sw_qc_run"
  )
  expect_error(
    sw_qc_plot(
      "nope",
      c("A", "B"),
      qc = list(),
      output_dir = withr::local_tempdir()
    ),
    "flowFrame"
  )
})

test_that("sw_qc_heatmap checks its input", {
  expect_error(sw_qc_heatmap(), "report_location")
  skip_if_not_installed("PeacoQC")
  expect_error(
    sw_qc_heatmap(withr::local_tempdir()),
    "No PeacoQC report found"
  )
})

test_that("sw_qc_plot writes one PNG and will not overwrite it", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("ggplot2")
  ff <- make_qc_ff()
  out <- withr::local_tempdir()
  res <- suppressMessages(sw_qc_run(ff, channels = c("A", "B")))
  expect_true(!is.null(res$PeacoQC))

  png <- suppressMessages(sw_qc_plot(ff, c("A", "B"), res, output_dir = out))
  expect_length(png, 1)
  expect_true(file.exists(png))
  expect_match(basename(png), "^PeacoQC_s1\\.png$")
  expect_error(
    sw_qc_plot(ff, c("A", "B"), res, output_dir = out),
    "not overwritten"
  )
})

test_that("sw_qc_heatmap reads the report written by a QC run", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("ComplexHeatmap")
  out <- withr::local_tempdir()
  for (nm in c("s1.fcs", "s2.fcs")) {
    suppressMessages(sw_qc_run(
      make_qc_ff(nm),
      channels = c("A", "B"),
      output_dir = out
    ))
  }
  expect_true(file.exists(file.path(
    out,
    "PeacoQC_results",
    "PeacoQC_report.txt"
  )))
  # From the folder and from the file itself.
  expect_false(is.null(suppressMessages(sw_qc_heatmap(out))))
  expect_false(is.null(suppressMessages(
    sw_qc_heatmap(file.path(out, "PeacoQC_results", "PeacoQC_report.txt"))
  )))
})
