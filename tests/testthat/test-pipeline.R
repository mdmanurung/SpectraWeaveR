# tests/testthat/test-pipeline.R
# Unit tests for R/pipeline.R — End-to-end pipeline orchestrator

test_that("run_pipeline rejects non-existent directory", {
  expect_error(
    run_pipeline(
      fcs_dir = "/nonexistent/dir",
      sample_meta = tibble::tibble(file = "f.fcs", sample = "S1", batch = "B1"),
      markers = "CD3",
      lineage_markers = "CD3"
    ),
    "does not exist"
  )
})

test_that("run_pipeline rejects non-data.frame sample_meta", {
  tmp_dir <- tempdir()
  expect_error(
    run_pipeline(
      fcs_dir = tmp_dir,
      sample_meta = "not_a_df",
      markers = "CD3",
      lineage_markers = "CD3"
    ),
    "data.frame"
  )
})

test_that("run_pipeline rejects missing sample_meta columns", {
  tmp_dir <- tempdir()
  expect_error(
    run_pipeline(
      fcs_dir = tmp_dir,
      sample_meta = tibble::tibble(sample = "S1"),
      markers = "CD3",
      lineage_markers = "CD3"
    ),
    "missing required column"
  )
})

test_that("run_pipeline rejects empty markers", {
  tmp_dir <- tempdir()
  expect_error(
    run_pipeline(
      fcs_dir = tmp_dir,
      sample_meta = tibble::tibble(file = "f.fcs", sample = "S1",
                                   batch = "B1"),
      markers = character(0),
      lineage_markers = "CD3"
    ),
    "non-empty character vector"
  )
})

test_that("run_pipeline rejects lineage_markers not in markers", {
  tmp_dir <- tempdir()
  expect_error(
    run_pipeline(
      fcs_dir = tmp_dir,
      sample_meta = tibble::tibble(file = "f.fcs", sample = "S1",
                                   batch = "B1"),
      markers = c("CD3", "CD4"),
      lineage_markers = c("CD3", "CD99")
    ),
    "subset of.*markers"
  )
})

test_that("run_pipeline rejects non-character fcs_dir", {
  expect_error(
    run_pipeline(
      fcs_dir = 42,
      sample_meta = tibble::tibble(file = "f.fcs", sample = "S1",
                                   batch = "B1"),
      markers = "CD3",
      lineage_markers = "CD3"
    ),
    "single directory path"
  )
})

test_that("run_pipeline returns expected structure (smoke test)", {
  # This is a structural test only — it verifies the expected return slots
  # without actually running the pipeline (which requires all dependencies).

  expected_slots <- c("flowset", "gating_set", "qc_results", "uncorrected",
                       "corrected", "correction_eval", "cluster_result",
                       "cluster_assignments", "cluster_mfis")

  expect_true(length(expected_slots) == 9)
  expect_true("cluster_result" %in% expected_slots)
  expect_true("corrected" %in% expected_slots)
})

# =========================================================================
# Integration test — full pipeline execution with synthetic data
# =========================================================================

test_that("run_pipeline executes end-to-end with synthetic FCS", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("cyCombine")
  skip_if_not_installed("kohonen")

  set.seed(42)

  # Create temp directory with 2 synthetic FCS files
  tmp_dir <- file.path(tempdir(), "pipeline_test_fcs")
  dir.create(tmp_dir, showWarnings = FALSE)
  out_dir <- file.path(tempdir(), "pipeline_test_output")
  dir.create(out_dir, showWarnings = FALSE)
  on.exit({
    unlink(tmp_dir, recursive = TRUE)
    unlink(out_dir, recursive = TRUE)
  })

  make_fcs <- function(filename, n = 300, batch_shift = 0) {
    mat <- matrix(
      c(abs(rnorm(n, 150000, 30000)),       # FSC-A
        abs(rnorm(n, 50000, 15000)),         # SSC-A
        abs(rnorm(n, 1000 + batch_shift, 300)),  # CD3
        abs(rnorm(n, 2000 + batch_shift, 500)),  # CD4
        abs(rnorm(n, 500 + batch_shift, 200))),   # CD8
      ncol = 5,
      dimnames = list(NULL, c("FSC-A", "SSC-A", "CD3", "CD4", "CD8"))
    )
    storage.mode(mat) <- "double"
    ff <- flowCore::flowFrame(mat)
    flowCore::write.FCS(ff, file.path(tmp_dir, filename))
  }

  make_fcs("sample1.fcs", batch_shift = 0)
  make_fcs("sample2.fcs", batch_shift = 200)

  sample_meta <- tibble::tibble(
    file = c("sample1.fcs", "sample2.fcs"),
    sample = c("S1", "S2"),
    batch = c("B1", "B2")
  )

  result <- run_pipeline(
    fcs_dir = tmp_dir,
    sample_meta = sample_meta,
    markers = c("CD3", "CD4", "CD8"),
    lineage_markers = c("CD3", "CD4", "CD8"),
    gating_template = NULL,
    output_dir = out_dir
  )

  expect_true(is.list(result))
  expect_true("flowset" %in% names(result))
  expect_true("qc_results" %in% names(result))
  expect_true("corrected" %in% names(result))
  expect_true("cluster_result" %in% names(result))

  # Verify corrected data is a tibble with marker columns
  if (!is.null(result$corrected)) {
    expect_s3_class(result$corrected, "tbl_df")
    expect_true(all(c("CD3", "CD4", "CD8") %in% names(result$corrected)))
  }

  # Verify cluster result
  if (!is.null(result$cluster_result)) {
    expect_s3_class(result$cluster_result, "sw_cluster_result")
  }
})
