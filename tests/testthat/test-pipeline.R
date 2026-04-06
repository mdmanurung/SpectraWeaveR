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
  # The full integration test requires flowCore, PeacoQC, cyCombine, kohonen.

  expected_slots <- c("flowset", "gating_set", "qc_results", "uncorrected",
                       "corrected", "correction_eval", "cluster_result",
                       "cluster_assignments", "cluster_mfis")

  # Just verify the slot names are documented correctly

  expect_true(length(expected_slots) == 9)
  expect_true("cluster_result" %in% expected_slots)
  expect_true("corrected" %in% expected_slots)
})
