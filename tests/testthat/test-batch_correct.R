# tests/testthat/test-batch_correct.R
# Unit tests for R/batch_correct.R — Batch correction (cyCombine)

# ===========================================================================
# sw_correct_prepare
# ===========================================================================

test_that("sw_correct_prepare creates correct tibble", {
  skip_if_not_installed("flowCore")

  # Create mock flowFrames
  mat1 <- matrix(c(1000, 2000, 3000, 4000, 5000, 6000), ncol = 2,
                 dimnames = list(NULL, c("CD3", "CD4")))
  mat2 <- matrix(c(1500, 2500, 3500, 4500, 5500, 6500), ncol = 2,
                 dimnames = list(NULL, c("CD3", "CD4")))

  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)

  ff_list <- list(S1 = ff1, S2 = ff2)

  sample_meta <- tibble::tibble(
    sample = c("S1", "S2"),
    batch = c("B1", "B2"),
    condition = c("Ctrl", "Stim")
  )

  result <- sw_correct_prepare(ff_list, sample_meta,
                                      markers = c("CD3", "CD4"),
                                      cofactor = 6000)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("CD3", "CD4", "sample", "batch", "condition", ".cell_id") %in%
                    names(result)))
  expect_equal(nrow(result), 6) # 3 cells per sample * 2 samples

  # Check arcsinh transformation was applied
  expect_equal(result$CD3[1], asinh(1000 / 6000), tolerance = 1e-10)
  expect_equal(result$CD4[1], asinh(4000 / 6000), tolerance = 1e-10)

  # Check metadata
  expect_true(all(result$sample[1:3] == "S1"))
  expect_true(all(result$batch[1:3] == "B1"))
})

test_that("sw_correct_prepare rejects unnamed list", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_correct_prepare(list(ff), tibble::tibble(sample = "S1",
                                                       batch = "B1"),
                              markers = c("CD3")),
    "named list"
  )
})

test_that("sw_correct_prepare rejects missing sample_meta columns", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_correct_prepare(list(S1 = ff),
                              tibble::tibble(sample = "S1"),
                              markers = c("CD3")),
    "missing required column"
  )
})

test_that("sw_correct_prepare rejects missing sample in metadata", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_correct_prepare(list(S1 = ff),
                              tibble::tibble(sample = "S2", batch = "B1"),
                              markers = c("CD3")),
    "not found in sample_meta"
  )
})

test_that("sw_correct_prepare validates cofactor", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_correct_prepare(list(S1 = ff),
                              tibble::tibble(sample = "S1", batch = "B1"),
                              markers = c("CD3"),
                              cofactor = -1),
    "positive number"
  )
})

# ===========================================================================
# sw_correct_run
# ===========================================================================

test_that("sw_correct_run rejects data without batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_run(df, markers = "CD3"),
               "batch")
})

test_that("sw_correct_run rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_run(df, markers = c("CD3", "CD99")),
               "not found")
})

test_that("sw_correct_run rejects non-data.frame", {
  expect_error(sw_correct_run("not_a_df", markers = "CD3"),
               "data.frame")
})

# ===========================================================================
# sw_correct_evaluate_quick
# ===========================================================================

test_that("sw_correct_evaluate_quick computes MAD metrics", {
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(50, mean = 1), rnorm(50, mean = 3)),
    CD4 = c(rnorm(50, mean = 2), rnorm(50, mean = 4)),
    batch = rep(c("B1", "B2"), each = 50)
  )

  # Simulate correction: reduce batch difference
  corrected <- tibble::tibble(
    CD3 = c(rnorm(50, mean = 2), rnorm(50, mean = 2.1)),
    CD4 = c(rnorm(50, mean = 3), rnorm(50, mean = 3.1)),
    batch = rep(c("B1", "B2"), each = 50)
  )

  result <- sw_correct_evaluate_quick(uncorrected, corrected,
                                   markers = c("CD3", "CD4"))

  expect_type(result, "list")
  expect_true(all(c("emd", "mad", "improved", "emd_reduction_pct") %in%
                    names(result)))
  expect_s3_class(result$mad, "tbl_df")
  expect_true("marker" %in% names(result$mad))
  expect_true("mad_before" %in% names(result$mad))
  expect_true("mad_after" %in% names(result$mad))

  # Correction should have improved things
  expect_true(result$improved)
})

test_that("sw_correct_evaluate_quick handles single batch", {
  df <- tibble::tibble(
    CD3 = rnorm(20),
    batch = "B1"
  )

  result <- sw_correct_evaluate_quick(df, df, markers = "CD3")
  expect_true(is.na(result$improved))
})

test_that("sw_correct_evaluate_quick rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_evaluate_quick(df, df, markers = "CD3"),
               "batch")
})

test_that("sw_correct_evaluate_quick rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_evaluate_quick(df, df, markers = "MISSING"),
               "not found")
})

# ===========================================================================
# sw_correct_normalize — modular workflow step 1
# ===========================================================================

test_that("sw_correct_normalize rejects non-data.frame", {
  expect_error(sw_correct_normalize("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_normalize rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_normalize(df, markers = "CD3"),
               "batch")
})

test_that("sw_correct_normalize rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_normalize(df, markers = character(0)),
               "non-empty")
})

test_that("sw_correct_normalize rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_normalize(df, markers = "MISSING"),
               "not found")
})

test_that("sw_correct_normalize validates norm_method", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_normalize(df, markers = "CD3", norm_method = "invalid"),
               "arg")
})

# ===========================================================================
# sw_correct_som — modular workflow step 2
# ===========================================================================

test_that("sw_correct_som rejects non-data.frame", {
  expect_error(sw_correct_som("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_som rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_som(df, markers = character(0)),
               "non-empty")
})

test_that("sw_correct_som rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_som(df, markers = "MISSING"),
               "not found")
})

# ===========================================================================
# sw_correct_apply — modular workflow step 3
# ===========================================================================

test_that("sw_correct_apply rejects non-data.frame", {
  expect_error(sw_correct_apply("bad", label = 1L, markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_apply rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_apply(df, label = rep(1L, 10), markers = "CD3"),
               "batch")
})

test_that("sw_correct_apply rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_apply(df, label = rep(1L, 10), markers = "MISSING"),
               "not found")
})

test_that("sw_correct_apply rejects non-numeric label", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_apply(df, label = letters[1:10], markers = "CD3"),
               "numeric")
})

test_that("sw_correct_apply rejects label length mismatch", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_apply(df, label = c(1L, 2L), markers = "CD3"),
               "must equal nrow")
})

# ===========================================================================
# sw_correct_detect_batch
# ===========================================================================

test_that("sw_correct_detect_batch rejects non-data.frame", {
  expect_error(sw_correct_detect_batch("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_detect_batch rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_correct_detect_batch(df, markers = "CD3"),
               "batch")
})

test_that("sw_correct_detect_batch rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_detect_batch(df, markers = "MISSING"),
               "not found")
})

# ===========================================================================
# sw_correct_emd
# ===========================================================================

test_that("sw_correct_emd rejects non-data.frame", {
  expect_error(sw_correct_emd("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_emd rejects missing label column", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_correct_emd(df, markers = "CD3"),
               "label")
})

test_that("sw_correct_emd rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  expect_error(sw_correct_emd(df, markers = "CD3"),
               "batch")
})

test_that("sw_correct_emd rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_emd(df, markers = "MISSING"),
               "not found")
})

test_that("sw_correct_emd rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_emd(df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_correct_evaluate_emd
# ===========================================================================

test_that("sw_correct_evaluate_emd rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_emd("bad", df, markers = "CD3"),
               "data.frame")
  expect_error(sw_correct_evaluate_emd(df, "bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_evaluate_emd rejects missing columns", {
  df_ok <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  df_no_batch <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  df_no_label <- tibble::tibble(CD3 = rnorm(10), batch = "B1")

  expect_error(sw_correct_evaluate_emd(df_no_batch, df_ok, markers = "CD3"),
               "batch")
  expect_error(sw_correct_evaluate_emd(df_ok, df_no_label, markers = "CD3"),
               "label")
})

test_that("sw_correct_evaluate_emd rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_emd(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_correct_evaluate_emd rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_emd(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_correct_evaluate_mad (cyCombine-backed)
# ===========================================================================

test_that("sw_correct_evaluate_mad rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_mad("bad", df, markers = "CD3"),
               "data.frame")
  expect_error(sw_correct_evaluate_mad(df, "bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_correct_evaluate_mad rejects missing columns", {
  df_ok <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  df_no_batch <- tibble::tibble(CD3 = rnorm(10), label = 1L)
  df_no_label <- tibble::tibble(CD3 = rnorm(10), batch = "B1")

  expect_error(sw_correct_evaluate_mad(df_no_batch, df_ok, markers = "CD3"),
               "batch")
  expect_error(sw_correct_evaluate_mad(df_ok, df_no_label, markers = "CD3"),
               "label")
})

test_that("sw_correct_evaluate_mad rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_mad(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_correct_evaluate_mad rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1", label = 1L)
  expect_error(sw_correct_evaluate_mad(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_plot_batch_densities
# ===========================================================================

test_that("sw_plot_batch_densities rejects non-data.frame", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_plot_batch_densities("bad", df, markers = "CD3"),
               "data.frame")
})

test_that("sw_plot_batch_densities rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_densities(df, df, markers = "MISSING"),
               "not found")
})

test_that("sw_plot_batch_densities rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_densities(df, df, markers = character(0)),
               "non-empty")
})

# ===========================================================================
# sw_plot_batch_dimred
# ===========================================================================

test_that("sw_plot_batch_dimred rejects non-data.frame", {
  expect_error(sw_plot_batch_dimred("bad", markers = "CD3"),
               "data.frame")
})

test_that("sw_plot_batch_dimred rejects missing markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = "MISSING"),
               "not found")
})

test_that("sw_plot_batch_dimred rejects missing batch column", {
  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_plot_batch_dimred(df, markers = "CD3"),
               "batch")
})

test_that("sw_plot_batch_dimred rejects empty markers", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = character(0)),
               "non-empty")
})

test_that("sw_plot_batch_dimred validates type argument", {
  df <- tibble::tibble(CD3 = rnorm(10), batch = "B1")
  expect_error(sw_plot_batch_dimred(df, markers = "CD3", type = "invalid"),
               "arg")
})

# ===========================================================================
# Integration tests — functional happy-path tests with cyCombine
# ===========================================================================

test_that("sw_correct_run runs end-to-end with cyCombine", {
  skip_if_not_installed("cyCombine")

  set.seed(42)
  n <- 100
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 1), rnorm(n, mean = 3)),
    CD4 = c(rnorm(n, mean = 2), rnorm(n, mean = 4)),
    batch = rep(c("B1", "B2"), each = n),
    sample = rep(c("S1", "S2"), each = n)
  )

  corrected <- sw_correct_run(uncorrected, markers = c("CD3", "CD4"),
                                xdim = 4, ydim = 4, rlen = 5, seed = 42,
                                design_audit = NULL)

  expect_identical(attr(corrected, "design_audit_check")$status, "skipped")
  expect_s3_class(corrected, "tbl_df")
  expect_true(all(c("CD3", "CD4") %in% names(corrected)))
  expect_equal(nrow(corrected), 2 * n)
  expect_true("batch" %in% names(corrected))
})

test_that("modular batch correction workflow (normalize -> create_som -> correct_data)", {
  skip_if_not_installed("cyCombine")

  set.seed(42)
  n <- 100
  df <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 1), rnorm(n, mean = 3)),
    CD4 = c(rnorm(n, mean = 2), rnorm(n, mean = 4)),
    batch = rep(c("B1", "B2"), each = n),
    sample = rep(c("S1", "S2"), each = n)
  )
  markers <- c("CD3", "CD4")

  # Step 1: Normalize
  normalized <- sw_correct_normalize(df, markers = markers, norm_method = "scale")
  expect_s3_class(normalized, "tbl_df")
  expect_equal(nrow(normalized), 2 * n)
  expect_true(all(markers %in% names(normalized)))

  # Step 2: Create SOM
  labels <- sw_correct_som(normalized, markers = markers,
                          xdim = 4, ydim = 4, rlen = 5, seed = 42)
  expect_true(is.numeric(labels) || is.integer(labels))
  expect_equal(length(labels), 2 * n)
  expect_true(all(labels >= 1))

  # Step 3: Correct data (using original, unnormalized df)
  corrected <- sw_correct_apply(df, label = labels, markers = markers)
  expect_s3_class(corrected, "tbl_df")
  expect_equal(nrow(corrected), 2 * n)
  expect_true(all(markers %in% names(corrected)))
})

test_that("sw_correct_prepare rejects duplicate sample names", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)

  expect_error(
    sw_correct_prepare(
      list(S1 = ff, S2 = ff),
      tibble::tibble(sample = c("S1", "S1"), batch = c("B1", "B2")),
      markers = "CD3"
    ),
    "duplicate sample name"
  )
})

test_that("sw_plot_batch_densities produces a ggplot object", {
  skip_if_not_installed("cyCombine")
  skip_if_not_installed("ggplot2")

  set.seed(42)
  n <- 50
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 1), rnorm(n, mean = 3)),
    batch = rep(c("B1", "B2"), each = n)
  )
  corrected <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 2), rnorm(n, mean = 2.1)),
    batch = rep(c("B1", "B2"), each = n)
  )

  p <- sw_plot_batch_densities(uncorrected, corrected, markers = "CD3")
  expect_s3_class(p, "gg")
})

# ===========================================================================
# Row id column required by cyCombine >= 0.3.0 evaluation functions
# ===========================================================================

test_that("sw_correct_prepare adds a unique row id column", {
  skip_if_not_installed("flowCore")

  mat <- matrix(c(1000, 2000, 3000, 4000), ncol = 2,
                dimnames = list(NULL, c("CD3", "CD4")))
  ff <- flowCore::flowFrame(mat)
  result <- sw_correct_prepare(
    list(S1 = ff, S2 = ff),
    tibble::tibble(sample = c("S1", "S2"), batch = c("B1", "B2")),
    markers = c("CD3", "CD4")
  )

  expect_identical(result$id, seq_len(nrow(result)))
})

test_that("sw_correct_evaluate_emd and _mad run on data without an id column", {
  skip_if_not_installed("cyCombine")
  skip_if_not_installed("emdist")

  set.seed(42)
  n <- 60
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 1), rnorm(n, mean = 3)),
    CD4 = c(rnorm(n, mean = 2), rnorm(n, mean = 4)),
    batch = rep(c("B1", "B2"), each = n),
    label = rep(c(1L, 2L), times = n)
  )
  corrected <- uncorrected
  corrected$CD3 <- corrected$CD3 - c(rep(0, n), rep(2, n))
  expect_false("id" %in% names(uncorrected))

  emd_eval <- suppressMessages(
    sw_correct_evaluate_emd(uncorrected, corrected, markers = c("CD3", "CD4"))
  )
  expect_true(all(c("reduction", "emd") %in% names(emd_eval)))

  mad_eval <- suppressMessages(
    sw_correct_evaluate_mad(uncorrected, corrected, markers = c("CD3", "CD4"))
  )
  expect_true(all(c("score", "mad") %in% names(mad_eval)))
})

test_that("no package function uses the magrittr pipe, which the namespace does not import", {
  # testthat attaches magrittr, which hides a missing %>% import at test time
  ns <- asNamespace("SpectraWeaveR")
  uses_magrittr_pipe <- vapply(
    ls(ns, all.names = TRUE),
    \(nm) {
      obj <- get(nm, envir = ns)
      is.function(obj) && "%>%" %in% all.names(body(obj))
    },
    logical(1)
  )
  expect_identical(names(uses_magrittr_pipe)[uses_magrittr_pipe], character(0))
})

# ===========================================================================
# Design-audit gate (A3)
# ===========================================================================

make_gate_audit <- function(batches) {
  samples <- data.frame(
    sample_id = paste0("s", seq_along(batches)),
    batch = batches,
    condition = rep(c("ctrl", "case"), length.out = length(batches)),
    reference_role = "none",
    reference_id = NA_character_
  )
  sw_design_audit(
    samples,
    minimum_batches = 1,
    minimum_batches_origin = "test fixture",
    require_shared_reference = FALSE,
    protected = "condition"
  )
}

gate_df <- tibble::tibble(
  CD3 = c(1, 2, 3, 4),
  batch = c("B1", "B1", "B2", "B2"),
  sample = c("S1", "S1", "S2", "S2")
)

test_that("sw_correct_run errors when design_audit is not given", {
  skip_if_not_installed("cyCombine")
  expect_error(
    sw_correct_run(gate_df, markers = "CD3"),
    "Required argument\\(s\\) with no default: 'design_audit'"
  )
})

test_that("sw_correct_run stops on an ineligible audit and names the row", {
  skip_if_not_installed("cyCombine")
  audit <- make_gate_audit(c("B1", "B1", "B1", "B1"))
  row <- audit$eligibility[audit$eligibility$name == "cycombine", ]
  expect_false(row$eligible)
  one_batch_df <- transform(gate_df, batch = "B1")
  expect_error(
    sw_correct_run(one_batch_df, markers = "CD3", design_audit = audit),
    "Design audit row 'cycombine' \\(sw_correct_run\\) is ineligible"
  )
})

test_that("the gate passes an eligible audit and records the row", {
  audit <- make_gate_audit(c("b1", "b1", "b2", "b2", "b3", "b3"))
  rec <- .sw_gate_design_audit(audit, method = "cycombine", batches = c("b1", "b2", "b3"))
  expect_identical(rec$status, "checked")
  expect_identical(rec$row$name, "cycombine")
  expect_true(rec$row$eligible)
})

test_that("the gate records NULL as an explicit skip", {
  rec <- .sw_gate_design_audit(NULL, method = "cycombine", batches = "B1")
  expect_identical(rec$status, "skipped")
  expect_null(rec$row)
})

test_that("the gate rejects a non-audit object and an unknown method", {
  expect_error(
    .sw_gate_design_audit(list(), method = "cycombine", batches = "B1"),
    "sw_design_audit object or NULL"
  )
  audit <- make_gate_audit(c("b1", "b1", "b2", "b2"))
  expect_error(
    .sw_gate_design_audit(audit, method = "nope", batches = c("b1", "b2")),
    "no eligibility row named 'nope'"
  )
})

test_that("the gate stops when the audit describes other batches than the data", {
  audit <- make_gate_audit(c("b1", "b1", "b2", "b2", "b3", "b3"))
  expect_error(
    .sw_gate_design_audit(audit, method = "cycombine", batches = c("x1", "x2")),
    "describes other batches than the data"
  )
})

test_that("the gate stops when covar is outside the audit's protected factors", {
  audit <- make_gate_audit(c("b1", "b1", "b2", "b2", "b3", "b3"))
  expect_error(
    .sw_gate_design_audit(
      audit,
      method = "cycombine",
      batches = c("b1", "b2", "b3"),
      covar = "sex"
    ),
    "not among the audit's protected factors"
  )
  rec <- .sw_gate_design_audit(
    audit,
    method = "cycombine",
    batches = c("b1", "b2", "b3"),
    covar = "condition"
  )
  expect_identical(rec$status, "checked")
})

test_that("sw_correct_run with an eligible audit records a checked gate", {
  skip_if_not_installed("cyCombine")
  audit <- make_gate_audit(c("B1", "B1", "B2", "B2", "B3", "B3"))
  set.seed(42)
  n <- 60
  uncorrected <- tibble::tibble(
    CD3 = c(rnorm(n, mean = 1), rnorm(n, mean = 3), rnorm(n, mean = 2)),
    CD4 = c(rnorm(n, mean = 2), rnorm(n, mean = 4), rnorm(n, mean = 3)),
    batch = rep(c("B1", "B2", "B3"), each = n),
    sample = rep(c("S1", "S2", "S3"), each = n)
  )
  corrected <- sw_correct_run(
    uncorrected,
    markers = c("CD3", "CD4"),
    xdim = 3,
    ydim = 3,
    rlen = 3,
    seed = 42,
    design_audit = audit
  )
  rec <- attr(corrected, "design_audit_check")
  expect_identical(rec$status, "checked")
  expect_true(rec$row$eligible)

  # an audit for other batches stops the run before correction
  expect_error(
    sw_correct_run(
      uncorrected,
      markers = c("CD3", "CD4"),
      design_audit = make_gate_audit(c("a", "a", "b", "b", "c", "c"))
    ),
    "describes other batches than the data"
  )
})
