# tests/testthat/test-gate.R
# Unit tests for R/gate.R — Automated gating (openCyto)

test_that("sw_gate_template creates valid CSV", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  result <- sw_gate_template(tmp_file)

  expect_true(file.exists(tmp_file))
  expect_equal(result, tmp_file)

  # Read and validate structure
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 3)

  # Required columns for openCyto gating template
  required_cols <- c("alias", "pop", "parent", "dims", "gating_method")
  expect_true(all(required_cols %in% names(gt)))

  # Check aliases
  expect_true("nonDebris" %in% gt$alias)
  expect_true("singlets" %in% gt$alias)
  expect_true("lymphocytes" %in% gt$alias)
})

test_that("sw_gate_template creates parent directory", {
  tmp_dir <- file.path(tempdir(), "nested", "gating_test")
  tmp_file <- file.path(tmp_dir, "template.csv")
  on.exit(unlink(file.path(tempdir(), "nested"), recursive = TRUE))

  sw_gate_template(tmp_file)
  expect_true(file.exists(tmp_file))
})

test_that("sw_gate_template rejects invalid template_type", {
  expect_error(sw_gate_template(tempfile(), template_type = "invalid"),
               "Unknown template_type")
})

test_that("sw_gate_template rejects non-character output_file", {
  expect_error(sw_gate_template(42),
               "single file path")
})

# =========================================================================
# New template types
# =========================================================================

test_that("sw_gate_template creates myeloid template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  sw_gate_template(tmp_file, template_type = "myeloid")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 4)
  expect_true("nonDebris" %in% gt$alias)
  expect_true("singlets" %in% gt$alias)
  expect_true("CD3neg" %in% gt$alias)
  expect_true("HLADRpos" %in% gt$alias)
  # Verify parent hierarchy: every parent must be root or an alias
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates nk template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  sw_gate_template(tmp_file, template_type = "nk")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 4)
  expect_true("CD56pos" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates treg template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  sw_gate_template(tmp_file, template_type = "treg")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 5)
  expect_true("CD3pos" %in% gt$alias)
  expect_true("CD4pos" %in% gt$alias)
  expect_true("Treg" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates full_pbmc template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  sw_gate_template(tmp_file, template_type = "full_pbmc")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 8)
  # Should have branching: CD3pos and CD3neg both under singlets
  expect_true("CD3pos" %in% gt$alias)
  expect_true("CD3neg" %in% gt$alias)
  expect_true("CD19pos" %in% gt$alias)
  expect_true("CD14pos" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_run rejects empty fcs_files", {
  expect_error(sw_gate_run(character(0), "template.csv"),
               "non-empty character vector")
})

test_that("sw_gate_run rejects missing template file", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  # Create a dummy FCS file path (won't actually be read due to template check)
  expect_error(sw_gate_run("dummy.fcs", "/nonexistent/template.csv"),
               "not found")
})

test_that("sw_gate_extract rejects non-GatingSet input", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_gate_extract("not_a_gs", "/singlets"),
               "GatingSet")
})

test_that("sw_gate_extract rejects non-character node", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_gate_extract("not_a_gs", 42),
               "single character string")
})

# =========================================================================
# Integration tests — actual gating execution
# =========================================================================

test_that("sw_gate_run runs end-to-end with synthetic FCS", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("openCyto")

  set.seed(42)
  n <- 500
  # Bimodal FSC-A: debris (low) + cells (high)
  fsc_a <- c(abs(rnorm(100, 20000, 5000)),
             abs(rnorm(400, 150000, 30000)))
  fsc_h <- fsc_a * 0.7 + rnorm(n, sd = 5000)
  ssc_a <- abs(rnorm(n, 50000, 15000))

  mat <- cbind("FSC-A" = fsc_a, "FSC-H" = fsc_h, "SSC-A" = ssc_a)
  storage.mode(mat) <- "double"
  ff <- flowCore::flowFrame(mat)

  # Write to temp FCS
  tmp_fcs <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(ff, tmp_fcs)
  on.exit(unlink(tmp_fcs))

  # Build template
  tmp_template <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_template), add = TRUE)
  sw_gate_template(tmp_template)

  gs <- sw_gate_run(tmp_fcs, tmp_template)
  expect_true(methods::is(gs, "GatingSet"))

  # Verify gate nodes exist
  nodes <- flowWorkspace::gs_get_pop_paths(gs)
  expect_true(length(nodes) >= 2)
})

test_that("sw_gate_extract returns named list of flowFrames", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("openCyto")

  set.seed(42)
  n <- 500
  fsc_a <- c(abs(rnorm(100, 20000, 5000)),
             abs(rnorm(400, 150000, 30000)))
  fsc_h <- fsc_a * 0.7 + rnorm(n, sd = 5000)
  ssc_a <- abs(rnorm(n, 50000, 15000))

  mat <- cbind("FSC-A" = fsc_a, "FSC-H" = fsc_h, "SSC-A" = ssc_a)
  storage.mode(mat) <- "double"
  ff <- flowCore::flowFrame(mat)

  tmp_fcs <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(ff, tmp_fcs)
  on.exit(unlink(tmp_fcs))

  tmp_template <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_template), add = TRUE)
  sw_gate_template(tmp_template)

  gs <- sw_gate_run(tmp_fcs, tmp_template)
  nodes <- flowWorkspace::gs_get_pop_paths(gs)

  # Extract from first non-root node
  leaf <- nodes[length(nodes)]
  ff_list <- sw_gate_extract(gs, leaf)
  expect_true(is.list(ff_list))
  expect_true(length(ff_list) >= 1)
  expect_true(methods::is(ff_list[[1]], "flowFrame"))
})
