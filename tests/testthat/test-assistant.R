# tests/testthat/test-assistant.R
# Unit tests for R/assistant.R — LLM-powered pipeline builder assistant

# ===========================================================================
# Tool function tests (no LLM / ellmer required)
# ===========================================================================

# ---- .tool_list_fcs_files ----

test_that(".tool_list_fcs_files returns error for non-existent directory", {
  result <- SpectraWeaveR:::.tool_list_fcs_files("/nonexistent/dir")
  expect_match(result, "does not exist")
})

test_that(".tool_list_fcs_files finds FCS files", {
  tmp <- tempdir()
  fcs_dir <- file.path(tmp, "test_fcs_dir")
  dir.create(fcs_dir, showWarnings = FALSE)
  on.exit(unlink(fcs_dir, recursive = TRUE), add = TRUE)

  file.create(file.path(fcs_dir, "sample1.fcs"))
  file.create(file.path(fcs_dir, "sample2.fcs"))
  file.create(file.path(fcs_dir, "not_fcs.csv"))

  result <- SpectraWeaveR:::.tool_list_fcs_files(fcs_dir)
  expect_match(result, "Found 2 FCS")
  expect_match(result, "sample1.fcs")
  expect_match(result, "sample2.fcs")
  expect_no_match(result, "not_fcs.csv")
})

test_that(".tool_list_fcs_files reports when no FCS files found", {
  tmp <- tempdir()
  empty_dir <- file.path(tmp, "test_empty_fcs")
  dir.create(empty_dir, showWarnings = FALSE)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  file.create(file.path(empty_dir, "data.csv"))

  result <- SpectraWeaveR:::.tool_list_fcs_files(empty_dir)
  expect_match(result, "No FCS files")
})

test_that(".tool_list_fcs_files truncates long lists", {
  tmp <- tempdir()
  many_dir <- file.path(tmp, "test_many_fcs")
  dir.create(many_dir, showWarnings = FALSE)
  on.exit(unlink(many_dir, recursive = TRUE), add = TRUE)

  for (i in seq_len(25)) {
    file.create(file.path(many_dir, paste0("sample_", i, ".fcs")))
  }

  result <- SpectraWeaveR:::.tool_list_fcs_files(many_dir)
  expect_match(result, "Found 25 FCS")
  expect_match(result, "and 5 more")
})

# ---- .tool_list_directory ----

test_that(".tool_list_directory returns error for non-existent path", {
  result <- SpectraWeaveR:::.tool_list_directory("/nonexistent/dir")
  expect_match(result, "does not exist")
})

test_that(".tool_list_directory lists files and directories", {
  tmp <- tempdir()
  test_dir <- file.path(tmp, "test_list_dir")
  dir.create(test_dir, showWarnings = FALSE)
  dir.create(file.path(test_dir, "subdir"), showWarnings = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  file.create(file.path(test_dir, "file.txt"))

  result <- SpectraWeaveR:::.tool_list_directory(test_dir)
  expect_match(result, "subdir")
  expect_match(result, "file.txt")
})

# ---- .tool_read_csv_columns ----

test_that(".tool_read_csv_columns returns error for missing file", {
  result <- SpectraWeaveR:::.tool_read_csv_columns("/nonexistent/file.csv")
  expect_match(result, "does not exist")
})

test_that(".tool_read_csv_columns reads CSV correctly", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_meta.csv")
  on.exit(unlink(csv_path), add = TRUE)

  df <- data.frame(
    file = c("s1.fcs", "s2.fcs"),
    sample = c("S1", "S2"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )
  write.csv(df, csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_read_csv_columns(csv_path, n_rows = 2L)
  expect_match(result, "Columns \\(3\\)")
  expect_match(result, "file")
  expect_match(result, "sample")
  expect_match(result, "batch")
})

# ---- .tool_read_xlsx_columns ----

test_that(".tool_read_xlsx_columns returns error for missing file", {
  result <- SpectraWeaveR:::.tool_read_xlsx_columns("/nonexistent/file.xlsx")
  expect_match(result, "does not exist")
})

test_that(".tool_read_xlsx_columns reports missing readxl", {
  skip_if(requireNamespace("readxl", quietly = TRUE),
          "readxl is installed, cannot test missing-package path")

  tmp <- tempdir()
  xlsx_path <- file.path(tmp, "test.xlsx")
  file.create(xlsx_path)
  on.exit(unlink(xlsx_path), add = TRUE)

  result <- SpectraWeaveR:::.tool_read_xlsx_columns(xlsx_path)
  expect_match(result, "readxl")
})

# ---- .tool_validate_sample_meta ----

test_that(".tool_validate_sample_meta returns error for missing metadata", {
  result <- SpectraWeaveR:::.tool_validate_sample_meta(
    "/nonexistent/meta.csv", "/tmp"
  )
  expect_match(result, "not found")
})

test_that(".tool_validate_sample_meta returns error for missing fcs_dir", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_validate_meta.csv")
  on.exit(unlink(csv_path), add = TRUE)

  write.csv(data.frame(file = "a.fcs", sample = "A", batch = "B1"),
            csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_validate_sample_meta(
    csv_path, "/nonexistent/fcs_dir"
  )
  expect_match(result, "not found")
})

test_that(".tool_validate_sample_meta identifies missing columns", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_validate_meta2.csv")
  fcs_dir <- file.path(tmp, "test_validate_fcs")
  dir.create(fcs_dir, showWarnings = FALSE)
  on.exit({
    unlink(csv_path)
    unlink(fcs_dir, recursive = TRUE)
  }, add = TRUE)

  # CSV with non-standard columns
  write.csv(data.frame(filename = "a.fcs", patient_id = "P1", run = "R1"),
            csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_validate_sample_meta(csv_path, fcs_dir)
  expect_match(result, "Missing expected columns")
})

test_that(".tool_validate_sample_meta validates successfully", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_validate_ok.csv")
  fcs_dir <- file.path(tmp, "test_validate_ok_fcs")
  dir.create(fcs_dir, showWarnings = FALSE)
  on.exit({
    unlink(csv_path)
    unlink(fcs_dir, recursive = TRUE)
  }, add = TRUE)

  # Create matching files
  file.create(file.path(fcs_dir, "s1.fcs"))
  file.create(file.path(fcs_dir, "s2.fcs"))
  write.csv(data.frame(
    file = c("s1.fcs", "s2.fcs"),
    sample = c("S1", "S2"),
    batch = c("B1", "B2")
  ), csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_validate_sample_meta(csv_path, fcs_dir)
  expect_match(result, "No issues found")
  expect_match(result, "Batch distribution")
})

# ---- .tool_check_batch_balance ----

test_that(".tool_check_batch_balance returns error for missing file", {
  result <- SpectraWeaveR:::.tool_check_batch_balance(
    "/nonexistent/meta.csv", "batch"
  )
  expect_match(result, "not found")
})

test_that(".tool_check_batch_balance returns error for missing column", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_balance.csv")
  on.exit(unlink(csv_path), add = TRUE)

  write.csv(data.frame(x = 1:3), csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_check_batch_balance(csv_path, "batch")
  expect_match(result, "not found")
})

test_that(".tool_check_batch_balance summarises batches correctly", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_balance2.csv")
  on.exit(unlink(csv_path), add = TRUE)

  write.csv(data.frame(
    sample = paste0("S", 1:6),
    batch = c("B1", "B1", "B1", "B2", "B2", "B2")
  ), csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_check_batch_balance(csv_path, "batch")
  expect_match(result, "Number of batches: 2")
  expect_match(result, "Total samples: 6")
  expect_match(result, "reasonably balanced")
})

test_that(".tool_check_batch_balance warns about imbalance", {
  tmp <- tempdir()
  csv_path <- file.path(tmp, "test_imbalance.csv")
  on.exit(unlink(csv_path), add = TRUE)

  write.csv(data.frame(
    sample = paste0("S", 1:7),
    batch = c("B1", "B1", "B1", "B1", "B1", "B2", "B2")
  ), csv_path, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_check_batch_balance(csv_path, "batch")
  expect_match(result, "WARNING")
  expect_match(result, "unbalanced")
})

# ---- .tool_read_fcs_header ----

test_that(".tool_read_fcs_header returns error for missing file", {
  result <- SpectraWeaveR:::.tool_read_fcs_header("/nonexistent/file.fcs")
  expect_match(result, "does not exist")
})

# ---- .tool_detect_channels ----

test_that(".tool_detect_channels returns error for missing file", {
  result <- SpectraWeaveR:::.tool_detect_channels("/nonexistent/file.fcs")
  expect_match(result, "does not exist")
})

# ---- .tool_validate_markers ----

test_that(".tool_validate_markers returns error for missing file", {
  result <- SpectraWeaveR:::.tool_validate_markers(
    c("CD3", "CD4"), "/nonexistent/file.fcs"
  )
  expect_match(result, "not found")
})

# ===========================================================================
# System prompt loading
# ===========================================================================

test_that(".load_prompt loads pipeline_builder prompt", {
  prompt <- SpectraWeaveR:::.load_prompt("pipeline_builder")
  expect_type(prompt, "character")
  expect_true(nzchar(prompt))
  expect_match(prompt, "SpectraWeaveR")
})

test_that(".load_prompt loads batch_correction_builder prompt", {
  prompt <- SpectraWeaveR:::.load_prompt("batch_correction_builder")
  expect_type(prompt, "character")
  expect_match(prompt, "batch")
})

test_that(".load_prompt loads gating_builder prompt", {
  prompt <- SpectraWeaveR:::.load_prompt("gating_builder")
  expect_type(prompt, "character")
  expect_match(prompt, "gating")
})

test_that(".load_prompt loads unmixing_builder prompt", {
  prompt <- SpectraWeaveR:::.load_prompt("unmixing_builder")
  expect_type(prompt, "character")
  expect_match(prompt, "unmixing")
})

test_that(".load_prompt fails for non-existent prompt", {
  expect_error(
    SpectraWeaveR:::.load_prompt("nonexistent_prompt"),
    "not found"
  )
})

# ===========================================================================
# Pipeline code generation
# ===========================================================================

test_that("sw_generate_pipeline_code requires mandatory fields", {
  expect_error(
    sw_generate_pipeline_code(list(fcs_dir = "/data")),
    "Missing required config"
  )
})

test_that("sw_generate_pipeline_code generates run_pipeline code", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3", "CD4", "CD8"),
    lineage_markers = c("CD3", "CD4"),
    batch_col = "batch",
    sample_col = "sample",
    file_col = "file",
    cofactor = 6000,
    n_metaclusters = 20,
    output_dir = "results",
    seed = 42
  )

  code <- sw_generate_pipeline_code(config, style = "run_pipeline",
                                     output = "string")
  expect_type(code, "character")
  expect_match(code, "library\\(SpectraWeaveR\\)")
  expect_match(code, "run_pipeline")
  expect_match(code, "/data/fcs")
  expect_match(code, "CD3")
  expect_match(code, "cofactor")
})

test_that("sw_generate_pipeline_code generates composable code", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3", "CD4"),
    lineage_markers = c("CD3"),
    cofactor = 6000
  )

  code <- sw_generate_pipeline_code(config, style = "composable",
                                     output = "string")
  expect_type(code, "character")
  expect_match(code, "sw_pipeline")
  expect_match(code, "sw_remove_margins")
  expect_match(code, "sw_signal_qc")
  expect_match(code, "sw_cluster")
})

test_that("sw_generate_pipeline_code handles column renaming", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3"),
    lineage_markers = c("CD3"),
    batch_col = "BatchID",
    sample_col = "PatientID",
    file_col = "filename"
  )

  code <- sw_generate_pipeline_code(config, output = "string")
  expect_match(code, 'BatchID.*batch')
  expect_match(code, 'PatientID.*sample')
  expect_match(code, 'filename.*file')
})

test_that("sw_generate_pipeline_code handles XLSX metadata", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.xlsx",
    markers = c("CD3"),
    lineage_markers = c("CD3")
  )

  code <- sw_generate_pipeline_code(config, output = "string")
  expect_match(code, "readxl")
  expect_match(code, "read_excel")
})

test_that("sw_generate_pipeline_code includes unmixing when specified", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3"),
    lineage_markers = c("CD3"),
    unmix_from_raw = TRUE,
    control_dir = "/data/controls",
    unstained_fcs = "/data/unstained.fcs",
    cytometer = "aurora"
  )

  code <- sw_generate_pipeline_code(config, output = "string")
  expect_match(code, "unmix_from_raw.*TRUE")
  expect_match(code, "control_dir")
})

test_that("sw_generate_pipeline_code includes gating template", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3"),
    lineage_markers = c("CD3"),
    gating_template = "/data/gates.csv"
  )

  code <- sw_generate_pipeline_code(config, output = "string")
  expect_match(code, "gating_template")
  expect_match(code, "gates.csv")
})

test_that("sw_generate_pipeline_code writes to file", {
  config <- list(
    fcs_dir = "/data/fcs",
    sample_meta_path = "/data/metadata.csv",
    markers = c("CD3"),
    lineage_markers = c("CD3")
  )

  tmp_file <- tempfile(fileext = ".R")
  on.exit(unlink(tmp_file), add = TRUE)

  sw_generate_pipeline_code(config, output = "file", path = tmp_file)
  expect_true(file.exists(tmp_file))
  contents <- paste(readLines(tmp_file), collapse = "\n")
  expect_match(contents, "library\\(SpectraWeaveR\\)")
})

test_that("sw_generate_pipeline_code errors when file output lacks path", {
  config <- list(
    fcs_dir = "/data",
    sample_meta_path = "/data/m.csv",
    markers = "CD3",
    lineage_markers = "CD3"
  )

  expect_error(
    sw_generate_pipeline_code(config, output = "file"),
    "path.*required"
  )
})

# ===========================================================================
# Dependency checks
# ===========================================================================

test_that(".check_ellmer provides clear install instructions", {
  # We can only test the message format if ellmer is not installed
  skip_if(requireNamespace("ellmer", quietly = TRUE),
          "ellmer is installed, cannot test missing-package path")

  expect_error(SpectraWeaveR:::.check_ellmer(), "ellmer")
})

test_that("sw_assistant requires ellmer", {
  skip_if(requireNamespace("ellmer", quietly = TRUE),
          "ellmer is installed, cannot test missing-package path")

  expect_error(sw_assistant(), "ellmer")
})

test_that("sw_quick_pipeline validates description argument", {
  skip_if(!requireNamespace("ellmer", quietly = TRUE),
          "ellmer required")

  expect_error(sw_quick_pipeline(""), "non-empty")
  expect_error(sw_quick_pipeline(42), "non-empty")
})

# ===========================================================================
# Integration tests (require ellmer + API key)
# ===========================================================================

test_that("sw_assistant returns chat object in non-interactive mode", {
  skip_if_not_installed("ellmer")
  skip_if(Sys.getenv("OPENAI_API_KEY") == "", "No OpenAI API key")

  chat <- sw_assistant(interactive = FALSE, privacy_notice = FALSE)
  expect_true(inherits(chat, "Chat"))
})

test_that("sw_assistant_configure returns chat object", {
  skip_if_not_installed("ellmer")
  skip_if(Sys.getenv("OPENAI_API_KEY") == "", "No OpenAI API key")

  chat <- sw_assistant_configure()
  expect_true(inherits(chat, "Chat"))
})

test_that("sw_pipeline_config_type returns ellmer type object", {
  skip_if_not_installed("ellmer")

  config_type <- sw_pipeline_config_type()
  expect_true(inherits(config_type, "TypeObject"))
})

# ===========================================================================
# MCP server and btw integration tests
# ===========================================================================

test_that("sw_mcp_server requires mcptools", {
  skip_if(requireNamespace("mcptools", quietly = TRUE),
          "mcptools is installed; cannot test missing-package error")

  expect_error(sw_mcp_server(), "mcptools")
})

test_that("sw_mcp_server has include_btw parameter", {
  params <- names(formals(sw_mcp_server))
  expect_true("include_btw" %in% params)
})

test_that(".sw_tool_list returns list of tools", {
  skip_if_not_installed("ellmer")

  tools <- SpectraWeaveR:::.sw_tool_list()
  expect_true(is.list(tools))
  expect_true(length(tools) == 9)  # 9 SpectraWeaveR custom tools
})

test_that("sw_assistant accepts include_btw parameter", {
  params <- names(formals(sw_assistant))
  expect_true("include_btw" %in% params)
})

test_that("sw_assistant_configure accepts include_btw parameter", {
  params <- names(formals(sw_assistant_configure))
  expect_true("include_btw" %in% params)
})

test_that("sw_assistant works with include_btw = FALSE", {
  skip_if_not_installed("ellmer")
  skip_if(Sys.getenv("OPENAI_API_KEY") == "", "No OpenAI API key")

  chat <- sw_assistant(interactive = FALSE, privacy_notice = FALSE,
                       include_btw = FALSE)
  expect_true(inherits(chat, "Chat"))
})

# ===========================================================================
# Privacy safeguard tests
# ===========================================================================

# ---- .detect_sensitive_columns ----

test_that(".detect_sensitive_columns detects PII column names", {
  sensitive <- SpectraWeaveR:::.detect_sensitive_columns(
    c("patient_id", "name", "dob", "mrn", "ssn", "diagnosis",
      "batch", "file", "sample", "CD3", "CD4")
  )
  expect_true("patient_id" %in% sensitive)
  expect_true("name" %in% sensitive)
  expect_true("dob" %in% sensitive)
  expect_true("mrn" %in% sensitive)
  expect_true("ssn" %in% sensitive)
  expect_true("diagnosis" %in% sensitive)
  expect_false("batch" %in% sensitive)
  expect_false("file" %in% sensitive)
  expect_false("sample" %in% sensitive)
  expect_false("CD3" %in% sensitive)
})

test_that(".detect_sensitive_columns is case-insensitive", {
  sensitive <- SpectraWeaveR:::.detect_sensitive_columns(
    c("PatientID", "PATIENT_ID", "Patient_Id", "Name", "DOB")
  )
  expect_equal(length(sensitive), 5)
})

test_that(".detect_sensitive_columns returns empty for safe columns", {
  sensitive <- SpectraWeaveR:::.detect_sensitive_columns(
    c("batch", "file", "sample", "condition", "CD3", "marker_count")
  )
  expect_equal(length(sensitive), 0)
})

test_that(".detect_sensitive_columns handles empty input", {
  expect_equal(length(SpectraWeaveR:::.detect_sensitive_columns(character(0))), 0)
})

# ---- .redact_preview ----

test_that(".redact_preview redacts sensitive columns", {
  df <- data.frame(
    patient_id = c("P001", "P002"),
    batch = c("B1", "B2"),
    diagnosis = c("Leukemia", "Lymphoma"),
    stringsAsFactors = FALSE
  )
  redacted <- SpectraWeaveR:::.redact_preview(df, c("patient_id", "diagnosis"))
  expect_equal(redacted$patient_id, c("[REDACTED]", "[REDACTED]"))
  expect_equal(redacted$diagnosis, c("[REDACTED]", "[REDACTED]"))
  expect_equal(redacted$batch, c("B1", "B2"))  # unchanged
})

test_that(".redact_preview with no sensitive columns returns unchanged df", {
  df <- data.frame(batch = c("B1", "B2"), file = c("a.fcs", "b.fcs"))
  result <- SpectraWeaveR:::.redact_preview(df, character(0))
  expect_equal(result, df)
})

# ---- .tool_read_csv_columns with privacy modes ----

test_that(".tool_read_csv_columns columns_only suppresses rows", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(patient_id = "P001", batch = "B1", CD3 = 1.5),
            tmp, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_read_csv_columns(tmp, columns_only = TRUE)
  expect_match(result, "patient_id")
  expect_match(result, "Row preview suppressed")
  expect_no_match(result, "P001")
})

test_that(".tool_read_csv_columns standard mode redacts PII", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(patient_id = "P001", batch = "B1", CD3 = 1.5),
            tmp, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_read_csv_columns(tmp, sensitivity = "standard")
  expect_match(result, "REDACTED")
  expect_match(result, "patient_id")
  expect_no_match(result, "P001")
  expect_match(result, "B1")  # non-sensitive values still shown
})

test_that(".tool_read_csv_columns strict mode forces columns_only", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(batch = "B1", CD3 = 1.5), tmp, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_read_csv_columns(tmp, sensitivity = "strict")
  expect_match(result, "strict privacy mode")
  expect_match(result, "batch")
  expect_no_match(result, "B1")
})

test_that(".tool_read_csv_columns standard mode passes through safe data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(batch = "B1", file = "s1.fcs", CD3 = 1.5),
            tmp, row.names = FALSE)

  result <- SpectraWeaveR:::.tool_read_csv_columns(tmp, sensitivity = "standard")
  expect_no_match(result, "REDACTED")
  expect_match(result, "B1")  # safe preview shown
})

# ---- sensitivity parameter acceptance ----

test_that("all assistant functions accept sensitivity parameter", {
  for (fn_name in c("sw_assistant", "sw_assistant_configure",
                    "sw_assistant_batch_correction", "sw_assistant_gating",
                    "sw_assistant_unmixing", "sw_mcp_server")) {
    fn <- get(fn_name, envir = asNamespace("SpectraWeaveR"))
    expect_true("sensitivity" %in% names(formals(fn)),
                info = paste(fn_name, "missing sensitivity param"))
  }
})

test_that("sensitivity rejects invalid values", {
  skip_if_not_installed("ellmer")
  skip_if(Sys.getenv("OPENAI_API_KEY") == "", "No OpenAI API key")

  expect_error(
    sw_assistant(interactive = FALSE, privacy_notice = FALSE,
                 sensitivity = "invalid"),
    "arg"
  )
})

# ---- system prompts contain Data Privacy ----

test_that("all system prompts contain Data Privacy section", {
  prompt_dir <- system.file("prompts", package = "SpectraWeaveR")
  skip_if(!nzchar(prompt_dir), "prompts directory not found")

  for (name in c("pipeline_builder", "batch_correction_builder",
                  "gating_builder", "unmixing_builder")) {
    path <- file.path(prompt_dir, paste0(name, ".md"))
    skip_if(!file.exists(path), paste("prompt file missing:", name))
    content <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(content, "Data Privacy",
                 info = paste(name, "missing Data Privacy section"))
  }
})
