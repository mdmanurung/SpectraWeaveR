#' @title LLM-Powered Pipeline Builder Assistant
#'
#' @description
#' Interactive assistant powered by large language models (via the ellmer
#' package) that guides users through configuring a SpectraWeaveR analysis
#' pipeline. The assistant asks structured questions, inspects user files
#' and directories via registered tools, and generates runnable R code.
#'
#' @name assistant
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------

.check_ellmer <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop(
      "Package 'ellmer' is required for the LLM assistant.\n",
      "Install it with: install.packages('ellmer')",
      call. = FALSE
    )
  }
}

#' Null-coalescing operator (backport for R < 4.4.0)
#' @noRd
.null_default <- function(x, default) {
  if (is.null(x)) default else x
}

# ---------------------------------------------------------------------------
# Privacy helpers
# ---------------------------------------------------------------------------

#' Regex patterns matching common PII / sensitive column names
#' @noRd
.PII_PATTERNS <- c(
  "^patient[_\\.]?id$", "^pat[_\\.]?id$", "^subject[_\\.]?id$",
  "^subj[_\\.]?id$", "^participant[_\\.]?id$",
  "^name$", "^first[_\\.]?name$", "^last[_\\.]?name$",
  "^surname$", "^full[_\\.]?name$", "^patient[_\\.]?name$",
  "^dob$", "^date[_\\.]?of[_\\.]?birth$", "^birth[_\\.]?date$",
  "^mrn$", "^medical[_\\.]?record", "^record[_\\.]?number$",
  "^ssn$", "^social[_\\.]?security",
  "^diagnosis$", "^dx$", "^icd", "^icd[_\\.]?code$",
  "^address$", "^zip$", "^zip[_\\.]?code$", "^postal",
  "^phone$", "^email$", "^telephone$",
  "^ethnicity$", "^race$"
)

#' Detect columns with names matching common PII patterns
#'
#' @param col_names Character vector of column names.
#' @return Character vector of column names that match PII patterns
#'   (may be empty).
#' @noRd
.detect_sensitive_columns <- function(col_names) {
  if (length(col_names) == 0) return(character(0))
  pattern <- paste(.PII_PATTERNS, collapse = "|")
  col_names[grepl(pattern, col_names, ignore.case = TRUE)]
}

#' Redact values in sensitive columns of a data.frame
#'
#' Replaces all values in columns matching \code{sensitive_cols} with
#' \code{"[REDACTED]"}.
#'
#' @param df A data.frame.
#' @param sensitive_cols Character vector of column names to redact.
#' @return The modified data.frame with sensitive column values replaced.
#' @noRd
.redact_preview <- function(df, sensitive_cols) {
  for (col in intersect(sensitive_cols, names(df))) {
    df[[col]] <- rep("[REDACTED]", nrow(df))
  }
  df
}

# ---------------------------------------------------------------------------
# System prompt helpers
# ---------------------------------------------------------------------------

#' Load a system prompt from inst/prompts
#'
#' @param name Base name of the prompt file (without .md extension).
#' @return Character string with the prompt content.
#' @noRd
.load_prompt <- function(name) {
  path <- system.file("prompts", paste0(name, ".md"),
                       package = "SpectraWeaveR")
  if (!nzchar(path) || !file.exists(path)) {
    stop("Prompt file '", name, ".md' not found in inst/prompts/.",
         call. = FALSE)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

# ---------------------------------------------------------------------------
# Tool implementations
# ---------------------------------------------------------------------------

#' List FCS files in a directory
#' @param dir_path Absolute path to a directory.
#' @return Character summary of FCS files found.
#' @noRd
.tool_list_fcs_files <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    return(paste0("Error: Directory does not exist: ", dir_path))
  }

  fcs <- list.files(dir_path, pattern = "\\.fcs$", ignore.case = TRUE)
  if (length(fcs) == 0) {
    return("No FCS files found in this directory.")
  }

  header <- paste0("Found ", length(fcs), " FCS file(s):\n")
  listed <- paste0("- ", utils::head(fcs, 20), collapse = "\n")
  more <- if (length(fcs) > 20) {
    paste0("\n... and ", length(fcs) - 20, " more")
  } else {
    ""
  }
  paste0(header, listed, more)
}

#' List files and subdirectories
#' @param dir_path Absolute path to a directory.
#' @return Character summary of directory contents.
#' @noRd
.tool_list_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    return(paste0("Error: Directory does not exist: ", dir_path))
  }

  entries <- list.files(dir_path, all.files = FALSE)
  if (length(entries) == 0) {
    return("Directory is empty.")
  }

  dirs <- entries[dir.exists(file.path(dir_path, entries))]
  files <- setdiff(entries, dirs)

  parts <- character(0)
  if (length(dirs) > 0) {
    parts <- c(parts,
               paste0("Directories (", length(dirs), "):"),
               paste0("  ", utils::head(dirs, 20)))
    if (length(dirs) > 20) {
      parts <- c(parts, paste0("  ... and ", length(dirs) - 20, " more"))
    }
  }
  if (length(files) > 0) {
    parts <- c(parts,
               paste0("Files (", length(files), "):"),
               paste0("  ", utils::head(files, 30)))
    if (length(files) > 30) {
      parts <- c(parts, paste0("  ... and ", length(files) - 30, " more"))
    }
  }
  paste(parts, collapse = "\n")
}

#' Read FCS file header to extract channel/marker information
#' @param file_path Path to a single FCS file.
#' @return Character summary of channels and markers.
#' @noRd
.tool_read_fcs_header <- function(file_path) {
  if (!file.exists(file_path)) {
    return(paste0("Error: File does not exist: ", file_path))
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    return("Error: Package 'flowCore' is not installed. Cannot read FCS files.")
  }

  tryCatch({
    header <- flowCore::read.FCSheader(file_path)
    kw <- header[[1]]

    # Extract parameter info
    par_count <- as.integer(kw[["$PAR"]])
    channels <- character(par_count)
    markers <- character(par_count)
    for (i in seq_len(par_count)) {
      channels[i] <- kw[[paste0("$P", i, "N")]]
      if (is.null(channels[i])) channels[i] <- ""
      markers[i] <- kw[[paste0("$P", i, "S")]]
      if (is.null(markers[i])) markers[i] <- ""
    }

    # Build formatted output
    lines <- paste0(
      "  ", channels,
      ifelse(nzchar(markers), paste0(" -> ", markers), "")
    )

    paste0(
      "FCS file: ", basename(file_path), "\n",
      "Parameters: ", par_count, "\n",
      "Channels (name -> marker):\n",
      paste(lines, collapse = "\n")
    )
  }, error = function(e) {
    paste0("Error reading FCS header: ", conditionMessage(e))
  })
}

#' Read CSV columns and preview rows
#' @param file_path Path to a CSV file.
#' @param n_rows Number of preview rows (default 5).
#' @param columns_only If TRUE, return column names and types only (no rows).
#' @param sensitivity Privacy level: "standard" (redact PII) or "strict"
#'   (columns only, no row values).
#' @return Character summary of columns and optionally preview data.
#' @noRd
.tool_read_csv_columns <- function(file_path, n_rows = 5L,
                                   columns_only = FALSE,
                                   sensitivity = "standard") {
  if (!file.exists(file_path)) {
    return(paste0("Error: File does not exist: ", file_path))
  }

  if (sensitivity == "strict") columns_only <- TRUE

  tryCatch({
    df <- utils::read.csv(file_path, nrows = n_rows, stringsAsFactors = FALSE)
    cols <- names(df)

    col_info <- paste0("Columns (", length(cols), "): ",
                       paste(cols, collapse = ", "))

    if (isTRUE(columns_only)) {
      col_types <- vapply(df, function(x) class(x)[1], character(1))
      type_info <- paste0("  ", cols, ": ", col_types, collapse = "\n")
      return(paste0(col_info, "\n\nColumn types:\n", type_info,
                    "\n\n(Row preview suppressed",
                    if (sensitivity == "strict") " — strict privacy mode" else "",
                    ")"))
    }

    sensitive <- .detect_sensitive_columns(cols)
    note <- ""
    if (length(sensitive) > 0) {
      df <- .redact_preview(df, sensitive)
      note <- paste0(
        "\n\nNOTE: Values redacted in column(s) ",
        paste(sensitive, collapse = ", "),
        " (potential PII detected)."
      )
    }

    preview <- utils::capture.output(print(utils::head(df, n_rows)))
    paste0(col_info, "\n\nPreview:\n", paste(preview, collapse = "\n"), note)
  }, error = function(e) {
    paste0("Error reading CSV: ", conditionMessage(e))
  })
}

#' Read XLSX columns and preview rows
#' @param file_path Path to an XLSX file.
#' @param sheet Sheet name or index (default 1).
#' @param n_rows Number of preview rows (default 5).
#' @param columns_only If TRUE, return column names and types only (no rows).
#' @param sensitivity Privacy level: "standard" (redact PII) or "strict"
#'   (columns only, no row values).
#' @return Character summary of columns and optionally preview data.
#' @noRd
.tool_read_xlsx_columns <- function(file_path, sheet = 1L, n_rows = 5L,
                                    columns_only = FALSE,
                                    sensitivity = "standard") {
  if (!file.exists(file_path)) {
    return(paste0("Error: File does not exist: ", file_path))
  }

  if (!requireNamespace("readxl", quietly = TRUE)) {
    return("Error: Package 'readxl' is not installed. Install with: install.packages('readxl')")
  }

  if (sensitivity == "strict") columns_only <- TRUE

  tryCatch({
    df <- readxl::read_excel(file_path, sheet = sheet, n_max = n_rows)
    cols <- names(df)

    col_info <- paste0("Columns (", length(cols), "): ",
                       paste(cols, collapse = ", "))

    if (isTRUE(columns_only)) {
      col_types <- vapply(df, function(x) class(x)[1], character(1))
      type_info <- paste0("  ", cols, ": ", col_types, collapse = "\n")
      return(paste0(col_info, "\n\nColumn types:\n", type_info,
                    "\n\n(Row preview suppressed",
                    if (sensitivity == "strict") " — strict privacy mode" else "",
                    ")"))
    }

    sensitive <- .detect_sensitive_columns(cols)
    note <- ""
    if (length(sensitive) > 0) {
      df <- .redact_preview(df, sensitive)
      note <- paste0(
        "\n\nNOTE: Values redacted in column(s) ",
        paste(sensitive, collapse = ", "),
        " (potential PII detected)."
      )
    }

    preview <- utils::capture.output(print(utils::head(df, n_rows)))
    paste0(col_info, "\n\nPreview:\n", paste(preview, collapse = "\n"), note)
  }, error = function(e) {
    paste0("Error reading XLSX: ", conditionMessage(e))
  })
}

#' Validate sample metadata against FCS files
#' @param meta_path Path to a CSV or XLSX metadata file.
#' @param fcs_dir Path to the FCS files directory.
#' @return Character summary of validation results.
#' @noRd
.tool_validate_sample_meta <- function(meta_path, fcs_dir) {
  if (!file.exists(meta_path)) {
    return(paste0("Error: Metadata file not found: ", meta_path))
  }
  if (!dir.exists(fcs_dir)) {
    return(paste0("Error: FCS directory not found: ", fcs_dir))
  }

  tryCatch({
    # Read metadata
    ext <- tolower(tools::file_ext(meta_path))
    if (ext == "xlsx" || ext == "xls") {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        return("Error: Package 'readxl' needed for XLSX files.")
      }
      meta <- as.data.frame(readxl::read_excel(meta_path))
    } else {
      meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
    }

    issues <- character(0)
    info <- character(0)

    # Check key columns
    cols <- names(meta)
    info <- c(info, paste0("Columns found: ", paste(cols, collapse = ", ")))
    info <- c(info, paste0("Rows: ", nrow(meta)))

    # Detect potentially sensitive columns
    sensitive <- .detect_sensitive_columns(cols)
    if (length(sensitive) > 0) {
      info <- c(info,
        paste0("PRIVACY NOTE: Potentially sensitive column(s) detected: ",
               paste(sensitive, collapse = ", "),
               ". Consider using coded identifiers."))
    }

    # Check for expected columns
    expected <- c("file", "sample", "batch")
    missing <- setdiff(expected, cols)
    if (length(missing) > 0) {
      issues <- c(issues,
        paste0("Missing expected columns: ", paste(missing, collapse = ", "),
               ". Available: ", paste(cols, collapse = ", ")))
    }

    # Check file column against FCS directory
    if ("file" %in% cols) {
      fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$",
                              ignore.case = TRUE)
      meta_files <- meta$file
      not_found <- setdiff(meta_files, fcs_files)
      if (length(not_found) > 0) {
        issues <- c(issues,
          paste0(length(not_found), " metadata file(s) not found in FCS dir: ",
                 paste(utils::head(not_found, 5), collapse = ", ")))
      }
      extra_fcs <- setdiff(fcs_files, meta_files)
      if (length(extra_fcs) > 0) {
        info <- c(info,
          paste0(length(extra_fcs),
                 " FCS file(s) in directory not listed in metadata"))
      }
    }

    # Check batch distribution
    if ("batch" %in% cols) {
      batch_tab <- table(meta$batch)
      info <- c(info,
        paste0("Batch distribution: ",
               paste(paste0(names(batch_tab), "=", batch_tab),
                     collapse = ", ")))
    }

    # Build result
    result <- paste0("Validation results:\n",
                     paste0("  ", info, collapse = "\n"))
    if (length(issues) > 0) {
      result <- paste0(result, "\n\nIssues:\n",
                       paste0("  WARNING: ", issues, collapse = "\n"))
    } else {
      result <- paste0(result, "\n\nNo issues found.")
    }
    result
  }, error = function(e) {
    paste0("Error validating metadata: ", conditionMessage(e))
  })
}

#' Detect and classify channels from an FCS file
#' @param file_path Path to a single FCS file.
#' @return Character summary of channel classification.
#' @noRd
.tool_detect_channels <- function(file_path) {
  if (!file.exists(file_path)) {
    return(paste0("Error: File does not exist: ", file_path))
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    return("Error: Package 'flowCore' is not installed.")
  }

  tryCatch({
    ff <- flowCore::read.FCS(file_path, truncate_max_range = FALSE,
                             which.lines = 1)
    all_channels <- flowCore::colnames(ff)

    # Classify channels
    scatter_pat <- "^(FSC|SSC)"
    time_pat <- "^Time$"

    scatter <- grep(scatter_pat, all_channels, value = TRUE,
                    ignore.case = TRUE)
    time_ch <- grep(time_pat, all_channels, value = TRUE, ignore.case = TRUE)
    fluor <- setdiff(all_channels, c(scatter, time_ch))

    # Get marker names if available
    params <- flowCore::pData(flowCore::parameters(ff))
    marker_map <- stats::setNames(params$desc, params$name)
    marker_map <- marker_map[!is.na(marker_map) & nzchar(marker_map)]

    fluor_with_markers <- vapply(fluor, function(ch) {
      m <- marker_map[ch]
      if (!is.na(m) && nzchar(m)) paste0(ch, " -> ", m) else ch
    }, character(1), USE.NAMES = FALSE)

    paste0(
      "Channel classification for: ", basename(file_path), "\n\n",
      "Scatter channels (", length(scatter), "):\n",
      paste0("  ", scatter, collapse = "\n"), "\n\n",
      "Time channel(s): ",
      if (length(time_ch) > 0) paste(time_ch, collapse = ", ") else "none",
      "\n\n",
      "Fluorochrome channels (", length(fluor), "):\n",
      paste0("  ", fluor_with_markers, collapse = "\n")
    )
  }, error = function(e) {
    paste0("Error reading FCS file: ", conditionMessage(e))
  })
}

#' Check batch balance in sample metadata
#' @param meta_path Path to a CSV or XLSX metadata file.
#' @param batch_col Name of the batch column.
#' @param sensitivity Privacy level: "standard" or "strict".
#' @return Character summary of batch distribution.
#' @noRd
.tool_check_batch_balance <- function(meta_path, batch_col,
                                      sensitivity = "standard") {
  if (!file.exists(meta_path)) {
    return(paste0("Error: File not found: ", meta_path))
  }

  tryCatch({
    ext <- tolower(tools::file_ext(meta_path))
    if (ext == "xlsx" || ext == "xls") {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        return("Error: Package 'readxl' needed for XLSX files.")
      }
      meta <- as.data.frame(readxl::read_excel(meta_path))
    } else {
      meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
    }

    if (!batch_col %in% names(meta)) {
      return(paste0("Error: Column '", batch_col,
                    "' not found. Available: ",
                    paste(names(meta), collapse = ", ")))
    }

    batches <- meta[[batch_col]]
    batch_tab <- table(batches)
    n_batches <- length(batch_tab)
    sizes <- as.integer(batch_tab)

    result <- paste0(
      "Batch balance summary (column: '", batch_col, "'):\n",
      "  Number of batches: ", n_batches, "\n",
      "  Total samples: ", nrow(meta), "\n",
      "  Samples per batch:\n",
      paste0("    ", names(batch_tab), ": ", batch_tab, collapse = "\n")
    )

    # Check balance
    if (n_batches < 2) {
      result <- paste0(result,
        "\n\n  WARNING: Only 1 batch found. Batch correction requires",
        "\n  at least 2 batches.")
    } else if (max(sizes) > 2 * min(sizes)) {
      result <- paste0(result,
        "\n\n  WARNING: Batches are unbalanced (largest is >2x smallest).",
        "\n  This may affect batch correction quality.")
    } else {
      result <- paste0(result,
        "\n\n  Batches appear reasonably balanced.")
    }

    # Check condition distribution if present
    cond_cols <- intersect(c("condition", "group", "treatment"), names(meta))
    if (length(cond_cols) > 0) {
      cond_col <- cond_cols[1]
      cross_tab <- table(meta[[batch_col]], meta[[cond_col]])
      result <- paste0(result,
        "\n\n  Condition ('", cond_col, "') × Batch cross-tabulation:\n",
        paste(utils::capture.output(print(cross_tab)), collapse = "\n"))
      if (sensitivity == "strict") {
        result <- paste0(result,
          "\n\n  NOTE: Condition labels are shown for design validation. ",
          "If these contain sensitive clinical information, ",
          "consider using coded labels.")
      }
    }

    # Report sensitive columns in metadata
    sensitive <- .detect_sensitive_columns(names(meta))
    if (length(sensitive) > 0) {
      result <- paste0(result,
        "\n\n  PRIVACY NOTE: Potentially sensitive column(s) detected: ",
        paste(sensitive, collapse = ", "), ".")
    }

    result
  }, error = function(e) {
    paste0("Error checking batch balance: ", conditionMessage(e))
  })
}

#' Validate that specified markers exist in FCS data
#' @param markers Character vector of marker names to validate.
#' @param fcs_file Path to a single FCS file.
#' @return Character summary of validation results.
#' @noRd
.tool_validate_markers <- function(markers, fcs_file) {
  if (!file.exists(fcs_file)) {
    return(paste0("Error: File not found: ", fcs_file))
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    return("Error: Package 'flowCore' is not installed.")
  }

  tryCatch({
    ff <- flowCore::read.FCS(fcs_file, truncate_max_range = FALSE,
                             which.lines = 1)
    all_channels <- flowCore::colnames(ff)

    # Also check marker descriptions
    params <- flowCore::pData(flowCore::parameters(ff))
    all_markers <- params$desc[!is.na(params$desc) & nzchar(params$desc)]
    all_names <- unique(c(all_channels, all_markers))

    found <- markers[markers %in% all_names]
    missing <- markers[!markers %in% all_names]

    result <- paste0("Marker validation against: ", basename(fcs_file), "\n")
    if (length(found) > 0) {
      result <- paste0(result,
        "  Found (", length(found), "): ",
        paste(found, collapse = ", "), "\n")
    }
    if (length(missing) > 0) {
      result <- paste0(result,
        "  NOT FOUND (", length(missing), "): ",
        paste(missing, collapse = ", "), "\n",
        "\n  Available channels/markers:\n",
        paste0("    ", all_names, collapse = "\n"))
    } else {
      result <- paste0(result, "  All markers found.")
    }
    result
  }, error = function(e) {
    paste0("Error validating markers: ", conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Tool registration
# ---------------------------------------------------------------------------

#' Register all SpectraWeaveR tools with a chat object
#'
#' Uses closure-based wrappers so that the \code{sensitivity} level is
#' baked into each tool function — the LLM cannot override it.
#'
#' @param chat An ellmer Chat object.
#' @param include_btw Logical; if \code{TRUE} (default), also register btw
#'   tools for R environment introspection when the btw package is available.
#' @param sensitivity Character; \code{"standard"} (default) auto-redacts PII
#'   in metadata previews; \code{"strict"} suppresses all row previews and
#'   disables btw env/files tools.
#' @return The chat object (modified in place) with tools registered.
#' @noRd
.register_tools <- function(chat, include_btw = TRUE,
                            sensitivity = "standard") {
  # --- btw tools for R environment introspection ---
  if (isTRUE(include_btw) && requireNamespace("btw", quietly = TRUE)) {
    btw_groups <- if (sensitivity == "strict") {
      c("docs", "session")
    } else {
      c("env", "docs", "session", "files")
    }
    tryCatch({
      chat$register_tools(btw::btw_tools(btw_groups))
      message("Registered btw tools (", paste(btw_groups, collapse = ", "),
              ") for R environment introspection.")
    }, error = function(e) {
      warning("Could not register btw tools: ", e$message, call. = FALSE)
    })
  }

  # --- Closure wrappers that capture `sensitivity` ---
  # The LLM sees columns_only as a tool parameter but cannot override
  # the sensitivity level, which is controlled by the user.
  local_read_csv <- function(file_path, n_rows = 5L, columns_only = FALSE) {
    .tool_read_csv_columns(file_path, n_rows = n_rows,
                           columns_only = columns_only,
                           sensitivity = sensitivity)
  }

  local_read_xlsx <- function(file_path, sheet = 1L, n_rows = 5L,
                              columns_only = FALSE) {
    .tool_read_xlsx_columns(file_path, sheet = sheet, n_rows = n_rows,
                            columns_only = columns_only,
                            sensitivity = sensitivity)
  }

  local_check_batch <- function(meta_path, batch_col) {
    .tool_check_batch_balance(meta_path, batch_col,
                              sensitivity = sensitivity)
  }

  # --- SpectraWeaveR-specific tools ---
  chat$register_tool(ellmer::tool(
    .tool_list_fcs_files,
    "List FCS files in a directory to help identify sample files.",
    dir_path = ellmer::type_string(
      "Absolute path to directory containing FCS files."
    )
  ))

  chat$register_tool(ellmer::tool(
    .tool_list_directory,
    "List files and subdirectories at a given path.",
    dir_path = ellmer::type_string("Absolute path to a directory.")
  ))

  chat$register_tool(ellmer::tool(
    .tool_read_fcs_header,
    "Read channel and marker names from the header of an FCS file (no expression data).",
    file_path = ellmer::type_string("Absolute path to an FCS file.")
  ))

  chat$register_tool(ellmer::tool(
    local_read_csv,
    "Preview column names and first few rows of a CSV file. Potentially sensitive columns are automatically redacted.",
    file_path = ellmer::type_string("Absolute path to a CSV file."),
    n_rows = ellmer::type_integer(
      "Number of rows to preview (default 5).",
      required = FALSE
    ),
    columns_only = ellmer::type_boolean(
      "If TRUE, return only column names and types without row values. Recommended for files that may contain sensitive information.",
      required = FALSE
    )
  ))

  chat$register_tool(ellmer::tool(
    local_read_xlsx,
    "Preview column names and first few rows of an Excel XLSX file. Potentially sensitive columns are automatically redacted.",
    file_path = ellmer::type_string("Absolute path to an XLSX file."),
    sheet = ellmer::type_integer(
      "Sheet index (default 1).",
      required = FALSE
    ),
    n_rows = ellmer::type_integer(
      "Number of rows to preview (default 5).",
      required = FALSE
    ),
    columns_only = ellmer::type_boolean(
      "If TRUE, return only column names and types without row values. Recommended for files that may contain sensitive information.",
      required = FALSE
    )
  ))

  chat$register_tool(ellmer::tool(
    .tool_validate_sample_meta,
    "Validate that sample metadata has required columns and matches FCS files.",
    meta_path = ellmer::type_string(
      "Absolute path to a CSV or XLSX metadata file."
    ),
    fcs_dir = ellmer::type_string(
      "Absolute path to the FCS files directory."
    )
  ))

  chat$register_tool(ellmer::tool(
    .tool_detect_channels,
    "Read an FCS file and classify channels into scatter, time, and fluorochrome categories.",
    file_path = ellmer::type_string("Absolute path to an FCS file.")
  ))

  chat$register_tool(ellmer::tool(
    local_check_batch,
    "Summarise the sample distribution across batches in a metadata file.",
    meta_path = ellmer::type_string(
      "Absolute path to a CSV or XLSX metadata file."
    ),
    batch_col = ellmer::type_string(
      "Name of the column containing batch labels."
    )
  ))

  chat$register_tool(ellmer::tool(
    .tool_validate_markers,
    "Check that specified marker names exist in an FCS file.",
    markers = ellmer::type_array(
      ellmer::type_string("A marker name."),
      "Character vector of marker names to validate."
    ),
    fcs_file = ellmer::type_string("Absolute path to an FCS file.")
  ))

  invisible(chat)
}

# ---------------------------------------------------------------------------
# Chat object creation
# ---------------------------------------------------------------------------

#' Create an LLM Chat Object for SpectraWeaveR
#'
#' Creates and configures a chat object from ellmer with the SpectraWeaveR
#' system prompt and tools pre-registered.
#'
#' @param provider Character string specifying the LLM provider. One of
#'   \code{"openai"} (default), \code{"anthropic"}, \code{"google"},
#'   \code{"ollama"}, or any provider supported by ellmer.
#' @param model Character string specifying the model name (optional).
#'   If \code{NULL}, the provider's default model is used.
#' @param system_prompt Character string or \code{NULL}. If \code{NULL},
#'   the default pipeline builder prompt is loaded.
#' @param include_btw Logical; if \code{TRUE} (default), register btw tools
#'   for R environment introspection alongside SpectraWeaveR tools.
#' @param sensitivity Character; \code{"standard"} (auto-redacts PII) or
#'   \code{"strict"} (columns only, no row previews, limited btw tools).
#' @param ... Additional arguments passed to the ellmer chat constructor
#'   (e.g. \code{api_key}).
#'
#' @return An ellmer \code{Chat} object with SpectraWeaveR tools registered.
#'
#' @noRd
.create_chat <- function(provider = "openai",
                         model = NULL,
                         system_prompt = NULL,
                         include_btw = TRUE,
                         sensitivity = "standard",
                         ...) {
  .check_ellmer()

  if (is.null(system_prompt)) {
    system_prompt <- .load_prompt("pipeline_builder")
  }

  sensitivity <- match.arg(sensitivity, c("standard", "strict"))

  # Append btw context hint to system prompt when btw is available
  if (isTRUE(include_btw) && requireNamespace("btw", quietly = TRUE)) {
    btw_hint <- if (sensitivity == "strict") {
      paste0(
        "\n\nYou have access to btw tools for R documentation and session ",
        "info. Environment inspection and file exploration tools are ",
        "disabled in strict privacy mode to protect sensitive data."
      )
    } else {
      paste0(
        "\n\nYou have access to btw tools for R environment introspection. ",
        "You can use them to describe data frames in the user's session, ",
        "read package documentation, inspect session info, and explore files. ",
        "Use these tools proactively when the user mentions data they have ",
        "loaded in R."
      )
    }
    system_prompt <- paste0(system_prompt, btw_hint)
  }

  # Build arguments for provider constructor
  args <- list(system_prompt = system_prompt, ...)
  if (!is.null(model)) {
    args$model <- model
  }

  # Create chat using the appropriate provider
  chat <- switch(provider,
    openai = do.call(ellmer::chat_openai, args),
    anthropic = do.call(ellmer::chat_anthropic, args),
    google = do.call(ellmer::chat_google_gemini, args),
    ollama = do.call(ellmer::chat_ollama, args),
    stop("Unknown provider '", provider,
         "'. Use 'openai', 'anthropic', 'google', or 'ollama'.",
         call. = FALSE)
  )

  .register_tools(chat, include_btw = include_btw, sensitivity = sensitivity)
  chat
}

# ---------------------------------------------------------------------------
# Main user-facing functions
# ---------------------------------------------------------------------------

#' Launch the SpectraWeaveR LLM Pipeline Builder
#'
#' Starts an interactive conversation with an LLM assistant that helps you
#' configure a SpectraWeaveR analysis pipeline. The assistant asks about
#' your experiment, inspects your files and directories, and generates
#' runnable R code.
#'
#' @section Prerequisites:
#' \itemize{
#'   \item Install ellmer: \code{install.packages("ellmer")}
#'   \item Set an API key for your chosen provider (e.g.
#'     \code{Sys.setenv(OPENAI_API_KEY = "sk-...")}).
#'   \item For local models, install and start Ollama.
#' }
#'
#' @section Privacy:
#' The assistant sends file paths, directory listings, column names, and
#' channel names to the LLM provider. No expression data (cell
#' measurements) is transmitted. For fully local operation, use
#' \code{provider = "ollama"}.
#'
#' @param provider Character string specifying the LLM provider. One of
#'   \code{"openai"} (default), \code{"anthropic"}, \code{"google"},
#'   \code{"ollama"}.
#' @param model Optional model name. If \code{NULL}, provider default is
#'   used.
#' @param interactive Logical; if \code{TRUE} (default), launches a live
#'   console session. If \code{FALSE}, returns the configured \code{Chat}
#'   object for programmatic use.
#' @param privacy_notice Logical; if \code{TRUE} (default), displays a
#'   privacy notice before starting.
#' @param include_btw Logical; if \code{TRUE} (default), also register
#'   \pkg{btw} tools for R environment introspection (describing data
#'   frames, reading documentation, inspecting session info). Requires the
#'   \pkg{btw} package to be installed; silently skipped if not available.
#' @param sensitivity Character; privacy sensitivity level.
#'   \code{"standard"} (default) auto-detects and redacts columns matching
#'   common PII patterns (patient ID, name, DOB, MRN, etc.) in metadata
#'   previews. \code{"strict"} suppresses all metadata row previews
#'   (column names only), disables btw env/files tools, and adds advisory
#'   notes. Use \code{"strict"} for clinical data with patient identifiers.
#' @param ... Additional arguments passed to the ellmer chat constructor.
#'
#' @return If \code{interactive = TRUE}, returns the \code{Chat} object
#'   invisibly after the session ends. If \code{interactive = FALSE},
#'   returns the \code{Chat} object directly.
#'
#' @examples
#' \dontrun{
#' # Interactive session with OpenAI
#' sw_assistant()
#'
#' # Use Anthropic Claude
#' sw_assistant(provider = "anthropic")
#'
#' # Strict privacy for clinical data
#' sw_assistant(sensitivity = "strict")
#'
#' # Maximum privacy: local model + strict mode
#' sw_assistant(provider = "ollama", model = "llama3",
#'              sensitivity = "strict")
#' }
#'
#' @export
sw_assistant <- function(provider = "openai",
                         model = NULL,
                         interactive = TRUE,
                         privacy_notice = TRUE,
                         include_btw = TRUE,
                         sensitivity = c("standard", "strict"),
                         ...) {
  .check_ellmer()
  sensitivity <- match.arg(sensitivity)

  if (isTRUE(privacy_notice) && isTRUE(interactive)) {
    btw_note <- ""
    if (isTRUE(include_btw) && requireNamespace("btw", quietly = TRUE)) {
      btw_note <- if (sensitivity == "strict") {
        "  btw tools: docs and session only (env/files disabled).\n"
      } else {
        paste0(
          "  btw tools can inspect your R environment (data frames, files).\n",
          "  Disable with: include_btw = FALSE\n"
        )
      }
    }
    message(
      "---\n",
      "SpectraWeaveR LLM Assistant\n",
      "---\n",
      "This assistant sends file paths, directory listings, and column names\n",
      "to the '", provider, "' LLM provider.\n\n",
      "PRIVACY:\n",
      "  No cell expression data (fluorescence intensities) is transmitted.\n",
      "  Metadata previews: columns matching PII patterns are auto-redacted.\n",
      if (sensitivity == "strict") {
        "  STRICT MODE: Row previews fully suppressed. btw env/files disabled.\n"
      } else {
        "  For clinical data, use: sensitivity = \"strict\"\n"
      },
      btw_note,
      "  For fully local operation: provider = \"ollama\"\n",
      "---"
    )
  }

  chat <- .create_chat(provider = provider, model = model,
                        include_btw = include_btw,
                        sensitivity = sensitivity, ...)

  if (isTRUE(interactive)) {
    ellmer::live_console(chat)
    invisible(chat)
  } else {
    chat
  }
}

#' Configure an LLM Provider for SpectraWeaveR
#'
#' Creates a configured chat object without starting an interactive session.
#' Useful for programmatic pipeline generation.
#'
#' @param provider Character; LLM provider name.
#' @param model Character; model name (optional).
#' @param include_btw Logical; if \code{TRUE} (default), register btw tools
#'   for R environment introspection.
#' @param ... Additional arguments passed to the ellmer chat constructor.
#'
#' @return An ellmer \code{Chat} object with SpectraWeaveR tools registered.
#'
#' @examples
#' \dontrun{
#' chat <- sw_assistant_configure(provider = "openai", model = "gpt-4o")
#' response <- chat$chat("Help me set up batch correction")
#' }
#'
#' @export
sw_assistant_configure <- function(provider = "openai",
                                   model = NULL,
                                   include_btw = TRUE,
                                   sensitivity = c("standard", "strict"),
                                   ...) {
  sensitivity <- match.arg(sensitivity)
  .create_chat(provider = provider, model = model,
               include_btw = include_btw, sensitivity = sensitivity, ...)
}

# ---------------------------------------------------------------------------
# Specialized assistants
# ---------------------------------------------------------------------------

#' Launch a Batch Correction Assistant
#'
#' A focused assistant that helps set up batch correction for spectral flow
#' cytometry data using SpectraWeaveR's cyCombine-based workflow.
#'
#' @inheritParams sw_assistant
#'
#' @return See \code{\link{sw_assistant}}.
#'
#' @examples
#' \dontrun{
#' sw_assistant_batch_correction()
#' }
#'
#' @export
sw_assistant_batch_correction <- function(provider = "openai",
                                          model = NULL,
                                          interactive = TRUE,
                                          privacy_notice = TRUE,
                                          include_btw = TRUE,
                                          sensitivity = c("standard", "strict"),
                                          ...) {
  .check_ellmer()
  sensitivity <- match.arg(sensitivity)

  prompt <- .load_prompt("batch_correction_builder")
  chat <- .create_chat(provider = provider, model = model,
                        system_prompt = prompt,
                        include_btw = include_btw,
                        sensitivity = sensitivity, ...)

  if (isTRUE(privacy_notice) && isTRUE(interactive)) {
    message(
      "---\n",
      "SpectraWeaveR Batch Correction Assistant\n",
      "---\n",
      "This assistant helps you configure batch correction.\n",
      "File paths and metadata column names are sent to '", provider, "'.\n",
      "---"
    )
  }

  if (isTRUE(interactive)) {
    ellmer::live_console(chat)
    invisible(chat)
  } else {
    chat
  }
}

#' Launch a Gating Assistant
#'
#' A focused assistant that helps create gating templates and apply
#' automated gating to spectral flow cytometry data.
#'
#' @inheritParams sw_assistant
#'
#' @return See \code{\link{sw_assistant}}.
#'
#' @examples
#' \dontrun{
#' sw_assistant_gating()
#' }
#'
#' @export
sw_assistant_gating <- function(provider = "openai",
                                model = NULL,
                                interactive = TRUE,
                                privacy_notice = TRUE,
                                include_btw = TRUE,
                                sensitivity = c("standard", "strict"),
                                ...) {
  .check_ellmer()
  sensitivity <- match.arg(sensitivity)

  prompt <- .load_prompt("gating_builder")
  chat <- .create_chat(provider = provider, model = model,
                        system_prompt = prompt,
                        include_btw = include_btw,
                        sensitivity = sensitivity, ...)

  if (isTRUE(privacy_notice) && isTRUE(interactive)) {
    message(
      "---\n",
      "SpectraWeaveR Gating Assistant\n",
      "---\n",
      "This assistant helps you set up automated gating.\n",
      "File paths and channel names are sent to '", provider, "'.\n",
      "---"
    )
  }

  if (isTRUE(interactive)) {
    ellmer::live_console(chat)
    invisible(chat)
  } else {
    chat
  }
}

#' Launch an Unmixing Assistant
#'
#' A focused assistant that helps configure spectral unmixing from raw
#' spectral flow cytometry data via the AutoSpectral workflow.
#'
#' @inheritParams sw_assistant
#'
#' @return See \code{\link{sw_assistant}}.
#'
#' @examples
#' \dontrun{
#' sw_assistant_unmixing()
#' }
#'
#' @export
sw_assistant_unmixing <- function(provider = "openai",
                                  model = NULL,
                                  interactive = TRUE,
                                  privacy_notice = TRUE,
                                  include_btw = TRUE,
                                  sensitivity = c("standard", "strict"),
                                  ...) {
  .check_ellmer()
  sensitivity <- match.arg(sensitivity)

  prompt <- .load_prompt("unmixing_builder")
  chat <- .create_chat(provider = provider, model = model,
                        system_prompt = prompt,
                        include_btw = include_btw,
                        sensitivity = sensitivity, ...)

  if (isTRUE(privacy_notice) && isTRUE(interactive)) {
    message(
      "---\n",
      "SpectraWeaveR Unmixing Assistant\n",
      "---\n",
      "This assistant helps you configure spectral unmixing.\n",
      "File paths and directory listings are sent to '", provider, "'.\n",
      "---"
    )
  }

  if (isTRUE(interactive)) {
    ellmer::live_console(chat)
    invisible(chat)
  } else {
    chat
  }
}

# ---------------------------------------------------------------------------
# Pipeline code generation
# ---------------------------------------------------------------------------

#' Generate Pipeline Code from Configuration
#'
#' Takes a pipeline configuration list and generates clean, runnable R code
#' for a SpectraWeaveR analysis pipeline. Supports two output modes:
#' \code{run_pipeline()} for standard workflows, or composable
#' \code{sw_pipeline()} for custom workflows.
#'
#' @param config A named list with pipeline configuration. Expected fields:
#'   \describe{
#'     \item{\code{fcs_dir}}{Character; path to FCS files directory.}
#'     \item{\code{sample_meta_path}}{Character; path to sample metadata file.}
#'     \item{\code{markers}}{Character vector; marker names.}
#'     \item{\code{lineage_markers}}{Character vector; lineage markers for
#'       clustering.}
#'     \item{\code{batch_col}}{Character; name of the batch column in
#'       metadata.}
#'     \item{\code{sample_col}}{Character; name of the sample ID column.}
#'     \item{\code{file_col}}{Character; name of the column with FCS
#'       filenames.}
#'     \item{\code{condition_col}}{Character; condition column (optional).}
#'     \item{\code{cofactor}}{Numeric; arcsinh cofactor (default 6000).}
#'     \item{\code{n_metaclusters}}{Integer; number of metaclusters (default
#'       20).}
#'     \item{\code{gating_template}}{Character; path to gating template CSV
#'       or \code{NULL}.}
#'     \item{\code{output_dir}}{Character; output directory.}
#'     \item{\code{seed}}{Integer; random seed (default 42).}
#'     \item{\code{unmix_from_raw}}{Logical; whether to unmix from raw data.}
#'     \item{\code{cytometer}}{Character; cytometer type.}
#'     \item{\code{control_dir}}{Character; path to controls (for
#'       unmixing).}
#'     \item{\code{unstained_fcs}}{Character; path to unstained FCS.}
#'   }
#' @param style Character; code generation style.
#'   \code{"run_pipeline"} (default) generates a single \code{run_pipeline()}
#'   call.
#'   \code{"composable"} generates a \code{sw_pipeline()} workflow.
#' @param output Character; where to send the generated code.
#'   \code{"console"} (default) prints to the console.
#'   \code{"string"} returns the code as a character string.
#'   \code{"file"} writes to \code{path}.
#' @param path Character; file path for \code{output = "file"}.
#'
#' @return If \code{output = "string"}, returns the generated code as a
#'   character string. If \code{output = "console"} or \code{"file"},
#'   returns the code invisibly.
#'
#' @examples
#' config <- list(
#'   fcs_dir = "/data/fcs",
#'   sample_meta_path = "/data/metadata.csv",
#'   markers = c("CD3", "CD4", "CD8", "CD19"),
#'   lineage_markers = c("CD3", "CD4", "CD8"),
#'   batch_col = "batch",
#'   sample_col = "sample",
#'   file_col = "file",
#'   cofactor = 6000,
#'   n_metaclusters = 20,
#'   output_dir = "results",
#'   seed = 42
#' )
#' code <- sw_generate_pipeline_code(config, output = "string")
#' cat(code)
#'
#' @export
sw_generate_pipeline_code <- function(config,
                                      style = c("run_pipeline", "composable"),
                                      output = c("console", "string", "file"),
                                      path = NULL) {
  style <- match.arg(style)
  output <- match.arg(output)

  if (output == "file" && is.null(path)) {
    stop("'path' is required when output = 'file'.", call. = FALSE)
  }

  # Validate required fields
  required <- c("fcs_dir", "sample_meta_path", "markers", "lineage_markers")
  missing <- setdiff(required, names(config))
  if (length(missing) > 0) {
    stop("Missing required config fields: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Apply defaults
  config$cofactor <- .null_default(config$cofactor, 6000)
  config$n_metaclusters <- .null_default(config$n_metaclusters, 20L)
  config$seed <- .null_default(config$seed, 42L)
  config$output_dir <- .null_default(config$output_dir, "SpectraWeaveR_output")
  config$batch_col <- .null_default(config$batch_col, "batch")
  config$sample_col <- .null_default(config$sample_col, "sample")
  config$file_col <- .null_default(config$file_col, "file")
  config$unmix_from_raw <- .null_default(config$unmix_from_raw, FALSE)

  # Generate code
  code <- if (style == "run_pipeline") {
    .generate_run_pipeline_code(config)
  } else {
    .generate_composable_code(config)
  }

  # Output
  switch(output,
    console = {
      cat(code)
      invisible(code)
    },
    string = code,
    file = {
      writeLines(code, path)
      message("Pipeline code written to: ", path)
      invisible(code)
    }
  )
}

#' Generate run_pipeline() code
#' @param config Pipeline configuration list.
#' @return Character string of R code.
#' @noRd
.generate_run_pipeline_code <- function(config) {
  lines <- character(0)

  lines <- c(lines,
    "# SpectraWeaveR Pipeline",
    "# Generated by sw_generate_pipeline_code()",
    paste0("# ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "library(SpectraWeaveR)",
    ""
  )

  # Metadata loading
  meta_ext <- tolower(tools::file_ext(config$sample_meta_path))
  if (meta_ext %in% c("xlsx", "xls")) {
    lines <- c(lines,
      "# Load sample metadata",
      "library(readxl)",
      paste0('sample_meta <- as.data.frame(read_excel("',
             config$sample_meta_path, '"))'),
      ""
    )
  } else {
    lines <- c(lines,
      "# Load sample metadata",
      paste0('sample_meta <- read.csv("', config$sample_meta_path, '")'),
      ""
    )
  }

  # Column renaming if needed
  renames <- character(0)
  if (config$file_col != "file") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$file_col, '"] <- "file"'))
  }
  if (config$sample_col != "sample") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$sample_col, '"] <- "sample"'))
  }
  if (config$batch_col != "batch") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$batch_col, '"] <- "batch"'))
  }
  if (!is.null(config$condition_col) && config$condition_col != "condition") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$condition_col, '"] <- "condition"'))
  }

  if (length(renames) > 0) {
    lines <- c(lines,
      "# Rename columns to match SpectraWeaveR conventions",
      renames,
      ""
    )
  }

  # Markers
  lines <- c(lines,
    "# Define markers",
    paste0("markers <- c(",
           paste0('"', config$markers, '"', collapse = ", "), ")"),
    paste0("lineage_markers <- c(",
           paste0('"', config$lineage_markers, '"', collapse = ", "), ")"),
    ""
  )

  # Pipeline call
  pipeline_args <- c(
    paste0('  fcs_dir         = "', config$fcs_dir, '"'),
    "  sample_meta     = sample_meta",
    "  markers         = markers",
    "  lineage_markers = lineage_markers",
    paste0("  cofactor        = ", config$cofactor),
    paste0("  n_metaclusters  = ", config$n_metaclusters),
    paste0('  output_dir      = "', config$output_dir, '"'),
    paste0("  seed            = ", config$seed)
  )

  if (!is.null(config$gating_template)) {
    pipeline_args <- c(pipeline_args,
      paste0('  gating_template = "', config$gating_template, '"'))
  }

  if (isTRUE(config$unmix_from_raw)) {
    pipeline_args <- c(pipeline_args,
      "  unmix_from_raw  = TRUE",
      paste0('  control_dir     = "', .null_default(config$control_dir, ""), '"'),
      paste0('  unstained_fcs   = "', .null_default(config$unstained_fcs, ""), '"'),
      paste0('  cytometer       = "', .null_default(config$cytometer, "aurora"), '"')
    )
  }

  lines <- c(lines,
    "# Run the full pipeline",
    "results <- run_pipeline(",
    paste(pipeline_args, collapse = ",\n"),
    ")",
    "",
    "# Access results",
    "results$corrected            # batch-corrected expression data",
    "results$cluster_assignments  # per-cell cluster labels",
    "results$cluster_mfis         # median fluorescence per cluster",
    ""
  )

  paste(lines, collapse = "\n")
}

#' Generate composable sw_pipeline() code
#' @param config Pipeline configuration list.
#' @return Character string of R code.
#' @noRd
.generate_composable_code <- function(config) {
  lines <- character(0)

  lines <- c(lines,
    "# SpectraWeaveR Composable Pipeline",
    "# Generated by sw_generate_pipeline_code()",
    paste0("# ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "library(SpectraWeaveR)",
    ""
  )

  # Metadata loading
  meta_ext <- tolower(tools::file_ext(config$sample_meta_path))
  if (meta_ext %in% c("xlsx", "xls")) {
    lines <- c(lines,
      "# Load sample metadata",
      "library(readxl)",
      paste0('sample_meta <- as.data.frame(read_excel("',
             config$sample_meta_path, '"))'),
      ""
    )
  } else {
    lines <- c(lines,
      "# Load sample metadata",
      paste0('sample_meta <- read.csv("', config$sample_meta_path, '")'),
      ""
    )
  }

  # Column renaming
  renames <- character(0)
  if (config$file_col != "file") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$file_col, '"] <- "file"'))
  }
  if (config$sample_col != "sample") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$sample_col, '"] <- "sample"'))
  }
  if (config$batch_col != "batch") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$batch_col, '"] <- "batch"'))
  }
  if (!is.null(config$condition_col) && config$condition_col != "condition") {
    renames <- c(renames,
      paste0('names(sample_meta)[names(sample_meta) == "',
             config$condition_col, '"] <- "condition"'))
  }

  if (length(renames) > 0) {
    lines <- c(lines,
      "# Rename columns to match SpectraWeaveR conventions",
      renames,
      ""
    )
  }

  # Markers
  lines <- c(lines,
    "# Define markers",
    paste0("markers <- c(",
           paste0('"', config$markers, '"', collapse = ", "), ")"),
    paste0("lineage_markers <- c(",
           paste0('"', config$lineage_markers, '"', collapse = ", "), ")"),
    ""
  )

  # FCS loading
  lines <- c(lines,
    "# Step 1: Load FCS files",
    paste0('fcs_files <- list.files("', config$fcs_dir,
           '", pattern = "\\\\.fcs$",'),
    "                        full.names = TRUE, ignore.case = TRUE)",
    "fs <- sw_read_fcs(fcs_files)",
    ""
  )

  # Build pipeline
  lines <- c(lines,
    "# Build the composable pipeline",
    'pip <- sw_pipeline("spectral_analysis")',
    "",
    "# Step 2: Remove margin events (must be before transformation)",
    'pip <- sw_pipeline_add(pip, sw_step("remove_margins", sw_remove_margins))',
    "",
    "# Step 3: Signal QC (PeacoQC)",
    'pip <- sw_pipeline_add(pip, sw_step("signal_qc", sw_signal_qc,',
    "                                     list(IT_limit = 0.55, MAD = 6)))",
    ""
  )

  # Gating (optional)
  if (!is.null(config$gating_template)) {
    lines <- c(lines,
      "# Step 4: Automated gating",
      paste0('gating_template <- "', config$gating_template, '"'),
      "# Note: Gating requires a separate workflow with sw_gate()",
      ""
    )
  }

  lines <- c(lines,
    "# Execute the pipeline on each sample",
    "# (Composable pipeline operates on individual flowFrames)",
    "processed <- lapply(seq_along(fs), function(i) {",
    "  sw_pipeline_run(pip, input = fs[[i]], trace = TRUE)",
    "})",
    "names(processed) <- flowCore::sampleNames(fs)",
    "",
    "# Step 5: Prepare for batch correction",
    "uncorrected <- sw_prepare_for_correction(",
    "  ff_list     = processed,",
    "  sample_meta = sample_meta,",
    "  markers     = markers,",
    paste0("  cofactor    = ", config$cofactor),
    ")",
    "",
    "# Step 6: Batch correction",
    "corrected <- sw_batch_correct(",
    "  uncorrected,",
    "  markers     = markers,",
    paste0("  seed        = ", config$seed),
    ")",
    "",
    "# Step 7: Clustering",
    "cluster_result <- sw_cluster(",
    "  corrected,",
    "  lineage_markers = lineage_markers,",
    paste0("  n_metaclusters  = ", config$n_metaclusters, ","),
    paste0("  seed            = ", config$seed),
    ")",
    "",
    "# Access results",
    "assignments <- sw_get_cluster_assignments(cluster_result)",
    "mfis <- sw_cluster_mfis(cluster_result)",
    'sw_plot_clusters(cluster_result, "clusters.pdf")',
    ""
  )

  paste(lines, collapse = "\n")
}

#' Generate a Pipeline via One-Shot LLM Request
#'
#' Provides a natural language description of your experiment and data, and
#' receives generated SpectraWeaveR pipeline code in return — without an
#' interactive conversation.
#'
#' @param description Character string describing the experiment, data
#'   location, and analysis goals.
#' @param provider Character; LLM provider (default \code{"openai"}).
#' @param model Character; model name (optional).
#' @param ... Additional arguments passed to the ellmer chat constructor.
#'
#' @return Character string containing generated R code.
#'
#' @examples
#' \dontrun{
#' code <- sw_quick_pipeline(
#'   description = "I have 20 Aurora 5L FCS files in /data/fcs/,
#'     metadata in /data/meta.csv with columns filename, patient_id,
#'     batch_number, treatment. Markers: CD3, CD4, CD8, CD19, CD56.",
#'   provider = "openai"
#' )
#' cat(code)
#' }
#'
#' @export
sw_quick_pipeline <- function(description,
                              provider = "openai",
                              model = NULL,
                              ...) {
  .check_ellmer()

  if (!is.character(description) || !nzchar(description)) {
    stop("'description' must be a non-empty character string.", call. = FALSE)
  }

  prompt <- paste0(
    .load_prompt("pipeline_builder"),
    "\n\n## Mode: One-Shot Pipeline Generation\n",
    "The user will provide a single description of their experiment.\n",
    "Generate a complete, runnable SpectraWeaveR pipeline script.\n",
    "Use tools to inspect any file paths mentioned.\n",
    "Output ONLY the R code (no markdown fences)."
  )

  chat <- .create_chat(provider = provider, model = model,
                        system_prompt = prompt, ...)

  response <- chat$chat(description, echo = "none")
  response
}

# ---------------------------------------------------------------------------
# Pipeline config type (for structured extraction)
# ---------------------------------------------------------------------------

#' Get the Pipeline Configuration Type Schema
#'
#' Returns an ellmer type specification for structured extraction of a
#' pipeline configuration from an LLM conversation. Used internally by
#' the assistant to produce validated configuration objects.
#'
#' @return An ellmer type object describing the pipeline config schema.
#'
#' @examples
#' \dontrun{
#' config_type <- sw_pipeline_config_type()
#' chat <- sw_assistant(interactive = FALSE)
#' config <- chat$chat_structured(
#'   "Generate config for my experiment...",
#'   type = config_type
#' )
#' }
#'
#' @export
sw_pipeline_config_type <- function() {
  .check_ellmer()

  ellmer::type_object(
    fcs_dir = ellmer::type_string("Path to FCS files directory."),
    sample_meta_path = ellmer::type_string(
      "Path to sample metadata CSV or XLSX file."
    ),
    sample_col = ellmer::type_string("Column name for sample IDs."),
    batch_col = ellmer::type_string("Column name for batch labels."),
    condition_col = ellmer::type_string(
      "Column for biological condition (optional).",
      required = FALSE
    ),
    file_col = ellmer::type_string("Column containing FCS filenames."),
    markers = ellmer::type_array(
      ellmer::type_string("Marker name."),
      "All marker names for the pipeline."
    ),
    lineage_markers = ellmer::type_array(
      ellmer::type_string("Lineage marker name."),
      "Lineage markers used for clustering."
    ),
    cofactor = ellmer::type_number(
      "Arcsinh cofactor (default 6000 for spectral flow)."
    ),
    n_metaclusters = ellmer::type_integer(
      "Number of metaclusters (default 20)."
    ),
    gating_template = ellmer::type_string(
      "Path to gating template CSV, or null to skip.",
      required = FALSE
    ),
    gate_node = ellmer::type_string(
      "Gate node to extract (default 'lymphocytes').",
      required = FALSE
    ),
    output_dir = ellmer::type_string("Output directory for results."),
    unmix_from_raw = ellmer::type_boolean(
      "Whether to unmix from raw spectral data."
    ),
    cytometer = ellmer::type_enum(
      c("aurora", "auroraNL", "id7000", "s8", "a8", "a5se",
        "opteon", "mosaic", "xenith"),
      "Cytometer type.",
      required = FALSE
    ),
    control_dir = ellmer::type_string(
      "Path to control FCS files (for unmixing).",
      required = FALSE
    ),
    unstained_fcs = ellmer::type_string(
      "Path to unstained FCS file(s) (for unmixing).",
      required = FALSE
    ),
    seed = ellmer::type_integer("Random seed for reproducibility."),
    clustering_method = ellmer::type_enum(
      c("som", "fastpg"),
      "Clustering method.",
      required = FALSE
    ),
    norm_method = ellmer::type_enum(
      c("scale", "rank", "none"),
      "Normalization method for batch correction.",
      required = FALSE
    )
  )
}


# ---------------------------------------------------------------------------
# MCP server
# ---------------------------------------------------------------------------

#' Build the list of SpectraWeaveR tools as ellmer tool objects
#'
#' Uses closure wrappers for tools that need the \code{sensitivity} level.
#'
#' @param sensitivity Character; \code{"standard"} or \code{"strict"}.
#' @return A list of ellmer tool objects.
#' @noRd
.sw_tool_list <- function(sensitivity = "standard") {
  .check_ellmer()

  # Closure wrappers capturing sensitivity
  local_read_csv <- function(file_path, n_rows = 5L, columns_only = FALSE) {
    .tool_read_csv_columns(file_path, n_rows = n_rows,
                           columns_only = columns_only,
                           sensitivity = sensitivity)
  }
  local_read_xlsx <- function(file_path, sheet = 1L, n_rows = 5L,
                              columns_only = FALSE) {
    .tool_read_xlsx_columns(file_path, sheet = sheet, n_rows = n_rows,
                            columns_only = columns_only,
                            sensitivity = sensitivity)
  }
  local_check_batch <- function(meta_path, batch_col) {
    .tool_check_batch_balance(meta_path, batch_col,
                              sensitivity = sensitivity)
  }

  list(
    ellmer::tool(
      .tool_list_fcs_files,
      "List FCS files in a directory to help identify sample files.",
      dir_path = ellmer::type_string(
        "Absolute path to directory containing FCS files."
      )
    ),
    ellmer::tool(
      .tool_list_directory,
      "List files and subdirectories at a given path.",
      dir_path = ellmer::type_string("Absolute path to a directory.")
    ),
    ellmer::tool(
      .tool_read_fcs_header,
      "Read channel and marker names from the header of an FCS file (no expression data).",
      file_path = ellmer::type_string("Absolute path to an FCS file.")
    ),
    ellmer::tool(
      local_read_csv,
      "Preview column names and first few rows of a CSV file. Potentially sensitive columns are automatically redacted.",
      file_path = ellmer::type_string("Absolute path to a CSV file."),
      n_rows = ellmer::type_integer(
        "Number of rows to preview (default 5).",
        required = FALSE
      ),
      columns_only = ellmer::type_boolean(
        "If TRUE, return only column names and types without row values.",
        required = FALSE
      )
    ),
    ellmer::tool(
      local_read_xlsx,
      "Preview column names and first few rows of an Excel XLSX file. Potentially sensitive columns are automatically redacted.",
      file_path = ellmer::type_string("Absolute path to an XLSX file."),
      sheet = ellmer::type_integer(
        "Sheet index (default 1).",
        required = FALSE
      ),
      n_rows = ellmer::type_integer(
        "Number of rows to preview (default 5).",
        required = FALSE
      ),
      columns_only = ellmer::type_boolean(
        "If TRUE, return only column names and types without row values.",
        required = FALSE
      )
    ),
    ellmer::tool(
      .tool_validate_sample_meta,
      "Validate that sample metadata has required columns and matches FCS files.",
      meta_path = ellmer::type_string(
        "Absolute path to a CSV or XLSX metadata file."
      ),
      fcs_dir = ellmer::type_string(
        "Absolute path to the FCS files directory."
      )
    ),
    ellmer::tool(
      .tool_detect_channels,
      "Read an FCS file and classify channels into scatter, time, and fluorochrome categories.",
      file_path = ellmer::type_string("Absolute path to an FCS file.")
    ),
    ellmer::tool(
      local_check_batch,
      "Summarise the sample distribution across batches in a metadata file.",
      meta_path = ellmer::type_string(
        "Absolute path to a CSV or XLSX metadata file."
      ),
      batch_col = ellmer::type_string(
        "Name of the column containing batch labels."
      )
    ),
    ellmer::tool(
      .tool_validate_markers,
      "Check that specified marker names exist in an FCS file.",
      markers = ellmer::type_array(
        ellmer::type_string("A marker name."),
        "Character vector of marker names to validate."
      ),
      fcs_file = ellmer::type_string("Absolute path to an FCS file.")
    )
  )
}

#' Start an MCP Server Exposing SpectraWeaveR Tools
#'
#' Launches a \href{https://modelcontextprotocol.io}{Model Context Protocol}
#' (MCP) server that exposes SpectraWeaveR's file-inspection tools and,
#' optionally, \pkg{btw} R-environment introspection tools.
#'
#' External AI applications — such as Claude Code, VS Code Copilot Chat,
#' or Cursor — can connect to this server to inspect FCS files, validate
#' sample metadata, detect channels, and check batch balance directly from
#' the user's R session.
#'
#' @section Setup for Claude Code:
#' \preformatted{
#' claude mcp add -s user spectraweaver -- \\
#'   Rscript -e "SpectraWeaveR::sw_mcp_server()"
#' }
#'
#' @section Setup for VS Code (\code{.vscode/mcp.json}):
#' \preformatted{
#' {
#'   "mcpServers": {
#'     "spectraweaver": {
#'       "command": "Rscript",
#'       "args": ["-e", "SpectraWeaveR::sw_mcp_server()"]
#'     }
#'   }
#' }
#' }
#'
#' @param include_btw Logical; if \code{TRUE} (default), also expose
#'   \pkg{btw} tools for R environment introspection (data-frame
#'   descriptions, package documentation, session info). Requires
#'   \pkg{btw}; silently skipped if not installed.
#' @param sensitivity Character; \code{"standard"} (default) auto-redacts
#'   PII in metadata previews; \code{"strict"} suppresses row previews and
#'   limits btw tools to docs and session only.
#'
#' @return This function does not return. It blocks the R process to serve
#'   MCP requests. Terminate with Ctrl-C or by stopping the process.
#'
#' @examples
#' \dontrun{
#' # Start MCP server (blocks the process)
#' sw_mcp_server()
#'
#' # Strict privacy for clinical data
#' sw_mcp_server(sensitivity = "strict")
#' }
#'
#' @export
sw_mcp_server <- function(include_btw = TRUE,
                          sensitivity = c("standard", "strict")) {
  if (!requireNamespace("mcptools", quietly = TRUE)) {
    stop(
      "Package 'mcptools' is required for MCP server functionality.\n",
      "Install it with: install.packages('mcptools')",
      call. = FALSE
    )
  }
  .check_ellmer()
  sensitivity <- match.arg(sensitivity)

  # Collect tools with sensitivity baked in via closures
  tools <- .sw_tool_list(sensitivity = sensitivity)

  if (isTRUE(include_btw) && requireNamespace("btw", quietly = TRUE)) {
    btw_groups <- if (sensitivity == "strict") {
      c("docs", "session")
    } else {
      c("env", "docs", "session", "files")
    }
    tryCatch({
      btw_tool_list <- btw::btw_tools(btw_groups)
      tools <- c(btw_tool_list, tools)
    }, error = function(e) {
      warning("Could not load btw tools: ", e$message, call. = FALSE)
    })
  }

  mcptools::mcp_server(tools = tools)
}
