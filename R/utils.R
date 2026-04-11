#' @title Format Conversion Utilities for SpectraWeaveR
#'
#' @description
#' Bridge functions that handle conversions between flowFrame, flowSet, tibble,
#' and matrix objects required by different tools in the spectral flow cytometry
#' pipeline.
#'
#' @name utils
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Shared internal validation helpers
# ---------------------------------------------------------------------------

#' @noRd
.validate_markers <- function(markers) {
  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }
}

#' @noRd
.validate_df <- function(df, arg_name = "df") {
  if (!is.data.frame(df)) {
    stop("'", arg_name, "' must be a data.frame or tibble.", call. = FALSE)
  }
}

#' @noRd
.validate_markers_in_df <- function(df, markers, arg_name = "df") {
  missing_markers <- setdiff(markers, names(df))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in '", arg_name, "': ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }
}

#' Read FCS Files with Aurora-Safe Defaults
#'
#' Wrapper around \code{\link[flowCore]{read.FCS}} or
#' \code{\link[flowCore]{read.flowSet}} with
#' \code{truncate_max_range = FALSE} enforced. This is required for Cytek Aurora
#' data because fluorescence intensities (~4e6) exceed flowCore's default range
#' cutoffs.
#'
#' @param files Character vector of FCS file paths. A single path returns a
#'   \code{flowFrame}; multiple paths return a \code{flowSet}.
#' @param ... Additional arguments passed to \code{flowCore::read.FCS()} or
#'   \code{flowCore::read.flowSet()}.
#'
#' @return A \code{flowFrame} (single file) or \code{flowSet} (multiple files).
#'
#' @export
sw_read_fcs <- function(files, ...) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required. Install it from Bioconductor.",
         call. = FALSE)
  }

  if (!is.character(files) || length(files) == 0) {
    stop("'files' must be a non-empty character vector of FCS file paths.",
         call. = FALSE)
  }

  # Verify files exist
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("FCS file(s) not found: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  if (length(files) == 1) {
    flowCore::read.FCS(files, truncate_max_range = FALSE, ...)
  } else {
    flowCore::read.flowSet(files, truncate_max_range = FALSE, ...)
  }
}

#' Convert flowFrame to Tibble
#'
#' Extracts the expression matrix from a \code{flowFrame} and returns it as
#' a tibble with optional metadata columns for sample ID, batch, and condition.
#'
#' @param ff A \code{flowCore::flowFrame} object.
#' @param sample Character scalar identifying the sample (default: \code{NULL}).
#' @param batch Character or integer scalar identifying the batch
#'   (default: \code{NULL}).
#' @param condition Character scalar identifying the experimental condition
#'   (default: \code{NULL}).
#'
#' @return A \code{tibble} with one row per cell and one column per channel,
#'   plus metadata columns if specified.
#'
#' @examples
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   mat <- matrix(rnorm(30), ncol = 3,
#'                 dimnames = list(NULL, c("FSC-A", "SSC-A", "CD3")))
#'   ff <- flowCore::flowFrame(mat)
#'   df <- sw_flowframe_to_tibble(ff, sample = "S1", batch = "B1")
#'   head(df)
#' }
#'
#' @export
sw_flowframe_to_tibble <- function(ff, sample = NULL, batch = NULL,
                                   condition = NULL) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  mat <- flowCore::exprs(ff)
  df <- tibble::as_tibble(as.data.frame(mat))

  if (!is.null(sample)) df$sample <- sample
  if (!is.null(batch)) df$batch <- batch
  if (!is.null(condition)) df$condition <- condition

  df
}

#' Convert Tibble to flowFrame
#'
#' Converts a tibble (or data.frame) into a \code{flowCore::flowFrame}.
#' Only the specified marker columns are included in the expression matrix;
#' metadata columns are dropped.
#'
#' @param df A \code{data.frame} or \code{tibble} containing marker columns.
#' @param markers Character vector of column names to include as channels in the
#'   \code{flowFrame}. If \code{NULL}, all numeric columns are used.
#'
#' @return A \code{flowCore::flowFrame}.
#'
#' @examples
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   df <- data.frame(CD3 = rnorm(10), CD4 = rnorm(10), label = letters[1:10])
#'   ff <- sw_tibble_to_flowframe(df, markers = c("CD3", "CD4"))
#'   flowCore::colnames(ff)
#' }
#'
#' @export
sw_tibble_to_flowframe <- function(df, markers = NULL) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  .validate_df(df, "df")

  if (is.null(markers)) {
    # Use all numeric columns
    markers <- names(df)[vapply(df, is.numeric, logical(1))]
  }

  .validate_markers_in_df(df, markers, "df")

  mat <- as.matrix(df[, markers, drop = FALSE])
  storage.mode(mat) <- "double"
  flowCore::flowFrame(mat)
}

#' Convert Expression Matrix to Tibble
#'
#' Converts a numeric matrix to a tibble with optional metadata columns.
#' Useful for converting output from \code{flowCore::exprs()} or similar
#' matrix-returning functions.
#'
#' @param mat A numeric matrix with cells as rows and channels as columns.
#' @param colnames Optional character vector of column names. If \code{NULL},
#'   uses \code{colnames(mat)}.
#' @param sample Character scalar identifying the sample (default: \code{NULL}).
#' @param batch Character or integer scalar identifying the batch
#'   (default: \code{NULL}).
#' @param condition Character scalar identifying the experimental condition
#'   (default: \code{NULL}).
#'
#' @return A \code{tibble}.
#'
#' @examples
#' mat <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
#' sw_exprs_to_tibble(mat, sample = "S1", batch = "B1")
#'
#' @export
sw_exprs_to_tibble <- function(mat, colnames = NULL, sample = NULL,
                               batch = NULL, condition = NULL) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("'mat' must be a matrix or data.frame.", call. = FALSE)
  }

  if (!is.null(colnames)) {
    if (length(colnames) != ncol(mat)) {
      stop("Length of 'colnames' must match ncol(mat).", call. = FALSE)
    }
    base::colnames(mat) <- colnames
  }

  df <- tibble::as_tibble(as.data.frame(mat))

  if (!is.null(sample)) df$sample <- sample
  if (!is.null(batch)) df$batch <- batch
  if (!is.null(condition)) df$condition <- condition

  df
}

#' Get Fluorescent Channel Names
#'
#' Returns the names of fluorescent/marker channels from a \code{flowFrame} or
#' \code{flowSet}, excluding scatter (FSC/SSC) and Time channels.
#'
#' @param ff_or_fs A \code{flowFrame} or \code{flowSet} object.
#'
#' @return Character vector of fluorescent channel names.
#'
#' @export
sw_get_fluor_channels <- function(ff_or_fs) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (methods::is(ff_or_fs, "flowSet")) {
    ff <- ff_or_fs[[1]]
  } else if (methods::is(ff_or_fs, "flowFrame")) {
    ff <- ff_or_fs
  } else {
    stop("'ff_or_fs' must be a flowFrame or flowSet object.", call. = FALSE)
  }

  all_channels <- flowCore::colnames(ff)

  # Exclude scatter and time channels (case-insensitive)
  exclude_pattern <- "^(FSC|SSC|Time)(-[A-Za-z])?$"
  fluor <- all_channels[!grepl(exclude_pattern, all_channels,
                               ignore.case = TRUE)]

  fluor
}

#' Set Marker Names on a flowFrame
#'
#' Renames channel descriptions (markers) on a \code{flowFrame} according to
#' a named character vector mapping channel names to marker names.
#'
#' @param ff A \code{flowFrame} object.
#' @param marker_map A named character vector where names are channel names
#'   (e.g., \code{"BV421-A"}) and values are marker names (e.g., \code{"CD3"}).
#'
#' @return The modified \code{flowFrame} with updated marker names.
#'
#' @export
sw_set_marker_names <- function(ff, marker_map) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  if (!is.character(marker_map) || is.null(names(marker_map))) {
    stop("'marker_map' must be a named character vector.", call. = FALSE)
  }

  pdata <- flowCore::pData(flowCore::parameters(ff))
  channel_names <- as.character(pdata$name)

  for (ch in names(marker_map)) {
    idx <- which(channel_names == ch)
    if (length(idx) == 1) {
      pdata$desc[idx] <- marker_map[[ch]]
    } else {
      warning("Channel '", ch, "' not found in flowFrame; skipping.",
              call. = FALSE)
    }
  }

  flowCore::pData(flowCore::parameters(ff)) <- pdata
  ff
}

# ---------------------------------------------------------------------------
# Channel classification functions (ported from CytoPipeline)
# ---------------------------------------------------------------------------

#' Identify Signal Columns
#'
#' Returns a logical vector indicating which columns of a \code{flowFrame}
#' or \code{flowSet} represent true signal (not metadata columns like Time,
#' Original_ID, File, SampleID).
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param x A \code{flowCore::flowFrame} or \code{flowCore::flowSet}.
#' @param exclude_patterns Character vector of patterns to exclude
#'   (case-insensitive grep). Default: \code{c("Time", "Original_ID",
#'   "File", "SampleID")}.
#'
#' @return A named logical vector with length equal to \code{ncol(x)}.
#'
#' @examples
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   mat <- matrix(rnorm(40), ncol = 4,
#'                 dimnames = list(NULL, c("FSC-A", "SSC-A", "CD3", "Time")))
#'   ff <- flowCore::flowFrame(mat)
#'   sw_are_signal_cols(ff)
#' }
#'
#' @export
sw_are_signal_cols <- function(x,
                               exclude_patterns = c(
                                 "Time", "Original_ID",
                                 "File", "SampleID"
                               )) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(x, "flowFrame") && !methods::is(x, "flowSet")) {
    stop("'x' must be a flowFrame or flowSet object.", call. = FALSE)
  }

  vapply(flowCore::colnames(x),
    FUN.VALUE = logical(1),
    FUN = function(ch, patterns) {
      res <- TRUE
      for (pat in patterns) {
        res <- res & !grepl(pat, ch, ignore.case = TRUE)
      }
      res
    },
    patterns = exclude_patterns
  )
}

#' Identify Fluorochrome Columns
#'
#' Returns a logical vector indicating which columns of a \code{flowFrame}
#' or \code{flowSet} represent fluorochrome channels (excludes scatter
#' channels, time, and metadata in addition to the patterns excluded by
#' \code{\link{sw_are_signal_cols}}).
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param x A \code{flowCore::flowFrame} or \code{flowCore::flowSet}.
#' @param exclude_patterns Character vector of patterns to exclude
#'   (case-insensitive grep). Default: \code{c("FSC", "SSC", "Time",
#'   "Original_ID", "File", "SampleID")}.
#'
#' @return A named logical vector with length equal to \code{ncol(x)}.
#'
#' @export
sw_are_fluor_cols <- function(x,
                              exclude_patterns = c(
                                "FSC", "SSC",
                                "Time", "Original_ID",
                                "File", "SampleID"
                              )) {
  sw_are_signal_cols(x, exclude_patterns = exclude_patterns)
}

# ---------------------------------------------------------------------------
# Aggregate and sample (ported from CytoPipeline)
# ---------------------------------------------------------------------------

#' Aggregate and Subsample Multiple Flow Frames
#'
#' Pools events from multiple \code{flowFrame} objects in a \code{flowSet}
#' and subsamples to a target total number of events. Two strategies are
#' available: balanced (equal events per file) or forced total count.
#'
#' Three tracking columns are added: \code{File} (integer index),
#' \code{File_scattered} (noisy version for plotting), and
#' \code{Original_ID} (pre-subsample event index).
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param fs A \code{flowCore::flowSet} or \code{flowCore::flowFrame}
#'   (auto-wrapped).
#' @param n_total_events Integer; target total number of events.
#' @param setup Character; strategy for distributing events.
#'   \code{"forceNEvent"} (default) tries to reach exactly
#'   \code{n_total_events}; \code{"forceBalance"} uses the minimum per-file
#'   count, yielding balanced but potentially fewer events.
#' @param seed Integer; seed for reproducibility. Default: \code{NULL}.
#' @param channels Character vector of channels to keep. Default: \code{NULL}
#'   (all channels from the first file).
#' @param keep_order Logical; if \code{TRUE}, preserves chronological order
#'   within each file. Default: \code{FALSE}.
#'
#' @return A single \code{flowCore::flowFrame} with aggregated events.
#'
#' @export
sw_aggregate_and_sample <- function(fs,
                                    n_total_events,
                                    setup = c("forceNEvent", "forceBalance"),
                                    seed = NULL,
                                    channels = NULL,
                                    keep_order = FALSE) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.numeric(n_total_events) || length(n_total_events) != 1 ||
      n_total_events < 1) {
    stop("'n_total_events' must be a positive integer.", call. = FALSE)
  }

  if (methods::is(fs, "flowFrame")) {
    fs <- flowCore::flowSet(fs)
  }
  if (!methods::is(fs, "flowSet")) {
    stop("'fs' must be a flowSet or flowFrame object.", call. = FALSE)
  }

  setup <- match.arg(setup)

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
              add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  n_frames <- length(fs)
  ff_nrows <- vapply(seq_len(n_frames), function(i) {
    flowCore::nrow(fs[[i]])
  }, integer(1))

  min_nrow <- min(ff_nrows)
  sum_nrow <- sum(ff_nrows)
  c_frame <- ceiling(n_total_events / n_frames)

  # Decide how many rows to select from each frame
  if (c_frame <= min_nrow) {
    n_rows_select <- rep(c_frame, n_frames)
  } else if (setup == "forceBalance") {
    n_rows_select <- rep(min_nrow, n_frames)
  } else {
    # forceNEvent: fill from larger frames
    n_rows_select <- rep(0L, n_frames)
    trial_rows <- c_frame
    allocated <- 0
    still_to_allocate <- n_total_events

    while (still_to_allocate > 0 && allocated < sum_nrow) {
      for (i in seq_len(n_frames)) {
        n_rows_select[i] <- min(ff_nrows[i], trial_rows)
      }
      allocated <- sum(n_rows_select)
      still_to_allocate <- n_total_events - allocated
      n_spare <- n_frames - sum(n_rows_select == ff_nrows)
      if (n_spare > 0) {
        trial_rows <- trial_rows + ceiling(still_to_allocate / n_spare)
      } else {
        break
      }
    }
  }

  # Handle excess
  n_excess <- sum(n_rows_select) - n_total_events
  if (n_excess < 0) {
    warning("Could not reach ", n_total_events,
            " events; sampled ", sum(n_rows_select), " events.",
            call. = FALSE)
  }
  while (n_excess > 0) {
    n_withdraw <- min(n_excess, n_frames)
    where_withdraw <- sample.int(n_frames, size = n_withdraw, replace = FALSE)
    for (j in seq_along(where_withdraw)) {
      n_rows_select[where_withdraw[j]] <-
        n_rows_select[where_withdraw[j]] - 1L
    }
    n_excess <- n_excess - n_withdraw
  }

  # Build aggregated flowFrame — collect matrices first, then rbind once
  # to avoid O(n^2) memory allocation from growing rbind inside a loop.
  mat_list <- vector("list", n_frames)
  use_cols <- NULL

  for (i in seq_len(n_frames)) {
    current_ff <- fs[[i]]
    ids <- sample(seq_len(nrow(current_ff)), n_rows_select[i])
    if (keep_order) ids <- sort(ids)

    # Add tracking columns
    file_ids <- rep(i, n_rows_select[i])
    m <- cbind(
      File = file_ids,
      File_scattered = file_ids + stats::rnorm(length(file_ids), 0, 0.1),
      Original_ID = ids
    )
    current_ff <- flowCore::fr_append_cols(current_ff[ids, ], m)

    if (i == 1L) {
      if (!is.null(channels)) {
        # Resolve marker names to channel names if needed
        all_cols <- flowCore::colnames(current_ff)
        resolved <- intersect(channels, all_cols)
        if (length(resolved) == 0) {
          stop("No matching channels found.", call. = FALSE)
        }
        tracking <- c("File", "File_scattered", "Original_ID")
        use_cols <- c(resolved, tracking)
      } else {
        use_cols <- flowCore::colnames(current_ff)
      }
    } else {
      use_cols <- intersect(use_cols, flowCore::colnames(current_ff))
      if (length(use_cols) == 0) {
        stop("No common channels between flow frames.", call. = FALSE)
      }
    }

    mat_list[[i]] <- flowCore::exprs(current_ff)
  }

  # Subset to common columns and combine in one operation
  mat_list <- lapply(mat_list, function(m) m[, use_cols, drop = FALSE])
  combined_mat <- do.call(rbind, mat_list)

  # Build result flowFrame from the combined matrix
  result_ff <- flowCore::flowFrame(combined_mat)

  result_ff
}

# ---------------------------------------------------------------------------
# Collect events retained (audit trail)
# ---------------------------------------------------------------------------

#' Collect Events Retained at Each Pipeline Step
#'
#' Given a named list of intermediate results from pipeline execution (each
#' being a \code{flowFrame}, \code{flowSet}, or a list of \code{flowFrame}s),
#' returns a summary tibble showing how many events were retained at each step.
#'
#' This provides a full audit trail of event filtering across the pipeline,
#' inspired by CytoPipeline's \code{collectNbOfRetainedEvents()}.
#'
#' @param intermediates A named list of intermediate results. Each element
#'   can be:
#'   \itemize{
#'     \item A \code{flowCore::flowFrame}: counts its rows
#'     \item A \code{flowCore::flowSet}: counts total rows across all frames
#'     \item A list of \code{flowFrame}s: counts total rows across all frames
#'     \item A numeric scalar: used directly as the event count
#'     \item A tibble/data.frame: counts its rows
#'   }
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{\code{step}}{Character; step name}
#'     \item{\code{n_events}}{Integer; number of events at this step}
#'     \item{\code{pct_of_initial}}{Numeric; percentage of initial events
#'       retained}
#'     \item{\code{pct_of_previous}}{Numeric; percentage of previous step's
#'       events retained}
#'   }
#'
#' @examples
#' intermediates <- list(
#'   raw = 10000,
#'   after_qc = 8500,
#'   after_gating = data.frame(x = rnorm(6000))
#' )
#' sw_collect_events_retained(intermediates)
#'
#' @export
sw_collect_events_retained <- function(intermediates) {
  if (!is.list(intermediates) || length(intermediates) == 0) {
    stop("'intermediates' must be a non-empty named list.", call. = FALSE)
  }

  step_names <- names(intermediates)
  if (is.null(step_names) || any(step_names == "")) {
    stop("'intermediates' must be a named list.", call. = FALSE)
  }

  count_events <- function(obj) {
    if (is.numeric(obj) && length(obj) == 1) {
      return(as.integer(obj))
    }
    if (is.data.frame(obj)) {
      return(nrow(obj))
    }
    if (requireNamespace("flowCore", quietly = TRUE)) {
      if (methods::is(obj, "flowFrame")) {
        return(flowCore::nrow(obj))
      }
      if (methods::is(obj, "flowSet")) {
        return(sum(flowCore::fsApply(obj, flowCore::nrow)))
      }
    }
    if (is.list(obj)) {
      # List of flowFrames or similar
      counts <- vapply(obj, function(o) {
        if (requireNamespace("flowCore", quietly = TRUE) &&
            methods::is(o, "flowFrame")) {
          flowCore::nrow(o)
        } else if (is.data.frame(o)) {
          nrow(o)
        } else {
          NA_integer_
        }
      }, integer(1))

      if (all(is.na(counts))) {
        warning("Could not count events for step with list elements.",
                call. = FALSE)
        return(NA_integer_)
      }
      return(sum(counts, na.rm = TRUE))
    }

    warning("Unsupported type for event counting; returning NA.",
            call. = FALSE)
    NA_integer_
  }

  n_events <- vapply(intermediates, count_events, integer(1))

  initial <- n_events[1]
  pct_of_initial <- round(100 * n_events / initial, 2)

  pct_of_previous <- c(100, round(
    100 * n_events[-1] / n_events[-length(n_events)], 2
  ))

  tibble::tibble(
    step = step_names,
    n_events = n_events,
    pct_of_initial = pct_of_initial,
    pct_of_previous = pct_of_previous
  )
}
