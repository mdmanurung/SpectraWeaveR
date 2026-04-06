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
#' @export
sw_tibble_to_flowframe <- function(df, markers = NULL) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame or tibble.", call. = FALSE)
  }

  if (is.null(markers)) {
    # Use all numeric columns
    markers <- names(df)[vapply(df, is.numeric, logical(1))]
  }

  missing_markers <- setdiff(markers, names(df))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in data: ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }

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
