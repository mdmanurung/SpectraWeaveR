#' @noRd
.check_flowstate <- function() {
  if (!requireNamespace("flowstate", quietly = TRUE)) {
    stop(
      "Package 'flowstate' is required. Install it from GitHub: ",
      "pak::pak('nlaniewski/flowstate').",
      call. = FALSE
    )
  }
}

# A flowstate is an S3 list of data.tables: data, parameters, keywords and
# (optional) spill. Only data.frame behaviour is used here, so the checks
# work on plain data.frames too.
#' @noRd
.validate_flowstate <- function(fs, arg_name = "fs") {
  if (!inherits(fs, "flowstate")) {
    stop(
      "'",
      arg_name,
      "' must be a flowstate object ",
      "(see sw_io_read_flowstate()).",
      call. = FALSE
    )
  }
  missing_parts <- setdiff(c("data", "parameters", "keywords"), names(fs))
  if (length(missing_parts) > 0) {
    stop(
      "'",
      arg_name,
      "' is missing element(s): ",
      paste(missing_parts, collapse = ", "),
      call. = FALSE
    )
  }
  not_df <- !vapply(
    fs[c("data", "parameters", "keywords")],
    is.data.frame,
    logical(1)
  )
  if (any(not_df)) {
    stop(
      "'",
      arg_name,
      "' elements must be data.tables or data.frames: ",
      paste(names(not_df)[not_df], collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(fs$data) == 0) {
    stop("'", arg_name, "' has no events in [['data']].", call. = FALSE)
  }
  if (!"sample.id" %in% names(fs$data)) {
    stop(
      "'",
      arg_name,
      "' [['data']] has no 'sample.id' column.",
      call. = FALSE
    )
  }
  invisible(fs)
}

#' Read FCS Files as flowstate Objects
#'
#' Wrapper around \code{flowstate::read.flowstate()}. The flowstate package
#' reads the FCS DATA segment as 32-bit floats and applies no range cutoff.
#' Aurora intensities (~4e6) are therefore kept as recorded. The parameter
#' range (\code{$PnR}) is stored in \code{fs$parameters$R} and is used only
#' when you call \code{flowstate::select_nonsaturating()}.
#'
#' @param files Character vector of FCS file paths. One path returns a
#'   \code{flowstate}. Several paths return a named list of flowstates, or one
#'   combined flowstate when \code{concatenate = TRUE}.
#' @param colnames_type Required. \code{"N"} names the data columns by the
#'   \code{$PnN} keyword (channel name). \code{"S"} names them by \code{$PnS}
#'   (marker or stain name). The choice sets every downstream column name, so
#'   it has no default.
#' @param concatenate Logical. If \code{TRUE}, combine several files into one
#'   flowstate. Default: \code{FALSE}.
#' @param ... Further arguments passed to \code{flowstate::read.flowstate()},
#'   for example \code{sample.id}.
#'
#' @return A \code{flowstate} object, or a named list of them. Each
#'   flowstate is a list of data.tables: \code{data}, \code{parameters},
#'   \code{keywords} and (if the file has one) \code{spill}.
#'
#' @details
#' flowstate (AGPL) is called through its public API only. It is in neither
#' Imports nor Suggests of this package and must be installed by the user.
#' Mixed \code{$CYT} values across files are refused by flowstate.
#'
#' @examples
#' \dontrun{
#' fs <- sw_io_read_flowstate(
#'   list.files("fcs", pattern = "fcs$", full.names = TRUE),
#'   colnames_type = "S",
#'   concatenate = TRUE
#' )
#' }
#'
#' @export
sw_io_read_flowstate <- function(
  files,
  colnames_type,
  concatenate = FALSE,
  ...
) {
  if (!is.character(files) || length(files) == 0) {
    stop(
      "'files' must be a non-empty character vector of FCS file paths.",
      call. = FALSE
    )
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(
      "FCS file(s) not found: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    missing(colnames_type) ||
      !is.character(colnames_type) ||
      length(colnames_type) != 1 ||
      !colnames_type %in% c("N", "S")
  ) {
    stop(
      "'colnames_type' is required and must be \"N\" or \"S\".",
      call. = FALSE
    )
  }
  if (
    !is.logical(concatenate) || length(concatenate) != 1 || is.na(concatenate)
  ) {
    stop("'concatenate' must be TRUE or FALSE.", call. = FALSE)
  }
  if (concatenate && length(files) < 2) {
    stop("'concatenate = TRUE' needs at least two files.", call. = FALSE)
  }

  .check_flowstate()

  flowstate::read.flowstate(
    files,
    colnames.type = colnames_type,
    concatenate = concatenate,
    ...
  )
}

#' Apply an arcsinh Transform to a flowstate
#'
#' Wrapper around \code{flowstate::flowstate.transform()} with
#' \code{transform.func = "asinh"}. The cofactor is required. flowstate
#' defaults to 5000 and this package's design note says 6000 for spectral
#' flow. Both change every transformed value, so the choice is yours.
#'
#' @param fs A \code{flowstate} object.
#' @param cofactor Required positive number. Values are divided by it before
#'   \code{asinh()}.
#' @param j Optional character vector of data columns to transform. Default
#'   \code{NULL} lets flowstate pick the fluorescence and scatter columns it
#'   treats as transformable.
#'
#' @return \code{fs}, invisibly. \strong{The object is changed in place}
#'   (data.table reference semantics). Each transformed column gets a
#'   \code{transformed} attribute that stores the function and cofactor.
#'
#' @examples
#' \dontrun{
#' sw_io_flowstate_transform(fs, cofactor = 6000)
#' }
#'
#' @export
sw_io_flowstate_transform <- function(fs, cofactor, j = NULL) {
  .validate_flowstate(fs, "fs")
  if (
    missing(cofactor) ||
      !is.numeric(cofactor) ||
      length(cofactor) != 1 ||
      is.na(cofactor) ||
      cofactor <= 0
  ) {
    stop(
      "'cofactor' is required and must be a single positive number.",
      call. = FALSE
    )
  }
  if (!is.null(j) && !is.character(j)) {
    stop(
      "'j' must be NULL or a character vector of column names.",
      call. = FALSE
    )
  }
  if (!is.null(j)) {
    absent <- setdiff(j, names(fs$data))
    if (length(absent) > 0) {
      stop(
        "Column(s) not found in fs$data: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
  }

  .check_flowstate()

  flowstate::flowstate.transform(
    fs,
    j = j,
    transform.func = "asinh",
    cofactor = cofactor
  )
  invisible(fs)
}

#' Convert a flowstate to a Tibble
#'
#' Pure conversion. No flowstate or data.table call is made, so it works on
#' any object with the flowstate shape.
#'
#' @param fs A \code{flowstate} object.
#'
#' @return A tibble with one row per event and one column per data column.
#'   The \code{sample.id} factor becomes a character column named
#'   \code{sample}, which matches \code{\link{sw_io_ff_to_tibble}}. Column
#'   attributes (\code{N}, \code{alias}, \code{transformed}) are dropped. If
#'   any column was transformed, the tibble carries the attribute
#'   \code{sw_transform}: a named list with the function and cofactor per
#'   column. No such attribute means the values are linear.
#'
#' @examples
#' \dontrun{
#' tbl <- sw_io_flowstate_to_tibble(fs)
#' }
#'
#' @export
sw_io_flowstate_to_tibble <- function(fs) {
  .validate_flowstate(fs, "fs")

  df <- as.data.frame(fs$data)

  # Record transforms before the attributes are dropped.
  transformed <- lapply(df, attr, which = "transformed")
  transformed <- transformed[!vapply(transformed, is.null, logical(1))]

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) {
      as.character(x)
    } else {
      attributes(x) <- NULL
      x
    }
  })
  names(df)[names(df) == "sample.id"] <- "sample"

  out <- tibble::as_tibble(df)
  if (length(transformed) > 0) {
    attr(out, "sw_transform") <- transformed
  }
  out
}

#' Convert a flowstate to a flowFrame or flowSet
#'
#' Builds a \code{flowCore::flowFrame} from the numeric data columns. Column
#' names are the flowstate names, so they follow the \code{colnames_type} used
#' at read time. flowCore derives each parameter range from the data, so no
#' value is truncated. Aurora intensities (~4e6) pass through unchanged.
#'
#' @param fs A \code{flowstate} object.
#'
#' @return A \code{flowFrame} if \code{fs} holds one sample. A named
#'   \code{flowSet} (one frame per \code{sample.id}) if it holds several.
#'   Only numeric columns are kept. Keywords and spillover are not carried
#'   over. Values stay as stored: transformed columns remain transformed.
#'
#' @examples
#' \dontrun{
#' ff <- sw_io_flowstate_to_ff(fs)
#' }
#'
#' @export
sw_io_flowstate_to_ff <- function(fs) {
  .validate_flowstate(fs, "fs")
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  df <- as.data.frame(fs$data)
  sample_id <- as.character(df$sample.id)
  numeric_cols <- setdiff(
    names(df)[vapply(df, is.numeric, logical(1))],
    "sample.id"
  )
  if (length(numeric_cols) == 0) {
    stop("fs$data has no numeric columns to convert.", call. = FALSE)
  }

  to_frame <- function(rows) {
    mat <- as.matrix(df[rows, numeric_cols, drop = FALSE])
    storage.mode(mat) <- "double"
    dimnames(mat) <- list(NULL, numeric_cols)
    flowCore::flowFrame(mat)
  }

  samples <- unique(sample_id)
  if (length(samples) == 1) {
    return(to_frame(seq_len(nrow(df))))
  }
  frames <- lapply(samples, \(s) to_frame(which(sample_id == s)))
  names(frames) <- samples
  flowCore::flowSet(frames)
}

#' Write a flowstate as FCS 3.1 Files
#'
#' Wrapper around \code{flowstate::write.flowstate()}. It writes one file per
#' \code{sample.id}. The file name is the original \code{$FIL} plus
#' \code{_flowstateMOD} and the optional suffix.
#'
#' @param fs A \code{flowstate} object.
#' @param file_dir Required. An existing directory for the output files.
#'   flowstate would default to \code{tempdir()}, which loses the files at
#'   the end of the session.
#' @param fil_suffix Optional single string appended to the file name.
#' @param endianness \code{"little"} (default in flowstate) or \code{"big"}.
#'
#' @return \code{file_dir}, invisibly.
#'
#' @details
#' \strong{\code{fs} is changed in place.} flowstate undoes any transform
#' stored on the data columns and turns factor columns into numeric derived
#' parameters before it writes. Pass a copy if you need the original
#' (\code{data.table::copy()} on each element). Per the flowstate
#' documentation, files open in flowCore and flowstate. Proprietary
#' software such as FlowJo may not read them.
#'
#' @examples
#' \dontrun{
#' sw_io_write_flowstate(fs, file_dir = "out/fcs", fil_suffix = "clean")
#' }
#'
#' @export
sw_io_write_flowstate <- function(
  fs,
  file_dir,
  fil_suffix = NULL,
  endianness = c("little", "big")
) {
  .validate_flowstate(fs, "fs")
  if (
    missing(file_dir) ||
      !is.character(file_dir) ||
      length(file_dir) != 1 ||
      is.na(file_dir)
  ) {
    stop(
      "'file_dir' is required and must be a single directory path.",
      call. = FALSE
    )
  }
  if (!dir.exists(file_dir)) {
    stop("'file_dir' does not exist: ", file_dir, call. = FALSE)
  }
  if (
    !is.null(fil_suffix) &&
      (!is.character(fil_suffix) ||
        length(fil_suffix) != 1 ||
        is.na(fil_suffix))
  ) {
    stop("'fil_suffix' must be NULL or a single string.", call. = FALSE)
  }
  endianness <- match.arg(endianness)

  .check_flowstate()

  flowstate::write.flowstate(
    fs,
    file.dir = file_dir,
    fil.suffix = fil_suffix,
    endianness = endianness
  )
  invisible(file_dir)
}
