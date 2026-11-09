#' @noRd
.check_spectracle <- function() {
  if (!requireNamespace("spectracle", quietly = TRUE)) {
    stop(
      "Package 'spectracle' is required. Install it from GitHub: ",
      "pak::pak('nlaniewski/spectracle').",
      call. = FALSE
    )
  }
}

#' @noRd
.validate_af_choice <- function(af) {
  choices <- c("drop", "keep_as_row", "separate")
  if (missing(af) || !is.character(af) || length(af) != 1 || !af %in% choices) {
    stop(
      "'af' is required and must be one of: ",
      paste(sprintf("\"%s\"", choices), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(af)
}

#' Derive Spectra from Raw Single-Color Controls with spectracle
#'
#' Wrapper around \code{spectracle::spectracle()} followed by
#' \code{\link{sw_io_spectracle_to_spectra}}. spectracle (AGPL) reads Cytek
#' Aurora (SpectroFlo) controls named like \code{"CD3 BUV615 (Cells)"} and
#' \code{"Unstained (Cells)"}. It is not in Imports or Suggests of this
#' package and must be installed by the user.
#'
#' @param raw_controls Either one directory that holds the control FCS files,
#'   or two or more FCS file paths. A single \code{.fcs} path is refused,
#'   because spectracle cannot handle it.
#' @param af Required. How to treat the autofluorescence (AF) row. See
#'   \code{\link{sw_io_spectracle_to_spectra}}.
#' @param unstained_ids,filter_top_expressing,top_expressing_override Passed
#'   to \code{spectracle::spectracle()} as \code{unstained.ids},
#'   \code{filter.top.expressing} and \code{top.expressing.override}. All
#'   three change which events define a spectrum, so they default to
#'   \code{NULL} (spectracle's own default) and you set them on purpose.
#' @param rename,tissue_type Passed to
#'   \code{\link{sw_io_spectracle_to_spectra}}.
#'
#' @return A spectra matrix, or a list when \code{af = "separate"}. See
#'   \code{\link{sw_io_spectracle_to_spectra}}.
#'
#' @examples
#' \dontrun{
#' spectra <- sw_unmix_spectracle("controls/", af = "drop")
#' }
#'
#' @export
sw_unmix_spectracle <- function(
  raw_controls,
  af,
  unstained_ids = NULL,
  filter_top_expressing = NULL,
  top_expressing_override = NULL,
  rename = NULL,
  tissue_type = NULL
) {
  if (
    !is.character(raw_controls) ||
      length(raw_controls) == 0 ||
      anyNA(raw_controls)
  ) {
    stop(
      "'raw_controls' must be a directory or a character vector of FCS ",
      "file paths.",
      call. = FALSE
    )
  }
  if (length(raw_controls) == 1) {
    if (!dir.exists(raw_controls)) {
      stop(
        "'raw_controls' must be an existing directory when it is a ",
        "single path (a single FCS file is not supported): ",
        raw_controls,
        call. = FALSE
      )
    }
  } else {
    absent <- raw_controls[!file.exists(raw_controls)]
    if (length(absent) > 0) {
      stop(
        "Control file(s) not found: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    if (!all(grepl("\\.fcs$", raw_controls, ignore.case = TRUE))) {
      stop("All 'raw_controls' paths must end in .fcs.", call. = FALSE)
    }
  }
  .validate_af_choice(af)
  if (
    !is.null(unstained_ids) &&
      (!is.character(unstained_ids) || is.null(names(unstained_ids)))
  ) {
    stop(
      "'unstained_ids' must be NULL or a named character vector.",
      call. = FALSE
    )
  }
  if (
    !is.null(filter_top_expressing) &&
      !is.character(filter_top_expressing)
  ) {
    stop(
      "'filter_top_expressing' must be NULL or a character vector.",
      call. = FALSE
    )
  }
  if (
    !is.null(top_expressing_override) &&
      (!is.numeric(top_expressing_override) ||
        is.null(names(top_expressing_override)))
  ) {
    stop(
      "'top_expressing_override' must be NULL or a named numeric vector.",
      call. = FALSE
    )
  }
  .validate_rename(rename)

  .check_spectracle()

  sp <- spectracle::spectracle(
    raw.reference.controls = raw_controls,
    unstained.ids = unstained_ids,
    filter.top.expressing = filter_top_expressing,
    top.expressing.override = top_expressing_override,
    print.timings = FALSE
  )
  sw_io_spectracle_to_spectra(
    sp,
    af = af,
    rename = rename,
    tissue_type = tissue_type
  )
}

#' @noRd
.validate_rename <- function(rename) {
  if (is.null(rename)) {
    return(invisible(NULL))
  }
  if (
    !is.character(rename) ||
      is.null(names(rename)) ||
      anyNA(rename) ||
      anyNA(names(rename)) ||
      anyDuplicated(names(rename)) > 0 ||
      anyDuplicated(rename) > 0
  ) {
    stop(
      "'rename' must be NULL or a named character vector ",
      "c(old = \"new\") with unique names and unique values.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Convert a spectracle Table to a Spectra Matrix
#'
#' Pure conversion. It needs neither spectracle nor data.table.
#'
#' @section Input shape:
#' \code{spectracle::spectracle()} returns one row per control and tissue
#' type. The columns are metadata (\code{deriving.function}, \code{sample.id},
#' \code{tissue.type}, \code{N} = fluorophore, \code{S} = marker, cytometer
#' keywords, \code{detector} = peak detector, \code{laser}), then one numeric
#' column per detector, then \code{hash.md5}. A generalized autofluorescence
#' row with \code{N == "AF"} is appended at the end. The detector columns are
#' the only numeric columns, which is how they are found when
#' \code{detectors} is \code{NULL}.
#'
#' @section Detector names:
#' spectracle names Aurora detector columns from the flowstate alias (for
#' example \code{"UV1"}), where FCS channel names carry a suffix (for example
#' \code{"UV1-A"}). Names must match the unmixing data. Pass
#' \code{rename = c(UV1 = "UV1-A", ...)} to map them. The mapping is by name,
#' so column order is never changed, and every detector must appear in it.
#'
#' @section Scale:
#' spectracle scales each spectrum by its maximum and sets negative values of
#' non-AF rows to 0, giving values in \eqn{[0, 1]} with peak 1. The AF row is
#' scaled the same way but not clamped. This package's spectra are
#' L-infinity normalized: each row is divided by its maximum absolute value,
#' so the peak is 1. The rule applied here is
#' \code{row / max(abs(row))}. For a spectracle row this leaves the values
#' unchanged except for floating-point error. The function reports how many
#' rows it had to rescale by more than 1e-6.
#'
#' @param x The data.table (or data.frame) returned by
#'   \code{spectracle::spectracle()}.
#' @param af Required. What to do with the AF row(s) (\code{N == "AF"}).
#'   \describe{
#'     \item{\code{"drop"}}{Remove them. Returns the fluorophore matrix only.}
#'     \item{\code{"keep_as_row"}}{Keep them in the matrix, with row name
#'       \code{"AF"}. Several AF rows (Cells and Beads) get
#'       \code{"AF_<tissue.type>"}. The unmixing step then treats AF as one
#'       more fluorophore.}
#'     \item{\code{"separate"}}{Return \code{list(spectra, af)}. \code{af} is
#'       a one-row (or one row per tissue type) matrix with the same detector
#'       columns. AutoSpectral's \code{af_spectra} needs at least two
#'       rows for variant extraction, which one generalized AF row does not
#'       meet.}
#'   }
#' @param detectors Optional character vector of detector column names. Default
#'   \code{NULL} takes all numeric columns, in table order.
#' @param rename Optional named character vector \code{c(old = "new")} for
#'   detector names. See the section on detector names.
#' @param tissue_type Optional character vector. Keeps only rows whose
#'   \code{tissue.type} is in it (for example \code{"Cells"}). Needed when a
#'   fluorophore appears for both Cells and Beads, because row names must be
#'   unique.
#'
#' @return A numeric matrix with fluorophores in rows (names from \code{N}) and
#'   detectors in columns, in the table's order. With \code{af = "separate"},
#'   a list with elements \code{spectra} and \code{af}.
#'
#' @examples
#' sp <- data.frame(
#'   N = c("BV421", "FITC", "AF"),
#'   tissue.type = "Cells",
#'   V1 = c(1, 0, 0.2), B1 = c(0, 1, 0.5),
#'   hash.md5 = c("a", "b", "c")
#' )
#' sw_io_spectracle_to_spectra(sp, af = "drop")
#'
#' @export
sw_io_spectracle_to_spectra <- function(
  x,
  af,
  detectors = NULL,
  rename = NULL,
  tissue_type = NULL
) {
  if (!is.data.frame(x)) {
    stop(
      "'x' must be the data.table or data.frame from spectracle().",
      call. = FALSE
    )
  }
  if (!"N" %in% names(x)) {
    stop("'x' has no 'N' (fluorophore name) column.", call. = FALSE)
  }
  .validate_af_choice(af)
  .validate_rename(rename)
  if (
    !is.null(detectors) &&
      (!is.character(detectors) ||
        length(detectors) == 0 ||
        anyDuplicated(detectors) > 0)
  ) {
    stop("'detectors' must be NULL or unique column names.", call. = FALSE)
  }
  if (!is.null(tissue_type)) {
    if (!is.character(tissue_type) || length(tissue_type) == 0) {
      stop("'tissue_type' must be NULL or a character vector.", call. = FALSE)
    }
    if (!"tissue.type" %in% names(x)) {
      stop(
        "'tissue_type' was given but 'x' has no 'tissue.type' column.",
        call. = FALSE
      )
    }
  }

  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.null(tissue_type)) {
    keep <- as.character(df$tissue.type) %in% tissue_type
    if (!any(keep)) {
      stop(
        "No rows have tissue.type in: ",
        paste(tissue_type, collapse = ", "),
        call. = FALSE
      )
    }
    df <- df[keep, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("'x' has no rows.", call. = FALSE)
  }

  # Detector columns: explicit, or every numeric column in table order.
  if (is.null(detectors)) {
    detectors <- names(df)[vapply(df, is.numeric, logical(1))]
  }
  absent <- setdiff(detectors, names(df))
  if (length(absent) > 0) {
    stop(
      "Detector column(s) not found in 'x': ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(detectors) == 0) {
    stop("No numeric detector columns found in 'x'.", call. = FALSE)
  }
  non_numeric <- detectors[!vapply(df[detectors], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "Detector column(s) are not numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  mat <- as.matrix(df[, detectors, drop = FALSE])
  storage.mode(mat) <- "double"
  if (!all(is.finite(mat))) {
    stop("Detector values must be finite (no NA, NaN or Inf).", call. = FALSE)
  }
  tol <- 1e-6
  n_out <- sum(rowSums(mat < -tol | mat > 1 + tol) > 0)
  if (n_out > 0) {
    warning(
      n_out, " spectrum row(s) have values outside [0, 1] (spectracle scale). ",
      "Is 'x' a spectracle() result? Rows are still rescaled by their peak.",
      call. = FALSE
    )
  }

  # L-infinity rule: row / max(abs(row)).
  row_max <- apply(abs(mat), 1, max)
  if (any(row_max == 0)) {
    stop(
      "Spectrum with all-zero values: ",
      paste(df$N[row_max == 0], collapse = ", "),
      call. = FALSE
    )
  }
  n_rescaled <- sum(abs(row_max - 1) > tol)
  mat <- mat / row_max
  if (n_rescaled > 0) {
    message(
      n_rescaled,
      " of ",
      nrow(mat),
      " spectra had a peak other than 1 and were rescaled to peak 1."
    )
  }

  # Detector names: map by name, never by position.
  if (!is.null(rename)) {
    unmapped <- setdiff(detectors, names(rename))
    if (length(unmapped) > 0) {
      stop(
        "'rename' has no entry for detector(s): ",
        paste(unmapped, collapse = ", "),
        call. = FALSE
      )
    }
    colnames(mat) <- unname(rename[detectors])
  } else {
    colnames(mat) <- detectors
  }
  if (anyDuplicated(colnames(mat)) > 0) {
    stop("Detector names are not unique after renaming.", call. = FALSE)
  }
  stopifnot(
    "detector order changed" = identical(
      unname(colnames(mat)),
      unname(if (is.null(rename)) detectors else rename[detectors])
    )
  )

  # Row names from the fluorophore column.
  fluor <- as.character(df$N)
  if (anyNA(fluor) || any(!nzchar(fluor))) {
    stop("'N' (fluorophore name) has missing or empty values.", call. = FALSE)
  }
  is_af <- fluor == "AF"
  if (af != "drop" && !any(is_af)) {
    stop(
      "af = \"",
      af,
      "\" but 'x' has no row with N == \"AF\".",
      call. = FALSE
    )
  }
  if (sum(is_af) > 1) {
    if (!"tissue.type" %in% names(df)) {
      stop(
        "Several AF rows and no 'tissue.type' column to tell them apart.",
        call. = FALSE
      )
    }
    fluor[is_af] <- paste0("AF_", as.character(df$tissue.type[is_af]))
  }
  rownames(mat) <- fluor
  if (anyDuplicated(fluor) > 0) {
    stop(
      "Duplicate row names: ",
      paste(unique(fluor[duplicated(fluor)]), collapse = ", "),
      ". Use 'tissue_type' to keep one tissue type.",
      call. = FALSE
    )
  }

  switch(
    af,
    drop = {
      out <- mat[!is_af, , drop = FALSE]
      if (nrow(out) == 0) {
        stop("No fluorophore rows remain after dropping AF.", call. = FALSE)
      }
      out
    },
    keep_as_row = mat,
    separate = {
      out <- mat[!is_af, , drop = FALSE]
      if (nrow(out) == 0) {
        stop("No fluorophore rows remain after separating AF.", call. = FALSE)
      }
      list(spectra = out, af = mat[is_af, , drop = FALSE])
    }
  )
}
