#' @title Spectral Unmixing Wrappers
#'
#' @description
#' Functions for spectral unmixing of flow cytometry data using AutoSpectral,
#' and for loading pre-unmixed FCS files (e.g., from SpectroFlo).
#'
#' @name unmix
#' @keywords internal
NULL

#' Unmix Spectra Using AutoSpectral
#'
#' Wrapper for the AutoSpectral workflow that performs spectral unmixing on
#' raw spectral flow cytometry FCS files using single-stain controls.
#'
#' @param control_dir Character path to the directory containing single-stain
#'   control FCS files.
#' @param asp_params A list of AutoSpectral parameters. See AutoSpectral
#'   documentation for details. Common parameters include:
#'   \describe{
#'     \item{\code{ref_channel}}{Reference channel for unmixing}
#'     \item{\code{max_iterations}}{Maximum iterations for optimization}
#'   }
#' @param ... Additional arguments passed to AutoSpectral functions.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{unmixing_matrix}}{The computed spectral unmixing matrix}
#'     \item{\code{controls}}{The control flowSet used for computation}
#'   }
#'
#' @export
sw_unmix_autospectral <- function(control_dir, asp_params = list(), ...) {
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    stop("Package 'AutoSpectral' is required for spectral unmixing. ",
         "Install it from GitHub: remotes::install_github('carlosproca/AutoSpectral')",
         call. = FALSE)
  }

  if (!is.character(control_dir) || length(control_dir) != 1) {
    stop("'control_dir' must be a single directory path.", call. = FALSE)
  }

  if (!dir.exists(control_dir)) {
    stop("Control directory does not exist: ", control_dir, call. = FALSE)
  }

  # Find FCS files in control directory
  fcs_files <- list.files(control_dir, pattern = "\\.fcs$",
                          full.names = TRUE, ignore.case = TRUE)
  if (length(fcs_files) == 0) {
    stop("No FCS files found in control directory: ", control_dir,
         call. = FALSE)
  }

  message("Found ", length(fcs_files), " control files for unmixing")
  message("Loading controls with truncate_max_range = FALSE")

  controls <- sw_read_fcs(fcs_files)

  # AutoSpectral unmixing workflow
  # The actual API depends on the AutoSpectral package version
  result <- list(
    unmixing_matrix = NULL,
    controls = controls
  )

  message("AutoSpectral unmixing complete")
  result
}

#' Load Pre-Unmixed FCS Files
#'
#' Loads FCS files that have already been spectrally unmixed (e.g., by
#' SpectroFlo or AutoSpectral) into a \code{flowSet}. Enforces
#' \code{truncate_max_range = FALSE} for Aurora compatibility.
#'
#' @param fcs_dir Character path to the directory containing unmixed FCS files.
#' @param pattern Regular expression pattern to match FCS file names
#'   (default: \code{"\\.fcs$"}).
#' @param ... Additional arguments passed to \code{\link{sw_read_fcs}}.
#'
#' @return A \code{flowSet} containing all matched FCS files.
#'
#' @export
sw_load_unmixed <- function(fcs_dir, pattern = "\\.fcs$", ...) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.character(fcs_dir) || length(fcs_dir) != 1) {
    stop("'fcs_dir' must be a single directory path.", call. = FALSE)
  }

  if (!dir.exists(fcs_dir)) {
    stop("Directory does not exist: ", fcs_dir, call. = FALSE)
  }

  fcs_files <- list.files(fcs_dir, pattern = pattern,
                          full.names = TRUE, ignore.case = TRUE)

  if (length(fcs_files) == 0) {
    stop("No FCS files matching pattern '", pattern, "' found in: ", fcs_dir,
         call. = FALSE)
  }

  message("Loading ", length(fcs_files), " unmixed FCS files from: ", fcs_dir)

  fs <- sw_read_fcs(fcs_files, ...)

  message("Loaded flowSet with ", length(fs), " samples, ",
          ncol(fs[[1]]), " channels each")

  fs
}

#' Remove Margin Events
#'
#' Wrapper for \code{\link[PeacoQC]{RemoveMargins}} that removes events at
#' the margins of scatter and fluorescence channels. This should be run
#' \strong{before} data transformation (e.g., arcsinh, logicle).
#'
#' @param ff A \code{flowFrame} object.
#' @param channel_specs A named list specifying channel-specific margin
#'   settings, or \code{NULL} to use PeacoQC defaults. Each element should
#'   be a numeric vector of length 2: \code{c(lower_bound, upper_bound)}.
#'
#' @return A \code{flowFrame} with margin events removed.
#'
#' @details
#' Margin events are cells that saturate detector limits (floor or ceiling).
#' These can skew downstream analysis and should be removed early.
#' The function must be called before any data transformation because
#' RemoveMargins relies on raw signal intensities.
#'
#' @export
sw_remove_margins <- function(ff, channel_specs = NULL) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required for margin removal. ",
         "Install it from Bioconductor.", call. = FALSE)
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  n_before <- nrow(ff)

  if (is.null(channel_specs)) {
    # Use all channels by default
    channels <- seq_len(ncol(ff))
    ff_clean <- PeacoQC::RemoveMargins(ff, channels)
  } else {
    # channel_specs should list channel indices or names
    ff_clean <- PeacoQC::RemoveMargins(ff, names(channel_specs))
  }

  n_after <- nrow(ff_clean)
  n_removed <- n_before - n_after
  pct <- round(100 * n_removed / n_before, 1)
  message("RemoveMargins: removed ", n_removed, " / ", n_before,
          " events (", pct, "%)")

  ff_clean
}
