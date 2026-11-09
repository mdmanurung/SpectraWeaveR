# R/qc_plots.R
# Thin wrappers around PeacoQC's two plot functions. PeacoQC is GPL (>= 3), so
# nothing is copied: both functions call PeacoQC::PlotPeacoQC() and
# PeacoQC::PeacoQCHeatmap() and only check arguments around them.

#' Plot where PeacoQC removed events in one file
#'
#' @description
#' Draw the event rate and the chosen channels over time, with the events that
#' PeacoQC removed marked. It wraps [PeacoQC::PlotPeacoQC()] and writes one
#' PNG file. The plot does not change the data.
#'
#' @details
#' `output_dir` has no default, so a plot never lands in the working folder by
#' accident. PeacoQC writes into `<output_dir>/PeacoQC_plots/` and names the
#' file `<prefix><FILENAME keyword without extension>.png`. It overwrites a
#' file with the same name; this function stops instead when that file exists.
#' Without a `FILENAME` keyword PeacoQC picks a random name, so no clash can
#' occur and no check is made.
#'
#' @param ff The `flowFrame` that was given to [sw_qc_run()].
#' @param channels Channel names or indices to draw. The same channels as the
#'   QC run, so the marks match what was tested.
#' @param qc A result from [sw_qc_run()], or the raw result of
#'   `PeacoQC::PeacoQC()`. It supplies which events were removed.
#' @param output_dir Folder for the PNG file. Created if missing.
#' @param ... Passed to [PeacoQC::PlotPeacoQC()], for example `display_cells`,
#'   `time_unit` or `prefix`. These change how the plot looks, not the QC.
#'
#' @return The path of the new PNG file, invisibly.
#'
#' @seealso [sw_qc_run()], [sw_qc_heatmap()]
#'
#' @examples
#' \dontrun{
#' res <- sw_qc_run(ff, channels = chans)
#' sw_qc_plot(ff, chans, res, output_dir = "qc_plots")
#' }
#'
#' @export
sw_qc_plot <- function(ff, channels, qc, output_dir, ...) {
  if (missing(output_dir)) {
    stop("`output_dir` is required and has no default.", call. = FALSE)
  }
  if (missing(ff) || missing(channels) || missing(qc)) {
    stop("`ff`, `channels` and `qc` are required.", call. = FALSE)
  }
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required for QC plots.", call. = FALSE)
  }
  if (!methods::is(ff, "flowFrame")) {
    stop("`ff` must be a flowFrame.", call. = FALSE)
  }
  if (
    !is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)
  ) {
    stop("`output_dir` must be one folder path.", call. = FALSE)
  }
  # Accept an sw_qc_run() result or the raw PeacoQC result.
  peaco <- if (is.list(qc) && !is.null(qc$PeacoQC)) qc$PeacoQC else qc
  if (!is.list(peaco) || is.null(peaco$GoodCells)) {
    stop(
      "`qc` must come from sw_qc_run() or PeacoQC::PeacoQC().",
      call. = FALSE
    )
  }
  # PeacoQC writes into <output_dir>/PeacoQC_plots/ and names the file from the
  # FILENAME keyword without its extension (a random name when it is missing).
  extra <- list(...)
  prefix <- if (is.null(extra$prefix)) "PeacoQC_" else extra$prefix
  plot_dir <- file.path(output_dir, "PeacoQC_plots")
  fname <- flowCore::keyword(ff)[["FILENAME"]]
  if (length(fname) == 1 && !is.na(fname)) {
    target <- file.path(
      plot_dir,
      paste0(prefix, tools::file_path_sans_ext(basename(fname)), ".png")
    )
    if (file.exists(target)) {
      stop("File exists, not overwritten: ", target, call. = FALSE)
    }
  }
  before <- list.files(plot_dir, pattern = "\\.png$", full.names = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  PeacoQC::PlotPeacoQC(
    ff = ff,
    channels = channels,
    output_directory = output_dir,
    display_peaks = peaco,
    ...
  )
  after <- list.files(plot_dir, pattern = "\\.png$", full.names = TRUE)
  invisible(setdiff(after, before))
}

#' Heatmap of PeacoQC results over many files
#'
#' @description
#' Draw one row per file with the share of events removed in total, by the
#' isolation tree step and by the MAD step, next to the settings used. It wraps
#' [PeacoQC::PeacoQCHeatmap()], which reads the report table that
#' [PeacoQC::PeacoQC()] writes into its output folder.
#'
#' @param report_location Path of the `PeacoQC_report.txt` file, or of the
#'   `output_dir` given to [sw_qc_run()] (the file sits in its
#'   `PeacoQC_results/` subfolder). Required.
#' @param ... Passed to [PeacoQC::PeacoQCHeatmap()], for example `show_values`,
#'   `show_row_names`, `latest_tests` or `title`.
#'
#' @return What [PeacoQC::PeacoQCHeatmap()] returns (a ComplexHeatmap object).
#'
#' @seealso [sw_qc_batch()], [sw_qc_plot()]
#'
#' @examples
#' \dontrun{
#' sw_qc_heatmap("qc_out/PeacoQC_report.txt")
#' }
#'
#' @export
sw_qc_heatmap <- function(report_location, ...) {
  if (missing(report_location)) {
    stop("`report_location` is required.", call. = FALSE)
  }
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required for QC plots.", call. = FALSE)
  }
  if (
    !is.character(report_location) ||
      length(report_location) != 1 ||
      is.na(report_location)
  ) {
    stop("`report_location` must be one path.", call. = FALSE)
  }
  path <- report_location
  if (dir.exists(path)) {
    # PeacoQC writes <output_dir>/PeacoQC_results/PeacoQC_report.txt.
    candidates <- file.path(
      path,
      c("PeacoQC_results/PeacoQC_report.txt", "PeacoQC_report.txt")
    )
    found <- candidates[file.exists(candidates)]
    path <- if (length(found) > 0) found[[1]] else candidates[[1]]
  }
  if (!file.exists(path)) {
    stop("No PeacoQC report found at: ", path, call. = FALSE)
  }
  PeacoQC::PeacoQCHeatmap(report_location = path, ...)
}
