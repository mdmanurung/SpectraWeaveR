#' @title Signal Quality Control with PeacoQC
#'
#' @description
#' Functions for signal stability quality control using the PeacoQC package.
#' PeacoQC uses Isolation Trees and Median Absolute Deviation (MAD) to
#' identify and remove aberrant signal events from flow cytometry data.
#'
#' @name qc
#' @keywords internal
NULL

#' Signal Quality Control on a Single flowFrame
#'
#' Wrapper for \code{\link[PeacoQC]{PeacoQC}} with sensible defaults for
#' spectral flow cytometry data. Identifies and removes events with aberrant
#' signal behavior using Isolation Trees and MAD-based methods.
#'
#' @param ff A \code{flowFrame} object (should be transformed before QC).
#' @param channels Character vector of channel names to evaluate. If
#'   \code{NULL}, all fluorescent channels (excluding FSC/SSC/Time) are used.
#' @param IT_limit Numeric; Isolation Tree contamination threshold
#'   (default: 0.55). Lower values are more aggressive in removing events.
#' @param MAD Numeric; Median Absolute Deviation threshold (default: 6).
#'   Lower values are more aggressive.
#' @param output_dir Character path for QC report output. If \code{NULL},
#'   no plots are generated.
#' @param ... Additional arguments passed to \code{PeacoQC::PeacoQC()}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{FinalFF}}{The cleaned \code{flowFrame} with aberrant events
#'       removed}
#'     \item{\code{GoodCells}}{Logical vector indicating retained (TRUE) vs
#'       removed (FALSE) events}
#'     \item{\code{n_removed}}{Number of events removed}
#'     \item{\code{pct_removed}}{Percentage of events removed}
#'   }
#'
#' @export
sw_signal_qc <- function(ff, channels = NULL, IT_limit = 0.55, MAD = 6,
                         output_dir = NULL, ...) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required for signal QC. ",
         "Install it from Bioconductor.", call. = FALSE)
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  if (!is.numeric(IT_limit) || IT_limit <= 0 || IT_limit >= 1) {
    stop("'IT_limit' must be a number between 0 and 1 (exclusive).",
         call. = FALSE)
  }

  if (!is.numeric(MAD) || MAD <= 0) {
    stop("'MAD' must be a positive number.", call. = FALSE)
  }

  # Determine channels to evaluate
  if (is.null(channels)) {
    channels <- sw_get_fluor_channels(ff)
  }

  n_before <- nrow(ff)

  # Run PeacoQC
  qc_args <- list(
    ff = ff,
    channels = channels,
    determine_good_cells = "all",
    IT_limit = IT_limit,
    MAD = MAD,
    plot = !is.null(output_dir),
    save_fcs = FALSE
  )

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    qc_args$output_directory <- output_dir
  }

  qc_args <- c(qc_args, list(...))

  qc_result <- do.call(PeacoQC::PeacoQC, qc_args)

  # Extract results
  good_cells <- qc_result$GoodCells
  ff_clean <- ff[good_cells, ]

  n_after <- nrow(ff_clean)
  n_removed <- n_before - n_after
  pct_removed <- round(100 * n_removed / n_before, 1)

  message("PeacoQC: removed ", n_removed, " / ", n_before,
          " events (", pct_removed, "%)")

  list(
    FinalFF = ff_clean,
    GoodCells = good_cells,
    n_removed = n_removed,
    pct_removed = pct_removed
  )
}

#' Batch Signal Quality Control
#'
#' Applies \code{\link{sw_signal_qc}} to a list of \code{flowFrame} objects
#' and returns cleaned frames plus summary statistics.
#'
#' @param ff_list A named list of \code{flowFrame} objects.
#' @param ... Additional arguments passed to \code{\link{sw_signal_qc}}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{cleaned}}{Named list of cleaned \code{flowFrame} objects}
#'     \item{\code{summary}}{A \code{tibble} with per-sample QC statistics}
#'   }
#'
#' @export
sw_signal_qc_batch <- function(ff_list, ...) {
  if (!is.list(ff_list) || length(ff_list) == 0) {
    stop("'ff_list' must be a non-empty list of flowFrame objects.",
         call. = FALSE)
  }

  sample_names <- names(ff_list)
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_along(ff_list))
    names(ff_list) <- sample_names
  }

  cleaned <- list()
  stats <- list()

  for (sn in sample_names) {
    message("Running QC on: ", sn)
    qc_result <- sw_signal_qc(ff_list[[sn]], ...)

    cleaned[[sn]] <- qc_result$FinalFF

    stats[[sn]] <- tibble::tibble(
      sample = sn,
      n_before = nrow(ff_list[[sn]]),
      n_after = nrow(qc_result$FinalFF),
      n_removed = qc_result$n_removed,
      pct_removed = qc_result$pct_removed
    )
  }

  summary_df <- dplyr::bind_rows(stats)

  list(
    cleaned = cleaned,
    summary = summary_df
  )
}

#' QC Summary Report
#'
#' Tabulates cells removed per sample and flags samples with high removal
#' rates (>30% by default).
#'
#' @param qc_results The output from \code{\link{sw_signal_qc_batch}}.
#' @param threshold Numeric; percentage threshold for flagging samples
#'   (default: 30). Samples with removal rate exceeding this value are flagged.
#'
#' @return A \code{tibble} with columns: sample, n_before, n_after, n_removed,
#'   pct_removed, and flagged (logical).
#'
#' @export
sw_qc_summary <- function(qc_results, threshold = 30) {
  if (!is.list(qc_results) || !"summary" %in% names(qc_results)) {
    stop("'qc_results' must be the output from sw_signal_qc_batch().",
         call. = FALSE)
  }

  if (!is.numeric(threshold) || threshold < 0 || threshold > 100) {
    stop("'threshold' must be a number between 0 and 100.", call. = FALSE)
  }

  summary_df <- qc_results$summary
  summary_df$flagged <- summary_df$pct_removed > threshold

  n_flagged <- sum(summary_df$flagged)
  if (n_flagged > 0) {
    warning(n_flagged, " sample(s) exceeded the ", threshold,
            "% removal threshold: ",
            paste(summary_df$sample[summary_df$flagged], collapse = ", "),
            call. = FALSE)
  }

  summary_df
}
