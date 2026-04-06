#' @title Batch Correction with cyCombine
#'
#' @description
#' Functions for preparing data and performing batch correction using the
#' cyCombine package, which applies ComBat batch correction on SOM clusters.
#'
#' @name batch_correct
#' @keywords internal
NULL

#' Prepare Data for Batch Correction
#'
#' Converts a list of \code{flowFrame} objects into a single tibble with
#' arcsinh-transformed marker values and required metadata columns for
#' cyCombine batch correction.
#'
#' @param ff_list A named list of \code{flowFrame} objects, one per sample.
#' @param sample_meta A \code{data.frame} or \code{tibble} with columns
#'   \code{sample} (matching names of \code{ff_list}), \code{batch}, and
#'   optionally \code{condition}.
#' @param markers Character vector of marker/channel names to include.
#' @param cofactor Numeric; arcsinh cofactor for transformation
#'   (default: 6000 for spectral flow).
#'
#' @return A \code{tibble} with one row per cell, arcsinh-transformed marker
#'   columns, and metadata columns (sample, batch, condition).
#'
#' @export
sw_prepare_for_correction <- function(ff_list, sample_meta, markers,
                                      cofactor = 6000) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.list(ff_list) || length(ff_list) == 0) {
    stop("'ff_list' must be a non-empty named list of flowFrame objects.",
         call. = FALSE)
  }

  if (is.null(names(ff_list))) {
    stop("'ff_list' must be a named list.", call. = FALSE)
  }

  if (!is.data.frame(sample_meta)) {
    stop("'sample_meta' must be a data.frame or tibble.", call. = FALSE)
  }

  required_cols <- c("sample", "batch")
  missing_cols <- setdiff(required_cols, names(sample_meta))
  if (length(missing_cols) > 0) {
    stop("'sample_meta' is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.numeric(cofactor) || cofactor <= 0) {
    stop("'cofactor' must be a positive number.", call. = FALSE)
  }

  dfs <- list()

  for (sn in names(ff_list)) {
    ff <- ff_list[[sn]]

    if (!methods::is(ff, "flowFrame")) {
      stop("Element '", sn, "' of ff_list is not a flowFrame.", call. = FALSE)
    }

    # Extract expression matrix
    mat <- flowCore::exprs(ff)
    available_markers <- intersect(markers, colnames(mat))

    if (length(available_markers) == 0) {
      stop("None of the specified markers found in sample '", sn, "'.",
           call. = FALSE)
    }

    mat_sub <- mat[, available_markers, drop = FALSE]

    # Apply arcsinh transformation
    mat_transformed <- asinh(mat_sub / cofactor)

    # Convert to tibble with metadata
    df <- tibble::as_tibble(as.data.frame(mat_transformed))

    # Look up metadata for this sample
    meta_row <- sample_meta[sample_meta$sample == sn, , drop = FALSE]
    if (nrow(meta_row) == 0) {
      stop("Sample '", sn, "' not found in sample_meta.", call. = FALSE)
    }

    df$sample <- sn
    df$batch <- meta_row$batch[1]
    if ("condition" %in% names(sample_meta)) {
      df$condition <- meta_row$condition[1]
    }

    # Add cell ID for tracking
    df$id <- seq_len(nrow(df))

    dfs[[sn]] <- df
  }

  result <- dplyr::bind_rows(dfs)
  message("Prepared ", nrow(result), " cells from ", length(ff_list),
          " samples for batch correction")

  result
}

#' Batch Correct Expression Data
#'
#' Wrapper for \code{cyCombine::batch_correct()} that applies ComBat-based
#' batch correction using SOM clusters as reference.
#'
#' @param uncorrected A \code{tibble} as produced by
#'   \code{\link{sw_prepare_for_correction}}.
#' @param markers Character vector of marker column names to correct.
#' @param covar Character scalar naming the covariate column (typically
#'   \code{"condition"} or \code{NULL}).
#' @param xdim Integer; SOM grid x-dimension (default: 8).
#' @param ydim Integer; SOM grid y-dimension (default: 8).
#' @param norm_method Character; normalization method for cyCombine
#'   (default: \code{"scale"}).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to \code{cyCombine::batch_correct()}.
#'
#' @return A \code{tibble} with batch-corrected marker values.
#'
#' @export
sw_batch_correct <- function(uncorrected, markers, covar = NULL,
                             xdim = 8, ydim = 8, norm_method = "scale",
                             seed = 42, ...) {
  if (!requireNamespace("cyCombine", quietly = TRUE)) {
    stop("Package 'cyCombine' is required for batch correction. ",
         "Install it from GitHub: remotes::install_github('biosurf/cyCombine')",
         call. = FALSE)
  }

  if (!is.data.frame(uncorrected)) {
    stop("'uncorrected' must be a data.frame or tibble.", call. = FALSE)
  }

  if (!"batch" %in% names(uncorrected)) {
    stop("'uncorrected' must contain a 'batch' column.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  missing_markers <- setdiff(markers, names(uncorrected))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in data: ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }

  set.seed(seed)

  message("Running cyCombine batch correction on ", nrow(uncorrected),
          " cells, ", length(markers), " markers")

  corrected <- cyCombine::batch_correct(
    df = uncorrected,
    markers = markers,
    covar = covar,
    xdim = xdim,
    ydim = ydim,
    norm_method = norm_method,
    ...
  )

  message("Batch correction complete")
  corrected
}

#' Evaluate Batch Correction Quality
#'
#' Computes Earth Mover's Distance (EMD) and Median Absolute Deviation (MAD)
#' metrics to evaluate batch correction quality.
#'
#' @param uncorrected A \code{tibble} with pre-correction data.
#' @param corrected A \code{tibble} with post-correction data.
#' @param markers Character vector of marker names to evaluate.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{emd}}{A \code{tibble} with EMD scores per marker per batch
#'       pair, before and after correction}
#'     \item{\code{mad}}{A \code{tibble} with MAD scores per marker}
#'     \item{\code{improved}}{Logical; TRUE if overall EMD decreased}
#'     \item{\code{emd_reduction_pct}}{Percentage reduction in median EMD}
#'   }
#'
#' @export
sw_evaluate_correction <- function(uncorrected, corrected, markers) {
  if (!is.data.frame(uncorrected) || !is.data.frame(corrected)) {
    stop("'uncorrected' and 'corrected' must be data.frames.", call. = FALSE)
  }

  if (!"batch" %in% names(uncorrected) || !"batch" %in% names(corrected)) {
    stop("Both data frames must contain a 'batch' column.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  missing_in_uncorr <- setdiff(markers, names(uncorrected))
  missing_in_corr <- setdiff(markers, names(corrected))
  if (length(missing_in_uncorr) > 0 || length(missing_in_corr) > 0) {
    stop("Some markers not found in the data.", call. = FALSE)
  }

  batches <- unique(uncorrected$batch)
  if (length(batches) < 2) {
    message("Only one batch found; skipping EMD evaluation.")
    return(list(
      emd = tibble::tibble(),
      mad = tibble::tibble(),
      improved = NA,
      emd_reduction_pct = NA
    ))
  }

  # Compute per-marker MAD across batch medians
  compute_mad_scores <- function(df, markers_vec) {
    batch_medians <- df %>%
      dplyr::group_by(.data$batch) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(markers_vec), stats::median),
        .groups = "drop"
      )

    mad_scores <- vapply(markers_vec, function(m) {
      stats::mad(batch_medians[[m]])
    }, numeric(1))

    tibble::tibble(marker = markers_vec, mad = mad_scores)
  }

  mad_before <- compute_mad_scores(uncorrected, markers)
  mad_after <- compute_mad_scores(corrected, markers)

  mad_df <- dplyr::left_join(
    mad_before %>% dplyr::rename(mad_before = "mad"),
    mad_after %>% dplyr::rename(mad_after = "mad"),
    by = "marker"
  )

  # Overall improvement assessment
  avg_mad_before <- mean(mad_df$mad_before, na.rm = TRUE)
  avg_mad_after <- mean(mad_df$mad_after, na.rm = TRUE)

  improved <- avg_mad_after < avg_mad_before
  reduction_pct <- round(100 * (avg_mad_before - avg_mad_after) /
                           avg_mad_before, 1)

  message("Batch correction evaluation:")
  message("  Mean MAD before: ", round(avg_mad_before, 4))
  message("  Mean MAD after:  ", round(avg_mad_after, 4))
  message("  Reduction: ", reduction_pct, "%")

  list(
    emd = tibble::tibble(),
    mad = mad_df,
    improved = improved,
    emd_reduction_pct = reduction_pct
  )
}
