#' @title Batch Correction with cyCombine
#'
#' @description
#' Functions for preparing data and performing batch correction using the
#' cyCombine package, which applies ComBat batch correction on SOM clusters.
#' Supports both an all-in-one workflow (\code{\link{sw_batch_correct}}) and a
#' modular workflow (\code{\link{sw_normalize}} \eqn{\to}
#' \code{\link{sw_create_som}} \eqn{\to} \code{\link{sw_correct_data}}).
#'
#' Diagnostic and evaluation helpers are provided for batch-effect detection
#' (\code{\link{sw_detect_batch_effect}}), quality metrics
#' (\code{\link{sw_compute_emd}}, \code{\link{sw_evaluate_emd}},
#' \code{\link{sw_evaluate_mad}}), and visualization
#' (\code{\link{sw_plot_batch_densities}}, \code{\link{sw_plot_batch_dimred}}).
#'
#' @name batch_correct
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Internal helper: check cyCombine availability
# ---------------------------------------------------------------------------
.check_cycombine <- function() {
  if (!requireNamespace("cyCombine", quietly = TRUE)) {
    stop("Package 'cyCombine' is required for batch correction. ",
         "Install it from GitHub: remotes::install_github('biosurf/cyCombine')",
         call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Internal helper: common input validation for batch-correction functions
# ---------------------------------------------------------------------------
.validate_batch_df <- function(df, markers = NULL, arg_name = "df") {

  if (!is.data.frame(df)) {
    stop("'", arg_name, "' must be a data.frame or tibble.", call. = FALSE)
  }

  if (!"batch" %in% names(df)) {
    stop("'", arg_name, "' must contain a 'batch' column.", call. = FALSE)
  }

  if (!is.null(markers)) {
    if (!is.character(markers) || length(markers) == 0) {
      stop("'markers' must be a non-empty character vector.", call. = FALSE)
    }
    missing_markers <- setdiff(markers, names(df))
    if (length(missing_markers) > 0) {
      stop("Marker(s) not found in data: ",
           paste(missing_markers, collapse = ", "), call. = FALSE)
    }
  }
}

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
  if (!requireNamespace("flowCore", quietly = TRUE)) { # nocov
    stop("Package 'flowCore' is required.", call. = FALSE) # nocov
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
#' All-in-one wrapper for \code{cyCombine::batch_correct()} that normalises,
#' clusters with a Self-Organising Map, and applies ComBat-based batch
#' correction.  For finer control use the modular functions
#' \code{\link{sw_normalize}}, \code{\link{sw_create_som}}, and
#' \code{\link{sw_correct_data}}.
#'
#' @param uncorrected A \code{tibble} as produced by
#'   \code{\link{sw_prepare_for_correction}}.
#' @param markers Character vector of marker column names to correct.
#' @param covar Character scalar naming the covariate column (typically
#'   \code{"condition"} or \code{NULL}).
#' @param label Optional pre-computed cluster labels (integer vector, same
#'   length as \code{nrow(uncorrected)}).  When supplied the internal SOM step
#'   is skipped.
#' @param xdim Integer; SOM grid x-dimension (default: 8).
#' @param ydim Integer; SOM grid y-dimension (default: 8).
#' @param rlen Integer; SOM training length (default: 10).
#' @param norm_method Character; normalization method for cyCombine
#'   (default: \code{"scale"}).  One of \code{"scale"}, \code{"rank"},
#'   \code{"none"}.
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to \code{cyCombine::batch_correct()}.
#'
#' @return A \code{tibble} with batch-corrected marker values.
#'
#' @export
sw_batch_correct <- function(uncorrected, markers, covar = NULL,
                             label = NULL,
                             xdim = 8, ydim = 8, rlen = 10,
                             norm_method = "scale",
                             seed = 42, ...) {
  .check_cycombine()
  .validate_batch_df(uncorrected, markers, arg_name = "uncorrected")

  set.seed(seed)

  message("Running cyCombine batch correction on ", nrow(uncorrected),
          " cells, ", length(markers), " markers")

  corrected <- cyCombine::batch_correct(
    df = uncorrected,
    label = label,
    markers = markers,
    covar = covar,
    xdim = xdim,
    ydim = ydim,
    rlen = rlen,
    norm_method = norm_method,
    ...
  )

  message("Batch correction complete")
  corrected
}

#' Evaluate Batch Correction Quality
#'
#' Quick evaluation of batch correction quality using per-marker MAD of batch
#' medians.  For richer diagnostics see \code{\link{sw_evaluate_emd}} and
#' \code{\link{sw_evaluate_mad}}, which delegate to cyCombine's cluster-aware
#' EMD and MAD computations.
#'
#' @param uncorrected A \code{tibble} with pre-correction data.
#' @param corrected A \code{tibble} with post-correction data.
#' @param markers Character vector of marker names to evaluate.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{emd}}{A \code{tibble} of per-marker, per-cluster EMD
#'       scores (before and after correction).  Populated when both data
#'       frames contain a \code{label} column and \pkg{cyCombine} is
#'       available; otherwise an empty tibble.}
#'     \item{\code{mad}}{A \code{tibble} with MAD scores per marker}
#'     \item{\code{improved}}{Logical; TRUE if overall MAD decreased}
#'     \item{\code{emd_reduction_pct}}{Percentage reduction in mean MAD}
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

  # Compute EMD if cluster labels are available and cyCombine is installed
  emd_result <- tibble::tibble()
  if ("label" %in% names(uncorrected) && "label" %in% names(corrected) &&
      requireNamespace("cyCombine", quietly = TRUE)) {
    tryCatch({
      emd_before <- cyCombine::compute_emd(
        uncorrected, markers = markers,
        cell_col = "label", batch_col = "batch"
      )
      emd_after <- cyCombine::compute_emd(
        corrected, markers = markers,
        cell_col = "label", batch_col = "batch"
      )
      join_cols <- setdiff(
        intersect(names(emd_before), names(emd_after)), "emd"
      )
      emd_result <- dplyr::left_join(
        emd_before %>% dplyr::rename(emd_before = "emd"),
        emd_after %>% dplyr::rename(emd_after = "emd"),
        by = join_cols
      )
      emd_median_before <- stats::median(emd_before$emd, na.rm = TRUE)
      emd_median_after <- stats::median(emd_after$emd, na.rm = TRUE)
      message("  EMD median before: ", round(emd_median_before, 4))
      message("  EMD median after:  ", round(emd_median_after, 4))
    }, error = function(e) {
      message("  EMD computation skipped: ", conditionMessage(e))
    })
  }

  list(
    emd = emd_result,
    mad = mad_df,
    improved = improved,
    emd_reduction_pct = reduction_pct
  )
}

# ===========================================================================
# Modular Workflow
# ===========================================================================

#' Normalise Data for Batch Correction
#'
#' Applies per-batch normalisation to marker columns, preparing the data for
#' SOM clustering.
#'
#' @param df A \code{tibble} with marker columns and a \code{batch} column,
#'   typically produced by \code{\link{sw_prepare_for_correction}}.
#' @param markers Character vector of marker column names to normalise.
#' @param norm_method Normalisation method.  One of \code{"scale"} (z-score,
#'   default for single-study data), \code{"rank"} (rank-based, recommended
#'   for multi-study merges), or \code{"none"}.
#' @param ... Additional arguments passed to \code{cyCombine::normalize()}.
#'
#' @return A \code{tibble} with normalised marker values, suitable for
#'   \code{\link{sw_create_som}}.
#'
#' @seealso \code{\link{sw_create_som}}, \code{\link{sw_correct_data}},
#'   \code{\link{sw_batch_correct}}
#'
#' @export
sw_normalize <- function(df, markers, norm_method = "scale", ...) {
  .check_cycombine()
  .validate_batch_df(df, markers)

  norm_method <- match.arg(norm_method,
                           c("scale", "rank", "none"))

  message("Normalising ", length(markers), " markers (method: ",
          norm_method, ")")

  cyCombine::normalize(df, markers = markers,
                       norm_method = norm_method, ...)
}


#' Create Self-Organising Map for Batch Correction
#'
#' Clusters cells using a Self-Organising Map (SOM).  The resulting labels are
#' passed to \code{\link{sw_correct_data}} so that ComBat is applied within
#' each cluster.
#'
#' @param df A \code{tibble} of normalised marker values, typically produced
#'   by \code{\link{sw_normalize}}.
#' @param markers Character vector of marker column names used for clustering.
#' @param xdim Integer; SOM grid x-dimension (default: 8).
#' @param ydim Integer; SOM grid y-dimension (default: 8).
#' @param rlen Integer; SOM training iterations (default: 10).
#' @param seed Integer; random seed (default: 42).
#' @param ... Additional arguments passed to \code{cyCombine::create_som()}.
#'
#' @return An integer vector of cluster labels, one per row in \code{df}.
#'
#' @seealso \code{\link{sw_normalize}}, \code{\link{sw_correct_data}},
#'   \code{\link{sw_batch_correct}}
#'
#' @export
sw_create_som <- function(df, markers, xdim = 8, ydim = 8, rlen = 10,
                          seed = 42, ...) {
  .check_cycombine()

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  missing_markers <- setdiff(markers, names(df))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in data: ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }

  set.seed(seed)

  message("Creating SOM (", xdim, "x", ydim, " grid, rlen=", rlen, ")")

  cyCombine::create_som(df, markers = markers,
                        xdim = xdim, ydim = ydim,
                        rlen = rlen, seed = seed, ...)
}


#' Apply ComBat Batch Correction with Pre-Computed Labels
#'
#' Runs ComBat batch correction within each cluster defined by the supplied
#' labels.  This is the final step of the modular workflow
#' (\code{\link{sw_normalize}} \eqn{\to} \code{\link{sw_create_som}} \eqn{\to}
#' \code{sw_correct_data}).
#'
#' @param df A \code{tibble} with the **original** (un-normalised) marker
#'   values.  Normalisation is only used for clustering; the correction itself
#'   operates on the original scale.
#' @param label Integer vector of cluster labels (from
#'   \code{\link{sw_create_som}} or any external clustering).
#' @param markers Character vector of marker column names to correct.
#' @param covar Character scalar naming the biological covariate column to
#'   preserve (e.g. \code{"condition"}), or \code{NULL}.
#' @param parametric Logical; use parametric ComBat (default: \code{TRUE}).
#' @param ... Additional arguments passed to \code{cyCombine::correct_data()}.
#'
#' @return A \code{tibble} with batch-corrected marker values.
#'
#' @seealso \code{\link{sw_normalize}}, \code{\link{sw_create_som}},
#'   \code{\link{sw_batch_correct}}
#'
#' @export
sw_correct_data <- function(df, label, markers, covar = NULL,
                            parametric = TRUE, ...) {
  .check_cycombine()
  .validate_batch_df(df, markers)

  if (!is.numeric(label) && !is.integer(label)) {
    stop("'label' must be a numeric/integer vector of cluster assignments.",
         call. = FALSE)
  }

  if (length(label) != nrow(df)) {
    stop("Length of 'label' (", length(label),
         ") must equal nrow(df) (", nrow(df), ").", call. = FALSE)
  }

  message("Applying ComBat correction across ", length(unique(label)),
          " clusters, ", length(markers), " markers")

  cyCombine::correct_data(
    df, label = label, markers = markers,
    covar = covar, parametric = parametric, ...
  )
}


# ===========================================================================
# Batch-Effect Detection
# ===========================================================================

#' Detect Batch Effects
#'
#' Produces diagnostic density plots, per-marker EMD scores, and an MDS plot
#' of median marker expression across batches.  Delegates to
#' \code{cyCombine::detect_batch_effect_express()}.
#'
#' @param df A \code{tibble} with marker columns, a \code{batch} column, and
#'   optionally a \code{sample} column.
#' @param markers Character vector of marker column names to evaluate.
#' @param out_dir Optional output directory.
#'   If \code{NULL} (default), plots are returned as a list of \pkg{ggplot2}
#'   objects.  If a path is given, plots are saved there as PNG files.
#' @param downsample Optional integer; subsample to this many cells before
#'   computing.  Useful for large datasets.
#' @param seed Integer; random seed (default: 42).
#' @param ... Additional arguments passed to
#'   \code{cyCombine::detect_batch_effect_express()}.
#'
#' @return A list of \code{ggplot} objects (when \code{out_dir = NULL}), or
#'   \code{NULL} invisibly after saving to \code{out_dir}.
#'
#' @seealso \code{\link{sw_evaluate_emd}}, \code{\link{sw_evaluate_mad}},
#'   \code{\link{sw_plot_batch_densities}}
#'
#' @export
sw_detect_batch_effect <- function(df, markers, out_dir = NULL,
                                   downsample = NULL, seed = 42, ...) {
  .check_cycombine()
  .validate_batch_df(df, markers)

  set.seed(seed)

  message("Detecting batch effects across ",
          length(unique(df$batch)), " batches, ",
          length(markers), " markers")

  cyCombine::detect_batch_effect_express(
    df, markers = markers, out_dir = out_dir,
    downsample = downsample, seed = seed, ...
  )
}


# ===========================================================================
# Evaluation — EMD
# ===========================================================================

#' Compute Earth Mover's Distance
#'
#' Computes Earth Mover's Distance (EMD) between batch distributions for each
#' marker within each SOM cluster.  The data must contain a \code{label}
#' column (cluster assignments); use \code{\link{sw_create_som}} to create
#' one if needed.
#'
#' @param df A \code{tibble} with marker columns, \code{batch}, and
#'   \code{label} columns.
#' @param markers Character vector of marker column names.
#' @param cell_col Name of the cluster label column (default: \code{"label"}).
#' @param batch_col Name of the batch column (default: \code{"batch"}).
#' @param binSize Numeric; bin size for EMD discretisation (default: 0.1).
#' @param ... Additional arguments passed to \code{cyCombine::compute_emd()}.
#'
#' @return A named list of EMD matrices (per cluster, per marker), as returned
#'   by \code{cyCombine::compute_emd()}.
#'
#' @seealso \code{\link{sw_evaluate_emd}}, \code{\link{sw_evaluate_mad}}
#'
#' @export
sw_compute_emd <- function(df, markers, cell_col = "label",
                           batch_col = "batch", binSize = 0.1, ...) {
  .check_cycombine()

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  if (!cell_col %in% names(df)) {
    stop("Column '", cell_col, "' not found in data. ",
         "Run sw_create_som() first to add cluster labels.", call. = FALSE)
  }

  if (!batch_col %in% names(df)) {
    stop("Column '", batch_col, "' not found in data.", call. = FALSE)
  }

  missing_markers <- setdiff(markers, names(df))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in data: ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }

  message("Computing EMD for ", length(markers), " markers across ",
          length(unique(df[[batch_col]])), " batches")

  cyCombine::compute_emd(
    df, markers = markers, cell_col = cell_col,
    batch_col = batch_col, binSize = binSize, ...
  )
}


#' Evaluate Batch Correction with EMD
#'
#' Compares Earth Mover's Distance distributions before and after batch
#' correction.  Both data frames must have a \code{label} column with
#' matching cluster assignments.
#'
#' @param uncorrected A \code{tibble} with pre-correction data (must include
#'   \code{batch} and \code{label} columns).
#' @param corrected A \code{tibble} with post-correction data (must include
#'   \code{batch} and \code{label} columns).
#' @param markers Character vector of marker column names to evaluate.
#' @param cell_col Name of the cluster label column (default: \code{"label"}).
#' @param batch_col Name of the batch column (default: \code{"batch"}).
#' @param binSize Numeric; bin size for EMD discretisation (default: 0.1).
#' @param ... Additional arguments passed to \code{cyCombine::evaluate_emd()}.
#'
#' @return A named list as returned by \code{cyCombine::evaluate_emd()},
#'   typically containing:
#'   \describe{
#'     \item{\code{violin}}{A \pkg{ggplot2} violin plot comparing EMDs}
#'     \item{\code{scatter}}{A \pkg{ggplot2} scatter plot of per-marker EMDs}
#'     \item{\code{reduction}}{Percentage reduction in median EMD}
#'     \item{\code{emd}}{A \code{tibble} of per-marker, per-cluster EMD
#'       values}
#'   }
#'
#' @seealso \code{\link{sw_compute_emd}}, \code{\link{sw_evaluate_mad}},
#'   \code{\link{sw_evaluate_correction}}
#'
#' @export
sw_evaluate_emd <- function(uncorrected, corrected, markers,
                            cell_col = "label", batch_col = "batch",
                            binSize = 0.1, ...) {
  .check_cycombine()

  for (nm in c("uncorrected", "corrected")) {
    d <- get(nm)
    if (!is.data.frame(d)) {
      stop("'", nm, "' must be a data.frame or tibble.", call. = FALSE)
    }
    for (col in c(batch_col, cell_col)) {
      if (!col %in% names(d)) {
        stop("Column '", col, "' not found in '", nm, "'.", call. = FALSE)
      }
    }
    missing_m <- setdiff(markers, names(d))
    if (length(missing_m) > 0) {
      stop("Marker(s) not found in '", nm, "': ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  message("Evaluating correction quality (EMD) for ",
          length(markers), " markers")

  cyCombine::evaluate_emd(
    uncorrected, corrected, markers = markers,
    cell_col = cell_col, batch_col = batch_col,
    binSize = binSize, ...
  )
}


# ===========================================================================
# Evaluation — MAD (cyCombine-backed)
# ===========================================================================

#' Evaluate Batch Correction with MAD
#'
#' Computes Median Absolute Deviation (MAD) of batch medians within SOM
#' clusters, before and after correction.  Both data frames must have a
#' \code{label} column with matching cluster assignments.
#'
#' @param uncorrected A \code{tibble} with pre-correction data (must include
#'   \code{batch} and \code{label} columns).
#' @param corrected A \code{tibble} with post-correction data (must include
#'   \code{batch} and \code{label} columns).
#' @param markers Character vector of marker column names to evaluate.
#' @param cell_col Name of the cluster label column (default: \code{"label"}).
#' @param batch_col Name of the batch column (default: \code{"batch"}).
#' @param ... Additional arguments passed to \code{cyCombine::evaluate_mad()}.
#'
#' @return A named list as returned by \code{cyCombine::evaluate_mad()},
#'   typically containing:
#'   \describe{
#'     \item{\code{score}}{Overall MAD score}
#'     \item{\code{mad}}{A \code{tibble} with per-marker MAD values}
#'   }
#'
#' @seealso \code{\link{sw_evaluate_emd}}, \code{\link{sw_evaluate_correction}}
#'
#' @export
sw_evaluate_mad <- function(uncorrected, corrected, markers,
                            cell_col = "label", batch_col = "batch", ...) {
  .check_cycombine()

  for (nm in c("uncorrected", "corrected")) {
    d <- get(nm)
    if (!is.data.frame(d)) {
      stop("'", nm, "' must be a data.frame or tibble.", call. = FALSE)
    }
    for (col in c(batch_col, cell_col)) {
      if (!col %in% names(d)) {
        stop("Column '", col, "' not found in '", nm, "'.", call. = FALSE)
      }
    }
    missing_m <- setdiff(markers, names(d))
    if (length(missing_m) > 0) {
      stop("Marker(s) not found in '", nm, "': ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  message("Evaluating correction quality (MAD) for ",
          length(markers), " markers")

  cyCombine::evaluate_mad(
    uncorrected, corrected, markers = markers,
    cell_col = cell_col, batch_col = batch_col, ...
  )
}


# ===========================================================================
# Visualisation
# ===========================================================================

#' Plot Marker Density Distributions by Batch
#'
#' Creates density ridgeline plots comparing marker distributions across
#' batches, before and after correction.  Delegates to
#' \code{cyCombine::plot_density()}.
#'
#' @param uncorrected A \code{tibble} with pre-correction data.
#' @param corrected A \code{tibble} with post-correction data.
#' @param markers Character vector of marker column names to plot.
#' @param filename Optional file path for saving the plot (e.g.
#'   \code{"density.pdf"}).  If \code{NULL} (default), the plot is returned
#'   but not saved.
#' @param ncol Integer; number of marker columns per page (default: 6).
#' @param ... Additional arguments passed to \code{cyCombine::plot_density()}.
#'
#' @return A \code{ggplot} object (or list of \code{ggplot} objects),
#'   invisibly when \code{filename} is given.
#'
#' @seealso \code{\link{sw_plot_batch_dimred}},
#'   \code{\link{sw_detect_batch_effect}}
#'
#' @export
sw_plot_batch_densities <- function(uncorrected, corrected, markers,
                                    filename = NULL, ncol = 6, ...) {
  .check_cycombine()

  if (!is.data.frame(uncorrected) || !is.data.frame(corrected)) {
    stop("'uncorrected' and 'corrected' must be data.frames.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  for (nm in c("uncorrected", "corrected")) {
    missing_m <- setdiff(markers, names(get(nm)))
    if (length(missing_m) > 0) {
      stop("Marker(s) not found in '", nm, "': ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
  }

  message("Creating density plots for ", length(markers), " markers")

  cyCombine::plot_density(
    uncorrected, corrected,
    markers = markers,
    filename = filename,
    ncol = ncol, ...
  )
}


#' Plot Dimensionality Reduction Coloured by Batch
#'
#' Creates a UMAP or PCA projection of the data, coloured by batch, to
#' visualise batch effects or the success of correction.  Delegates to
#' \code{cyCombine::plot_dimred()}.
#'
#' @param df A \code{tibble} with marker columns and a \code{batch} column.
#' @param markers Character vector of marker column names to include in the
#'   dimensionality reduction.
#' @param type Character; type of dimensionality reduction.  One of
#'   \code{"umap"} (default) or \code{"pca"}.
#' @param name Character; descriptive label for the plot title (default:
#'   \code{"Batch"}).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to \code{cyCombine::plot_dimred()}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{sw_plot_batch_densities}},
#'   \code{\link{sw_detect_batch_effect}}
#'
#' @export
sw_plot_batch_dimred <- function(df, markers, type = "umap",
                                 name = "Batch", seed = 42, ...) {
  .check_cycombine()

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  missing_markers <- setdiff(markers, names(df))
  if (length(missing_markers) > 0) {
    stop("Marker(s) not found in data: ",
         paste(missing_markers, collapse = ", "), call. = FALSE)
  }

  if (!"batch" %in% names(df)) {
    stop("'df' must contain a 'batch' column.", call. = FALSE)
  }

  type <- match.arg(type, c("umap", "pca"))

  set.seed(seed)

  message("Computing ", toupper(type), " projection for ",
          nrow(df), " cells")

  cyCombine::plot_dimred(
    df, markers = markers, type = type,
    name = name, seed = seed, ...
  )
}
