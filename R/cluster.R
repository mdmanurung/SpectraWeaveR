#' @title Clustering with FlowSOM
#'
#' @description
#' Functions for cell clustering, visualization, and projection using the
#' FlowSOM package.
#'
#' @name cluster
#' @keywords internal
NULL

#' Cluster Cells Using FlowSOM
#'
#' Converts input data to a \code{flowFrame}, runs FlowSOM self-organizing
#' map clustering, and applies consensus metaclustering.
#'
#' @param corrected A \code{tibble}, \code{data.frame}, or numeric matrix
#'   containing the (typically batch-corrected) expression data.
#' @param lineage_markers Character vector of marker names to use for
#'   clustering. These should be lineage markers (e.g., CD3, CD4, CD8), not
#'   functional markers.
#' @param xdim Integer; SOM grid x-dimension (default: 10).
#' @param ydim Integer; SOM grid y-dimension (default: 10).
#' @param n_metaclusters Integer; number of metaclusters to identify
#'   (default: 20).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to \code{FlowSOM::FlowSOM()}.
#'
#' @return A FlowSOM object containing the trained SOM and metacluster
#'   assignments.
#'
#' @export
sw_cluster <- function(corrected, lineage_markers, xdim = 10, ydim = 10,
                       n_metaclusters = 20, seed = 42, ...) {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("Package 'FlowSOM' is required for clustering. ",
         "Install it from Bioconductor.", call. = FALSE)
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.character(lineage_markers) || length(lineage_markers) == 0) {
    stop("'lineage_markers' must be a non-empty character vector.",
         call. = FALSE)
  }

  if (!is.numeric(n_metaclusters) || n_metaclusters < 2) {
    stop("'n_metaclusters' must be an integer >= 2.", call. = FALSE)
  }

  # Convert input to flowFrame
  if (is.matrix(corrected)) {
    missing_m <- setdiff(lineage_markers, colnames(corrected))
    if (length(missing_m) > 0) {
      stop("Markers not found in matrix columns: ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
    mat <- corrected[, lineage_markers, drop = FALSE]
    storage.mode(mat) <- "double"
    ff <- flowCore::flowFrame(mat)
  } else if (is.data.frame(corrected)) {
    missing_m <- setdiff(lineage_markers, names(corrected))
    if (length(missing_m) > 0) {
      stop("Markers not found in data columns: ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
    mat <- as.matrix(corrected[, lineage_markers, drop = FALSE])
    storage.mode(mat) <- "double"
    ff <- flowCore::flowFrame(mat)
  } else {
    stop("'corrected' must be a matrix, data.frame, or tibble.", call. = FALSE)
  }

  set.seed(seed)

  message("Running FlowSOM on ", nrow(ff), " cells with ",
          length(lineage_markers), " markers (", xdim, "x", ydim,
          " grid, ", n_metaclusters, " metaclusters)")

  fsom <- FlowSOM::FlowSOM(
    ff,
    colsToUse = lineage_markers,
    xdim = xdim,
    ydim = ydim,
    nClus = n_metaclusters,
    seed = seed,
    ...
  )

  message("FlowSOM clustering complete")
  fsom
}

#' Get Metacluster Assignments
#'
#' Extracts per-cell metacluster assignments from a FlowSOM result.
#'
#' @param fsom_result A FlowSOM object (output from \code{\link{sw_cluster}}).
#'
#' @return An integer vector of metacluster assignments, one per cell.
#'
#' @export
sw_get_cluster_assignments <- function(fsom_result) {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("Package 'FlowSOM' is required.", call. = FALSE)
  }

  if (!methods::is(fsom_result, "FlowSOM")) {
    stop("'fsom_result' must be a FlowSOM object.", call. = FALSE)
  }

  FlowSOM::GetMetaclusters(fsom_result)
}

#' Compute Median Fluorescence Intensities per Metacluster
#'
#' Calculates the median expression of each marker within each metacluster.
#'
#' @param fsom_result A FlowSOM object (output from \code{\link{sw_cluster}}).
#'
#' @return A \code{tibble} with one row per metacluster and one column per
#'   marker, plus a \code{metacluster} column.
#'
#' @export
sw_cluster_mfis <- function(fsom_result) {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("Package 'FlowSOM' is required.", call. = FALSE)
  }

  if (!methods::is(fsom_result, "FlowSOM")) {
    stop("'fsom_result' must be a FlowSOM object.", call. = FALSE)
  }

  mfi_mat <- FlowSOM::GetClusterMFIs(fsom_result, prettyColnames = FALSE)

  df <- tibble::as_tibble(as.data.frame(mfi_mat))
  df$metacluster <- seq_len(nrow(mfi_mat))

  # Reorder: metacluster first
  df <- df[, c("metacluster", setdiff(names(df), "metacluster"))]
  df
}

#' Plot FlowSOM Clusters
#'
#' Generates FlowSOM summary and star plots and saves them to a file.
#'
#' @param fsom_result A FlowSOM object (output from \code{\link{sw_cluster}}).
#' @param plot_file Character path for the output PDF file
#'   (default: \code{"FlowSOM_clusters.pdf"}).
#'
#' @return The file path (invisibly).
#'
#' @export
sw_plot_clusters <- function(fsom_result, plot_file = "FlowSOM_clusters.pdf") {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("Package 'FlowSOM' is required.", call. = FALSE)
  }

  if (!methods::is(fsom_result, "FlowSOM")) {
    stop("'fsom_result' must be a FlowSOM object.", call. = FALSE)
  }

  if (!is.character(plot_file) || length(plot_file) != 1) {
    stop("'plot_file' must be a single file path.", call. = FALSE)
  }

  # Ensure output directory exists
  out_dir <- dirname(plot_file)
  if (!dir.exists(out_dir) && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE)
  }

  grDevices::pdf(plot_file, width = 12, height = 8)
  tryCatch({
    FlowSOM::FlowSOMmary(fsom_result)
    bg_values <- FlowSOM::GetMetaclusters(fsom_result)
    FlowSOM::PlotStars(fsom_result, backgroundValues = bg_values)
  }, finally = {
    grDevices::dev.off()
  })

  message("FlowSOM plots saved to: ", plot_file)
  invisible(plot_file)
}

#' Map New Data onto Trained FlowSOM
#'
#' Projects new samples onto an existing trained FlowSOM model using
#' \code{FlowSOM::NewData()}.
#'
#' @param fsom_result A FlowSOM object (the trained model).
#' @param new_ff A \code{flowFrame} with new data to project.
#'
#' @return A FlowSOM object with the new data mapped onto the existing SOM.
#'
#' @export
sw_map_new_data <- function(fsom_result, new_ff) {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("Package 'FlowSOM' is required.", call. = FALSE)
  }

  if (!methods::is(fsom_result, "FlowSOM")) {
    stop("'fsom_result' must be a FlowSOM object.", call. = FALSE)
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(new_ff, "flowFrame")) {
    stop("'new_ff' must be a flowFrame object.", call. = FALSE)
  }

  message("Mapping ", nrow(new_ff), " new cells onto trained FlowSOM")

  fsom_new <- FlowSOM::NewData(fsom_result, new_ff)

  message("Mapping complete")
  fsom_new
}
