#' @title Clustering with kohonen SOM or FastPG
#'
#' @description
#' Functions for cell clustering, visualization, and projection using the
#' kohonen R package (self-organizing maps) or FastPG (fast graph-based
#' community detection).
#'
#' @name cluster
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Extract the expression matrix from heterogeneous input
#'
#' Accepts a matrix, data.frame, or tibble and returns a numeric matrix
#' restricted to the requested markers.
#'
#' @param corrected Input data.
#' @param lineage_markers Character vector of marker names.
#' @return A numeric matrix with rows = cells and columns = lineage_markers.
#' @noRd
.extract_marker_matrix <- function(corrected, lineage_markers) {
  if (is.matrix(corrected)) {
    missing_m <- setdiff(lineage_markers, colnames(corrected))
    if (length(missing_m) > 0) {
      stop("Markers not found in matrix columns: ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
    mat <- corrected[, lineage_markers, drop = FALSE]
  } else if (is.data.frame(corrected)) {
    missing_m <- setdiff(lineage_markers, names(corrected))
    if (length(missing_m) > 0) {
      stop("Markers not found in data columns: ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
    mat <- as.matrix(corrected[, lineage_markers, drop = FALSE])
  } else {
    stop("'corrected' must be a matrix, data.frame, or tibble.", call. = FALSE)
  }
  storage.mode(mat) <- "double"
  mat
}

# ---------------------------------------------------------------------------
# Main clustering entry point
# ---------------------------------------------------------------------------

#' Cluster Cells Using kohonen SOM or FastPG
#'
#' Clusters cells from expression data using either a self-organizing map
#' (SOM) via the \pkg{kohonen} package with hierarchical metaclustering, or
#' fast graph-based community detection via the \pkg{FastPG} package.
#'
#' @param corrected A \code{tibble}, \code{data.frame}, or numeric matrix
#'   containing the (typically batch-corrected) expression data.
#' @param lineage_markers Character vector of marker names to use for
#'   clustering. These should be lineage markers (e.g., CD3, CD4, CD8), not
#'   functional markers.
#' @param method Character; clustering method to use. One of \code{"som"}
#'   (default) for kohonen SOM with hierarchical metaclustering, or
#'   \code{"fastpg"} for FastPG graph-based clustering.
#' @param xdim Integer; SOM grid x-dimension (default: 10). Used only when
#'   \code{method = "som"}.
#' @param ydim Integer; SOM grid y-dimension (default: 10). Used only when
#'   \code{method = "som"}.
#' @param rlen Integer; number of times the complete dataset is presented to
#'   the SOM network (default: 10). Used only when \code{method = "som"}.
#' @param n_metaclusters Integer; number of metaclusters to identify
#'   (default: 20). Used only when \code{method = "som"}.
#' @param k Integer; number of nearest neighbours for the FastPG shared
#'   nearest-neighbour graph (default: 30). Used only when
#'   \code{method = "fastpg"}.
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to the underlying clustering
#'   function (\code{kohonen::som} or \code{FastPG::fastCluster}).
#'
#' @return An \code{sw_cluster_result} list with components:
#'   \describe{
#'     \item{\code{assignments}}{Integer vector of cluster/metacluster labels
#'       (one per cell).}
#'     \item{\code{data}}{The numeric matrix of expression values used for
#'       clustering.}
#'     \item{\code{lineage_markers}}{Character vector of marker names.}
#'     \item{\code{method}}{The clustering method used.}
#'     \item{\code{model}}{The underlying model object (\code{kohonen} SOM
#'       object for \code{"som"}, or \code{NULL} for \code{"fastpg"}).}
#'     \item{\code{n_clusters}}{Integer; total number of clusters found.}
#'   }
#'
#' @export
sw_cluster <- function(corrected, lineage_markers,
                       method = c("som", "fastpg"),
                       xdim = 10, ydim = 10, rlen = 10,
                       n_metaclusters = 20,
                       k = 30,
                       seed = 42, ...) {
  # --- Input validation (runs before dependency checks) ---
  if (!is.character(lineage_markers) || length(lineage_markers) == 0) {
    stop("'lineage_markers' must be a non-empty character vector.",
         call. = FALSE)
  }

  method <- match.arg(method)

  if (method == "som") {
    if (!is.numeric(n_metaclusters) || n_metaclusters < 2) {
      stop("'n_metaclusters' must be an integer >= 2.", call. = FALSE)
    }
  }

  mat <- .extract_marker_matrix(corrected, lineage_markers)

  if (method == "som") {
    .cluster_som(mat, lineage_markers, xdim = xdim, ydim = ydim,
                 rlen = rlen, n_metaclusters = n_metaclusters,
                 seed = seed, ...)
  } else {
    .cluster_fastpg(mat, lineage_markers, k = k, seed = seed, ...)
  }
}

# ---------------------------------------------------------------------------
# kohonen SOM method
# ---------------------------------------------------------------------------

#' @noRd
.cluster_som <- function(mat, lineage_markers, xdim, ydim, rlen,
                         n_metaclusters, seed, ...) {
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required for SOM clustering. ",
         "Install it with install.packages('kohonen').", call. = FALSE)
  }

  set.seed(seed)

  grid <- kohonen::somgrid(xdim = xdim, ydim = ydim, topo = "hexagonal")

  message("Running kohonen SOM on ", nrow(mat), " cells with ",
          length(lineage_markers), " markers (", xdim, "x", ydim,
          " grid, ", n_metaclusters, " metaclusters)")

  som_model <- kohonen::som(mat, grid = grid, rlen = rlen, ...)

  # Hierarchical metaclustering of SOM codes
  codes <- som_model$codes[[1]]
  d <- stats::dist(codes)
  hc <- stats::hclust(d, method = "ward.D2")
  node_metaclusters <- stats::cutree(hc, k = n_metaclusters)

  # Map node metaclusters to cell-level assignments
  cell_assignments <- node_metaclusters[som_model$unit.classif]

  message("kohonen SOM clustering complete")


  result <- list(
    assignments = as.integer(cell_assignments),
    data = mat,
    lineage_markers = lineage_markers,
    method = "som",
    model = som_model,
    n_clusters = as.integer(n_metaclusters),
    node_metaclusters = node_metaclusters
  )
  class(result) <- "sw_cluster_result"
  result
}

# ---------------------------------------------------------------------------
# FastPG method
# ---------------------------------------------------------------------------

#' @noRd
.cluster_fastpg <- function(mat, lineage_markers, k, seed, ...) {
  if (!requireNamespace("FastPG", quietly = TRUE)) {
    stop("Package 'FastPG' is required for FastPG clustering. ",
         "Install it from GitHub: remotes::install_github('sararselitsky/FastPG').",
         call. = FALSE)
  }

  set.seed(seed)

  message("Running FastPG on ", nrow(mat), " cells with ",
          length(lineage_markers), " markers (k=", k, ")")

  fpg <- FastPG::fastCluster(mat, k = k, ...)

  communities <- as.integer(fpg$communities)
  n_clusters <- length(unique(communities))

  message("FastPG clustering complete (", n_clusters, " communities found)")

  result <- list(
    assignments = communities,
    data = mat,
    lineage_markers = lineage_markers,
    method = "fastpg",
    model = NULL,
    n_clusters = n_clusters
  )
  class(result) <- "sw_cluster_result"
  result
}

# ---------------------------------------------------------------------------
# Accessor functions
# ---------------------------------------------------------------------------

#' Get Cluster Assignments
#'
#' Extracts per-cell cluster assignments from a clustering result.
#'
#' @param cluster_result An \code{sw_cluster_result} object (output from
#'   \code{\link{sw_cluster}}).
#'
#' @return An integer vector of cluster assignments, one per cell.
#'
#' @export
sw_get_cluster_assignments <- function(cluster_result) {
  if (!inherits(cluster_result, "sw_cluster_result")) {
    stop("'cluster_result' must be an sw_cluster_result object.",
         call. = FALSE)
  }

  cluster_result$assignments
}

#' Compute Median Fluorescence Intensities per Cluster
#'
#' Calculates the median expression of each marker within each cluster.
#'
#' @param cluster_result An \code{sw_cluster_result} object (output from
#'   \code{\link{sw_cluster}}).
#'
#' @return A \code{tibble} with one row per cluster and one column per
#'   marker, plus a \code{cluster} column.
#'
#' @export
sw_cluster_mfis <- function(cluster_result) {
  if (!inherits(cluster_result, "sw_cluster_result")) {
    stop("'cluster_result' must be an sw_cluster_result object.",
         call. = FALSE)
  }

  mat <- cluster_result$data
  assignments <- cluster_result$assignments
  markers <- cluster_result$lineage_markers

  clusters <- sort(unique(assignments))
  mfi_list <- lapply(clusters, function(cl) {
    idx <- which(assignments == cl)
    vapply(markers, function(m) stats::median(mat[idx, m]), numeric(1))
  })
  mfi_mat <- do.call(rbind, mfi_list)

  df <- tibble::as_tibble(as.data.frame(mfi_mat))
  df$cluster <- clusters

  # Reorder: cluster first
  df <- df[, c("cluster", setdiff(names(df), "cluster"))]
  df
}

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

#' Plot Cluster Heatmap
#'
#' Generates a heatmap of median fluorescence intensities per cluster and
#' saves it to a PDF file. Uses \pkg{pheatmap} if available, otherwise
#' falls back to \code{\link[stats]{heatmap}}.
#'
#' @param cluster_result An \code{sw_cluster_result} object (output from
#'   \code{\link{sw_cluster}}).
#' @param plot_file Character path for the output PDF file
#'   (default: \code{"cluster_heatmap.pdf"}).
#'
#' @return The file path (invisibly).
#'
#' @export
sw_plot_clusters <- function(cluster_result,
                             plot_file = "cluster_heatmap.pdf") {
  if (!inherits(cluster_result, "sw_cluster_result")) {
    stop("'cluster_result' must be an sw_cluster_result object.",
         call. = FALSE)
  }

  if (!is.character(plot_file) || length(plot_file) != 1) {
    stop("'plot_file' must be a single file path.", call. = FALSE)
  }

  # Ensure output directory exists
  out_dir <- dirname(plot_file)
  if (!dir.exists(out_dir) && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE)
  }

  mfi <- sw_cluster_mfis(cluster_result)
  mat_plot <- as.matrix(mfi[, setdiff(names(mfi), "cluster")])
  rownames(mat_plot) <- paste0("C", mfi$cluster)

  grDevices::pdf(plot_file, width = 12, height = 8)
  tryCatch({
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(mat_plot, scale = "column",
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         main = paste0("Cluster MFI Heatmap (",
                                       cluster_result$method, ")"))
    } else {
      stats::heatmap(mat_plot, scale = "column",
                     main = paste0("Cluster MFI Heatmap (",
                                   cluster_result$method, ")"))
    }
  }, finally = {
    grDevices::dev.off()
  })

  message("Cluster heatmap saved to: ", plot_file)
  invisible(plot_file)
}

# ---------------------------------------------------------------------------
# Prediction / mapping
# ---------------------------------------------------------------------------

#' Map New Data onto Trained SOM
#'
#' Projects new samples onto an existing trained kohonen SOM model using
#' \code{kohonen::map.kohonen()}. Only available when the clustering was
#' performed with \code{method = "som"}.
#'
#' @param cluster_result An \code{sw_cluster_result} object (output from
#'   \code{\link{sw_cluster}} with \code{method = "som"}).
#' @param newdata A \code{tibble}, \code{data.frame}, or numeric matrix
#'   containing the new expression data. Must include all lineage markers
#'   used in the original clustering.
#'
#' @return An \code{sw_cluster_result} object with the new data mapped onto
#'   the existing SOM and metacluster assignments.
#'
#' @export
sw_predict_clusters <- function(cluster_result, newdata) {
  if (!inherits(cluster_result, "sw_cluster_result")) {
    stop("'cluster_result' must be an sw_cluster_result object.",
         call. = FALSE)
  }

  if (cluster_result$method != "som") {
    stop("Prediction is only supported for method = 'som'.", call. = FALSE)
  }

  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required.", call. = FALSE)
  }

  markers <- cluster_result$lineage_markers
  new_mat <- .extract_marker_matrix(newdata, markers)

  message("Mapping ", nrow(new_mat), " new cells onto trained SOM")

  mapped <- kohonen::map.kohonen(cluster_result$model, newdata = new_mat)
  node_metaclusters <- cluster_result$node_metaclusters
  cell_assignments <- node_metaclusters[mapped$unit.classif]

  message("Mapping complete")

  result <- list(
    assignments = as.integer(cell_assignments),
    data = new_mat,
    lineage_markers = markers,
    method = "som",
    model = cluster_result$model,
    n_clusters = cluster_result$n_clusters,
    node_metaclusters = node_metaclusters
  )
  class(result) <- "sw_cluster_result"
  result
}
