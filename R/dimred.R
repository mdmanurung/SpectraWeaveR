#' @title Dimensionality Reduction
#'
#' @description
#' Standalone dimensionality reduction functions for visualising
#' high-dimensional spectral flow cytometry data.  Supports UMAP (via
#' \pkg{uwot}) and PCA (via \code{stats::prcomp}).
#'
#' @name dimred
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# sw_dimred_run
# ---------------------------------------------------------------------------

#' Run Dimensionality Reduction
#'
#' Computes a 2D embedding of marker expression data using UMAP or PCA.
#' Returns a tibble with coordinates and all non-marker metadata columns
#' from the input preserved.
#'
#' @param df A \code{data.frame} or \code{tibble} containing marker columns
#'   and optional metadata columns (e.g., batch, sample, cluster).
#' @param markers Character vector of marker column names to include in
#'   the reduction.
#' @param method Character; dimensionality reduction method.
#'   One of \code{"umap"} (default) or \code{"pca"}.
#' @param n_dims Integer; number of output dimensions (default: 2).
#' @param n_neighbors Integer; number of neighbours for UMAP
#'   (default: 15).  Ignored for PCA.
#' @param min_dist Numeric; minimum distance parameter for UMAP
#'   (default: 0.2).  Ignored for PCA.
#' @param max_cells Integer or \code{NULL}; if the input has more rows
#'   than \code{max_cells}, a random subsample is drawn before computing
#'   the reduction (default: 50000).  Set \code{NULL} to disable.
#' @param scale Logical; whether to centre and scale columns before
#'   reduction (default: \code{TRUE}).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments passed to \code{uwot::umap()} or
#'   \code{stats::prcomp()}.
#'
#' @return A \code{tibble} with columns \code{dim1}, \code{dim2}
#'   (and \code{dim3}, etc. when \code{n_dims > 2}), plus all non-marker
#'   columns from the input.  An attribute \code{"method"} records the
#'   method used, and for PCA an attribute \code{"variance_explained"}
#'   stores the proportion of variance explained per component.
#'
#' @seealso \code{\link{sw_plot_dimred}},
#'   \code{\link{sw_plot_batch_dimred}}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(CD3 = rnorm(200), CD4 = rnorm(200), batch = "B1")
#' # PCA (no extra dependency)
#' result <- sw_dimred_run(df, markers = c("CD3", "CD4"), method = "pca")
#' head(result)
#' # UMAP (requires uwot)
#' result <- sw_dimred_run(df, markers = c("CD3", "CD4"), method = "umap")
#' }
#'
#' @export
sw_dimred_run <- function(df, markers,
                          method = c("umap", "pca"),
                          n_dims = 2L,
                          n_neighbors = 15L,
                          min_dist = 0.2,
                          max_cells = 50000L,
                          scale = TRUE,
                          seed = 42L,
                          ...) {
  # --- Input validation ---
  .validate_df(df, "df")
  .validate_markers(markers)
  .validate_markers_in_df(df, markers, "df")

  method <- match.arg(method)

  if (!is.numeric(n_dims) || length(n_dims) != 1 || n_dims < 1) {
    stop("'n_dims' must be a positive integer.", call. = FALSE)
  }
  n_dims <- as.integer(n_dims)

  # --- Subsample if needed ---
  set.seed(seed)

  idx <- seq_len(nrow(df))
  if (!is.null(max_cells) && nrow(df) > max_cells) {
    message("Subsampling from ", nrow(df), " to ", max_cells, " cells")
    idx <- sort(sample.int(nrow(df), max_cells))
    df <- df[idx, , drop = FALSE]
  }

  # --- Extract marker matrix ---
  mat <- as.matrix(df[, markers, drop = FALSE])
  storage.mode(mat) <- "double"

  # Identify non-marker (metadata) columns
  meta_cols <- setdiff(names(df), markers)
  meta_df <- df[, meta_cols, drop = FALSE]

  # --- Run reduction ---
  if (method == "umap") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Package 'uwot' is required for UMAP. ",
           "Install it with: install.packages('uwot')", call. = FALSE)
    }

    message("Computing UMAP (", nrow(mat), " cells, ",
            length(markers), " markers)")

    coords <- uwot::umap(mat,
                          n_components = n_dims,
                          n_neighbors = n_neighbors,
                          min_dist = min_dist,
                          scale = scale,
                          ret_model = FALSE,
                          ...)

    var_explained <- NULL

  } else {
    # PCA
    message("Computing PCA (", nrow(mat), " cells, ",
            length(markers), " markers)")

    pca <- stats::prcomp(mat, center = scale, scale. = scale,
                         rank. = n_dims, ...)
    coords <- pca$x[, seq_len(n_dims), drop = FALSE]

    total_var <- sum(pca$sdev^2)
    var_explained <- (pca$sdev[seq_len(n_dims)]^2) / total_var
  }

  # --- Build output tibble ---
  coord_names <- paste0("dim", seq_len(n_dims))
  colnames(coords) <- coord_names
  result <- tibble::as_tibble(as.data.frame(coords))

  # Append metadata columns
  if (ncol(meta_df) > 0) {
    result <- dplyr::bind_cols(result, meta_df)
  }

  attr(result, "method") <- method
  if (!is.null(var_explained)) {
    attr(result, "variance_explained") <- var_explained
  }

  result
}


# ---------------------------------------------------------------------------
# sw_plot_dimred
# ---------------------------------------------------------------------------

#' Plot Dimensionality Reduction
#'
#' Creates a scatter plot of a dimensionality reduction result, coloured by
#' a user-specified column.  Numeric columns use a continuous viridis scale;
#' categorical columns use discrete colours.
#'
#' @param dimred_result A \code{tibble} returned by
#'   \code{\link{sw_dimred_run}}, containing \code{dim1} and \code{dim2}
#'   columns.
#' @param color_by Character; name of the column to map to point colour.
#'   Can be a marker name, \code{"cluster"}, \code{"batch"}, etc.
#' @param point_size Numeric; size of points (default: 0.5).
#' @param alpha Numeric; point transparency (default: 0.6).
#' @param title Character or \code{NULL}; plot title.
#' @param ... Additional arguments passed to
#'   \code{ggplot2::geom_point()}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{sw_dimred_run}}
#'
#' @examples
#' \dontrun{
#' dr <- data.frame(dim1 = rnorm(100), dim2 = rnorm(100),
#'                  cluster = sample(1:3, 100, replace = TRUE))
#' attr(dr, "method") <- "pca"
#' sw_plot_dimred(dr, color_by = "cluster")
#' }
#'
#' @export
sw_plot_dimred <- function(dimred_result, color_by = NULL,
                           point_size = 0.5, alpha = 0.6,
                           title = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Install it with: install.packages('ggplot2')", call. = FALSE)
  }

  .validate_df(dimred_result, "dimred_result")

  if (!all(c("dim1", "dim2") %in% names(dimred_result))) {
    stop("'dimred_result' must contain 'dim1' and 'dim2' columns. ",
         "Use sw_dimred_run() to compute coordinates.", call. = FALSE)
  }

  if (!is.null(color_by) && !color_by %in% names(dimred_result)) {
    stop("Column '", color_by, "' not found in dimred_result.", call. = FALSE)
  }

  method <- attr(dimred_result, "method")
  if (is.null(method)) method <- "DR"
  axis_label <- toupper(method)

  if (is.null(color_by)) {
    p <- ggplot2::ggplot(dimred_result,
                         ggplot2::aes(x = .data$dim1, y = .data$dim2)) +
      ggplot2::geom_point(size = point_size, alpha = alpha, ...)
  } else {
    p <- ggplot2::ggplot(
      dimred_result,
      ggplot2::aes(x = .data$dim1, y = .data$dim2,
                   color = .data[[color_by]])
    ) +
      ggplot2::geom_point(size = point_size, alpha = alpha, ...)

    # Use continuous scale for numeric, discrete for factor/character
    if (is.numeric(dimred_result[[color_by]])) {
      p <- p + ggplot2::scale_color_viridis_c()
    }
  }

  p <- p +
    ggplot2::labs(
      x = paste0(axis_label, " 1"),
      y = paste0(axis_label, " 2"),
      title = title,
      color = color_by
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  p
}
