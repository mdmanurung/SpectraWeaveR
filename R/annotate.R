#' @title Cell Type Annotation
#'
#' @description
#' Functions for annotating cell clusters from spectral flow cytometry data.
#' Supports automated annotation via cosine similarity against a built-in
#' PBMC reference panel (\code{\link{sw_annotate_run}}), manual
#' annotation via a user-supplied cluster-to-label mapping
#' (\code{\link{sw_annotate_manual}}), and annotation visualization
#' (\code{\link{sw_plot_annotation}}).
#'
#' The built-in reference covers 35 populations across 5 lineages (PBMC core,
#' myeloid subtypes, T-cell subtypes, B-cell subtypes, and NK subtypes) and
#' 30 common spectral flow cytometry markers.  A binary marker-relevance mask
#' focuses the cosine similarity computation on the markers that are
#' biologically informative for each population, reducing noise from markers
#' that do not discriminate between cell types.
#'
#' @name annotate
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Validate and normalise annotation input
#'
#' Accepts either an \code{sw_cluster_result} object (from
#' \code{\link{sw_cluster_run}}) or a data.frame/tibble produced by
#' \code{\link{sw_cluster_mfi}}.  Returns a tibble with a \code{cluster}
#' column and numeric marker columns.
#'
#' @param x Input object.
#' @return A data.frame with at least a \code{cluster} column.
#' @noRd
.check_annotation_input <- function(x) {
  if (inherits(x, "sw_cluster_result")) {
    return(sw_cluster_mfi(x))
  }
  if (!is.data.frame(x)) {
    stop(
      "'x' must be an 'sw_cluster_result' object or a data.frame/tibble ",
      "from sw_cluster_mfi().",
      call. = FALSE
    )
  }
  if (!"cluster" %in% names(x)) {
    stop(
      "'x' must contain a 'cluster' column when passing a data.frame/tibble.",
      call. = FALSE
    )
  }
  numeric_cols <- vapply(x, is.numeric, logical(1))
  if (sum(numeric_cols) < 2) {
    stop(
      "'x' must have at least one numeric marker column in addition to 'cluster'.",
      call. = FALSE
    )
  }
  x
}

#' Cosine similarity between two numeric vectors
#'
#' Returns 0 when either vector has zero L2-norm (avoids division by zero).
#'
#' @param a,b Numeric vectors of the same length.
#' @return A single numeric in \eqn{[-1, 1]}.
#' @noRd
.cosine_sim <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (denom < .Machine$double.eps) return(0)
  sum(a * b) / denom
}

#' Score one cluster vector against all reference populations
#'
#' For each reference population, restricts to the informative (mask == 1)
#' markers that are also present in \code{available_markers}, then computes
#' the cosine similarity between the min-max-normalised cluster MFI vector
#' and the reference vector (also normalised to the 0--3 range).
#'
#' @param mfi_vec Named numeric vector of MFI values for one cluster.
#' @param ref_matrix Numeric matrix (populations \eqn{\times} markers).
#' @param ref_mask Integer matrix (populations \eqn{\times} markers, 0/1).
#' @param available_markers Character vector of markers present in both the
#'   data and the reference.
#' @return A tibble with columns \code{population} and \code{score}.
#' @noRd
.score_cluster <- function(mfi_vec, ref_matrix, ref_mask, available_markers) {
  pops <- rownames(ref_matrix)

  scores <- vapply(seq_len(nrow(ref_matrix)), function(i) {
    mask_i     <- as.integer(ref_mask[i, available_markers, drop = TRUE])
    informative <- available_markers[mask_i == 1L]

    if (length(informative) < 2L) return(0)

    a <- mfi_vec[informative]
    b <- ref_matrix[i, informative, drop = TRUE]

    # Min-max normalise the cluster vector within the informative markers
    a_rng <- range(a, na.rm = TRUE)
    if (diff(a_rng) < .Machine$double.eps) return(0)
    a_norm <- (a - a_rng[1]) / diff(a_rng)

    # Normalise the reference (scale already 0--3)
    b_norm <- b / 3

    .cosine_sim(a_norm, b_norm)
  }, numeric(1))

  tibble::tibble(population = pops, score = scores)
}

# ---------------------------------------------------------------------------
# sw_annotate_load_ref
# ---------------------------------------------------------------------------

#' Load a Built-in Cell Type Reference Panel
#'
#' Reads the built-in reference expression matrix and marker-relevance mask
#' shipped with SpectraWeaveR.  The reference covers 35 PBMC, myeloid, T, B,
#' and NK sub-populations across 30 common spectral flow cytometry markers.
#'
#' @param name Character; name of the reference panel to load.  Currently
#'   only \code{"pbmc"} is supported (default).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{matrix}}{A numeric matrix (populations \eqn{\times} markers)
#'       with expression levels on a 0--3 scale (0 = negative, 1 = dim,
#'       2 = positive, 3 = bright).}
#'     \item{\code{mask}}{A binary integer matrix (populations \eqn{\times}
#'       markers) where 1 indicates the marker is informative for identifying
#'       that population.}
#'   }
#'   Row names of both matrices are population labels; column names are
#'   marker names.
#'
#' @details
#' The reference is stored as two CSV files in
#' \code{inst/extdata/pbmc_reference_matrix.csv} and
#' \code{inst/extdata/pbmc_reference_mask.csv}.  Users can build a custom
#' reference by constructing a list with the same structure:
#' \code{list(matrix = my_mat, mask = my_mask)}.
#'
#' @examples
#' ref <- sw_annotate_load_ref("pbmc")
#' dim(ref$matrix)   # 35 populations × 30 markers
#' rownames(ref$matrix)[1:5]
#'
#' @export
sw_annotate_load_ref <- function(name = "pbmc") {
  if (!is.character(name) || length(name) != 1L || nchar(name) == 0L) {
    stop("'name' must be a single non-empty character string.", call. = FALSE)
  }
  name <- match.arg(name, "pbmc")

  mat_path  <- system.file("extdata", "pbmc_reference_matrix.csv",
                            package = "SpectraWeaveR")
  mask_path <- system.file("extdata", "pbmc_reference_mask.csv",
                            package = "SpectraWeaveR")

  if (nchar(mat_path) == 0L || nchar(mask_path) == 0L) {
    stop(
      "Built-in reference data files not found. ",
      "Please reinstall SpectraWeaveR.",
      call. = FALSE
    )
  }

  mat_df  <- utils::read.csv(mat_path,  row.names = 1L,
                              check.names = FALSE, stringsAsFactors = FALSE)
  mask_df <- utils::read.csv(mask_path, row.names = 1L,
                              check.names = FALSE, stringsAsFactors = FALSE)

  mat  <- as.matrix(mat_df)
  mask <- as.matrix(mask_df)
  storage.mode(mat)  <- "double"
  storage.mode(mask) <- "integer"

  list(matrix = mat, mask = mask)
}

# ---------------------------------------------------------------------------
# sw_annotate_run
# ---------------------------------------------------------------------------

#' Automatically Annotate Clusters by Cell Type
#'
#' Assigns cell type labels to each cluster by comparing cluster median
#' fluorescence intensity (MFI) profiles against a reference panel using
#' cosine similarity.  By default the built-in 35-population PBMC reference
#' is used; users can supply a custom reference via the \code{reference}
#' argument.
#'
#' @param x An \code{sw_cluster_result} object (returned by
#'   \code{\link{sw_cluster_run}}) or a \code{data.frame}/\code{tibble} with a
#'   \code{cluster} column and numeric marker columns (e.g., the output of
#'   \code{\link{sw_cluster_mfi}}).
#' @param reference A named list with elements \code{matrix} and \code{mask}
#'   (as returned by \code{\link{sw_annotate_load_ref}}).  If \code{NULL}
#'   (default) the built-in PBMC reference is loaded automatically.
#' @param markers Character vector of marker names to use for annotation.
#'   If \code{NULL} (default) all numeric columns in \code{x} that are also
#'   present in the reference are used.
#' @param min_score Numeric threshold in \eqn{[0, 1]}.  Clusters whose
#'   top cosine similarity score is below this value are labelled
#'   \code{NA} rather than assigned a cell type.  Default: 0.3.
#' @param ... Reserved for future use.
#'
#' @return A \code{tibble} with one row per cluster and columns:
#'   \describe{
#'     \item{\code{cluster}}{Cluster identifier.}
#'     \item{\code{cell_type}}{Top annotation label, or \code{NA} if the
#'       best score is below \code{min_score}.}
#'     \item{\code{score}}{Cosine similarity score for the top annotation.}
#'     \item{\code{second_type}}{Second-best annotation label.}
#'     \item{\code{second_score}}{Cosine similarity score for the second
#'       annotation.}
#'     \item{\code{n_markers_used}}{Number of markers in the overlap between
#'       the data and the reference.}
#'   }
#'
#' @details
#' The annotation algorithm:
#' \enumerate{
#'   \item Extract the MFI tibble (via \code{\link{sw_cluster_mfi}} when
#'     \code{x} is an \code{sw_cluster_result}).
#'   \item Identify the intersection of available marker names and reference
#'     column names.
#'   \item For each cluster, min-max-normalise the MFI vector within the
#'     informative markers for each reference population (those where
#'     \code{mask == 1}), then compute the cosine similarity against the
#'     corresponding reference vector (normalised to \eqn{[0, 1]} by
#'     dividing by 3).
#'   \item Return the top-1 and top-2 results per cluster.
#' }
#'
#' @seealso \code{\link{sw_annotate_manual}}, \code{\link{sw_plot_annotation}},
#'   \code{\link{sw_cluster_mfi}}, \code{\link{sw_annotate_load_ref}}
#'
#' @examples
#' \dontrun{
#' # After clustering:
#' result <- sw_cluster_run(corrected, lineage_markers = markers)
#'
#' # Auto-annotate using built-in PBMC reference
#' annotation <- sw_annotate_run(result)
#' head(annotation)
#'
#' # Use only specific markers for annotation
#' annotation <- sw_annotate_run(result, markers = c("CD3", "CD4", "CD8a"))
#'
#' # Use a custom reference
#' my_ref <- list(
#'   matrix = my_matrix,  # populations × markers, values 0–3
#'   mask   = my_mask     # populations × markers, values 0/1
#' )
#' annotation <- sw_annotate_run(result, reference = my_ref)
#' }
#'
#' @export
sw_annotate_run <- function(x, reference = NULL, markers = NULL,
                                  min_score = 0.3, ...) {

  mfi_df <- .check_annotation_input(x)

  if (is.null(reference)) {
    reference <- sw_annotate_load_ref("pbmc")
  }

  if (!is.list(reference) || !all(c("matrix", "mask") %in% names(reference))) {
    stop(
      "'reference' must be a list with 'matrix' and 'mask' elements ",
      "(as returned by sw_annotate_load_ref()).",
      call. = FALSE
    )
  }

  ref_matrix <- reference$matrix
  ref_mask   <- reference$mask

  if (!is.matrix(ref_matrix) || !is.matrix(ref_mask)) {
    stop("'reference$matrix' and 'reference$mask' must be matrices.",
         call. = FALSE)
  }
  if (!identical(dim(ref_matrix), dim(ref_mask))) {
    stop(
      "'reference$matrix' and 'reference$mask' must have the same dimensions.",
      call. = FALSE
    )
  }
  if (!identical(dimnames(ref_matrix), dimnames(ref_mask))) {
    stop(
      "'reference$matrix' and 'reference$mask' must have identical ",
      "row and column names.",
      call. = FALSE
    )
  }

  if (!is.numeric(min_score) || length(min_score) != 1L ||
      min_score < 0 || min_score > 1) {
    stop("'min_score' must be a single number in [0, 1].", call. = FALSE)
  }

  # Determine marker columns
  meta_cols   <- c("cluster")
  all_numeric <- setdiff(names(mfi_df), meta_cols)
  all_numeric <- all_numeric[vapply(mfi_df[all_numeric], is.numeric, logical(1))]

  if (!is.null(markers)) {
    if (!is.character(markers) || length(markers) == 0L) {
      stop("'markers' must be a non-empty character vector or NULL.", call. = FALSE)
    }
    missing_m <- setdiff(markers, all_numeric)
    if (length(missing_m) > 0L) {
      stop("Markers not found in MFI data: ",
           paste(missing_m, collapse = ", "), call. = FALSE)
    }
    marker_cols <- markers
  } else {
    marker_cols <- all_numeric
  }

  # Intersect with reference columns
  available_markers <- intersect(marker_cols, colnames(ref_matrix))
  if (length(available_markers) < 2L) {
    stop(
      "Fewer than 2 markers overlap between data and reference. ",
      "Check that marker names match the reference ",
      "(e.g., 'CD4', 'CD8a', 'HLA_DR', 'GranzymeB').",
      call. = FALSE
    )
  }

  # Score each cluster
  n_clusters <- nrow(mfi_df)
  results    <- vector("list", n_clusters)

  for (i in seq_len(n_clusters)) {
    mfi_vec        <- as.numeric(mfi_df[i, available_markers, drop = TRUE])
    names(mfi_vec) <- available_markers

    scores_df <- .score_cluster(mfi_vec, ref_matrix, ref_mask, available_markers)
    scores_df <- scores_df[order(scores_df$score, decreasing = TRUE), ]

    top1_type  <- as.character(scores_df$population[1L])
    top1_score <- scores_df$score[1L]
    top2_type  <- if (nrow(scores_df) >= 2L)
                    as.character(scores_df$population[2L])
                  else NA_character_
    top2_score <- if (nrow(scores_df) >= 2L)
                    scores_df$score[2L]
                  else NA_real_

    # Apply minimum score threshold
    if (!is.na(top1_score) && top1_score < min_score) {
      top1_type <- NA_character_
    }

    results[[i]] <- tibble::tibble(
      cluster        = mfi_df[["cluster"]][i],
      cell_type      = top1_type,
      score          = top1_score,
      second_type    = top2_type,
      second_score   = top2_score,
      n_markers_used = length(available_markers)
    )
  }

  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# sw_annotate_manual
# ---------------------------------------------------------------------------

#' Manually Annotate Clusters with User-supplied Labels
#'
#' Appends a \code{cell_type} column to the cluster MFI tibble using a
#' user-supplied named character vector that maps cluster identifiers to
#' cell type labels.
#'
#' @param x An \code{sw_cluster_result} object or a \code{data.frame}/
#'   \code{tibble} from \code{\link{sw_cluster_mfi}} (must have a
#'   \code{cluster} column).
#' @param annotation_map A named character vector where names are cluster
#'   identifiers (coerced to character) and values are cell type labels.
#'   Example: \code{c("1" = "CD4+ T cell", "2" = "B cell")}.
#'
#' @return A \code{tibble} identical to the MFI tibble from
#'   \code{x} with a \code{cell_type} column inserted immediately after
#'   \code{cluster}.  Clusters absent from \code{annotation_map} receive
#'   \code{NA} and a warning is issued.
#'
#' @seealso \code{\link{sw_annotate_run}}, \code{\link{sw_plot_annotation}}
#'
#' @examples
#' \dontrun{
#' result <- sw_cluster_run(corrected, lineage_markers = markers)
#' mfis   <- sw_cluster_mfi(result)
#'
#' annotation <- sw_annotate_manual(mfis, c(
#'   "1"  = "CD4+ T cell",
#'   "2"  = "CD8+ T cell",
#'   "3"  = "B cell",
#'   "4"  = "NK cell",
#'   "5"  = "CD14+ Monocyte"
#' ))
#' head(annotation)
#' }
#'
#' @export
sw_annotate_manual <- function(x, annotation_map) {
  mfi_df <- .check_annotation_input(x)

  if (!is.character(annotation_map)) {
    stop(
      "'annotation_map' must be a named character vector, e.g., ",
      "c('1' = 'CD4+ T cell', '2' = 'B cell').",
      call. = FALSE
    )
  }
  if (is.null(names(annotation_map)) || any(nchar(names(annotation_map)) == 0L)) {
    stop(
      "'annotation_map' must be fully named (all elements have non-empty names).",
      call. = FALSE
    )
  }

  # Coerce cluster IDs to character for matching
  cluster_ids <- as.character(mfi_df[["cluster"]])
  map_keys    <- names(annotation_map)

  missing_ids <- setdiff(cluster_ids, map_keys)
  if (length(missing_ids) > 0L) {
    warning(
      "Cluster ID(s) not found in 'annotation_map': ",
      paste(missing_ids, collapse = ", "),
      ". Setting cell_type to NA for these clusters.",
      call. = FALSE
    )
  }

  cell_types        <- annotation_map[cluster_ids]
  names(cell_types) <- NULL

  tibble::add_column(mfi_df, cell_type = unname(cell_types), .after = "cluster")
}

# ---------------------------------------------------------------------------
# sw_plot_annotation
# ---------------------------------------------------------------------------

#' Visualize Cluster Cell Type Annotations
#'
#' Produces either a heatmap of cluster MFI profiles annotated with cell type
#' labels, or a UMAP scatter plot coloured by cell type.
#'
#' @param annotation A \code{data.frame}/\code{tibble} containing at least
#'   a \code{cluster} column (and a \code{cell_type} column for UMAP plots).
#'   Typically the output of \code{\link{sw_annotate_run}} or
#'   \code{\link{sw_annotate_manual}}.
#' @param cluster_result An \code{sw_cluster_result} object.  Required when
#'   \code{type = "heatmap"} and \code{annotation} does not already contain
#'   marker MFI columns.
#' @param dimred A \code{tibble} from \code{\link{sw_dimred_run}} with columns
#'   \code{dim1}, \code{dim2}, and \code{cluster}.  Required when
#'   \code{type = "umap"}.
#' @param type Character; plot type.  One of \code{"heatmap"} (default) or
#'   \code{"umap"}.
#'
#' @return A \pkg{ggplot2} object (for UMAP) or a \pkg{pheatmap} object
#'   (for heatmap, when \pkg{pheatmap} is installed; otherwise a
#'   \pkg{ggplot2} tile plot).
#'
#' @seealso \code{\link{sw_annotate_run}}, \code{\link{sw_dimred_run}}
#'
#' @examples
#' \dontrun{
#' annotation <- sw_annotate_run(result)
#'
#' # Heatmap
#' sw_plot_annotation(annotation, cluster_result = result)
#'
#' # UMAP overlay (requires prior dimensionality reduction with cluster column)
#' dr <- sw_dimred_run(dplyr::mutate(corrected,
#'                     cluster = sw_cluster_assignments(result)),
#'                     markers = lineage_markers)
#' sw_plot_annotation(annotation, dimred = dr, type = "umap")
#' }
#'
#' @export
sw_plot_annotation <- function(annotation,
                                cluster_result = NULL,
                                dimred         = NULL,
                                type           = c("heatmap", "umap")) {
  type <- match.arg(type)

  if (!is.data.frame(annotation)) {
    stop(
      "'annotation' must be a data.frame or tibble with a 'cluster' column.",
      call. = FALSE
    )
  }
  if (!"cluster" %in% names(annotation)) {
    stop("'annotation' must contain a 'cluster' column.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for sw_plot_annotation().\n",
      "Install with: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  # ---- UMAP plot ---------------------------------------------------------
  if (type == "umap") {
    if (is.null(dimred)) {
      stop(
        "'dimred' must be provided for type = 'umap'. ",
        "Run sw_dimred_run() first and pass the result.",
        call. = FALSE
      )
    }
    if (!is.data.frame(dimred) ||
        !all(c("dim1", "dim2") %in% names(dimred))) {
      stop(
        "'dimred' must be a tibble with 'dim1' and 'dim2' columns ",
        "(output of sw_dimred_run()).",
        call. = FALSE
      )
    }
    if (!"cluster" %in% names(dimred)) {
      stop(
        "'dimred' must contain a 'cluster' column for annotation mapping.",
        call. = FALSE
      )
    }
    if (!"cell_type" %in% names(annotation)) {
      stop(
        "'annotation' must contain a 'cell_type' column for type = 'umap'.",
        call. = FALSE
      )
    }

    # Merge annotation into dimred
    ann_sub  <- annotation[, c("cluster", "cell_type"), drop = FALSE]
    plot_df  <- merge(dimred, ann_sub, by = "cluster", all.x = TRUE)

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x      = .data[["dim1"]],
          y      = .data[["dim2"]],
          colour = .data[["cell_type"]]
        )
      ) +
      ggplot2::geom_point(size = 0.5, alpha = 0.6) +
      ggplot2::labs(x = "Dim 1", y = "Dim 2", colour = "Cell type") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right")

    return(p)
  }

  # ---- Heatmap -----------------------------------------------------------
  # Determine which columns hold marker MFI values
  non_marker_cols <- c("cluster", "cell_type", "score", "second_type",
                       "second_score", "n_markers_used")
  marker_cols <- setdiff(names(annotation), non_marker_cols)
  numeric_mask <- vapply(annotation[marker_cols], is.numeric, logical(1))
  marker_cols  <- marker_cols[numeric_mask]

  if (length(marker_cols) == 0L) {
    if (is.null(cluster_result)) {
      stop(
        "For type = 'heatmap', either 'annotation' must include numeric ",
        "marker MFI columns (from sw_annotate_manual()) or 'cluster_result' ",
        "must be provided.",
        call. = FALSE
      )
    }
    if (!inherits(cluster_result, "sw_cluster_result")) {
      stop("'cluster_result' must be an 'sw_cluster_result' object.",
           call. = FALSE)
    }
    mfi_df      <- sw_cluster_mfi(cluster_result)
    marker_cols <- setdiff(names(mfi_df), "cluster")
  } else {
    mfi_df <- annotation[, c("cluster", marker_cols), drop = FALSE]
  }

  # Build row labels: cell type (if available) or "Cluster N"
  has_labels <- "cell_type" %in% names(annotation)
  if (has_labels) {
    label_map <- stats::setNames(
      as.character(annotation[["cell_type"]]),
      as.character(annotation[["cluster"]])
    )
    row_labels <- label_map[as.character(mfi_df[["cluster"]])]
    row_labels[is.na(row_labels)] <- paste0(
      "Cluster ", mfi_df[["cluster"]][is.na(row_labels)]
    )
  } else {
    row_labels <- paste0("Cluster ", mfi_df[["cluster"]])
  }

  mat            <- as.matrix(mfi_df[, marker_cols, drop = FALSE])
  rownames(mat)  <- row_labels

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    p <- pheatmap::pheatmap(
      mat,
      scale        = "column",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      fontsize_row = 8L,
      fontsize_col = 8L,
      main         = "Cluster MFI Heatmap",
      silent       = TRUE
    )
    return(p)
  }

  # ggplot2 tile fallback
  n_clusters  <- nrow(mat)
  n_markers   <- length(marker_cols)
  long_df <- data.frame(
    cluster_label = rep(row_labels, times = n_markers),
    marker        = rep(marker_cols, each  = n_clusters),
    value         = as.vector(mat),
    stringsAsFactors = FALSE
  )

  # Column-scale: z-score per marker
  for (m in marker_cols) {
    idx <- long_df[["marker"]] == m
    v   <- long_df[["value"]][idx]
    mu  <- mean(v, na.rm = TRUE)
    sig <- stats::sd(v, na.rm = TRUE)
    long_df[["scaled"]][idx] <- if (!is.na(sig) && sig > 0) (v - mu) / sig else 0
  }

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x    = .data[["marker"]],
      y    = .data[["cluster_label"]],
      fill = .data[["scaled"]]
    )
  ) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Cluster MFI Heatmap",
    x     = "Marker",
    y     = "Cluster",
    fill  = "Scaled\nexpression"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45L, hjust = 1)
  )
}
