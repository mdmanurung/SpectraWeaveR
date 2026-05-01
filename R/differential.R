#' @title Differential Abundance and Expression Analysis
#'
#' @description
#' Functions for testing differential cluster abundance
#' (\code{\link{sw_differential_abundance}}) and differential marker
#' expression (\code{\link{sw_differential_expression}}) between experimental
#' conditions in spectral flow cytometry data.
#'
#' Differential abundance uses a negative-binomial GLM (\code{"edgeR"}) or a
#' voom-based linear model (\code{"voom"}) on per-sample cluster counts via the
#' \pkg{edgeR} and \pkg{limma} Bioconductor packages.
#'
#' Differential expression applies \pkg{limma}'s \code{lmFit} / \code{eBayes}
#' pipeline to per-sample-per-cluster median arcsinh expressions, treating
#' each marker as a separate feature.
#'
#' Both modules require Bioconductor packages in \code{Suggests} and guard
#' their usage with \code{requireNamespace()} checks.
#'
#' @name differential
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Check diffcyt / edgeR availability
#' @noRd
.check_diffcyt_deps <- function() {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "Package 'edgeR' is required for differential abundance analysis.\n",
      "Install with: BiocManager::install('edgeR')",
      call. = FALSE
    )
  }
}

#' Check limma availability
#' @noRd
.check_limma <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop(
      "Package 'limma' is required for differential expression analysis.\n",
      "Install with: BiocManager::install('limma')",
      call. = FALSE
    )
  }
}

#' Validate a sample-level metadata frame
#' @noRd
.validate_meta_da <- function(meta, required_cols, arg_name = "meta") {
  if (!is.data.frame(meta)) {
    stop("'", arg_name, "' must be a data.frame or tibble.", call. = FALSE)
  }
  missing_c <- setdiff(required_cols, names(meta))
  if (length(missing_c) > 0L) {
    stop(
      "'", arg_name, "' is missing required column(s): ",
      paste(missing_c, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Build a cluster-count matrix from a cell-level tibble
#'
#' Rows = clusters (sorted), columns = samples (sorted).  Returns an
#' integer matrix plus the sorted sample and cluster vectors.
#'
#' @noRd
.build_count_matrix <- function(x, cluster_col, sample_col) {
  clusters  <- sort(unique(x[[cluster_col]]))
  samples   <- sort(unique(x[[sample_col]]))

  counts <- matrix(
    0L,
    nrow     = length(clusters),
    ncol     = length(samples),
    dimnames = list(as.character(clusters), as.character(samples))
  )

  for (sn in samples) {
    sub <- x[x[[sample_col]] == sn, , drop = FALSE]
    for (cl in clusters) {
      counts[as.character(cl), as.character(sn)] <-
        sum(sub[[cluster_col]] == cl)
    }
  }

  list(counts = counts, clusters = clusters, samples = as.character(samples))
}

#' Compute per-sample-per-cluster medians
#'
#' Returns a named list (one element per cluster) of numeric matrices with
#' rows = samples, columns = markers.  Samples with fewer than
#' \code{min_cells} events for a cluster get NA rows.
#'
#' @noRd
.compute_cluster_medians <- function(corrected, cluster_col, markers,
                                     sample_col, min_cells = 1L) {
  clusters <- sort(unique(corrected[[cluster_col]]))
  samples  <- sort(unique(corrected[[sample_col]]))

  medians_list <- lapply(clusters, function(cl) {
    sub_cl <- corrected[corrected[[cluster_col]] == cl, , drop = FALSE]
    mat    <- matrix(
      NA_real_,
      nrow     = length(samples),
      ncol     = length(markers),
      dimnames = list(as.character(samples), markers)
    )
    for (sn in samples) {
      sub_s <- sub_cl[sub_cl[[sample_col]] == sn, markers, drop = FALSE]
      if (nrow(sub_s) >= min_cells) {
        mat[as.character(sn), ] <-
          apply(sub_s, 2L, stats::median, na.rm = TRUE)
      }
    }
    mat
  })
  names(medians_list) <- as.character(clusters)
  medians_list
}

# ---------------------------------------------------------------------------
# sw_differential_abundance
# ---------------------------------------------------------------------------

#' Differential Cluster Abundance Analysis
#'
#' Tests whether cluster proportions differ between experimental groups using
#' a negative-binomial GLM (\code{method = "edgeR"}) or a voom-based linear
#' model (\code{method = "voom"}).  Both methods require the \pkg{edgeR}
#' Bioconductor package; the voom method additionally requires \pkg{limma}.
#'
#' @param x A \code{data.frame}/\code{tibble} with at least a
#'   \code{"cluster"} column and a column named by \code{sample_col}.
#'   Typically the corrected tibble from \code{\link{sw_batch_correct}} or
#'   \code{\link{sw_correct_data}} with cluster assignments attached (e.g.,
#'   via \code{dplyr::mutate(corrected,
#'   cluster = sw_get_cluster_assignments(result))}).
#' @param meta A sample-level \code{data.frame}/\code{tibble} with columns
#'   \code{sample_col} and \code{group_col}.  One row per sample.
#' @param group_col Character; name of the column in \code{meta} that holds
#'   the experimental group labels (e.g., \code{"condition"}).  Must have
#'   exactly two levels for a pairwise test.
#' @param sample_col Character; name of the sample identifier column.
#'   Default: \code{"sample"}.
#' @param method Character; DA test method.  One of \code{"edgeR"} (default,
#'   negative-binomial QL F-test) or \code{"voom"} (limma-voom on library-
#'   size-normalised counts).
#' @param min_prop Numeric; minimum mean cluster proportion across all samples.
#'   Clusters below this threshold are excluded before testing.
#'   Default: \code{0.01}.
#' @param fdr_threshold Numeric; FDR significance threshold stored as an
#'   attribute on the returned tibble (used by \code{\link{sw_plot_abundance}}).
#'   Default: \code{0.05}.
#' @param ... Reserved for future use.
#'
#' @return A \code{tibble} of class \code{c("sw_da_result", "tbl_df")} with
#'   one row per (tested) cluster and columns:
#'   \describe{
#'     \item{\code{cluster}}{Cluster identifier.}
#'     \item{\code{logFC}}{log2 fold-change (group B vs group A).}
#'     \item{\code{pvalue}}{Unadjusted p-value.}
#'     \item{\code{adj_pvalue}}{BH-adjusted p-value.}
#'     \item{\code{mean_prop_A}}{Mean proportion of this cluster in group A.}
#'     \item{\code{mean_prop_B}}{Mean proportion of this cluster in group B.}
#'   }
#'   Attributes: \code{groups} (character[2]), \code{group_col} (character),
#'   \code{fdr_threshold} (numeric).
#'
#' @seealso \code{\link{sw_differential_expression}},
#'   \code{\link{sw_plot_abundance}}
#'
#' @examples
#' \dontrun{
#' # corrected: tibble from sw_batch_correct() with cluster column added
#' corrected_cl <- dplyr::mutate(
#'   corrected,
#'   cluster = sw_get_cluster_assignments(cluster_result)
#' )
#' meta <- data.frame(sample = unique(corrected_cl$sample),
#'                    condition = c("ctrl","ctrl","trt","trt","ctrl","trt"))
#' da <- sw_differential_abundance(corrected_cl, meta, group_col = "condition")
#' da[da$adj_pvalue < 0.05, ]
#' }
#'
#' @export
sw_differential_abundance <- function(x, meta, group_col,
                                       sample_col    = "sample",
                                       method        = c("edgeR", "voom"),
                                       min_prop      = 0.01,
                                       fdr_threshold = 0.05,
                                       ...) {
  method <- match.arg(method)
  .check_diffcyt_deps()
  if (method == "voom") .check_limma()

  # --- Input validation ---------------------------------------------------
  if (!is.data.frame(x)) {
    stop(
      "'x' must be a data.frame or tibble with 'cluster' and sample columns.",
      call. = FALSE
    )
  }
  if (!"cluster" %in% names(x)) {
    stop("'x' must contain a 'cluster' column.", call. = FALSE)
  }
  if (!sample_col %in% names(x)) {
    stop("'x' must contain a '", sample_col, "' column.", call. = FALSE)
  }
  .validate_meta_da(meta, c(sample_col, group_col), "meta")

  if (!is.numeric(min_prop) || length(min_prop) != 1L ||
      min_prop < 0 || min_prop >= 1) {
    stop("'min_prop' must be a number in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("'fdr_threshold' must be a number in (0, 1).", call. = FALSE)
  }

  # --- Build count matrix -------------------------------------------------
  cm       <- .build_count_matrix(x, cluster_col = "cluster",
                                   sample_col = sample_col)
  counts   <- cm$counts
  clusters <- cm$clusters
  samples  <- cm$samples

  totals <- colSums(counts)
  if (any(totals == 0L)) {
    zero_sn <- samples[totals == 0L]
    stop(
      "Sample(s) with zero events: ",
      paste(zero_sn, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Filter low-abundance clusters --------------------------------------
  props      <- sweep(counts, 2L, totals, "/")
  mean_props <- rowMeans(props)
  keep       <- mean_props >= min_prop

  if (sum(keep) == 0L) {
    stop("No clusters pass the min_prop = ", min_prop, " filter.", call. = FALSE)
  }
  counts_f <- counts[keep, , drop = FALSE]

  # --- Align metadata to column order of counts ---------------------------
  meta_ord  <- meta[match(samples, meta[[sample_col]]), , drop = FALSE]
  group_vec <- as.factor(meta_ord[[group_col]])

  if (nlevels(group_vec) != 2L) {
    stop(
      "'group_col' must have exactly 2 levels for a pairwise test; found ",
      nlevels(group_vec), ".",
      call. = FALSE
    )
  }
  groups <- levels(group_vec)

  design <- stats::model.matrix(~ group_vec)

  # --- Run test -----------------------------------------------------------
  if (method == "edgeR") {
    dge <- edgeR::DGEList(counts = counts_f, lib.size = totals)
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmQLFit(dge, design = design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2L)
    tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table

    cl_names <- rownames(tab)
    pvals    <- tab[["PValue"]]
    lfc      <- tab[["logFC"]]
  } else {
    # voom
    dge  <- edgeR::DGEList(counts = counts_f, lib.size = totals)
    dge  <- edgeR::calcNormFactors(dge)
    v    <- limma::voom(dge, design = design)
    fit  <- limma::lmFit(v, design = design)
    fit  <- limma::eBayes(fit)
    tab  <- limma::topTable(fit, coef = 2L, n = Inf, sort.by = "none")

    cl_names <- rownames(tab)
    pvals    <- tab[["P.Value"]]
    lfc      <- tab[["logFC"]]
  }

  # Compute group mean proportions (using filtered clusters only)
  prop_f   <- sweep(counts_f[cl_names, , drop = FALSE], 2L, totals, "/")
  grp_A    <- group_vec == groups[1L]
  grp_B    <- group_vec == groups[2L]
  mean_A   <- rowMeans(prop_f[, grp_A, drop = FALSE])
  mean_B   <- rowMeans(prop_f[, grp_B, drop = FALSE])

  result <- tibble::tibble(
    cluster    = cl_names,
    logFC      = lfc,
    pvalue     = pvals,
    adj_pvalue = stats::p.adjust(pvals, method = "BH"),
    mean_prop_A = mean_A,
    mean_prop_B = mean_B
  )

  attr(result, "groups")        <- groups
  attr(result, "group_col")     <- group_col
  attr(result, "fdr_threshold") <- fdr_threshold
  class(result) <- c("sw_da_result", class(result))
  result
}

# ---------------------------------------------------------------------------
# sw_differential_expression
# ---------------------------------------------------------------------------

#' Differential Marker Expression Analysis
#'
#' Tests whether per-cluster marker expression differs between experimental
#' groups by applying \pkg{limma}'s \code{lmFit}/\code{eBayes} pipeline to
#' per-sample-per-cluster median arcsinh expressions.
#'
#' Each cluster is tested independently.  For each cluster, a per-sample
#' median expression matrix (samples \eqn{\times} markers) is constructed
#' using only samples with at least \code{min_cells} events.  The matrix is
#' then passed directly to limma without further transformation (the arcsinh
#' scale with cofactor 6000 is retained).
#'
#' @param corrected A \code{data.frame}/\code{tibble} containing per-cell
#'   marker expressions, a cluster identifier column (\code{cluster_col}),
#'   and a sample identifier column (\code{sample_col}).  Typically the output
#'   of \code{\link{sw_batch_correct}} or \code{\link{sw_correct_data}} with
#'   cluster assignments added.
#' @param meta A sample-level \code{data.frame}/\code{tibble} with columns
#'   \code{sample_col} and \code{group_col}.
#' @param group_col Character; name of the group column in \code{meta}.
#' @param sample_col Character; name of the sample column.  Default:
#'   \code{"sample"}.
#' @param cluster_col Character; name of the cluster column.  Default:
#'   \code{"cluster"}.
#' @param markers Character vector of marker names to test.  If \code{NULL}
#'   (default), all numeric columns in \code{corrected} that are not metadata
#'   columns (\code{batch}, \code{condition}, \code{.cell_id}, \code{label})
#'   are used.
#' @param min_cells Integer; minimum number of events a sample must contribute
#'   to a cluster to be included in the per-cluster test.  Default: \code{3}.
#' @param fdr_threshold Numeric; FDR threshold stored as an attribute on the
#'   result (used by \code{\link{sw_plot_volcano}}).  Default: \code{0.05}.
#' @param ... Reserved for future use.
#'
#' @return A \code{tibble} of class \code{c("sw_de_result", "tbl_df")} in
#'   long format with one row per cluster-marker pair and columns:
#'   \describe{
#'     \item{\code{cluster}}{Cluster identifier.}
#'     \item{\code{marker}}{Marker name.}
#'     \item{\code{logFC}}{log2 fold-change (group B vs group A).}
#'     \item{\code{pvalue}}{Unadjusted p-value (from limma eBayes).}
#'     \item{\code{adj_pvalue}}{BH-adjusted p-value (within cluster).}
#'     \item{\code{mean_A}}{Mean arcsinh expression in group A.}
#'     \item{\code{mean_B}}{Mean arcsinh expression in group B.}
#'   }
#'   Attributes: \code{groups} (character[2]), \code{group_col} (character),
#'   \code{fdr_threshold} (numeric).
#'
#' @seealso \code{\link{sw_differential_abundance}},
#'   \code{\link{sw_plot_volcano}}, \code{\link{sw_plot_boxplots}}
#'
#' @examples
#' \dontrun{
#' corrected_cl <- dplyr::mutate(
#'   corrected,
#'   cluster = sw_get_cluster_assignments(cluster_result)
#' )
#' meta <- data.frame(sample = unique(corrected_cl$sample),
#'                    condition = c("ctrl","ctrl","trt","trt","ctrl","trt"))
#' de <- sw_differential_expression(corrected_cl, meta,
#'                                   group_col = "condition",
#'                                   markers   = lineage_markers)
#' de[de$adj_pvalue < 0.05, ]
#' }
#'
#' @export
sw_differential_expression <- function(corrected, meta, group_col,
                                        sample_col  = "sample",
                                        cluster_col = "cluster",
                                        markers     = NULL,
                                        min_cells   = 3L,
                                        fdr_threshold = 0.05,
                                        ...) {
  .check_limma()

  # --- Input validation ---------------------------------------------------
  if (!is.data.frame(corrected)) {
    stop("'corrected' must be a data.frame or tibble.", call. = FALSE)
  }
  if (!cluster_col %in% names(corrected)) {
    stop("'corrected' must contain a '", cluster_col, "' column.", call. = FALSE)
  }
  if (!sample_col %in% names(corrected)) {
    stop("'corrected' must contain a '", sample_col, "' column.", call. = FALSE)
  }
  .validate_meta_da(meta, c(sample_col, group_col), "meta")

  if (!is.numeric(min_cells) || length(min_cells) != 1L || min_cells < 1) {
    stop("'min_cells' must be a positive integer.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)

  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("'fdr_threshold' must be a number in (0, 1).", call. = FALSE)
  }

  # --- Determine markers --------------------------------------------------
  if (is.null(markers)) {
    skip_cols <- c(cluster_col, sample_col, "batch", "condition",
                   ".cell_id", "label")
    candidates <- setdiff(names(corrected), skip_cols)
    markers <- candidates[
      vapply(corrected[candidates], is.numeric, logical(1))
    ]
  } else {
    missing_m <- setdiff(markers, names(corrected))
    if (length(missing_m) > 0L) {
      stop(
        "Markers not found in 'corrected': ",
        paste(missing_m, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (length(markers) == 0L) {
    stop("No numeric marker columns found in 'corrected'.", call. = FALSE)
  }

  # --- Align meta ---------------------------------------------------------
  samples  <- sort(unique(corrected[[sample_col]]))
  clusters <- sort(unique(corrected[[cluster_col]]))

  meta_ord  <- meta[match(samples, meta[[sample_col]]), , drop = FALSE]
  group_vec <- as.factor(meta_ord[[group_col]])

  if (nlevels(group_vec) != 2L) {
    stop(
      "'group_col' must have exactly 2 levels; found ",
      nlevels(group_vec), ".",
      call. = FALSE
    )
  }
  groups <- levels(group_vec)

  # --- Per-cluster limma tests --------------------------------------------
  all_results <- list()

  for (cl in clusters) {
    sub_cl <- corrected[corrected[[cluster_col]] == cl, , drop = FALSE]

    # Build per-sample median matrix
    med_mat <- matrix(
      NA_real_,
      nrow     = length(samples),
      ncol     = length(markers),
      dimnames = list(as.character(samples), markers)
    )
    valid_samples <- character(0L)

    for (sn in as.character(samples)) {
      sub_s <- sub_cl[sub_cl[[sample_col]] == sn, markers, drop = FALSE]
      if (nrow(sub_s) >= min_cells) {
        med_mat[sn, ]  <- apply(sub_s, 2L, stats::median, na.rm = TRUE)
        valid_samples  <- c(valid_samples, sn)
      }
    }

    # Need at least (n_groups + 1) valid samples for a rank-sufficient design
    if (length(valid_samples) < (nlevels(group_vec) + 1L)) next

    med_sub    <- med_mat[valid_samples, , drop = FALSE]
    grp_sub    <- group_vec[match(valid_samples, as.character(samples))]
    design_sub <- stats::model.matrix(~ grp_sub)

    # Skip if design is rank-deficient (e.g., all valid samples are one group)
    if (qr(design_sub)$rank < ncol(design_sub)) next

    # limma: rows = markers (features), cols = samples
    fit <- limma::lmFit(t(med_sub), design_sub)
    fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef = 2L, n = Inf,
                            sort.by = "none", adjust.method = "BH")

    # Group means
    grp_A_idx <- grp_sub == groups[1L]
    grp_B_idx <- grp_sub == groups[2L]

    mean_A <- if (sum(grp_A_idx) > 0L)
      colMeans(med_sub[grp_A_idx, , drop = FALSE], na.rm = TRUE)
    else
      rep(NA_real_, length(markers))

    mean_B <- if (sum(grp_B_idx) > 0L)
      colMeans(med_sub[grp_B_idx, , drop = FALSE], na.rm = TRUE)
    else
      rep(NA_real_, length(markers))

    all_results[[as.character(cl)]] <- tibble::tibble(
      cluster    = cl,
      marker     = rownames(tab),
      logFC      = tab[["logFC"]],
      pvalue     = tab[["P.Value"]],
      adj_pvalue = tab[["adj.P.Val"]],
      mean_A     = mean_A[rownames(tab)],
      mean_B     = mean_B[rownames(tab)]
    )
  }

  if (length(all_results) == 0L) {
    stop(
      "No clusters had enough samples with sufficient cells ",
      "(min_cells = ", min_cells, ") for testing.",
      call. = FALSE
    )
  }

  result <- do.call(rbind, all_results)

  attr(result, "groups")        <- groups
  attr(result, "group_col")     <- group_col
  attr(result, "fdr_threshold") <- fdr_threshold
  class(result) <- c("sw_de_result", class(result))
  result
}

# ---------------------------------------------------------------------------
# sw_plot_volcano
# ---------------------------------------------------------------------------

#' Volcano Plot for Differential Expression
#'
#' Visualizes the results of \code{\link{sw_differential_expression}} as a
#' volcano plot (log2 fold-change on the x-axis, \eqn{-\log_{10}} adjusted
#' p-value on the y-axis).  Significant markers are coloured; top hits are
#' labelled.  When multiple clusters are present the plot is faceted.
#'
#' @param de_result Output of \code{\link{sw_differential_expression}}: a
#'   \code{data.frame}/\code{tibble} with columns \code{cluster},
#'   \code{marker}, \code{logFC}, and \code{adj_pvalue}.
#' @param cluster Character or integer vector of cluster identifier(s) to
#'   plot.  \code{NULL} (default) plots all clusters, faceted.
#' @param fdr_threshold Numeric; significance threshold for colouring points.
#'   Default: \code{0.05}.
#' @param fc_threshold Numeric; absolute log2 FC threshold.  Default:
#'   \code{0.5}.
#' @param label_top Integer; number of top significant markers to label.
#'   Default: \code{10}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @seealso \code{\link{sw_differential_expression}},
#'   \code{\link{sw_plot_boxplots}}
#'
#' @examples
#' \dontrun{
#' de <- sw_differential_expression(corrected_cl, meta,
#'                                   group_col = "condition")
#' sw_plot_volcano(de)
#' # Single cluster only:
#' sw_plot_volcano(de, cluster = 1)
#' }
#'
#' @export
sw_plot_volcano <- function(de_result, cluster = NULL,
                             fdr_threshold = 0.05,
                             fc_threshold  = 0.5,
                             label_top     = 10L) {
  if (!is.data.frame(de_result)) {
    stop(
      "'de_result' must be the output of sw_differential_expression().",
      call. = FALSE
    )
  }
  req_cols <- c("cluster", "marker", "logFC", "pvalue", "adj_pvalue")
  missing_c <- setdiff(req_cols, names(de_result))
  if (length(missing_c) > 0L) {
    stop(
      "'de_result' must have columns: ",
      paste(req_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for sw_plot_volcano().\n",
      "Install with: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  if (!is.null(cluster)) {
    de_result <- de_result[de_result[["cluster"]] %in% cluster, , drop = FALSE]
    if (nrow(de_result) == 0L) {
      stop(
        "No rows found for cluster(s): ",
        paste(cluster, collapse = ", "),
        call. = FALSE
      )
    }
  }

  de_result[["neg_log10_p"]] <-
    -log10(de_result[["adj_pvalue"]] + .Machine$double.eps)
  de_result[["significant"]]  <-
    de_result[["adj_pvalue"]] < fdr_threshold &
    abs(de_result[["logFC"]]) > fc_threshold

  de_result[["direction"]] <- ifelse(
    de_result[["significant"]] & de_result[["logFC"]] > 0, "up",
    ifelse(
      de_result[["significant"]] & de_result[["logFC"]] < 0, "down",
      "ns"
    )
  )

  # Build label set: top N significant markers by adjusted p-value
  sig_df   <- de_result[de_result[["significant"]], , drop = FALSE]
  sig_df   <- sig_df[order(sig_df[["neg_log10_p"]], decreasing = TRUE), ]
  n_top    <- min(as.integer(label_top), nrow(sig_df))
  top_labs <- if (n_top > 0L) utils::head(sig_df, n_top) else sig_df[0L, ]

  colour_vals <- c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "grey70")

  n_clusters <- length(unique(de_result[["cluster"]]))
  ncol_wrap  <- min(n_clusters, 3L)

  p <- ggplot2::ggplot(
    de_result,
    ggplot2::aes(
      x      = .data[["logFC"]],
      y      = .data[["neg_log10_p"]],
      colour = .data[["direction"]]
    )
  ) +
  ggplot2::geom_point(size = 2L, alpha = 0.8) +
  ggplot2::scale_colour_manual(
    values = colour_vals,
    labels = c("up" = "Up", "down" = "Down", "ns" = "Not sig."),
    name   = NULL
  ) +
  ggplot2::geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    linetype   = "dashed",
    colour     = "grey50"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype   = "dashed",
    colour     = "grey50"
  ) +
  ggplot2::labs(
    x = expression(log[2] ~ FC),
    y = expression(-log[10] ~ adj. ~ p)
  ) +
  ggplot2::theme_bw()

  if (nrow(top_labs) > 0L) {
    p <- p + ggplot2::geom_text(
      data    = top_labs,
      mapping = ggplot2::aes(label = .data[["marker"]]),
      vjust   = -0.5,
      size    = 3L,
      inherit.aes = TRUE
    )
  }

  if (n_clusters > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~ cluster,
      ncol       = ncol_wrap,
      labeller   = ggplot2::label_both
    )
  } else {
    p <- p + ggplot2::ggtitle(paste("Cluster", unique(de_result[["cluster"]])))
  }

  p
}

# ---------------------------------------------------------------------------
# sw_plot_boxplots
# ---------------------------------------------------------------------------

#' Boxplots of Marker Expression by Group
#'
#' Produces violin + boxplots for one or more markers, faceted by marker and
#' coloured by experimental group.  Useful for visualising expression
#' differences identified by \code{\link{sw_differential_expression}}.
#'
#' @param corrected A \code{data.frame}/\code{tibble} with per-cell marker
#'   expressions, a group column, and optionally a cluster column.
#' @param markers Character vector of marker names to plot.
#' @param group_col Character; name of the group column in \code{corrected}.
#' @param sample_col Character; name of the sample column.  Default:
#'   \code{"sample"}.
#' @param cluster_col Character; name of the cluster column.  Default:
#'   \code{"cluster"}.
#' @param cluster Optional cluster identifier(s) to subset before plotting.
#'   \code{NULL} (default) plots all cells.
#' @param ... Reserved for future use.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @seealso \code{\link{sw_differential_expression}},
#'   \code{\link{sw_plot_volcano}}
#'
#' @examples
#' \dontrun{
#' # Plot CD4 and CD8a expression by condition in cluster 1
#' sw_plot_boxplots(corrected_cl, markers = c("CD4", "CD8a"),
#'                  group_col = "condition", cluster = 1)
#' }
#'
#' @export
sw_plot_boxplots <- function(corrected, markers, group_col,
                              sample_col  = "sample",
                              cluster_col = "cluster",
                              cluster     = NULL,
                              ...) {
  if (!is.data.frame(corrected)) {
    stop("'corrected' must be a data.frame or tibble.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for sw_plot_boxplots().\n",
      "Install with: install.packages('ggplot2')",
      call. = FALSE
    )
  }
  if (!is.character(markers) || length(markers) == 0L) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }
  missing_m <- setdiff(markers, names(corrected))
  if (length(missing_m) > 0L) {
    stop(
      "Markers not found in 'corrected': ",
      paste(missing_m, collapse = ", "),
      call. = FALSE
    )
  }
  if (!group_col %in% names(corrected)) {
    stop("Column '", group_col, "' not found in 'corrected'.", call. = FALSE)
  }
  if (!cluster_col %in% names(corrected)) {
    stop("Column '", cluster_col, "' not found in 'corrected'.", call. = FALSE)
  }

  plot_df <- corrected
  if (!is.null(cluster)) {
    plot_df <- plot_df[plot_df[[cluster_col]] %in% cluster, , drop = FALSE]
    if (nrow(plot_df) == 0L) {
      stop(
        "No rows found for cluster(s): ",
        paste(cluster, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # Pivot to long format without tidyr dependency
  n_rows <- nrow(plot_df)
  long_df <- do.call(rbind, lapply(markers, function(m) {
    data.frame(
      group_var = plot_df[[group_col]],
      marker    = m,
      value     = plot_df[[m]],
      stringsAsFactors = FALSE
    )
  }))

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x    = .data[["group_var"]],
      y    = .data[["value"]],
      fill = .data[["group_var"]]
    )
  ) +
  ggplot2::geom_violin(alpha = 0.6, trim = TRUE) +
  ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  ggplot2::facet_wrap(
    ~ marker,
    scales = "free_y",
    ncol   = min(length(markers), 4L)
  ) +
  ggplot2::labs(x = group_col, y = "Expression", fill = group_col) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(angle = 45L, hjust = 1),
    legend.position = "none"
  )
}

# ---------------------------------------------------------------------------
# sw_plot_abundance
# ---------------------------------------------------------------------------

#' Bar Chart of Differential Cluster Abundance
#'
#' Visualizes the output of \code{\link{sw_differential_abundance}} as either
#' a grouped bar chart of mean cluster proportions or a small-multiples box
#' plot (one panel per cluster) showing mean proportions by group.
#'
#' @param da_result Output of \code{\link{sw_differential_abundance}}: a
#'   \code{data.frame}/\code{tibble} with columns \code{cluster},
#'   \code{mean_prop_A}, and \code{mean_prop_B}.
#' @param type Character; \code{"bar"} (default) for a grouped bar chart or
#'   \code{"box"} for faceted small multiples.
#' @param ... Reserved for future use.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @seealso \code{\link{sw_differential_abundance}}
#'
#' @examples
#' \dontrun{
#' da <- sw_differential_abundance(corrected_cl, meta, group_col = "condition")
#' sw_plot_abundance(da)
#' sw_plot_abundance(da, type = "box")
#' }
#'
#' @export
sw_plot_abundance <- function(da_result, type = c("bar", "box"), ...) {
  type <- match.arg(type)

  if (!is.data.frame(da_result)) {
    stop(
      "'da_result' must be the output of sw_differential_abundance().",
      call. = FALSE
    )
  }
  req_cols <- c("cluster", "mean_prop_A", "mean_prop_B")
  missing_c <- setdiff(req_cols, names(da_result))
  if (length(missing_c) > 0L) {
    stop(
      "'da_result' must have columns: ",
      paste(req_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for sw_plot_abundance().\n",
      "Install with: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  groups  <- attr(da_result, "groups")
  if (is.null(groups) || length(groups) < 2L) groups <- c("Group A", "Group B")
  fdr_thr <- attr(da_result, "fdr_threshold")
  if (is.null(fdr_thr)) fdr_thr <- 0.05

  # Pivot to long
  long_df <- rbind(
    data.frame(
      cluster = da_result[["cluster"]],
      group   = groups[1L],
      prop    = da_result[["mean_prop_A"]],
      sig     = da_result[["adj_pvalue"]] < fdr_thr,
      stringsAsFactors = FALSE
    ),
    data.frame(
      cluster = da_result[["cluster"]],
      group   = groups[2L],
      prop    = da_result[["mean_prop_B"]],
      sig     = da_result[["adj_pvalue"]] < fdr_thr,
      stringsAsFactors = FALSE
    )
  )
  long_df[["group"]] <- factor(long_df[["group"]], levels = groups)

  if (type == "bar") {
    p <- ggplot2::ggplot(
      long_df,
      ggplot2::aes(
        x    = .data[["cluster"]],
        y    = .data[["prop"]],
        fill = .data[["group"]]
      )
    ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(x = "Cluster", y = "Mean proportion", fill = "Group") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45L, hjust = 1)
    )
  } else {
    # type == "box": faceted, one panel per cluster
    p <- ggplot2::ggplot(
      long_df,
      ggplot2::aes(
        x    = .data[["group"]],
        y    = .data[["prop"]],
        fill = .data[["group"]]
      )
    ) +
    ggplot2::geom_col(alpha = 0.7) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::labs(x = "Group", y = "Mean proportion", fill = "Group") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 45L, hjust = 1)
    )
  }

  p
}
