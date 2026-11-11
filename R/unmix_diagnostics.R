#' @title Unmixing Quality Diagnostics
#'
#' @description
#' Functions for evaluating spectral unmixing quality.  The spillover
#' spreading matrix (SSM) quantifies how much uncertainty each fluorochrome
#' introduces into other channels due to spectral overlap.  Per-channel
#' quality summaries flag markers with high coefficient of variation.
#'
#' @name unmix_diagnostics
#' @keywords internal
NULL


# ---------------------------------------------------------------------------
# Spillover Spreading Matrix
# ---------------------------------------------------------------------------

#' Cosine Similarity Matrix of Reference Spectra
#'
#' Computes a fluorophore-by-fluorophore matrix of absolute cosine similarity
#' between the reference spectra used for unmixing.  High values mark pairs
#' with similar spectra, which are harder to unmix.
#'
#' This is \strong{not} a spillover spreading matrix (SSM).  A true SSM
#' measures the spread in a channel caused by another fluorophore, and it is
#' computed from unmixed single-stain controls.  This function uses only the
#' reference spectra, so it cannot show spread.  For a real SSM, use the one
#' AutoSpectral writes during unmixing.  The unmix report and dashboard
#' (\code{\link{sw_unmix_report}()}, \code{\link{sw_unmix_dashboard}()}) read
#' it from the run folder.  The name \code{sw_unmix_spillover_matrix} is kept
#' for backward compatibility.
#'
#' @param reference_spectra A numeric matrix of reference spectra with
#'   fluorophore names as row names and detector names as column names.
#'   Each row is the emission spectrum of one fluorophore (as returned by
#'   \code{\link{sw_unmix_prepare}()$spectra}).
#' @param unmixed_ff An optional \code{flowFrame} of unmixed data.  When
#'   provided, the cosine similarity is averaged with the absolute Pearson
#'   correlation of the unmixed fluorophore channels.  When \code{NULL}
#'   (default), only the cosine similarity of the spectra is returned.
#'
#' @return A named list with class \code{"sw_ssm"} (the class name and the
#'   \code{spreading} column names say "SSM" and "spreading" but hold cosine
#'   similarity, not spillover spreading) containing:
#'   \describe{
#'     \item{\code{matrix}}{Numeric matrix (fluorophores x fluorophores) of
#'       absolute cosine similarity.  Diagonal entries are zero.}
#'     \item{\code{summary}}{A \code{tibble} with per-fluorophore
#'       \code{max_spreading}, \code{mean_spreading}, and
#'       \code{worst_partner} (the fluorophore with the highest similarity).
#'       The first two are the row maximum and mean of \code{matrix}.}
#'   }
#'
#' @details
#' The spectra-based matrix is the absolute normalised dot product between
#' reference spectra:
#'
#' \deqn{S_{ij} = \frac{|R_i \cdot R_j|}{\|R_i\| \, \|R_j\|}}
#'
#' where \eqn{R_i} is the reference spectrum (row vector) of fluorophore
#' \eqn{i}.
#' This is equivalent to the absolute cosine similarity and ranges from 0
#' (orthogonal spectra) to 1 (identical spectra).
#'
#' @seealso \code{\link{sw_plot_spillover_matrix}}, \code{\link{sw_unmix_quality}},
#'   \code{\link{sw_unmix_prepare}}
#'
#' @examples
#' # Create toy reference spectra (3 fluorophores, 5 detectors)
#' spectra <- matrix(c(1, 0.1, 0, 0, 0,
#'                     0, 1, 0.2, 0, 0,
#'                     0, 0, 0.1, 1, 0.3),
#'                   nrow = 3, byrow = TRUE,
#'                   dimnames = list(c("BV421", "PE", "APC"),
#'                                   paste0("Det", 1:5)))
#' ssm <- sw_unmix_spillover_matrix(spectra)
#' ssm$summary
#'
#' @export
sw_unmix_spillover_matrix <- function(reference_spectra,
                                          unmixed_ff = NULL) {
  # --- Validate reference_spectra ---
  if (!is.matrix(reference_spectra) || !is.numeric(reference_spectra)) {
    stop("'reference_spectra' must be a numeric matrix.", call. = FALSE)
  }

  if (is.null(rownames(reference_spectra))) {
    stop("'reference_spectra' must have row names (fluorophore names).",
         call. = FALSE)
  }

  if (nrow(reference_spectra) < 2) {
    stop("'reference_spectra' must have at least 2 rows (fluorophores).",
         call. = FALSE)
  }

  # --- Validate optional unmixed_ff ---
  if (!is.null(unmixed_ff)) {
    if (!requireNamespace("flowCore", quietly = TRUE)) {
      stop("Package 'flowCore' is required when 'unmixed_ff' is provided.",
           call. = FALSE)
    }
    if (!methods::is(unmixed_ff, "flowFrame")) {
      stop("'unmixed_ff' must be a flowFrame object.", call. = FALSE)
    }
  }

  fluors <- rownames(reference_spectra)
  n_fluors <- length(fluors)

  # --- Compute spectra-based SSM (absolute cosine similarity) ---
  # Normalise each spectrum to unit length
  norms <- sqrt(rowSums(reference_spectra^2))
  norms[norms == 0] <- 1  # Avoid division by zero

  R_norm <- reference_spectra / norms

  # Cosine similarity matrix: R_norm %*% t(R_norm)
  ssm <- abs(R_norm %*% t(R_norm))

  # Zero the diagonal (no self-spreading)
  diag(ssm) <- 0

  rownames(ssm) <- fluors
  colnames(ssm) <- fluors

  # --- Empirical adjustment from unmixed data (optional) ---
  if (!is.null(unmixed_ff)) {
    expr_mat <- flowCore::exprs(unmixed_ff)
    # Find matching fluorophore columns
    matched <- intersect(fluors, colnames(expr_mat))
    if (length(matched) >= 2) {
      sub_mat <- expr_mat[, matched, drop = FALSE]
      # Compute correlation-based spreading on unmixed data
      cor_mat <- stats::cor(sub_mat, use = "pairwise.complete.obs")
      cor_mat[is.na(cor_mat)] <- 0
      emp_ssm <- abs(cor_mat)
      diag(emp_ssm) <- 0

      # Average spectra-based and empirical for matched fluorophores
      for (fi in matched) {
        for (fj in matched) {
          if (fi != fj) {
            ssm[fi, fj] <- (ssm[fi, fj] + emp_ssm[fi, fj]) / 2
          }
        }
      }
    }
  }

  # --- Summary tibble ---
  summary_df <- tibble::tibble(
    fluorophore = fluors,
    max_spreading = apply(ssm, 1, max),
    mean_spreading = rowMeans(ssm),
    worst_partner = fluors[apply(ssm, 1, which.max)]
  )

  structure(
    list(
      matrix = ssm,
      summary = summary_df
    ),
    class = "sw_ssm"
  )
}


# ---------------------------------------------------------------------------
# Plot SSM
# ---------------------------------------------------------------------------

#' Plot Reference-Spectra Similarity Matrix
#'
#' Produces a heatmap of the \code{"sw_ssm"} object.  Its values are cosine
#' similarity, not spillover spreading (see
#' \code{\link{sw_unmix_spillover_matrix}}).  Uses
#' \pkg{pheatmap} if available, otherwise falls back to
#' \code{stats::heatmap()}.
#'
#' @param ssm An object of class \code{"sw_ssm"} as returned by
#'   \code{\link{sw_unmix_spillover_matrix}}.
#' @param plot_file Optional file path to save the plot as PDF.
#' @param ... Additional arguments passed to the underlying heatmap
#'   function.
#'
#' @return Invisible \code{NULL}; called for its side effect (plot).
#'
#' @seealso \code{\link{sw_unmix_spillover_matrix}}
#'
#' @export
sw_plot_spillover_matrix <- function(ssm, plot_file = NULL, ...) {
  if (!inherits(ssm, "sw_ssm")) {
    stop("'ssm' must be an object of class 'sw_ssm' ",
         "(from sw_unmix_spillover_matrix).", call. = FALSE)
  }

  mat <- ssm$matrix

  # Open PDF device if requested
  if (!is.null(plot_file)) {
    grDevices::pdf(plot_file, width = 10, height = 9)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
      mat,
      color = grDevices::colorRampPalette(
        c("white", "khaki1", "orange", "red3")
      )(100),
      main = "Spillover Spreading Matrix",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize_number = 7,
      ...
    )
  } else {
    stats::heatmap(mat, Rowv = NA, Colv = NA, scale = "none",
                   col = grDevices::heat.colors(50),
                   main = "Spillover Spreading Matrix", ...)
  }

  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Unmixing quality summary
# ---------------------------------------------------------------------------

#' Unmixing Quality Summary
#'
#' Computes per-channel quality metrics for unmixed flow cytometry data.
#' Reports the mean, standard deviation, coefficient of variation (CV),
#' and flags channels exceeding a CV threshold.
#'
#' @param unmixed_ff A \code{flowFrame} of unmixed data.
#' @param channels Character vector of channel names to evaluate.  If
#'   \code{NULL} (default), all fluorochrome channels (excluding scatter
#'   and metadata) are used.
#' @param cv_threshold Numeric; channels with CV exceeding this value
#'   are flagged (default: 0.5).
#'
#' @return A \code{tibble} with columns: \code{channel}, \code{mean},
#'   \code{sd}, \code{cv}, \code{pct_negative}, and \code{flagged}.
#'
#' @seealso \code{\link{sw_unmix_spillover_matrix}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   mat <- matrix(abs(rnorm(300, 5000, 2000)), ncol = 3,
#'                 dimnames = list(NULL, c("CD3", "CD4", "CD8")))
#'   ff <- flowCore::flowFrame(mat)
#'   quality <- sw_unmix_quality(ff, channels = c("CD3", "CD4", "CD8"))
#'   quality
#' }
#' }
#'
#' @export
sw_unmix_quality <- function(unmixed_ff, channels = NULL,
                                cv_threshold = 0.5) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(unmixed_ff, "flowFrame")) {
    stop("'unmixed_ff' must be a flowFrame object.", call. = FALSE)
  }

  if (!is.numeric(cv_threshold) || cv_threshold <= 0) {
    stop("'cv_threshold' must be a positive number.", call. = FALSE)
  }

  expr_mat <- flowCore::exprs(unmixed_ff)

  if (is.null(channels)) {
    fluor_mask <- sw_channel_is_fluor(unmixed_ff)
    channels <- flowCore::colnames(unmixed_ff)[fluor_mask]
  }

  missing_ch <- setdiff(channels, colnames(expr_mat))
  if (length(missing_ch) > 0) {
    stop("Channel(s) not found: ",
         paste(missing_ch, collapse = ", "), call. = FALSE)
  }

  stats_list <- lapply(channels, function(ch) {
    vals <- expr_mat[, ch]
    ch_mean <- mean(vals, na.rm = TRUE)
    ch_sd <- stats::sd(vals, na.rm = TRUE)
    ch_cv <- if (abs(ch_mean) > .Machine$double.eps) {
      ch_sd / abs(ch_mean)
    } else {
      NA_real_
    }
    pct_neg <- 100 * sum(vals < 0, na.rm = TRUE) / length(vals)

    tibble::tibble(
      channel = ch,
      mean = ch_mean,
      sd = ch_sd,
      cv = ch_cv,
      pct_negative = pct_neg,
      flagged = !is.na(ch_cv) && ch_cv > cv_threshold
    )
  })

  result <- dplyr::bind_rows(stats_list)

  n_flagged <- sum(result$flagged)
  if (n_flagged > 0) {
    message("Unmixing quality: ", n_flagged, " channel(s) flagged (CV > ",
            cv_threshold, "): ",
            paste(result$channel[result$flagged], collapse = ", "))
  } else {
    message("Unmixing quality: all channels within CV threshold (",
            cv_threshold, ")")
  }

  result
}
