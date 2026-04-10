#' @title Scale Transformation Utilities
#'
#' @description
#' Functions for estimating and applying scale transformations to flow
#' cytometry data. Ported from CytoPipeline (UCLouvain-CBIO) and adapted
#' for SpectraWeaveR's spectral flow cytometry workflows.
#'
#' @name transforms
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# estimateScaleTransforms
# ---------------------------------------------------------------------------

#' Estimate Scale Transformations
#'
#' Estimates per-channel scale transformations to obtain well-separated
#' positive and negative populations. Fluorochrome channels can be transformed
#' with logicle (via \code{flowCore::estimateLogicle}), while scatter channels
#' can be linearly scaled to match the dynamic range of a reference
#' fluorochrome channel.
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param ff A \code{flowCore::flowFrame}.
#' @param fluo_method Character; method for fluorochrome channels.
#'   One of \code{"estimateLogicle"} (default) or \code{"none"}.
#' @param scatter_method Character; method for scatter channels.
#'   One of \code{"none"} (default) or \code{"linearQuantile"}.
#' @param scatter_ref_marker Character; reference fluorochrome channel or marker
#'   name whose 5th/95th percentiles define the target range for linear scatter
#'   scaling. Required when \code{scatter_method = "linearQuantile"}.
#' @param specific_scatter_channels Character vector of scatter channel names
#'   that should receive the fluorochrome method instead of the scatter method.
#'   Default: \code{NULL}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'   Default: \code{FALSE}.
#' @param ... Additional arguments passed to
#'   \code{flowCore::estimateLogicle()}.
#'
#' @return A \code{flowCore::transformList} object.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   mat <- matrix(abs(rnorm(500, 50000, 15000)), ncol = 5,
#'                 dimnames = list(NULL, c("FSC-A", "SSC-A", "CD3", "CD4", "CD8")))
#'   ff <- flowCore::flowFrame(mat)
#'   trans <- sw_estimate_scale_transforms(ff, fluo_method = "estimateLogicle")
#' }
#' }
#'
#' @export
sw_estimate_scale_transforms <- function(
    ff,
    fluo_method = c("estimateLogicle", "none"),
    scatter_method = c("none", "linearQuantile"),
    scatter_ref_marker = NULL,
    specific_scatter_channels = NULL,
    verbose = FALSE,
    ...) {

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  fluo_method <- match.arg(fluo_method)
  scatter_method <- match.arg(scatter_method)

  transList <- NULL

  # --- Fluorochrome channels ---
  if (fluo_method == "estimateLogicle") {
    if (verbose) {
      message("Estimating logicle transformations for fluorochrome channels...")
    }
    fluoCols <- flowCore::colnames(ff)[sw_are_fluor_cols(ff)]
    if (length(fluoCols) > 0) {
      transList <- flowCore::estimateLogicle(ff, fluoCols, ...)
    }
  }

  # --- Scatter channels ---
  if (scatter_method == "linearQuantile") {
    if (is.null(scatter_ref_marker)) {
      stop("'scatter_ref_marker' is required when ",
           "scatter_method = 'linearQuantile'.", call. = FALSE)
    }

    if (verbose) {
      message("Estimating linear transformation for scatter channels; ",
              "reference marker = ", scatter_ref_marker, "...")
    }

    transList <- .compute_scatter_linear_scale(
      ff,
      transList = transList,
      reference_channel = scatter_ref_marker,
      verbose = verbose
    )
  }

  # --- Specific scatter channels to treat like fluo ---
  if (!is.null(specific_scatter_channels) &&
      fluo_method == "estimateLogicle") {
    scatter_channels <- flowCore::colnames(ff)[
      !sw_are_fluor_cols(ff) & sw_are_signal_cols(ff)
    ]
    effective <- intersect(specific_scatter_channels, scatter_channels)

    if (length(effective) > 0) {
      # Remove existing scatter transforms for these channels
      for (ch in effective) {
        if (!is.null(transList)) {
          transList@transforms[[ch]] <- NULL
        }
      }
      # Add logicle transforms instead
      logicle_trans <- flowCore::estimateLogicle(ff, effective, ...)
      if (is.null(transList)) {
        transList <- logicle_trans
      } else {
        transList <- c(transList, logicle_trans)
      }

      if (verbose) {
        message("Applied logicle to specific scatter channels: ",
                paste(effective, collapse = ", "))
      }
    }
  }

  if (is.null(transList)) {
    stop("No transformations were estimated. Check fluo_method and ",
         "scatter_method parameters.", call. = FALSE)
  }

  transList
}

# ---------------------------------------------------------------------------
# Internal: compute linear scale for scatter channels
# ---------------------------------------------------------------------------

#' @keywords internal
.compute_scatter_linear_scale <- function(
    ff,
    transList = NULL,
    reference_channel,
    verbose = FALSE) {

  # Resolve reference channel (could be marker name)
  ref_info <- flowCore::getChannelMarker(ff, reference_channel)
  reference_channel <- ref_info$name

  fluoCols <- flowCore::colnames(ff)[sw_are_fluor_cols(ff)]
  if (!(reference_channel %in% fluoCols)) {
    stop("'scatter_ref_marker' must be a fluorochrome channel.", call. = FALSE)
  }

  scatter_channels <- flowCore::colnames(ff)[
    !sw_are_fluor_cols(ff) & sw_are_signal_cols(ff)
  ]
  if (length(scatter_channels) == 0) {
    if (verbose) message("No scatter channels found to scale.")
    return(transList)
  }

  # Get reference quantiles (with optional transform applied)
  if (!is.null(transList) &&
      !is.null(transList@transforms[[reference_channel]])) {
    transfo <- flowCore::transformList(
      from = reference_channel,
      tfun = transList@transforms[[reference_channel]]@f
    )
    ff_t <- flowCore::transform(ff[, reference_channel], transfo)
  } else {
    ff_t <- ff[, reference_channel]
  }

  q5_goal <- stats::quantile(
    flowCore::exprs(ff_t)[, reference_channel], 0.05)
  q95_goal <- stats::quantile(
    flowCore::exprs(ff_t)[, reference_channel], 0.95)

  # Helper to add a linear transform for a scatter channel
  add_linear <- function(ch, a, b) {
    tf <- flowCore::linearTransform(a = a, b = b)
    if (is.null(transList)) {
      transList <<- flowCore::transformList(from = ch, tfun = tf)
    } else {
      transList@transforms[[ch]] <<- NULL
      transList <<- c(transList, flowCore::transformList(ch, tf))
    }
  }

  # Process area scatter channels and propagate to -H and -W
  for (prefix in c("FSC", "SSC")) {
    area_ch <- paste0(prefix, "-A")
    if (area_ch %in% scatter_channels) {
      q5 <- stats::quantile(flowCore::exprs(ff)[, area_ch], 0.05)
      q95 <- stats::quantile(flowCore::exprs(ff)[, area_ch], 0.95)
      a <- (q95_goal - q5_goal) / (q95 - q5)
      b <- q5_goal - q5 * a

      if (verbose) {
        message(sprintf("  %s: a=%.4e, b=%.4e (q5=%.1f->%.1f, q95=%.1f->%.1f)",
                        area_ch, a, b, q5, q5_goal, q95, q95_goal))
      }

      add_linear(area_ch, a, b)

      # Propagate same transform to -H and -W
      for (suffix in c("-H", "-W")) {
        related_ch <- paste0(prefix, suffix)
        if (related_ch %in% scatter_channels) {
          add_linear(related_ch, a, b)
          if (verbose) {
            message(sprintf("  %s: same transform as %s",
                            related_ch, area_ch))
          }
        }
      }
    }
  }

  transList
}

# ---------------------------------------------------------------------------
# applyScaleTransforms
# ---------------------------------------------------------------------------

#' Apply Scale Transformations
#'
#' Applies a pre-computed \code{flowCore::transformList} to a
#' \code{flowFrame} or \code{flowSet}.
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param x A \code{flowCore::flowFrame} or \code{flowCore::flowSet}.
#' @param trans_list A \code{flowCore::transformList} object, typically
#'   obtained from \code{\link{sw_estimate_scale_transforms}}.
#'
#' @return The transformed \code{flowFrame} or \code{flowSet}.
#'
#' @examples
#' \dontrun{
#' trans <- sw_estimate_scale_transforms(ff, fluo_method = "estimateLogicle")
#' ff_transformed <- sw_apply_scale_transforms(ff, trans)
#' }
#'
#' @export
sw_apply_scale_transforms <- function(x, trans_list) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(x, "flowFrame") && !methods::is(x, "flowSet")) {
    stop("'x' must be a flowFrame or flowSet object.", call. = FALSE)
  }

  if (!methods::is(trans_list, "transformList")) {
    stop("'trans_list' must be a flowCore::transformList object.",
         call. = FALSE)
  }

  flowCore::transform(x, trans_list)
}
