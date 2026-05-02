#' @title Gating and Event Filtering Utilities
#'
#' @description
#' Standalone gating functions for doublet removal, debris removal, and margin
#' event removal. Ported from CytoPipeline and CytoPipelineUtils
#' (UCLouvain-CBIO) and adapted for SpectraWeaveR.
#'
#' These functions complement the openCyto-based gating in \code{R/gate.R}
#' by providing lightweight, template-free alternatives for common
#' preprocessing gates.
#'
#' @name gating_utils
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# singletsGate — parallelogram singlet gate
# ---------------------------------------------------------------------------

#' Create a Singlets Gate
#'
#' Creates a parallelogram-shaped polygon gate for doublet removal based on
#' the ratio between two channels (typically FSC-A and FSC-H). Cells are
#' retained if their channel ratio is within \code{nmad} median absolute
#' deviations of the median ratio.
#'
#' Unlike PeacoQC's semi-conic gate, this produces a parallelogram with
#' vertical edges (constant channel-2 range at any channel-1 value), which
#' is more appropriate for spectral flow cytometry data.
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param ff A \code{flowCore::flowFrame}.
#' @param filter_id Character; name for the gate filter.
#'   Default: \code{"Singlets"}.
#' @param channel1 Character; the area channel (x-axis).
#'   Default: \code{"FSC-A"}.
#' @param channel2 Character; the height channel (y-axis).
#'   Default: \code{"FSC-H"}.
#' @param nmad Numeric; number of MADs above the median ratio to allow.
#'   Default: \code{4}.
#' @param verbose Logical; if \code{TRUE}, prints gate statistics.
#'   Default: \code{FALSE}.
#'
#' @return A \code{flowCore::polygonGate} object.
#'
#' @export
sw_gate_singlets <- function(ff,
                             filter_id = "Singlets",
                             channel1 = "FSC-A",
                             channel2 = "FSC-H",
                             nmad = 4,
                             verbose = FALSE) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  if (!is.character(channel1) || length(channel1) != 1) {
    stop("'channel1' must be a single channel name.", call. = FALSE)
  }
  if (!is.character(channel2) || length(channel2) != 1) {
    stop("'channel2' must be a single channel name.", call. = FALSE)
  }
  if (!is.numeric(nmad) || length(nmad) != 1 || nmad <= 0) {
    stop("'nmad' must be a single positive number.", call. = FALSE)
  }

  all_ch <- flowCore::colnames(ff)
  if (!(channel1 %in% all_ch)) {
    stop("Channel '", channel1, "' not found in flowFrame.", call. = FALSE)
  }
  if (!(channel2 %in% all_ch)) {
    stop("Channel '", channel2, "' not found in flowFrame.", call. = FALSE)
  }

  expr_data <- flowCore::exprs(ff)
  ch1_vals <- expr_data[, channel1]
  ch2_vals <- expr_data[, channel2]

  ratio <- ch1_vals / (1 + ch2_vals)
  ratio_median <- stats::median(ratio)
  ratio_mad <- stats::mad(ratio)

  ch1_median <- stats::median(ch1_vals)
  ch1_min <- min(ch1_vals)
  ch1_max <- max(ch1_vals)

  ch2_ref <- ch1_median / ratio_median - 1
  ch2_min <- ch1_median / (ratio_median + nmad * ratio_mad) - 1
  ch2_range <- ch2_ref - ch2_min

  if (verbose) {
    message(sprintf(
      "Singlets gate: median ratio=%.4f, width=%.4f, ch1 median=%.1f, ch2 range=%.1f",
      ratio_median, nmad * ratio_mad, ch1_median, ch2_range
    ))
  }

  ch2_target_min <- ch1_min / ratio_median - 1
  ch2_target_max <- ch1_max / ratio_median - 1

  # Parallelogram gate with vertical edges
  gate_matrix <- matrix(
    data = c(
      ch1_min, ch1_min, ch1_max, ch1_max,
      ch2_target_min - ch2_range,
      ch2_target_min + ch2_range,
      ch2_target_max + ch2_range,
      ch2_target_max - ch2_range
    ),
    ncol = 2,
    dimnames = list(NULL, c(channel1, channel2))
  )

  flowCore::polygonGate(filterId = filter_id, .gate = gate_matrix)
}


# ---------------------------------------------------------------------------
# removeDoublets — CytoPipeline algorithm
# ---------------------------------------------------------------------------

#' Remove Doublets Using CytoPipeline Algorithm
#'
#' Removes doublet events using the parallelogram singlet gate from
#' \code{\link{sw_gate_singlets}}. Can apply gates on multiple channel pairs
#' (e.g., FSC-A/FSC-H and SSC-A/SSC-H) and combine them.
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param ff A \code{flowCore::flowFrame}.
#' @param area_channels Character vector of area-type channel names
#'   (e.g., \code{c("FSC-A", "SSC-A")}). Length 1 or 2.
#' @param height_channels Character vector of height-type channel names
#'   (e.g., \code{c("FSC-H", "SSC-H")}). Must match length of
#'   \code{area_channels}.
#' @param nmads Numeric vector; MAD bandwidth per channel pair.
#'   Default: \code{rep(4, length(area_channels))}.
#' @param verbose Logical; if \code{TRUE}, prints gate statistics.
#'   Default: \code{FALSE}.
#'
#' @return A \code{flowFrame} with doublets removed.
#'
#' @export
sw_filter_doublets <- function(ff,
                               area_channels = c("FSC-A", "SSC-A"),
                               height_channels = c("FSC-H", "SSC-H"),
                               nmads = rep(4, length(area_channels)),
                               verbose = FALSE) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  n_filters <- length(area_channels)
  if (n_filters < 1 || n_filters > 2) {
    stop("'area_channels' must have length 1 or 2.", call. = FALSE)
  }
  if (length(height_channels) != n_filters) {
    stop("'height_channels' must have the same length as 'area_channels'.",
         call. = FALSE)
  }
  if (length(nmads) != n_filters) {
    stop("'nmads' must have the same length as 'area_channels'.",
         call. = FALSE)
  }

  for (i in seq_len(n_filters)) {
    current_gate <- sw_gate_singlets(
      ff,
      filter_id = paste0("Singlets_", area_channels[i]),
      channel1 = area_channels[i],
      channel2 = height_channels[i],
      nmad = nmads[i],
      verbose = verbose
    )

    if (i == 1) {
      combined_gate <- current_gate
    } else {
      combined_gate <- combined_gate & current_gate
    }
  }

  flt <- flowCore::filter(ff, combined_gate)
  flowCore::Subset(ff, flt)
}


# ---------------------------------------------------------------------------
# removeDoubletsPeacoQC — PeacoQC-based doublet removal
# ---------------------------------------------------------------------------

#' Remove Doublets Using PeacoQC
#'
#' Removes doublet events using \code{PeacoQC::RemoveDoublets()}.
#' Can apply on multiple channel pairs sequentially.
#'
#' Ported from CytoPipelineUtils (UCLouvain-CBIO).
#'
#' @param ff A \code{flowCore::flowFrame}.
#' @param area_channels Character vector of area-type channel names.
#' @param height_channels Character vector of height-type channel names.
#' @param nmads Numeric vector; MAD bandwidth per channel pair.
#'   Default: \code{rep(4, length(area_channels))}.
#' @param verbose Logical; if \code{TRUE}, prints removal statistics.
#'   Default: \code{TRUE}.
#' @param ... Additional arguments passed to
#'   \code{PeacoQC::RemoveDoublets()}.
#'
#' @return A \code{flowFrame} with doublets removed.
#'
#' @export
sw_filter_doublets_peacoqc <- function(ff,
                                       area_channels = c("FSC-A", "SSC-A"),
                                       height_channels = c("FSC-H", "SSC-H"),
                                       nmads = rep(4, length(area_channels)),
                                       verbose = TRUE,
                                       ...) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required. Install it from Bioconductor.",
         call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  if (flowCore::nrow(ff) == 0) {
    return(ff)
  }

  n_filters <- length(area_channels)
  if (n_filters < 1 || n_filters > 2) {
    stop("'area_channels' must have length 1 or 2.", call. = FALSE)
  }
  if (length(height_channels) != n_filters) {
    stop("'height_channels' must have the same length as 'area_channels'.",
         call. = FALSE)
  }
  if (length(nmads) != n_filters) {
    stop("'nmads' must have the same length as 'area_channels'.",
         call. = FALSE)
  }

  for (i in seq_len(n_filters)) {
    ff <- PeacoQC::RemoveDoublets(
      ff,
      channel1 = area_channels[i],
      channel2 = height_channels[i],
      nmad = nmads[i],
      verbose = verbose,
      ...
    )
  }

  ff
}


# ---------------------------------------------------------------------------
# removeDebrisManualGate
# ---------------------------------------------------------------------------

#' Remove Debris Using Manual Polygon Gate
#'
#' Removes debris events using a manual polygon gate in the FSC-A vs SSC-A
#' 2D representation. Coordinates are provided as a numeric vector:
#' first all x-coordinates (FSC), then all y-coordinates (SSC).
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param ff A \code{flowCore::flowFrame}.
#' @param fsc_channel Character; name of the forward scatter channel.
#'   Default: \code{"FSC-A"}.
#' @param ssc_channel Character; name of the side scatter channel.
#'   Default: \code{"SSC-A"}.
#' @param gate_data Numeric vector of polygon gate coordinates. The first
#'   half contains x-coordinates (FSC values), the second half contains
#'   y-coordinates (SSC values). Must have even length (at least 6 for a
#'   triangle).
#'
#' @return A \code{flowFrame} with debris events removed.
#'
#' @examples
#' \dontrun{
#' # Define a polygon gate in FSC-A vs SSC-A space
#' gate_coords <- c(
#'   73615, 110174, 213000, 201000, 126000,  # FSC-A x-coords
#'   47679, 260500, 260500, 113000, 35000     # SSC-A y-coords
#' )
#' ff_clean <- sw_filter_debris(ff, gate_data = gate_coords)
#' }
#'
#' @export
sw_filter_debris <- function(ff,
                                  fsc_channel = "FSC-A",
                                  ssc_channel = "SSC-A",
                                  gate_data) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  if (!is.character(fsc_channel) || length(fsc_channel) != 1) {
    stop("'fsc_channel' must be a single channel name.", call. = FALSE)
  }
  if (!is.character(ssc_channel) || length(ssc_channel) != 1) {
    stop("'ssc_channel' must be a single channel name.", call. = FALSE)
  }

  if (!is.numeric(gate_data) || length(gate_data) < 6 ||
      length(gate_data) %% 2 != 0) {
    stop("'gate_data' must be a numeric vector with even length >= 6.",
         call. = FALSE)
  }

  all_ch <- flowCore::colnames(ff)
  if (!(fsc_channel %in% all_ch)) {
    stop("Channel '", fsc_channel, "' not found in flowFrame.", call. = FALSE)
  }
  if (!(ssc_channel %in% all_ch)) {
    stop("Channel '", ssc_channel, "' not found in flowFrame.", call. = FALSE)
  }

  gate_matrix <- matrix(
    data = gate_data, ncol = 2,
    dimnames = list(NULL, c(fsc_channel, ssc_channel))
  )

  cells_gate <- flowCore::polygonGate(
    filterId = "Cells",
    .gate = gate_matrix
  )

  selected <- flowCore::filter(ff, cells_gate)
  flowCore::Subset(ff, selected)
}


# ---------------------------------------------------------------------------
# removeMarginsPeacoQC — enhanced margin removal
# ---------------------------------------------------------------------------

#' Remove Margin Events Using PeacoQC
#'
#' Enhanced wrapper around \code{PeacoQC::RemoveMargins()} that automatically
#' selects signal channels and supports per-channel range specifications.
#' Works on both \code{flowFrame} and \code{flowSet} inputs.
#'
#' Ported from CytoPipeline (UCLouvain-CBIO).
#'
#' @param x A \code{flowCore::flowFrame} or \code{flowCore::flowSet}.
#' @param channel_specifications A named list of lists, where each inner list
#'   contains \code{minRange} and \code{maxRange} values for a specific
#'   channel. Use channel names or marker names as list names.
#'   Special name \code{"AllFluoChannels"} sets defaults for all fluorochrome
#'   channels. Default: \code{NULL} (use PeacoQC internal defaults).
#' @param ... Additional arguments passed to
#'   \code{PeacoQC::RemoveMargins()}.
#'
#' @return A \code{flowFrame} or \code{flowSet} with margin events removed.
#'
#' @export
sw_filter_margins_peacoqc <- function(x, channel_specifications = NULL, ...) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("Package 'PeacoQC' is required. Install it from Bioconductor.",
         call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.null(channel_specifications)) {
    if (!is.list(channel_specifications)) {
      stop("'channel_specifications' must be a list of lists.", call. = FALSE)
    }
    if (!all(vapply(channel_specifications, length, integer(1)) == 2)) {
      stop("Each element of 'channel_specifications' must contain ",
           "exactly 2 values (minRange, maxRange).", call. = FALSE)
    }
  }

  process_one <- function(ff, channel_specifications) {
    channels_for_margins <- flowCore::colnames(ff)[sw_channel_is_signal(ff)]

    markers_for_margins <- flowCore::pData(
      flowCore::parameters(ff)
    )$desc[sw_channel_is_signal(ff)]

    pqc_specs <- channel_specifications

    if (!is.null(pqc_specs)) {
      # Handle AllFluoChannels default
      default_fluo <- pqc_specs[["AllFluoChannels"]]
      if (!is.null(default_fluo)) {
        pqc_specs[["AllFluoChannels"]] <- NULL
      }

      # Resolve marker names to channel names
      spec_names <- names(pqc_specs)
      for (i in seq_along(spec_names)) {
        ch_name <- spec_names[i]
        if (!(ch_name %in% flowCore::colnames(ff))) {
          marker_idx <- which(markers_for_margins == ch_name)
          if (length(marker_idx) == 0) {
            stop("channel_specifications: could not find '", ch_name,
                 "' as channel or marker.", call. = FALSE)
          }
          spec_names[i] <- channels_for_margins[marker_idx[1]]
        }
      }
      names(pqc_specs) <- spec_names

      # Apply default fluo specs to unlisted fluor channels
      if (!is.null(default_fluo)) {
        fluor_chs <- flowCore::colnames(ff)[sw_channel_is_fluor(ff)]
        for (ch in fluor_chs) {
          if (!(ch %in% spec_names)) {
            pqc_specs[[ch]] <- default_fluo
          }
        }
      }
    }

    PeacoQC::RemoveMargins(
      ff,
      channels = channels_for_margins,
      channel_specifications = pqc_specs,
      ...
    )
  }

  if (methods::is(x, "flowFrame")) {
    process_one(x, channel_specifications)
  } else if (methods::is(x, "flowSet")) {
    flowCore::fsApply(x, FUN = process_one, simplify = TRUE,
                      channel_specifications = channel_specifications)
  } else {
    stop("'x' must be a flowFrame or flowSet object.", call. = FALSE)
  }
}
