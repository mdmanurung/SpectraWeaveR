#' @title Cleaning Stage: Margins, Transform, Signal QC, Gates
#'
#' @description
#' One call for the cleaning stage of the validated pipeline order: margin
#' events, arcsinh transformation, PeacoQC signal QC, then non-debris,
#' singlet and live gates. Nothing after the live gate.
#'
#' @name clean_gate
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Required arguments: one table drives the missing-argument message, and the
# roxygen @param text and the ROADMAP register rows repeat its wording.
# ---------------------------------------------------------------------------

.sw_clean_args <- list(
  margin_range = list(
    controls = paste0(
      "Lower and upper limit (raw scale) of the margin step, applied to every ",
      "fluorescence channel. PeacoQC removes an event at or beyond a limit in ",
      "any fluorescence channel, and the event at the observed minimum or ",
      "maximum of a channel. Scatter channels keep the range stored in the ",
      "FCS file. A lower limit above 0 is read as 0."
    ),
    how = paste0(
      "Look at the raw fluorescence densities and the share of events flagged ",
      "per channel. ",
      "In a failed smoke run on 5 BMV donors with PeacoQC's default ranges, ",
      "PeacoQC warned that more than 99.34 % (A3), 99.81 % (A6), 99.48 % ",
      "(B11), 99.16 % (B2) and 99.86 % (D1) of the events were margin events ",
      "(experiments/runs/2026-11-11_clean-run-smoke/warnings.txt, a run ",
      "folder that git does not track). The comment on MAX_MARGIN_FRAC in ",
      "experiments/benchmarks/05_bmv-unmix-compare/12_cleanup.R notes that ",
      "the defaults removed 99.9%."
    ),
    reference = paste0(
      "c(-21761, 1e7) on fluorescence channels (arcsinh -2 to the ",
      "instrument top) as a test value"
    ),
    example = "margin_range = c(<min>, <max>)"
  ),
  debris_window = list(
    controls = paste0(
      "Lower and upper FSC-A limit (raw scale) of the events the nonDebris ",
      "density is built from. openCyto drops events outside [min, max] ",
      "before it looks for peaks and valleys, so the cut is a valley of the ",
      "density of the events inside the window and lies inside it. Events ",
      "below the cut are debris. Events above max are kept. If the window ",
      "holds one peak or no valley, the cut falls on the lowest event in the ",
      "window."
    ),
    how = paste0(
      "Plot raw FSC-A as a density and FSC-A against SSC-A (ggplot2::geom_hex() ",
      "on sw_io_ff_to_tibble(ff)). Debris is the low-FSC cloud beside the main ",
      "cells. Look at one sample per batch, because the valley can move. The ",
      "window must hold the valley you want. Use Inf as the upper bound to ",
      "keep every event above min in the density."
    ),
    reference = "c(2e5, 1.5e6) with mindensity",
    example = "debris_window = c(<min>, <max>)"
  ),
  singlet_wider_gate = list(
    controls = paste0(
      "Width of the openCyto singletGate on FSC-A against FSC-H. TRUE ",
      "widens it and keeps more events near the edge."
    ),
    how = paste0(
      "Plot FSC-A against FSC-H of non-debris events (raw scale) and look at ",
      "the width of the diagonal. Run both settings on a trial sample and ",
      "compare the singlets rows of the audit table."
    ),
    reference = "TRUE",
    example = "singlet_wider_gate = <TRUE or FALSE>"
  ),
  pool_samples = list(
    controls = paste0(
      "TRUE fits each gate once on the pooled events of all samples and ",
      "applies it to every sample. FALSE fits one gate per sample."
    ),
    how = paste0(
      "Plot per-sample densities of FSC-A and of the viability channel on ",
      "one panel. If the valleys line up across donors and batches, pooling ",
      "is defensible. If they shift by batch, do not pool. Compare the ",
      "$gates table of a trial run with each setting."
    ),
    reference = "TRUE (gates pooled over donors)",
    example = "pool_samples = <TRUE or FALSE>"
  ),
  viability_channel = list(
    controls = "Channel of the viability dye. Dead cells are the high tail.",
    how = paste0(
      "List channels with sw_channel_get_fluor(ff) and check the panel sheet ",
      "for the amine-reactive dye. Plot its density after arcsinh and look ",
      "for a dead peak or tail."
    ),
    reference = "\"LIVE DEAD NIR-A\"",
    example = "viability_channel = \"<channel name>\""
  ),
  viability_method = list(
    controls = paste0(
      "How the live cut is found: \"mindensity\" takes the valley of the ",
      "density, \"boundary\" takes the fixed cut in viability_cut."
    ),
    how = paste0(
      "Plot the density of the arcsinh-transformed viability channel per ",
      "sample (sw_gate_density_cut() shows one valley cut per vector). Use ",
      "\"mindensity\" when every sample shows two peaks. Use \"boundary\" ",
      "when a sample has no valley."
    ),
    reference = "\"mindensity\"",
    example = "viability_method = \"<mindensity or boundary>\""
  ),
  viability_cut = list(
    controls = paste0(
      "Live cut on the arcsinh scale for viability_method = \"boundary\" ",
      "(events below it are live). Must be NULL for \"mindensity\"."
    ),
    how = paste0(
      "For \"boundary\", read the cut from the same density plot as ",
      "viability_method, or from an unstained or fully live control. For ",
      "\"mindensity\", pass NULL on purpose."
    ),
    reference = "NULL (mindensity finds the valley)",
    example = "viability_cut = <number or NULL>"
  ),
  IT_limit = list(
    controls = paste0(
      "Isolation Tree limit of PeacoQC. It trades stringency against the ",
      "number of events kept."
    ),
    how = paste0(
      "Run sw_qc_run() on one trial sample per batch with a few values, read ",
      "pct_removed, and inspect the QC plot with sw_qc_plot(). Summarise ",
      "several samples with sw_qc_summary() on sw_qc_batch() output."
    ),
    reference = "0.55",
    example = "IT_limit = <number between 0 and 1>"
  ),
  MAD = list(
    controls = paste0(
      "Median Absolute Deviation limit of PeacoQC. It trades stringency ",
      "against the number of events kept."
    ),
    how = paste0(
      "Same trial runs as IT_limit: vary it with sw_qc_run(), read ",
      "pct_removed, inspect with sw_qc_plot(), compare samples with ",
      "sw_qc_summary()."
    ),
    reference = "6",
    example = "MAD = <positive number>"
  ),
  cofactor = list(
    controls = paste0(
      "Divisor of the arcsinh transformation of all fluorescence channels, ",
      "asinh(x / cofactor). Scatter stays raw. It shifts the QC and viability ",
      "cuts."
    ),
    how = paste0(
      "Transform a few channels, including the viability channel, with ",
      "several cofactors (flowCore::arcsinhTransform(a = 0, b = 1 / cofactor)) ",
      "and plot their densities. Pick the value where the negative peak is ",
      "compact and the positive peak is not squeezed."
    ),
    reference = "6000",
    example = "cofactor = <positive number>"
  )
)

# Message for the arguments the caller left out; one stop, never a warning.
.sw_clean_missing_message <- function(absent) {
  block <- vapply(
    absent,
    function(nm) {
      a <- .sw_clean_args[[nm]]
      paste0(
        nm,
        "\n",
        "  Controls: ",
        a$controls,
        "\n",
        "  How to choose: ",
        a$how,
        "\n",
        "  Reference: benchmark 05 used ",
        a$reference,
        " on 40-colour Aurora PBMC; not a recommendation.\n",
        "  Example: ",
        a$example
      )
    },
    character(1)
  )
  paste0(
    "sw_clean_run() needs ",
    length(absent),
    " required argument(s) with no ",
    "default: ",
    paste(absent, collapse = ", "),
    ".\n",
    "Choose each value yourself. The function fills in nothing.\n\n",
    paste(block, collapse = "\n\n")
  )
}

#' Clean Spectral Flow Data: Margins, Transform, Signal QC, Gates
#'
#' Runs the cleaning stage of the validated pipeline order on FCS files,
#' flowFrames or a flowSet, in this order:
#' \enumerate{
#'   \item margin events removed on the raw scale, with \code{margin_range} on
#'     the fluorescence channels (\code{\link{sw_filter_margins_peacoqc}});
#'   \item arcsinh transformation of all fluorescence channels
#'     (\code{\link{sw_transform_apply}}); scatter channels stay raw;
#'   \item PeacoQC signal QC (\code{\link{sw_qc_batch}} and
#'     \code{\link{sw_qc_run}});
#'   \item openCyto gates on the QC-cleaned events: \code{nonDebris},
#'     \code{singlets}, \code{live} (\code{\link{sw_gate_template}} and
#'     \code{\link{sw_gate_run}}).
#' }
#' Nothing after the live gate: no CD3 or lymphocyte gate.
#'
#' Every impactful value is a required argument. If any are missing, the
#' function stops once and its message names all of them, with what each
#' controls, how to choose it, and the benchmark 05 reference value (not a
#' recommendation). No step runs before that check.
#'
#' @param x Character vector of FCS file paths (read with
#'   \code{truncate_max_range = FALSE}), a \code{flowCore::flowFrame}, a
#'   \code{flowSet}, or a list of \code{flowFrame}s. A list or flowSet keeps its
#'   names; otherwise samples are named after the files or \code{sample_<i>}.
#' @param margin_range Numeric length 2, \code{c(min, max)}; required, no
#'   default. Lower and upper limit (raw scale) of the margin step, applied to
#'   every fluorescence channel; two finite numbers with \code{min} below
#'   \code{max}. Scatter channels keep the range stored in the FCS file. It
#'   is passed to \code{PeacoQC::RemoveMargins()} as the \code{minRange} and
#'   \code{maxRange} of every fluorescence channel (see Details). The value
#'   sets how many events every later step sees. The user approved it as a
#'   required argument in the 2026-11-11 session, because PeacoQC's own default
#'   ranges (the FCS parameter ranges) flagged almost all events in a smoke
#'   run (sources in Details).
#' @param debris_window Numeric length 2; required, no default. Lower and
#'   upper FSC-A limit (raw scale) of the events the \code{nonDebris} density
#'   is built from. It is written as \code{min} and \code{max} of the openCyto
#'   \code{mindensity} method, which drops events outside the window before
#'   it looks for peaks and valleys (see Details). The cut is therefore a
#'   valley of the density of the events inside the window and lies inside
#'   it. Events below the cut are debris; events above \code{max} are kept.
#'   Use \code{Inf} for no upper limit. It sets the debris boundary, so it
#'   changes every later count.
#' @param singlet_wider_gate Logical; required, no default. \code{wider_gate}
#'   of the openCyto \code{singletGate}. \code{TRUE} keeps more events near the
#'   edge of the singlet diagonal.
#' @param pool_samples Logical; required, no default. \code{TRUE}: one gate per
#'   node is fitted on the pooled events of all samples and applied to each,
#'   so every sample has the same cut. \code{FALSE}: one gate per sample. The
#'   template records this as \code{collapseDataForGating} and, when pooling,
#'   \code{groupBy} set to the number of samples (openCyto pools only when
#'   \code{groupBy} is not empty). It changes the cut of every sample.
#' @param viability_channel Character; required, no default. Channel of the
#'   viability dye, used by the \code{live} gate. A wrong channel gates on the
#'   wrong cells.
#' @param viability_method Character; required, no default. \code{"mindensity"}
#'   (valley of the density) or \code{"boundary"} (fixed cut in
#'   \code{viability_cut}).
#' @param viability_cut Required, no default. A single number on the arcsinh
#'   scale for \code{viability_method = "boundary"} (events below it are
#'   live), or \code{NULL} for \code{"mindensity"}. Pass \code{NULL}
#'   explicitly; leaving it out is an error.
#' @param IT_limit Numeric in (0, 1); required, no default. Isolation Tree
#'   limit of PeacoQC. It trades stringency against events kept.
#' @param MAD Numeric > 0; required, no default. MAD limit of PeacoQC. It
#'   trades stringency against events kept.
#' @param cofactor Numeric > 0; required, no default. Arcsinh cofactor:
#'   \code{asinh(x / cofactor)} on every fluorescence channel. It shifts the
#'   QC and viability cuts.
#'
#' @details
#' The margin step calls \code{\link{sw_filter_margins_peacoqc}} with
#' \code{channel_specifications = list(AllFluoChannels = margin_range)}. The
#' wrapper turns \code{AllFluoChannels} into one \code{c(minRange, maxRange)}
#' entry for each fluorescence channel (\code{sw_channel_is_fluor()}) and
#' hands them to \code{PeacoQC::RemoveMargins()}, which writes them over the
#' \code{minRange} and \code{maxRange} of those channels. Scatter channels
#' keep the \code{minRange} and \code{maxRange} of the FCS parameters. The
#' \code{debris_window} and the singlet gate act on raw scatter; the
#' viability cut acts on the arcsinh scale.
#'
#' \strong{Why \code{margin_range} has no default.} In a failed smoke run on 5
#' BMV donors with PeacoQC's default ranges, PeacoQC warned that more than
#' 99.34 % (A3), 99.81 % (A6), 99.48 % (B11), 99.16 % (B2) and 99.86 % (D1)
#' of the events were margin events
#' (\code{experiments/runs/2026-11-11_clean-run-smoke/warnings.txt}, a run
#' folder that git does not track). The comment on \code{MAX_MARGIN_FRAC} in
#' \code{experiments/benchmarks/05_bmv-unmix-compare/12_cleanup.R} notes that
#' the defaults removed 99.9%.
#'
#' \strong{The margin step removes more than events on the range limit.}
#' PeacoQC's \code{RemoveMargins()} keeps an event only if, in every scatter
#' and fluorescence channel, \code{e > max(min(minRange, 0), min(e))} and
#' \code{e < min(maxRange, max(e))}. A lower limit above 0 therefore acts as 0.
#' The package wrapper \code{sw_filter_margins_peacoqc()} picks these channels
#' with \code{sw_channel_is_signal()}, which leaves out \code{Time},
#' \code{Original_ID}, \code{File}, \code{SampleID} and \code{AF Index} by
#' default, so \code{Time} is not checked. Both comparisons are strict, so in
#' each channel of each sample the events at the observed minimum and the
#' observed maximum are removed, even when no event reaches \code{margin_range}
#' or the instrument range. The \code{margins} row of the \code{audit} table
#' (\code{n_removed}) shows how many events this step removed per sample.
#'
#' \strong{What \code{debris_window} does.} The \code{nonDebris} gate is the
#' openCyto \code{mindensity} method on FSC-A, called with \code{min} and
#' \code{max} (only \code{min} when \code{max} is \code{Inf}). In the
#' installed openCyto (2.18.0), \code{gate_mindensity()} starts with
#' \code{fr <- .truncate_flowframe(fr, channels = channel, min = min, max = max)},
#' and \code{.truncate_flowframe()} ends with
#' \code{Subset(flow_frame, truncate_filter)}, a rectangle gate on
#' \code{c(min, max)}. So the window sets the data range the density is built
#' from; it is not a separate bound applied to the cut. The peaks and valleys
#' come from that truncated vector \code{x}, with
#' \code{gate_range <- c(min(x), max(x))}. With two or more peaks and
#' valleys the cut is the deepest valley between the two highest peaks, or the
#' deepest valley overall if none lies between them. The gate
#' is \code{c(cutpoint, Inf)} on FSC-A, which keeps events at or above the
#' cut, so events below \code{min} are always debris and events above
#' \code{max} are never removed by this gate. With one peak, or no valley
#' inside the window, openCyto sets
#' \code{cutpoint <- ifelse(positive, gate_range[1], gate_range[2])}: the cut
#' falls on the lowest event inside the window and removes nothing at or
#' above \code{min}. \code{sw_clean_run()} reports this case: it warns and sets
#' \code{cut_at_window_floor} to \code{TRUE} in the \code{gates} table. With
#' \code{pool_samples = TRUE} the one pooled cut is compared with the lowest
#' in-window event of all samples together; with \code{FALSE} each sample's cut
#' is compared with the lowest in-window event of that sample. A
#' window that leaves out the real valley can also put the cut at a different
#' valley inside it, so check the \code{nonDebris} \code{definition} in the
#' \code{gates} table against the FSC-A density.
#'
#' A sample with no events left after any step stops the run, with the
#' sample and the step named. After the margin step the message also gives the
#' events before and after, the percentage removed and \code{margin_range}. If
#' PeacoQC raises an error on a sample (for example when too few events are
#' left for its analysis), the run stops with the sample name, the number of
#' events handed to PeacoQC and the original PeacoQC message. The run never
#' skips a sample and sets no minimum event count.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{events}}{Named list of \code{flowFrame}s: the live, singlet,
#'       non-debris, QC-cleaned events, fluorescence on the arcsinh scale.}
#'     \item{\code{audit}}{Tibble, one row per sample and step, with columns
#'       \code{sample}, \code{step} (\code{raw}, \code{margins}, \code{qc},
#'       \code{non_debris}, \code{singlets}, \code{live}), \code{n_events},
#'       \code{n_removed} (events lost at this step; \code{NA} for
#'       \code{raw}), \code{pct_of_initial}, \code{pct_of_previous}.}
#'     \item{\code{audit_total}}{The same, over all samples
#'       (\code{\link{sw_io_event_audit}} columns plus \code{n_removed}).}
#'     \item{\code{gates}}{Tibble, one row per sample and gate: \code{sample},
#'       \code{gate}, \code{class}, \code{channels}, \code{definition}
#'       (rectangle limits or polygon vertices) and \code{keeps} (\code{"inside"}
#'       or \code{"outside"} the shape; the live gate keeps events outside) and
#'       \code{cut_at_window_floor} (\code{nonDebris} rows only, \code{NA}
#'       otherwise: \code{TRUE} when the cut sits on the lowest event inside
#'       \code{debris_window} because it held one peak or no valley).}
#'     \item{\code{qc}}{The per-sample PeacoQC summary tibble of
#'       \code{\link{sw_qc_batch}} (\code{sample}, \code{n_before},
#'       \code{n_after}, \code{n_removed}, \code{pct_removed}), with the
#'       \code{IT_limit} and \code{MAD} used.}
#'     \item{\code{template}}{The gating template as a data frame.}
#'     \item{\code{parameters}}{The values of all required arguments.}
#'   }
#'
#' @examples
#' if (all(vapply(c("flowCore", "flowWorkspace", "openCyto", "PeacoQC"),
#'                requireNamespace, logical(1), quietly = TRUE))) {
#'   # A small synthetic frame: cells, debris (low FSC-A), doublets (low
#'   # FSC-H for their FSC-A) and dead cells (high viability dye).
#'   set.seed(1)
#'   n_cells <- 8000
#'   n_debris <- 600
#'   n_doublet <- 600
#'   n_dead <- 800
#'   fsc_a <- c(
#'     rnorm(n_cells, 6e5, 3e4), rnorm(n_debris, 5e4, 1.5e4),
#'     rnorm(n_doublet, 6e5, 3e4), rnorm(n_dead, 6e5, 3e4)
#'   )
#'   n <- length(fsc_a)
#'   fsc_h <- fsc_a * 0.8 + rnorm(n, 0, 8000)
#'   is_doublet <- seq_len(n) %in% (n_cells + n_debris + seq_len(n_doublet))
#'   fsc_h[is_doublet] <- fsc_a[is_doublet] * 0.45
#'   viability <- c(rnorm(n - n_dead, 1000, 600), rnorm(n_dead, 3e5, 5e4))
#'   mat <- cbind(
#'     "FSC-A" = fsc_a,
#'     "FSC-H" = fsc_h,
#'     "SSC-A" = abs(rnorm(n, 1e5, 2e4)),
#'     "Time" = seq_len(n) / 100,
#'     "BV421-A" = rnorm(n, 2000, 800),
#'     "PE-A" = rnorm(n, 3000, 900),
#'     "LIVE DEAD NIR-A" = viability
#'   )
#'   mat <- mat[sample(n), ] # mix the classes along the run
#'   mat[, "Time"] <- seq_len(n) / 100
#'   ff <- flowCore::flowFrame(mat)
#'   pd <- flowCore::pData(flowCore::parameters(ff))
#'   pd$minRange <- -1e6
#'   pd$maxRange <- 4e6
#'   flowCore::pData(flowCore::parameters(ff)) <- pd
#'
#'   # Every value below is a PLACEHOLDER for this synthetic frame. Choose
#'   # each one for your own data (see the message of a call without them).
#'   # PeacoQC warns that this frame has too few bins for its isolation
#'   # tree analysis; real files have more events.
#'   res <- sw_clean_run(
#'     ff,
#'     margin_range = c(-21761, 1e7), # PLACEHOLDER for this synthetic frame
#'     debris_window = c(1e4, 1.5e6), # PLACEHOLDER: must hold the debris valley
#'     singlet_wider_gate = TRUE, # PLACEHOLDER: TRUE or FALSE
#'     pool_samples = TRUE, # PLACEHOLDER: TRUE or FALSE
#'     viability_channel = "LIVE DEAD NIR-A", # PLACEHOLDER: your dye channel
#'     viability_method = "mindensity", # PLACEHOLDER: "mindensity" or "boundary"
#'     viability_cut = NULL, # PLACEHOLDER: NULL for "mindensity"
#'     IT_limit = 0.55, # PLACEHOLDER: between 0 and 1
#'     MAD = 6, # PLACEHOLDER: positive number
#'     cofactor = 6000 # PLACEHOLDER: positive number
#'   )
#'   res$audit
#'   res$gates
#' }
#'
#' @seealso \code{\link{sw_gate_template}}, \code{\link{sw_gate_run}},
#'   \code{\link{sw_qc_run}}, \code{\link{sw_io_event_audit}}
#' @export
sw_clean_run <- function(
  x,
  margin_range,
  debris_window,
  singlet_wider_gate,
  pool_samples,
  viability_channel,
  viability_method,
  viability_cut,
  IT_limit,
  MAD,
  cofactor
) {
  # --- Required arguments: collect all, stop once -------------------------
  is_missing <- c(
    margin_range = missing(margin_range),
    debris_window = missing(debris_window),
    singlet_wider_gate = missing(singlet_wider_gate),
    pool_samples = missing(pool_samples),
    viability_channel = missing(viability_channel),
    viability_method = missing(viability_method),
    viability_cut = missing(viability_cut),
    IT_limit = missing(IT_limit),
    MAD = missing(MAD),
    cofactor = missing(cofactor)
  )
  if (any(is_missing)) {
    stop(
      .sw_clean_missing_message(names(is_missing)[is_missing]),
      call. = FALSE
    )
  }

  # --- Input validation, before any dependency check ----------------------
  is_frame <- function(o) methods::is(o, "flowFrame")
  x_ok <- is.character(x) &&
    length(x) > 0 ||
    is_frame(x) ||
    methods::is(x, "flowSet") ||
    is.list(x) && length(x) > 0 && all(vapply(x, is_frame, logical(1)))
  if (!x_ok) {
    stop(
      "'x' must be FCS file paths, a flowFrame, a flowSet or a list of ",
      "flowFrames.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(margin_range) ||
      length(margin_range) != 2 ||
      !all(is.finite(margin_range)) ||
      margin_range[1] >= margin_range[2]
  ) {
    stop(
      "'margin_range' must be c(min, max): two finite numbers with min ",
      "below max.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(debris_window) ||
      length(debris_window) != 2 ||
      anyNA(debris_window) ||
      !is.finite(debris_window[1]) ||
      debris_window[1] >= debris_window[2]
  ) {
    stop(
      "'debris_window' must be c(min, max) with a finite min below max ",
      "(max may be Inf).",
      call. = FALSE
    )
  }
  for (nm in c("singlet_wider_gate", "pool_samples")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1 || is.na(v)) {
      stop("'", nm, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (
    !is.character(viability_channel) ||
      length(viability_channel) != 1 ||
      is.na(viability_channel) ||
      !nzchar(viability_channel)
  ) {
    stop("'viability_channel' must be a single channel name.", call. = FALSE)
  }
  if (
    !is.character(viability_method) ||
      length(viability_method) != 1 ||
      !viability_method %in% c("mindensity", "boundary")
  ) {
    stop(
      "'viability_method' must be \"mindensity\" or \"boundary\".",
      call. = FALSE
    )
  }
  if (viability_method == "boundary") {
    if (
      !is.numeric(viability_cut) ||
        length(viability_cut) != 1 ||
        !is.finite(viability_cut)
    ) {
      stop(
        "'viability_cut' must be a single finite number when ",
        "viability_method = \"boundary\".",
        call. = FALSE
      )
    }
  } else if (!is.null(viability_cut)) {
    stop(
      "'viability_cut' must be NULL when viability_method = ",
      "\"mindensity\".",
      call. = FALSE
    )
  }
  if (
    !is.numeric(IT_limit) ||
      length(IT_limit) != 1 ||
      is.na(IT_limit) ||
      IT_limit <= 0 ||
      IT_limit >= 1
  ) {
    stop(
      "'IT_limit' must be a number between 0 and 1 (exclusive).",
      call. = FALSE
    )
  }
  if (!is.numeric(MAD) || length(MAD) != 1 || is.na(MAD) || MAD <= 0) {
    stop("'MAD' must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(cofactor) ||
      length(cofactor) != 1 ||
      !is.finite(cofactor) ||
      cofactor <= 0
  ) {
    stop("'cofactor' must be a single positive number.", call. = FALSE)
  }

  if (is.character(x) && !all(file.exists(x))) {
    stop(
      "FCS file(s) not found: ",
      paste(x[!file.exists(x)], collapse = ", "),
      call. = FALSE
    )
  }

  # --- Dependency checks --------------------------------------------------
  for (pkg in c("flowCore", "flowWorkspace", "openCyto", "PeacoQC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }

  # --- Samples ------------------------------------------------------------
  if (is.character(x)) {
    ffs <- lapply(x, function(f) sw_io_read_fcs(f))
    names(ffs) <- tools::file_path_sans_ext(basename(x))
  } else if (is_frame(x)) {
    ffs <- list(x)
    id <- flowCore::identifier(x)
    names(ffs) <- if (nzchar(id) && id != "anonymous") id else "sample_1"
  } else if (methods::is(x, "flowSet")) {
    ffs <- stats::setNames(
      lapply(flowCore::sampleNames(x), function(sn) x[[sn]]),
      flowCore::sampleNames(x)
    )
  } else {
    ffs <- x
    if (is.null(names(ffs))) names(ffs) <- paste0("sample_", seq_along(ffs))
  }
  if (anyDuplicated(names(ffs)) || any(!nzchar(names(ffs)))) {
    stop("Sample names must be unique and non-empty.", call. = FALSE)
  }
  sample_names <- names(ffs)

  needed <- c("FSC-A", "FSC-H", viability_channel)
  for (sn in sample_names) {
    absent <- setdiff(needed, flowCore::colnames(ffs[[sn]]))
    if (length(absent) > 0) {
      stop(
        "Sample '",
        sn,
        "' lacks channel(s): ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
  }
  fluor <- flowCore::colnames(ffs[[1]])[sw_channel_is_fluor(ffs[[1]])]
  if (!viability_channel %in% fluor) {
    stop(
      "'viability_channel' (",
      viability_channel,
      ") is not a ",
      "fluorescence channel.",
      call. = FALSE
    )
  }
  for (sn in sample_names) {
    if (
      !setequal(
        flowCore::colnames(ffs[[sn]])[sw_channel_is_fluor(ffs[[sn]])],
        fluor
      )
    ) {
      stop(
        "Sample '",
        sn,
        "' has different fluorescence channels from '",
        sample_names[1],
        "'.",
        call. = FALSE
      )
    }
  }

  # Counts per sample and step; a step that empties a sample stops the run.
  counts <- list(raw = vapply(ffs, nrow, integer(1)))
  record <- function(step, frames) {
    n <- vapply(frames, nrow, integer(1))
    empty <- names(n)[n == 0]
    if (length(empty) > 0) {
      stop(
        "No events left after step '",
        step,
        "' in sample(s): ",
        paste(empty, collapse = ", "),
        call. = FALSE
      )
    }
    counts[[step]] <<- n
  }
  record("raw", ffs)

  # --- 1. Margins (raw scale) --------------------------------------------
  # margin_range reaches PeacoQC::RemoveMargins() as channel_specifications:
  # "AllFluoChannels" gives every fluorescence channel c(minRange, maxRange);
  # scatter keeps the range of the FCS parameters.
  n_before_margins <- counts[["raw"]]
  ffs <- lapply(
    ffs,
    sw_filter_margins_peacoqc,
    channel_specifications = list(AllFluoChannels = margin_range)
  )
  n_after_margins <- vapply(ffs, nrow, integer(1))
  if (any(n_after_margins == 0)) {
    stop(
      .sw_clean_margins_empty_message(
        n_before_margins,
        n_after_margins,
        margin_range
      ),
      call. = FALSE
    )
  }
  record("margins", ffs)

  # --- 2. Arcsinh on fluorescence channels --------------------------------
  trans_list <- flowCore::transformList(
    fluor,
    lapply(fluor, function(ch) {
      flowCore::arcsinhTransform(
        transformationId = paste0("asinh_", ch),
        a = 0,
        b = 1 / cofactor,
        c = 0
      )
    })
  )
  ffs <- lapply(ffs, sw_transform_apply, trans_list = trans_list)

  # --- 3. PeacoQC ---------------------------------------------------------
  # QC runs on the transformed channels only: the Original_ID column that
  # the margin step appends is a counter, not a signal.
  # One sample per call, so a PeacoQC error can name its sample.
  qc_each <- lapply(sample_names, function(sn) {
    n_to_qc <- nrow(ffs[[sn]])
    tryCatch(
      sw_qc_batch(
        ffs[sn],
        channels = fluor,
        IT_limit = IT_limit,
        MAD = MAD
      ),
      error = function(e) {
        stop(
          "PeacoQC failed on sample '",
          sn,
          "' with ",
          n_to_qc,
          " events handed to it (after the margin step). PeacoQC said: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  })
  ffs <- stats::setNames(
    lapply(qc_each, function(q) q$cleaned[[1]]),
    sample_names
  )
  record("qc", ffs)
  qc_summary <- dplyr::bind_rows(lapply(qc_each, function(q) q$summary))
  qc_summary$IT_limit <- IT_limit
  qc_summary$MAD <- MAD

  # --- 4. Gates: nonDebris, singlets, live --------------------------------
  template <- .sw_clean_template(
    debris_window,
    singlet_wider_gate,
    pool_samples,
    viability_channel,
    viability_method,
    viability_cut,
    n_samples = length(sample_names)
  )
  tmpl_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmpl_file), add = TRUE)
  utils::write.csv(template, tmpl_file, row.names = FALSE)

  fs <- flowCore::flowSet(ffs)
  flowCore::sampleNames(fs) <- sample_names
  gs <- tryCatch(
    sw_gate_run(fs, gating_template = tmpl_file, transform_channels = NULL),
    error = function(e) {
      stop(
        "Gating failed (check debris_window and the viability settings): ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  nodes <- c(
    non_debris = "/nonDebris",
    singlets = "/nonDebris/singlets",
    live = "/nonDebris/singlets/live"
  )
  for (step in names(nodes)) {
    counts[[step]] <- vapply(
      sample_names,
      function(sn) {
        as.integer(flowWorkspace::gh_pop_get_count(gs[[sn]], nodes[[step]]))
      },
      integer(1)
    )
    if (any(counts[[step]] == 0)) {
      stop(
        "No events left after step '",
        step,
        "' in sample(s): ",
        paste(sample_names[counts[[step]] == 0], collapse = ", "),
        call. = FALSE
      )
    }
  }

  events <- sw_gate_extract(gs, nodes[["live"]])
  events <- lapply(events, function(f) {
    if (methods::is(f, "cytoframe")) {
      flowWorkspace::cytoframe_to_flowFrame(f)
    } else {
      f
    }
  })

  # --- Audit and gate tables ----------------------------------------------
  audit_one <- function(n) {
    a <- sw_io_event_audit(as.list(n))
    a$n_events <- unname(a$n_events)
    a$n_removed <- c(NA_integer_, -diff(a$n_events))
    a
  }
  audit <- dplyr::bind_rows(lapply(sample_names, function(sn) {
    a <- audit_one(vapply(counts, function(v) v[[sn]], integer(1)))
    a$sample <- sn
    a[, c("sample", setdiff(names(a), "sample"))]
  }))
  audit_total <- audit_one(vapply(counts, sum, integer(1)))

  gates <- dplyr::bind_rows(lapply(sample_names, function(sn) {
    dplyr::bind_rows(lapply(
      c("nonDebris", "singlets", "live"),
      function(g) .sw_clean_gate_row(sn, g, gs[[sn]], nodes)
    ))
  }))
  # "+" keeps events inside the gate shape, "-" keeps events outside it.
  gates$keeps <- ifelse(
    template$pop[match(gates$gate, template$alias)] == "+",
    "inside",
    "outside"
  )

  # mindensity puts the cut on the lowest event inside debris_window when the
  # window holds one peak or no valley. Flag that case per sample.
  window_floor <- .sw_clean_window_floor(ffs, debris_window, pool_samples)
  debris_cut <- vapply(
    sample_names,
    function(sn) {
      g <- flowWorkspace::gh_pop_get_gate(gs[[sn]], nodes[["non_debris"]])
      g@min[["FSC-A"]]
    },
    numeric(1)
  )
  at_floor <- unname(debris_cut <= window_floor)
  gates$cut_at_window_floor <- NA
  gates$cut_at_window_floor[gates$gate == "nonDebris"] <-
    at_floor[match(gates$sample[gates$gate == "nonDebris"], sample_names)]
  if (any(at_floor)) {
    warning(
      "The nonDebris cut sits on the lowest event inside debris_window ",
      "(FSC-A = ",
      paste(signif(debris_cut[at_floor], 6), collapse = ", "),
      ") for sample(s): ",
      paste(sample_names[at_floor], collapse = ", "),
      ". openCyto found one density peak or no density valley inside the ",
      "window, so it removes nothing at or above debris_window[1]. Check the ",
      "window against the FSC-A density.",
      call. = FALSE
    )
  }

  list(
    events = events,
    audit = audit,
    audit_total = audit_total,
    gates = gates,
    qc = qc_summary,
    template = template,
    parameters = list(
      margin_range = margin_range,
      debris_window = debris_window,
      singlet_wider_gate = singlet_wider_gate,
      pool_samples = pool_samples,
      viability_channel = viability_channel,
      viability_method = viability_method,
      viability_cut = viability_cut,
      IT_limit = IT_limit,
      MAD = MAD,
      cofactor = cofactor
    )
  )
}


# Gating template: the nonDebris and singlets rows come from
# sw_gate_template(); the live row is appended here.
.sw_clean_template <- function(
  debris_window,
  singlet_wider_gate,
  pool_samples,
  viability_channel,
  viability_method,
  viability_cut,
  n_samples
) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  # "myeloid" needs no flowClust arguments; only its first two rows are kept.
  suppressMessages(sw_gate_template(
    tmp,
    template_type = "myeloid",
    debris_min = debris_window[1],
    singlet_wider_gate = singlet_wider_gate,
    pool_samples = pool_samples
  ))
  tmpl <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  tmpl <- tmpl[tmpl$alias %in% c("nonDebris", "singlets"), , drop = FALSE]

  # openCyto pools samples only when collapseDataForGating is TRUE and groupBy
  # is not empty: with an empty groupBy it fits one gate per sample. A number
  # N in groupBy puts N samples in one group, so the number of samples gives
  # one group for the whole run.
  group_by <- if (pool_samples) as.integer(n_samples) else NA_integer_
  tmpl$groupBy <- group_by

  if (is.finite(debris_window[2])) {
    tmpl$gating_args[tmpl$alias == "nonDebris"] <- paste0(
      "min = ",
      debris_window[1],
      ", max = ",
      debris_window[2]
    )
  }

  live <- if (viability_method == "mindensity") {
    .make_gate_row_raw(
      "live",
      "-",
      "singlets",
      viability_channel,
      "mindensity",
      collapseDataForGating = pool_samples,
      groupBy = group_by
    )
  } else {
    .make_gate_row_raw(
      "live",
      "+",
      "singlets",
      viability_channel,
      "boundary",
      paste0("min = -Inf, max = ", viability_cut),
      collapseDataForGating = pool_samples,
      groupBy = group_by
    )
  }
  rbind(tmpl, live)
}

# Lowest FSC-A inside debris_window, per sample. The cut sits there when the
# window holds one peak or no valley. Pooled: the gate is fitted once on the
# events of all samples, so every sample gets the one pooled floor. Not
# pooled: each sample gets the floor of its own events.
.sw_clean_window_floor <- function(frames, debris_window, pool_samples) {
  in_window <- lapply(frames, function(f) {
    v <- flowCore::exprs(f)[, "FSC-A"]
    v[v >= debris_window[1] & v <= debris_window[2]]
  })
  if (pool_samples) {
    rep(min(unlist(in_window)), length(frames))
  } else {
    vapply(in_window, min, numeric(1))
  }
}

# Error text for samples the margin step emptied: counts, share removed and
# the range that did it.
.sw_clean_margins_empty_message <- function(n_before, n_after, margin_range) {
  emptied <- names(n_after)[n_after == 0]
  lines <- vapply(
    emptied,
    function(sn) {
      sprintf(
        "  %s: %d events before, %d after, %s%% removed",
        sn,
        n_before[[sn]],
        n_after[[sn]],
        format(round(
          100 * (n_before[[sn]] - n_after[[sn]]) / n_before[[sn]],
          1
        ))
      )
    },
    character(1)
  )
  paste0(
    "No events left after step 'margins' in sample(s): ",
    paste(emptied, collapse = ", "),
    ".\n",
    paste(lines, collapse = "\n"),
    "\nmargin_range = c(",
    paste(
      vapply(margin_range, format, character(1), digits = 15),
      collapse = ", "
    ),
    "). Check it against the raw fluorescence densities of these samples."
  )
}

# One row of the gates table for one sample and one gate.
.sw_clean_gate_row <- function(sample, gate, gh, nodes) {
  node <- nodes[[switch(
    gate,
    nonDebris = "non_debris",
    singlets = "singlets",
    live = "live"
  )]]
  g <- flowWorkspace::gh_pop_get_gate(gh, node)
  if (methods::is(g, "rectangleGate")) {
    channels <- names(g@min)
    definition <- paste0(
      channels,
      ": [",
      g@min,
      ", ",
      g@max,
      "]",
      collapse = "; "
    )
  } else if (methods::is(g, "polygonGate")) {
    channels <- colnames(g@boundaries)
    definition <- paste(
      apply(g@boundaries, 1, function(r) {
        paste0("(", paste(signif(r, 6), collapse = ", "), ")")
      }),
      collapse = " "
    )
  } else {
    channels <- NA_character_
    definition <- NA_character_
  }
  tibble::tibble(
    sample = sample,
    gate = gate,
    class = class(g)[1],
    channels = paste(channels, collapse = ","),
    definition = definition
  )
}
