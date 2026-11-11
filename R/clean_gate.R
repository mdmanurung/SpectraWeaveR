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
  debris_window = list(
    controls = paste0(
      "Lower and upper FSC-A bound (raw scale) of the cut point of the ",
      "nonDebris gate. Events below the cut are debris."
    ),
    how = paste0(
      "Plot raw FSC-A as a density and FSC-A against SSC-A (ggplot2::geom_hex() ",
      "on sw_io_ff_to_tibble(ff)). Debris is the low-FSC cloud beside the main ",
      "cells. Look at one sample per batch, because the valley can move. Use ",
      "Inf as the upper bound for no cap."
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
#'   \item margin events removed on the raw scale
#'     (\code{\link{sw_filter_margins_peacoqc}});
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
#' @param debris_window Numeric length 2; required, no default. Lower and
#'   upper FSC-A bound (raw scale) of the cut point of the \code{nonDebris}
#'   gate. It is written as \code{min} and \code{max} of the openCyto
#'   \code{mindensity} method; use \code{Inf} for no upper bound. It sets
#'   the debris boundary, so it changes every later count.
#' @param singlet_wider_gate Logical; required, no default. \code{wider_gate}
#'   of the openCyto \code{singletGate}. \code{TRUE} keeps more events near the
#'   edge of the singlet diagonal.
#' @param pool_samples Logical; required, no default. \code{TRUE}: one gate per
#'   node is fitted on the pooled events of all samples and applied to each.
#'   \code{FALSE}: one gate per sample. It changes the cut of every sample.
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
#' Margins use PeacoQC's per-channel ranges from the FCS parameters (the
#' package's existing \code{\link{sw_filter_margins_peacoqc}} with no
#' channel specifications). The \code{debris_window} and the singlet gate act
#' on raw scatter; the viability cut acts on the arcsinh scale.
#'
#' A sample with no events left after any step stops the run, with the
#' sample and the step named. The run never skips a sample.
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
#'       or \code{"outside"} the shape; the live gate keeps events outside).}
#'     \item{\code{qc}}{The per-sample PeacoQC summary tibble of
#'       \code{\link{sw_qc_batch}} (\code{sample}, \code{n_before},
#'       \code{n_after}, \code{n_removed}, \code{pct_removed}), with the
#'       \code{IT_limit} and \code{MAD} used.}
#'     \item{\code{template}}{The gating template as a data frame.}
#'     \item{\code{parameters}}{The values of all required arguments.}
#'   }
#'
#' @examples
#' \donttest{
#' # Benchmark 05 values, shown for reference only; they are not defaults.
#' if (all(vapply(c("flowCore", "flowWorkspace", "openCyto", "PeacoQC"),
#'                requireNamespace, logical(1), quietly = TRUE))) {
#'   res <- sw_clean_run(
#'     "sample1.fcs",
#'     debris_window = c(2e5, 1.5e6),
#'     singlet_wider_gate = TRUE,
#'     pool_samples = TRUE,
#'     viability_channel = "LIVE DEAD NIR-A",
#'     viability_method = "mindensity",
#'     viability_cut = NULL,
#'     IT_limit = 0.55,
#'     MAD = 6,
#'     cofactor = 6000
#'   )
#'   res$audit
#' }
#' }
#'
#' @seealso \code{\link{sw_gate_template}}, \code{\link{sw_gate_run}},
#'   \code{\link{sw_qc_run}}, \code{\link{sw_io_event_audit}}
#' @export
sw_clean_run <- function(
  x,
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
  ffs <- lapply(ffs, sw_filter_margins_peacoqc)
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
  qc <- sw_qc_batch(ffs, channels = fluor, IT_limit = IT_limit, MAD = MAD)
  ffs <- qc$cleaned[sample_names]
  record("qc", ffs)
  qc_summary <- qc$summary
  qc_summary$IT_limit <- IT_limit
  qc_summary$MAD <- MAD

  # --- 4. Gates: nonDebris, singlets, live --------------------------------
  template <- .sw_clean_template(
    debris_window,
    singlet_wider_gate,
    pool_samples,
    viability_channel,
    viability_method,
    viability_cut
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

  list(
    events = events,
    audit = audit,
    audit_total = audit_total,
    gates = gates,
    qc = qc_summary,
    template = template,
    parameters = list(
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
  viability_cut
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
      collapseDataForGating = pool_samples
    )
  } else {
    .make_gate_row_raw(
      "live",
      "+",
      "singlets",
      viability_channel,
      "boundary",
      paste0("min = -Inf, max = ", viability_cut),
      collapseDataForGating = pool_samples
    )
  }
  rbind(tmpl, live)
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
