#' @title Spectral Unmixing Wrappers
#'
#' @description
#' Functions for spectral unmixing of flow cytometry data using AutoSpectral,
#' and for loading pre-unmixed FCS files (e.g., from SpectroFlo).
#'
#' The wrappers target AutoSpectral 1.8.8 and follow its automated workflow:
#' \enumerate{
#'   \item \code{\link{sw_unmix_setup}}: \code{get.autospectral.param()},
#'     \code{create.control.file()} or \code{check.control.file()}
#'   \item \code{\link{sw_unmix_prepare}}: \code{get.spectra.automated()}
#'     and \code{reload.flow.control()} (or the legacy gating path)
#'   \item \code{\link{sw_unmix_extract_af}}: \code{get.af.spectra()}
#'   \item \code{\link{sw_unmix_extract_variants}}:
#'     \code{get.spectral.variants()}
#'   \item \code{\link{sw_unmix_run}}: \code{unmix.fcs()} or
#'     \code{unmix.folder()}
#'   \item \code{\link{sw_unmix_pipeline}}: one-call orchestrator
#' }
#'
#' Each AutoSpectral call is built by an internal \code{.sw_args_*()}
#' function. It returns the argument list under AutoSpectral's own names and
#' leaves out every argument the caller did not set, so AutoSpectral's
#' defaults apply. A parity test checks these names against the 1.8.8
#' formals.
#'
#' @name unmix
#' @keywords internal
NULL

# ---- Helper: check AutoSpectral availability ----
.check_autospectral <- function() {
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    stop(
      "Package 'AutoSpectral' is required. ",
      "Install it from GitHub: ",
      "remotes::install_github('DrCytometer/AutoSpectral')",
      call. = FALSE
    )
  }
}

# ---- Helper: supported cytometers (get.autospectral.param 1.8.8) ----
.valid_cytometers <- function() {
  c(
    "aurora",
    "auroraNL",
    "id7000",
    "a5se",
    "opteon",
    "mosaic",
    "xenith",
    "cytostellar",
    "a8",
    "s8",
    "discover"
  )
}

# ---- Helper: valid unmixing methods (unmix.fcs 1.8.8) ----
.valid_unmix_methods <- function() {
  c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson")
}

# ---- Helper: override asp output directories ----
.override_asp_dirs <- function(asp, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  asp$figure.gate.dir <- file.path(output_dir, "figure_gate")
  asp$figure.af.dir <- file.path(output_dir, "figure_autofluorescence")
  asp$figure.clean.control.dir <- file.path(output_dir, "figure_clean_controls")
  asp$figure.spectral.ribbon.dir <- file.path(
    output_dir,
    "figure_spectral_ribbon"
  )
  asp$figure.spectra.dir <- file.path(output_dir, "figure_spectra")
  asp$figure.similarity.heatmap.dir <- file.path(
    output_dir,
    "figure_similarity_heatmap"
  )
  asp$figure.scatter.dir.base <- file.path(output_dir, "figure_scatter")
  asp$table.spectra.dir <- file.path(output_dir, "table_spectra")
  asp$unmixed.fcs.dir <- file.path(output_dir, "unmixed_fcs")
  asp$variant.dir <- file.path(output_dir, "figure_spectral_variants")
  asp$fix.unmixing.dir <- file.path(output_dir, "figure_fix_my_unmix")

  asp
}

# ---- Helpers: validate optional arguments ----
.check_flag <- function(x, name) {
  if (!is.null(x) && (!is.logical(x) || length(x) != 1 || is.na(x))) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
}

.check_threads <- function(threads) {
  if (
    !is.null(threads) &&
      (!is.numeric(threads) || length(threads) != 1 || threads < 1)
  ) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }
}

.check_som_dim <- function(som_dim) {
  if (
    !is.null(som_dim) &&
      (!is.numeric(som_dim) || length(som_dim) != 1 || som_dim < 2)
  ) {
    stop("'som_dim' must be a single integer >= 2.", call. = FALSE)
  }
}

.check_chunk_size <- function(chunk_size) {
  if (
    !is.null(chunk_size) &&
      (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 1)
  ) {
    stop("'chunk_size' must be a positive number.", call. = FALSE)
  }
}

# ============================================================================
# Argument builders: one per AutoSpectral call
# ============================================================================

# Drop NULL entries so AutoSpectral's own defaults apply
.sw_drop_null <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

# Append user arguments from `...`; they keep AutoSpectral's names
.sw_add_dots <- function(args, dots, fn) {
  if (length(dots) == 0) {
    return(args)
  }
  dot_names <- names(dots)
  if (is.null(dot_names) || any(dot_names == "")) {
    stop(
      "Arguments passed through '...' to AutoSpectral::",
      fn,
      "() must be named.",
      call. = FALSE
    )
  }
  clash <- intersect(dot_names, names(args))
  if (length(clash) > 0) {
    stop(
      "Argument(s) ",
      paste(clash, collapse = ", "),
      " for AutoSpectral::",
      fn,
      "() are set by the wrapper; use the wrapper argument instead of '...'.",
      call. = FALSE
    )
  }
  c(args, dots)
}

# Pass figures = FALSE only when the setup turned figures off
.sw_figures_arg <- function(asp, dots) {
  if ("figures" %in% names(dots) || isTRUE(asp$figures)) NULL else FALSE
}

.sw_args_get_autospectral_param <- function(cytometer, figures) {
  list(cytometer = cytometer, figures = figures)
}

.sw_args_create_control_file <- function(control_dir, asp, output_dir, legacy) {
  .sw_drop_null(list(
    control.dir = control_dir,
    asp = asp,
    output.dir = output_dir,
    legacy = if (legacy) TRUE
  ))
}

.sw_args_check_control_file <- function(
  control_dir,
  control_file,
  asp,
  legacy
) {
  .sw_drop_null(list(
    control.dir = control_dir,
    control.def.file = control_file,
    asp = asp,
    legacy = if (legacy) TRUE
  ))
}

.sw_args_get_spectra_automated <- function(setup, dots = list()) {
  args <- .sw_drop_null(list(
    control.dir = setup$control_dir,
    control.def.file = setup$control_file,
    asp = setup$asp,
    figures = .sw_figures_arg(setup$asp, dots)
  ))
  .sw_add_dots(args, dots, "get.spectra.automated")
}

.sw_args_reload_flow_control <- function(setup) {
  list(
    control.dir = setup$control_dir,
    control.def.file = setup$control_file,
    asp = setup$asp
  )
}

.sw_args_define_flow_control <- function(
  setup,
  gating_system,
  gate_list,
  parallel,
  threads
) {
  .sw_drop_null(list(
    control.dir = setup$control_dir,
    control.def.file = setup$control_file,
    asp = setup$asp,
    gate = TRUE,
    gating.system = gating_system,
    gate.list = gate_list,
    parallel = parallel,
    verbose = TRUE,
    threads = threads
  ))
}

.sw_args_clean_controls <- function(
  flow_control,
  asp,
  af_remove,
  parallel,
  threads
) {
  .sw_drop_null(list(
    flow.control = flow_control,
    asp = asp,
    af.remove = af_remove,
    universal.negative = TRUE,
    downsample = TRUE,
    parallel = parallel,
    verbose = TRUE,
    threads = threads
  ))
}

.sw_args_get_fluorophore_spectra <- function(flow_control, asp, clean) {
  .sw_drop_null(list(
    flow.control = flow_control,
    asp = asp,
    use.clean.expr = clean,
    figures = asp$figures
  ))
}

.sw_args_get_af_spectra <- function(
  unstained_fcs,
  asp,
  spectra,
  title = NULL,
  som_dim = NULL,
  refine = NULL,
  parallel = NULL,
  threads = NULL,
  dots = list()
) {
  # A user title in `...` wins over the per-tissue title
  if ("title" %in% names(dots)) {
    title <- NULL
  }
  args <- .sw_drop_null(list(
    unstained.sample = unstained_fcs,
    asp = asp,
    spectra = spectra,
    som.dim = if (!is.null(som_dim)) as.integer(som_dim),
    figures = .sw_figures_arg(asp, dots),
    title = title,
    refine = refine,
    parallel = parallel,
    threads = if (!is.null(threads)) as.integer(threads)
  ))
  .sw_add_dots(args, dots, "get.af.spectra")
}

.sw_args_get_spectral_variants <- function(
  setup,
  spectra,
  n_cells = NULL,
  som_dim = NULL,
  parallel = NULL,
  threads = NULL,
  dots = list()
) {
  args <- .sw_drop_null(list(
    control.dir = setup$control_dir,
    control.def.file = setup$control_file,
    asp = setup$asp,
    spectra = spectra,
    figures = .sw_figures_arg(setup$asp, dots),
    parallel = parallel,
    threads = if (!is.null(threads)) as.integer(threads),
    n.cells = if (!is.null(n_cells)) as.integer(n_cells),
    som.dim = if (!is.null(som_dim)) as.integer(som_dim)
  ))
  .sw_add_dots(args, dots, "get.spectral.variants")
}

# Shared part of unmix.fcs() and unmix.folder() calls
.sw_args_unmix_common <- function(
  spectra,
  asp,
  flow_control,
  method,
  af_spectra = NULL,
  spectra_variants = NULL,
  output_dir = NULL,
  file_suffix = NULL,
  speed = NULL,
  parallel = NULL,
  threads = NULL,
  chunk_size = NULL
) {
  .sw_drop_null(list(
    spectra = spectra,
    asp = asp,
    flow.control = flow_control,
    method = method,
    af.spectra = af_spectra,
    spectra.variants = spectra_variants,
    output.dir = output_dir,
    file.suffix = file_suffix,
    speed = speed,
    parallel = parallel,
    threads = if (!is.null(threads)) as.integer(threads),
    chunk.size = if (!is.null(chunk_size)) as.integer(chunk_size)
  ))
}

.sw_args_unmix_fcs <- function(fcs_file, ..., dots = list()) {
  args <- c(list(fcs.file = fcs_file), .sw_args_unmix_common(...))
  .sw_add_dots(args, dots, "unmix.fcs")
}

.sw_args_unmix_folder <- function(fcs_dir, ..., dots = list()) {
  args <- c(list(fcs.dir = fcs_dir), .sw_args_unmix_common(...))
  .sw_add_dots(args, dots, "unmix.folder")
}

# Call an AutoSpectral function; stop on argument names it does not accept.
# Several AutoSpectral functions take `...`, which would swallow a misspelt
# or removed argument without an error.
.sw_call_autospectral <- function(fn, args, ns = "AutoSpectral") {
  f <- getExportedValue(ns, fn)
  accepted <- setdiff(names(formals(f)), "...")
  unknown <- setdiff(names(args), accepted)
  if (length(unknown) > 0) {
    stop(
      ns,
      "::",
      fn,
      "() has no argument(s): ",
      paste(unknown, collapse = ", "),
      ". The wrappers target AutoSpectral 1.8.8; installed version: ",
      as.character(utils::packageVersion(ns)),
      ".",
      call. = FALSE
    )
  }
  do.call(f, args)
}

# ============================================================================
# Function 1: sw_unmix_setup
# ============================================================================

#' Initialize AutoSpectral Parameters and Control File
#'
#' Sets up the AutoSpectral parameter list for a given cytometer, and either
#' creates a draft control file CSV from the FCS filenames in
#' \code{control_dir} or checks an existing one.
#'
#' @param control_dir Character path to the directory containing single-stain
#'   control FCS files.
#' @param cytometer Character string specifying the cytometer type. One of
#'   \code{"aurora"} (default), \code{"auroraNL"}, \code{"id7000"},
#'   \code{"a5se"}, \code{"opteon"}, \code{"mosaic"}, \code{"xenith"},
#'   \code{"cytostellar"}, \code{"a8"}, \code{"s8"} or \code{"discover"}
#'   (the values \code{AutoSpectral::get.autospectral.param()} accepts in
#'   1.8.8).
#' @param control_file Character path to an existing control file CSV, or
#'   \code{NULL} (default) to auto-generate one from the FCS file names.
#' @param output_dir Character path for output files (figures, tables).
#'   Default: \code{"SpectraWeaveR_unmix"}.
#' @param figures Logical; whether to generate QC plots throughout the
#'   workflow. Default: \code{TRUE}.
#' @param workflow Character; \code{"automated"} (default) or
#'   \code{"legacy"}. Sets AutoSpectral's \code{legacy} argument of
#'   \code{create.control.file()} and \code{check.control.file()}. A
#'   control file for the legacy gating path (\code{sw_unmix_prepare(
#'   workflow = "legacy")}) needs the gating columns that only
#'   \code{legacy = TRUE} writes.
#'
#' @return A named list (\code{sw_setup}) with components:
#'   \describe{
#'     \item{\code{asp}}{The AutoSpectral parameter list, with output
#'       directories routed under \code{output_dir}.}
#'     \item{\code{control_file}}{The path to the control CSV.}
#'     \item{\code{control_dir}}{The control directory path.}
#'     \item{\code{workflow}}{The workflow the control file was made or
#'       checked for.}
#'     \item{\code{control_issues}}{What \code{check.control.file()}
#'       returned (\code{TRUE} or a data frame of issues), or \code{NULL}
#'       for a new draft file.}
#'   }
#'
#' @details
#' This is Step 1 of the SpectraWeaveR unmixing workflow. If
#' \code{control_file = NULL}, a draft CSV is written to \code{output_dir}
#' and a message asks the user to review it before
#' \code{\link{sw_unmix_prepare}}.
#'
#' \code{check.control.file()} runs in AutoSpectral's default non-strict
#' mode: it returns a table of issues instead of stopping. The wrapper
#' warns when that table holds errors and keeps it in
#' \code{control_issues}.
#'
#' @examples
#' \dontrun{
#' setup <- sw_unmix_setup("path/to/controls", cytometer = "aurora")
#' }
#'
#' @export
sw_unmix_setup <- function(
  control_dir,
  cytometer = "aurora",
  control_file = NULL,
  output_dir = "SpectraWeaveR_unmix",
  figures = TRUE,
  workflow = c("automated", "legacy")
) {
  # Validate inputs before checking dependency
  if (!is.character(control_dir) || length(control_dir) != 1) {
    stop("'control_dir' must be a single directory path.", call. = FALSE)
  }
  if (!dir.exists(control_dir)) {
    stop("Control directory does not exist: ", control_dir, call. = FALSE)
  }

  if (!is.character(cytometer) || length(cytometer) != 1) {
    stop("'cytometer' must be a single character string.", call. = FALSE)
  }
  valid <- .valid_cytometers()
  if (!cytometer %in% valid) {
    stop(
      "'cytometer' must be one of: ",
      paste(valid, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a single directory path.", call. = FALSE)
  }

  if (!is.logical(figures) || length(figures) != 1) {
    stop("'figures' must be TRUE or FALSE.", call. = FALSE)
  }

  workflow <- match.arg(workflow)
  legacy <- workflow == "legacy"

  if (!is.null(control_file)) {
    if (!is.character(control_file) || length(control_file) != 1) {
      stop("'control_file' must be a single file path.", call. = FALSE)
    }
    if (!file.exists(control_file)) {
      stop("Control file does not exist: ", control_file, call. = FALSE)
    }
  }

  .check_autospectral()

  # Initialize AutoSpectral parameters
  message("Initializing AutoSpectral parameters for cytometer: ", cytometer)
  asp <- .sw_call_autospectral(
    "get.autospectral.param",
    .sw_args_get_autospectral_param(cytometer, figures)
  )

  # Override output directories to route under output_dir
  asp <- .override_asp_dirs(asp, output_dir)

  control_issues <- NULL
  if (is.null(control_file)) {
    message("Creating draft control file from FCS files in: ", control_dir)
    # create.control.file() never overwrites; the new CSV is the one added
    csv_before <- list.files(
      output_dir,
      pattern = "^fcs_control_file.*\\.csv$",
      full.names = TRUE
    )
    .sw_call_autospectral(
      "create.control.file",
      .sw_args_create_control_file(control_dir, asp, output_dir, legacy)
    )
    csv_after <- list.files(
      output_dir,
      pattern = "^fcs_control_file.*\\.csv$",
      full.names = TRUE
    )
    control_file <- setdiff(csv_after, csv_before)
    if (length(control_file) != 1) {
      stop("Failed to create control file in: ", output_dir, call. = FALSE)
    }

    message("\n*** IMPORTANT ***")
    message("A draft control file has been created at: ", control_file)
    message(
      "Please review and edit it before proceeding to sw_unmix_prepare()."
    )
    message("See AutoSpectral documentation for control file details:")
    message(
      "  https://drcytometer.github.io/AutoSpectral/articles/",
      "02_Control_File_example.html"
    )
  } else {
    message("Checking control file: ", control_file)
    control_issues <- .sw_call_autospectral(
      "check.control.file",
      .sw_args_check_control_file(control_dir, control_file, asp, legacy)
    )
    # Non-strict mode returns issues instead of stopping
    if (
      is.data.frame(control_issues) && any(control_issues$severity == "error")
    ) {
      warning(
        "check.control.file() found errors: ",
        paste(
          unique(control_issues$rule[control_issues$severity == "error"]),
          collapse = ", "
        ),
        ". See setup$control_issues.",
        call. = FALSE
      )
    } else if (isTRUE(control_issues)) {
      message("Control file check passed.")
    }
  }

  structure(
    list(
      asp = asp,
      control_file = control_file,
      control_dir = control_dir,
      workflow = workflow,
      control_issues = control_issues
    ),
    class = "sw_setup"
  )
}

# ============================================================================
# Function 2: sw_unmix_prepare
# ============================================================================

#' Prepare Controls and Extract Fluorophore Spectra
#'
#' Extracts fluorophore spectral signatures from the single-stain controls
#' and builds the \code{flow.control} structure needed for unmixing.
#'
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param workflow Character; \code{"automated"} (default) or
#'   \code{"legacy"}. \code{"automated"} is AutoSpectral's recommended path:
#'   \code{get.spectra.automated()} then \code{reload.flow.control()}.
#'   \code{"legacy"} stays available: \code{define.flow.control()},
#'   \code{clean.controls()} and \code{get.fluorophore.spectra()}. The
#'   legacy path needs a control file with gating columns (see
#'   \code{workflow} in \code{\link{sw_unmix_setup}}).
#' @param gating_system Legacy only. Character; \code{"landmarks"} (default)
#'   or \code{"density"}. Note that AutoSpectral's own default is
#'   \code{"density"}.
#' @param gate_list Legacy only. Optional named list of pre-defined gates
#'   (from \code{AutoSpectral::define.gate.landmarks()} or
#'   \code{AutoSpectral::define.gate.density()}). Default: \code{NULL}.
#' @param clean Legacy only. Logical; whether to run
#'   \code{clean.controls()}. Default: \code{TRUE}.
#' @param af_remove Legacy only. Logical; whether to remove autofluorescence
#'   intrusions during cleaning. Default: \code{TRUE}.
#' @param parallel Legacy only. Logical; parallel processing. Default:
#'   \code{FALSE}.
#' @param threads Legacy only. Integer or \code{NULL}; number of threads.
#'   Default: \code{NULL} (AutoSpectral's default).
#' @param ... Automated only. Further arguments for
#'   \code{AutoSpectral::get.spectra.automated()}, under AutoSpectral's own
#'   names (e.g. \code{n.candidates}, \code{refine}). Nothing is set for
#'   you, except \code{figures = FALSE} when the setup turned figures off.
#'
#' @return A named list (\code{sw_controls}) with components:
#'   \describe{
#'     \item{\code{flow_control}}{The full \code{flow.control} data
#'       structure from AutoSpectral.}
#'     \item{\code{spectra}}{Fluorophore spectral signature matrix
#'       (fluorophores in rows, detectors in columns), L-infinity
#'       normalized.}
#'     \item{\code{setup}}{Pass-through of the \code{sw_setup} object.}
#'     \item{\code{workflow}}{The workflow used.}
#'   }
#'
#' @details
#' This is Step 2 of the SpectraWeaveR unmixing workflow. The legacy-only
#' arguments are ignored, with a warning, when \code{workflow =
#' "automated"}.
#'
#' @export
sw_unmix_prepare <- function(
  setup,
  workflow = c("automated", "legacy"),
  gating_system = c("landmarks", "density"),
  gate_list = NULL,
  clean = TRUE,
  af_remove = TRUE,
  parallel = FALSE,
  threads = NULL,
  ...
) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop(
      "'setup' must be an sw_setup object from sw_unmix_setup().",
      call. = FALSE
    )
  }

  workflow <- match.arg(workflow)
  gating_given <- !missing(gating_system) # match.arg() below makes missing() FALSE
  gating_system <- match.arg(gating_system)

  if (!is.logical(clean) || length(clean) != 1) {
    stop("'clean' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(af_remove) || length(af_remove) != 1) {
    stop("'af_remove' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("'parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  .check_threads(threads)

  dots <- list(...)
  if (workflow == "automated") {
    # Name the legacy-only arguments the caller set; they have no effect here
    legacy_set <- c(
      gating_system = gating_given,
      gate_list = !missing(gate_list),
      clean = !missing(clean),
      af_remove = !missing(af_remove),
      parallel = !missing(parallel),
      threads = !missing(threads)
    )
    if (any(legacy_set)) {
      warning(
        "Ignored with workflow = 'automated': ",
        paste(names(legacy_set)[legacy_set], collapse = ", "),
        ". These apply to workflow = 'legacy' only.",
        call. = FALSE
      )
    }
  } else {
    if (length(dots) > 0) {
      stop(
        "'...' is forwarded to get.spectra.automated() and is not used ",
        "with workflow = 'legacy'.",
        call. = FALSE
      )
    }
    if (identical(setup$workflow, "automated")) {
      stop(
        "The control file was made for the automated workflow and lacks ",
        "gating columns. Run sw_unmix_setup(workflow = 'legacy') first.",
        call. = FALSE
      )
    }
  }

  .check_autospectral()

  asp <- setup$asp

  if (workflow == "automated") {
    message("\n=== Extracting fluorophore spectra (get.spectra.automated) ===")
    spectra <- .sw_call_autospectral(
      "get.spectra.automated",
      .sw_args_get_spectra_automated(setup, dots)
    )
    message("\n=== Building flow.control (reload.flow.control) ===")
    flow_control <- .sw_call_autospectral(
      "reload.flow.control",
      .sw_args_reload_flow_control(setup)
    )
  } else {
    # Step 2a: Load controls, apply gating, organize into flow.control
    message("\n=== Loading and gating controls ===")
    flow_control <- .sw_call_autospectral(
      "define.flow.control",
      .sw_args_define_flow_control(
        setup,
        gating_system,
        gate_list,
        parallel,
        threads
      )
    )

    # Step 2b: Clean controls (optional)
    if (clean) {
      message("\n=== Cleaning controls ===")
      flow_control <- .sw_call_autospectral(
        "clean.controls",
        .sw_args_clean_controls(flow_control, asp, af_remove, parallel, threads)
      )
    }

    # Step 2c: Extract fluorophore spectra
    message("\n=== Extracting fluorophore spectra ===")
    spectra <- .sw_call_autospectral(
      "get.fluorophore.spectra",
      .sw_args_get_fluorophore_spectra(flow_control, asp, clean)
    )
  }

  message(
    "Extracted spectra for ",
    nrow(spectra),
    " fluorophores across ",
    ncol(spectra),
    " detectors"
  )

  structure(
    list(
      flow_control = flow_control,
      spectra = spectra,
      setup = setup,
      workflow = workflow
    ),
    class = "sw_controls"
  )
}

# ============================================================================
# Function 3: sw_unmix_extract_af
# ============================================================================

#' Extract Per-Cell Autofluorescence Spectra
#'
#' Extracts autofluorescence (AF) spectral variation from unstained samples
#' with \code{AutoSpectral::get.af.spectra()}. For multi-tissue experiments,
#' provide a named list of unstained files to get tissue-matched AF spectra.
#'
#' @param unstained_fcs Character path to an unstained FCS file, or a named
#'   list of paths for tissue-specific AF (e.g.,
#'   \code{list(spleen = "spleen_unstained.fcs",
#'              lung = "lung_unstained.fcs")}).
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param refine Logical or \code{NULL}; second-pass refinement of
#'   high-error cells. \code{NULL} (default) uses AutoSpectral's default
#'   (\code{TRUE} in 1.8.8).
#' @param som_dim Integer or \code{NULL}; SOM grid dimension. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{10} in 1.8.8).
#' @param parallel Logical or \code{NULL}; parallel processing. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{TRUE} in 1.8.8).
#' @param threads Integer or \code{NULL}; number of threads. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{if (parallel) 0 else 1},
#'   where 0 means all cores).
#' @param ... Further arguments for \code{AutoSpectral::get.af.spectra()},
#'   under AutoSpectral's own names (e.g. \code{af.assign.method}).
#'
#' @return If \code{unstained_fcs} is a single path: a matrix of AF spectra
#'   (spectra in rows, detectors in columns). If a named list: a named list
#'   of such matrices.
#'
#' @details
#' This is Step 3 of the SpectraWeaveR unmixing workflow. For a named list,
#' each tissue's plots and table get the title \code{"AF <tissue>"} unless
#' \code{title} is passed in \code{...}.
#'
#' Installation of \code{AutoSpectralRcpp} is strongly recommended for
#' faster processing.
#'
#' @export
sw_unmix_extract_af <- function(
  unstained_fcs,
  setup,
  spectra,
  refine = NULL,
  som_dim = NULL,
  parallel = NULL,
  threads = NULL,
  ...
) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop(
      "'setup' must be an sw_setup object from sw_unmix_setup().",
      call. = FALSE
    )
  }

  if (!is.matrix(spectra)) {
    stop(
      "'spectra' must be a matrix of fluorophore spectral signatures.",
      call. = FALSE
    )
  }

  .check_flag(refine, "refine")
  .check_som_dim(som_dim)
  .check_flag(parallel, "parallel")
  .check_threads(threads)

  # Validate unstained_fcs type early (before AutoSpectral check)
  if (
    !(is.character(unstained_fcs) && length(unstained_fcs) == 1) &&
      !is.list(unstained_fcs)
  ) {
    stop(
      "'unstained_fcs' must be a single file path or a named list of file paths.",
      call. = FALSE
    )
  }

  .check_autospectral()

  # Warn if AutoSpectralRcpp is not available
  if (!requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    message(
      "Note: Install AutoSpectralRcpp for faster AF extraction: ",
      "remotes::install_github('DrCytometer/AutoSpectralRcpp')"
    )
  }

  asp <- setup$asp
  dots <- list(...)

  if (is.character(unstained_fcs)) {
    if (!file.exists(unstained_fcs)) {
      stop("Unstained FCS file does not exist: ", unstained_fcs, call. = FALSE)
    }
    message("Extracting AF spectra from: ", basename(unstained_fcs))
    af_spectra <- .sw_call_autospectral(
      "get.af.spectra",
      .sw_args_get_af_spectra(
        unstained_fcs,
        asp,
        spectra,
        som_dim = som_dim,
        refine = refine,
        parallel = parallel,
        threads = threads,
        dots = dots
      )
    )
    message("Extracted ", nrow(af_spectra), " AF spectral signatures")
    return(af_spectra)
  }

  # Named list of unstained files (multi-tissue)
  if (is.null(names(unstained_fcs)) || any(names(unstained_fcs) == "")) {
    stop(
      "'unstained_fcs' list must have non-empty names (e.g., tissue types).",
      call. = FALSE
    )
  }

  af_list <- list()
  for (tissue in names(unstained_fcs)) {
    fcs_path <- unstained_fcs[[tissue]]
    if (!is.character(fcs_path) || length(fcs_path) != 1) {
      stop(
        "Each element of 'unstained_fcs' must be a single file path. Problem with: ",
        tissue,
        call. = FALSE
      )
    }
    if (!file.exists(fcs_path)) {
      stop(
        "Unstained FCS file does not exist for '",
        tissue,
        "': ",
        fcs_path,
        call. = FALSE
      )
    }

    message(
      "Extracting AF spectra for '",
      tissue,
      "' from: ",
      basename(fcs_path)
    )
    af_list[[tissue]] <- .sw_call_autospectral(
      "get.af.spectra",
      .sw_args_get_af_spectra(
        fcs_path,
        asp,
        spectra,
        title = paste("AF", tissue),
        som_dim = som_dim,
        refine = refine,
        parallel = parallel,
        threads = threads,
        dots = dots
      )
    )
    message(
      "  -> ",
      nrow(af_list[[tissue]]),
      " AF signatures for '",
      tissue,
      "'"
    )
  }
  af_list
}

# ============================================================================
# Function 4: sw_unmix_extract_variants
# ============================================================================

#' Extract Spectral Variants for Per-Cell Fluorophore Optimization
#'
#' Maps fluorophore emission variability from single-stain controls with
#' \code{AutoSpectral::get.spectral.variants()}. The resulting spectral
#' variants are used during unmixing to optimize fluorophore signatures on a
#' per-cell basis, reducing spillover spread.
#'
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param n_cells Integer or \code{NULL}; maximum positive events per
#'   fluorophore. \code{NULL} (default) uses AutoSpectral's default
#'   (\code{10000} in 1.8.8).
#' @param som_dim Integer or \code{NULL}; SOM grid dimension. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{5} in 1.8.8).
#' @param parallel Logical or \code{NULL}; parallel processing. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{TRUE} in 1.8.8).
#' @param threads Integer or \code{NULL}; number of threads. \code{NULL}
#'   (default) uses AutoSpectral's default (\code{asp$worker.process.n}).
#' @param ... Further arguments for
#'   \code{AutoSpectral::get.spectral.variants()}, under AutoSpectral's own
#'   names. Use \code{unstained.sample = "path.fcs"} when the control file
#'   has no \code{"AF"} row.
#'
#' @return The list from \code{get.spectral.variants()} with class
#'   \code{sw_variants} (elements include \code{thresholds} and
#'   \code{variants}).
#'
#' @details
#' This is Step 4 of the SpectraWeaveR unmixing workflow. It is optional
#' but recommended: per-cell fluorophore optimization can reduce spillover
#' spread, especially for tandem dyes.
#'
#' The wrapper takes no AF spectra. Since AutoSpectral 1.6.0,
#' \code{get.spectral.variants()} derives AF in situ from the unstained
#' files in the control table and ignores \code{af.spectra} and
#' \code{refine}, so the earlier \code{af_spectra} and \code{refine}
#' arguments were removed.
#'
#' @export
sw_unmix_extract_variants <- function(
  setup,
  spectra,
  n_cells = NULL,
  som_dim = NULL,
  parallel = NULL,
  threads = NULL,
  ...
) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop(
      "'setup' must be an sw_setup object from sw_unmix_setup().",
      call. = FALSE
    )
  }
  if (!is.matrix(spectra)) {
    stop(
      "'spectra' must be a matrix of fluorophore spectral signatures.",
      call. = FALSE
    )
  }
  dots <- list(...)
  # Catch the removed wrapper arguments by name
  removed <- intersect(names(dots), c("af_spectra", "refine", "af.spectra"))
  if (length(removed) > 0) {
    stop(
      "Removed argument(s): ",
      paste(removed, collapse = ", "),
      ". AutoSpectral 1.8.8 derives AF in situ for spectral variants and ",
      "has no second-pass refinement here.",
      call. = FALSE
    )
  }
  .check_som_dim(som_dim)
  if (
    !is.null(n_cells) &&
      (!is.numeric(n_cells) || length(n_cells) != 1 || n_cells < 1)
  ) {
    stop("'n_cells' must be a positive integer.", call. = FALSE)
  }
  .check_flag(parallel, "parallel")
  .check_threads(threads)

  .check_autospectral()

  message("\n=== Extracting spectral variants ===")
  variants <- .sw_call_autospectral(
    "get.spectral.variants",
    .sw_args_get_spectral_variants(
      setup,
      spectra,
      n_cells = n_cells,
      som_dim = som_dim,
      parallel = parallel,
      threads = threads,
      dots = dots
    )
  )

  message(
    "Extracted spectral variants for ",
    length(variants$variants),
    " fluorophores"
  )

  structure(variants, class = "sw_variants")
}

# ============================================================================
# Function 5: sw_unmix_run
# ============================================================================

#' Unmix Spectral Flow Cytometry Data
#'
#' Performs spectral unmixing on fully-stained sample FCS files with
#' \code{AutoSpectral::unmix.fcs()} (files) or
#' \code{AutoSpectral::unmix.folder()} (a directory).
#'
#' @param input Character path to a single FCS file, a character vector of
#'   FCS file paths, or a directory containing FCS files.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param flow_control The \code{flow.control} list from
#'   \code{\link{sw_unmix_prepare}}.
#' @param method Character; unmixing algorithm, as in \code{unmix.fcs()}
#'   1.8.8. One of \code{"AutoSpectral"} (default), \code{"OLS"},
#'   \code{"WLS"}, \code{"Poisson"} or \code{"FastPoisson"}.
#' @param af_spectra A matrix of AF spectra with at least 2 rows (from
#'   \code{\link{sw_unmix_extract_af}}). Required for
#'   \code{method = "AutoSpectral"}.
#' @param spectra_variants A spectral variants list (from
#'   \code{\link{sw_unmix_extract_variants}}). Optional; enables per-cell
#'   fluorophore optimization.
#' @param speed Character; precision-speed trade-off for per-cell fluorophore
#'   optimization. \code{"fast"} (1 variant/cell, default), \code{"medium"}
#'   (3), or \code{"slow"} (10).
#' @param output_dir Character path for output unmixed FCS files. Default:
#'   \code{NULL} (uses \code{setup$asp$unmixed.fcs.dir}).
#' @param file_suffix Character; optional suffix appended to output file
#'   names. Default: \code{NULL}.
#' @param parallel Logical or \code{NULL}; parallel processing. \code{NULL}
#'   (default) uses AutoSpectral's default, which is \code{TRUE} for
#'   \code{unmix.fcs()} and \code{FALSE} for \code{unmix.folder()} in 1.8.8.
#' @param threads Integer or \code{NULL}; number of threads. \code{NULL}
#'   (default) uses AutoSpectral's default.
#' @param chunk_size Number or \code{NULL}; events per processing chunk.
#'   \code{NULL} (default) uses AutoSpectral's default (\code{2e6}).
#' @param ... Further arguments for \code{unmix.fcs()} or
#'   \code{unmix.folder()}, under AutoSpectral's own names (e.g.
#'   \code{weighted}, \code{include.raw}).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{output_dir}}{Path to the directory containing unmixed
#'       FCS files.}
#'     \item{\code{output_files}}{Character vector of FCS files this call
#'       wrote or replaced in \code{output_dir}.}
#'     \item{\code{spectra}}{The spectra matrix used for unmixing.}
#'     \item{\code{method}}{The unmixing method used.}
#'   }
#'
#' @details
#' This is Step 5 of the SpectraWeaveR unmixing workflow. The function
#' auto-detects whether \code{input} is a single file, multiple files, or
#' a directory.
#'
#' Installation of \code{AutoSpectralRcpp} is strongly recommended for
#' faster processing; \code{"FastPoisson"} uses it.
#'
#' @export
sw_unmix_run <- function(
  input,
  spectra,
  setup,
  flow_control,
  method = c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson"),
  af_spectra = NULL,
  spectra_variants = NULL,
  speed = c("fast", "medium", "slow"),
  output_dir = NULL,
  file_suffix = NULL,
  parallel = NULL,
  threads = NULL,
  chunk_size = NULL,
  ...
) {
  # Validate inputs before checking dependency
  if (!is.character(input) || length(input) == 0) {
    stop(
      "'input' must be a non-empty character vector (file path(s) or directory).",
      call. = FALSE
    )
  }
  if (!is.matrix(spectra)) {
    stop(
      "'spectra' must be a matrix of fluorophore spectral signatures.",
      call. = FALSE
    )
  }
  if (!inherits(setup, "sw_setup")) {
    stop(
      "'setup' must be an sw_setup object from sw_unmix_setup().",
      call. = FALSE
    )
  }
  if (!is.list(flow_control)) {
    stop(
      "'flow_control' must be the flow.control list from sw_unmix_prepare().",
      call. = FALSE
    )
  }
  .check_flag(parallel, "parallel")
  .check_threads(threads)
  .check_chunk_size(chunk_size)
  if (
    !is.null(file_suffix) &&
      (!is.character(file_suffix) || length(file_suffix) != 1)
  ) {
    stop(
      "'file_suffix' must be NULL or a single character string.",
      call. = FALSE
    )
  }

  method <- match.arg(method)
  speed <- match.arg(speed)

  # unmix.fcs() needs an AF matrix with >= 2 rows for AutoSpectral
  if (method == "AutoSpectral") {
    if (is.null(af_spectra)) {
      stop(
        "'af_spectra' is required for method = 'AutoSpectral'. ",
        "Use sw_unmix_extract_af() to obtain it, or choose ",
        "method = 'OLS' or 'WLS'.",
        call. = FALSE
      )
    }
    if (!is.matrix(af_spectra) || nrow(af_spectra) < 2) {
      stop(
        "'af_spectra' must be a matrix with at least 2 rows for ",
        "method = 'AutoSpectral'. For multi-tissue AF, pass one tissue's matrix.",
        call. = FALSE
      )
    }
  }

  .check_autospectral()

  if (method == "AutoSpectral" && is.null(spectra_variants)) {
    message(
      "Note: Providing 'spectra_variants' via ",
      "sw_unmix_extract_variants() may improve unmixing quality."
    )
  }

  # Warn if AutoSpectralRcpp is not available
  if (
    method %in%
      c("AutoSpectral", "FastPoisson") &&
      !requireNamespace("AutoSpectralRcpp", quietly = TRUE)
  ) {
    message(
      "Note: Install AutoSpectralRcpp for faster unmixing: ",
      "remotes::install_github('DrCytometer/AutoSpectralRcpp')"
    )
  }

  asp <- setup$asp
  dots <- list(...)

  if (is.null(output_dir)) {
    output_dir <- asp$unmixed.fcs.dir
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Snapshot output_dir to report only files this call writes
  list_fcs <- function() {
    list.files(
      output_dir,
      pattern = "\\.fcs$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  files_before <- list_fcs()
  mtime_before <- stats::setNames(file.mtime(files_before), files_before)

  common <- list(
    spectra = spectra,
    asp = asp,
    flow_control = flow_control,
    method = method,
    af_spectra = af_spectra,
    spectra_variants = spectra_variants,
    output_dir = output_dir,
    file_suffix = file_suffix,
    speed = speed,
    parallel = parallel,
    threads = threads,
    chunk_size = chunk_size
  )

  message("\n=== Unmixing with method: ", method, " ===")

  is_dir <- length(input) == 1 && dir.exists(input)
  if (is_dir) {
    fcs_files <- list.files(
      input,
      pattern = "\\.fcs$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(fcs_files) == 0) {
      stop("No FCS files found in directory: ", input, call. = FALSE)
    }
    message("Found ", length(fcs_files), " FCS files in: ", input)

    .sw_call_autospectral(
      "unmix.folder",
      do.call(
        .sw_args_unmix_folder,
        c(list(fcs_dir = input), common, list(dots = dots))
      )
    )
  } else {
    missing_files <- input[!file.exists(input)]
    if (length(missing_files) > 0) {
      stop(
        "FCS file(s) not found: ",
        paste(missing_files, collapse = ", "),
        call. = FALSE
      )
    }

    for (fcs_file in input) {
      message("Unmixing: ", basename(fcs_file))
      .sw_call_autospectral(
        "unmix.fcs",
        do.call(
          .sw_args_unmix_fcs,
          c(list(fcs_file = fcs_file), common, list(dots = dots))
        )
      )
    }
  }

  # New files, or files whose modification time changed
  files_after <- list_fcs()
  changed <- is.na(mtime_before[files_after]) |
    file.mtime(files_after) > mtime_before[files_after]
  output_files <- files_after[changed]

  message(
    "\nUnmixing complete. ",
    length(output_files),
    " file(s) written to: ",
    output_dir
  )

  list(
    output_dir = output_dir,
    output_files = output_files,
    spectra = spectra,
    method = method
  )
}

# ============================================================================
# Function 6: sw_unmix_pipeline (convenience orchestrator)
# ============================================================================

#' End-to-End Unmixing Pipeline
#'
#' One-call orchestrator that runs the AutoSpectral unmixing workflow:
#' setup, control preparation, AF extraction, spectral variant mapping, and
#' unmixing.
#'
#' @param control_dir Character path to single-stain control FCS files.
#' @param sample_input Character path to fully-stained FCS file(s) or a
#'   directory.
#' @param unstained_fcs Character path to an unstained FCS file, or a named
#'   list for tissue-specific AF (see \code{\link{sw_unmix_extract_af}}).
#'   Required when \code{method = "AutoSpectral"}.
#' @param cytometer Character; cytometer type. Default: \code{"aurora"}.
#' @param control_file Character path to an existing control CSV, or
#'   \code{NULL} to auto-generate. Default: \code{NULL}.
#' @param method Character; unmixing method, one of \code{"AutoSpectral"}
#'   (default), \code{"OLS"}, \code{"WLS"}, \code{"Poisson"} or
#'   \code{"FastPoisson"}.
#' @param speed Character; variant testing speed. Default: \code{"fast"}.
#' @param refine Logical or \code{NULL}; AF refinement in
#'   \code{get.af.spectra()}. \code{NULL} (default) uses AutoSpectral's
#'   default.
#' @param output_dir Character path for all output files. Default:
#'   \code{"SpectraWeaveR_unmix"}.
#' @param parallel Logical or \code{NULL}; parallel processing for AF,
#'   variants and unmixing (and for the legacy control path). \code{NULL}
#'   (default) uses each AutoSpectral function's default.
#' @param threads Integer or \code{NULL}; number of threads. Default:
#'   \code{NULL}.
#' @param workflow Character; \code{"automated"} (default) or
#'   \code{"legacy"}, passed to \code{\link{sw_unmix_setup}} and
#'   \code{\link{sw_unmix_prepare}}.
#' @param gating_system Character; gating algorithm for the legacy
#'   workflow. Default: \code{"landmarks"}.
#' @param som_dim Integer or \code{NULL}; SOM grid dimension for both AF
#'   and variant extraction. \code{NULL} (default) keeps AutoSpectral's
#'   separate defaults (10 for AF, 5 for variants in 1.8.8).
#' @param chunk_size Number or \code{NULL}; events per unmixing chunk.
#'   \code{NULL} (default) uses AutoSpectral's default.
#'
#' @return A named list with all intermediate results:
#'   \describe{
#'     \item{\code{setup}}{The \code{sw_setup} object.}
#'     \item{\code{flow_control}}{The \code{flow.control} structure.}
#'     \item{\code{spectra}}{The fluorophore spectral matrix.}
#'     \item{\code{af_spectra}}{The AF spectra (if extracted).}
#'     \item{\code{spectra_variants}}{The spectral variants (if
#'       extracted).}
#'     \item{\code{unmixed}}{The unmixing result from
#'       \code{\link{sw_unmix_run}}.}
#'   }
#'
#' @details
#' This convenience function calls \code{\link{sw_unmix_setup}},
#' \code{\link{sw_unmix_prepare}}, \code{\link{sw_unmix_extract_af}},
#' \code{\link{sw_unmix_extract_variants}}, and \code{\link{sw_unmix_run}}
#' in sequence. Steps 3--4 run only for \code{method = "AutoSpectral"}.
#'
#' \strong{Multi-tissue experiments}: When \code{unstained_fcs} is a named
#' list, the pipeline extracts AF spectra for each tissue but unmixes with
#' the first tissue's AF spectra. The pipeline cannot match sample files to
#' tissues. For tissue-matched unmixing, call \code{\link{sw_unmix_run}}
#' once per tissue group with the matching matrix.
#'
#' @export
sw_unmix_pipeline <- function(
  control_dir,
  sample_input,
  unstained_fcs = NULL,
  cytometer = "aurora",
  control_file = NULL,
  method = c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson"),
  speed = c("fast", "medium", "slow"),
  refine = NULL,
  output_dir = "SpectraWeaveR_unmix",
  parallel = NULL,
  threads = NULL,
  workflow = c("automated", "legacy"),
  gating_system = c("landmarks", "density"),
  som_dim = NULL,
  chunk_size = NULL
) {
  method <- match.arg(method)
  speed <- match.arg(speed)
  workflow <- match.arg(workflow)
  gating_system <- match.arg(gating_system)

  # Validate inputs
  if (!is.character(control_dir) || length(control_dir) != 1) {
    stop("'control_dir' must be a single directory path.", call. = FALSE)
  }
  if (!is.character(sample_input) || length(sample_input) == 0) {
    stop(
      "'sample_input' must be a non-empty character vector (file path(s) or directory).",
      call. = FALSE
    )
  }
  .check_flag(refine, "refine")
  .check_flag(parallel, "parallel")
  .check_threads(threads)
  .check_som_dim(som_dim)
  .check_chunk_size(chunk_size)

  # Only AutoSpectral unmixing uses AF spectra and variants
  needs_percell <- method == "AutoSpectral"
  if (needs_percell && is.null(unstained_fcs)) {
    stop(
      "'unstained_fcs' is required for method = '",
      method,
      "'. Provide path(s) to unstained FCS file(s).",
      call. = FALSE
    )
  }

  results <- list()

  # --- Step 1: Setup ---
  message("\n", strrep("=", 60))
  message("Step 1/5: AutoSpectral Setup")
  message(strrep("=", 60))
  setup <- sw_unmix_setup(
    control_dir = control_dir,
    cytometer = cytometer,
    control_file = control_file,
    output_dir = output_dir,
    figures = TRUE,
    workflow = workflow
  )
  results$setup <- setup

  # --- Step 2: Prepare controls ---
  message("\n", strrep("=", 60))
  message("Step 2/5: Prepare Controls")
  message(strrep("=", 60))
  # Legacy-only arguments go in only for the legacy path
  prepare_args <- list(setup = setup, workflow = workflow)
  if (workflow == "legacy") {
    prepare_args <- c(
      prepare_args,
      .sw_drop_null(list(
        gating_system = gating_system,
        parallel = parallel,
        threads = threads
      ))
    )
  }
  controls <- do.call(sw_unmix_prepare, prepare_args)
  results$flow_control <- controls$flow_control
  results$spectra <- controls$spectra

  # --- Steps 3-4: Per-cell optimization (conditional) ---
  af_spectra <- NULL
  spectra_variants <- NULL

  if (needs_percell) {
    message("\n", strrep("=", 60))
    message("Step 3/5: Extract AF Spectra")
    message(strrep("=", 60))
    af_spectra <- sw_unmix_extract_af(
      unstained_fcs = unstained_fcs,
      setup = setup,
      spectra = controls$spectra,
      refine = refine,
      som_dim = som_dim,
      parallel = parallel,
      threads = threads
    )
    results$af_spectra <- af_spectra

    # Variants derive AF in situ from the control table
    message("\n", strrep("=", 60))
    message("Step 4/5: Extract Spectral Variants")
    message(strrep("=", 60))
    spectra_variants <- sw_unmix_extract_variants(
      setup = setup,
      spectra = controls$spectra,
      som_dim = som_dim,
      parallel = parallel,
      threads = threads
    )
    results$spectra_variants <- spectra_variants
  } else {
    message("\n", strrep("=", 60))
    message("Steps 3-4/5: Skipped (not needed for method = '", method, "')")
    message(strrep("=", 60))
    results$af_spectra <- NULL
    results$spectra_variants <- NULL
  }

  # --- Step 5: Unmix ---
  message("\n", strrep("=", 60))
  message("Step 5/5: Unmix Samples")
  message(strrep("=", 60))

  # For multi-tissue AF, unmix with the first tissue's AF
  unmix_af <- if (is.list(af_spectra) && !is.matrix(af_spectra)) {
    message(
      "Note: Multiple tissue AF spectra provided. Using '",
      names(af_spectra)[1],
      "' for unmixing. ",
      "For tissue-specific AF, call sw_unmix_run() separately per tissue."
    )
    af_spectra[[1]]
  } else {
    af_spectra
  }

  unmixed <- sw_unmix_run(
    input = sample_input,
    spectra = controls$spectra,
    setup = setup,
    flow_control = controls$flow_control,
    method = method,
    af_spectra = unmix_af,
    spectra_variants = spectra_variants,
    speed = speed,
    parallel = parallel,
    threads = threads,
    chunk_size = chunk_size
  )
  results$unmixed <- unmixed

  message("\n", strrep("=", 60))
  message("Pipeline complete!")
  message(strrep("=", 60))

  results
}

#' Load Pre-Unmixed FCS Files
#'
#' Loads FCS files that have already been spectrally unmixed (e.g., by
#' SpectroFlo or AutoSpectral) into a \code{flowSet}. Enforces
#' \code{truncate_max_range = FALSE} for Aurora compatibility.
#'
#' @param fcs_dir Character path to the directory containing unmixed FCS files.
#' @param pattern Regular expression pattern to match FCS file names
#'   (default: \code{"\\.fcs$"}).
#' @param ... Additional arguments passed to \code{\link{sw_io_read_fcs}}.
#'
#' @return A \code{flowSet} containing all matched FCS files.
#'
#' @examples
#' \dontrun{
#' fs <- sw_io_load_unmixed("path/to/unmixed_fcs/")
#' flowCore::sampleNames(fs)
#' }
#'
#' @export
sw_io_load_unmixed <- function(fcs_dir, pattern = "\\.fcs$", ...) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!is.character(fcs_dir) || length(fcs_dir) != 1) {
    stop("'fcs_dir' must be a single directory path.", call. = FALSE)
  }

  if (!dir.exists(fcs_dir)) {
    stop("Directory does not exist: ", fcs_dir, call. = FALSE)
  }

  fcs_files <- list.files(
    fcs_dir,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(fcs_files) == 0) {
    stop(
      "No FCS files matching pattern '",
      pattern,
      "' found in: ",
      fcs_dir,
      call. = FALSE
    )
  }

  message("Loading ", length(fcs_files), " unmixed FCS files from: ", fcs_dir)

  fs <- sw_io_read_fcs(fcs_files, ...)

  message(
    "Loaded flowSet with ",
    length(fs),
    " samples, ",
    ncol(fs[[1]]),
    " channels each"
  )

  fs
}

#' Remove Margin Events
#'
#' Wrapper for \code{\link[PeacoQC]{RemoveMargins}} that removes events at
#' the margins of scatter and fluorescence channels. This should be run
#' \strong{before} data transformation (e.g., arcsinh, logicle).
#'
#' @param ff A \code{flowFrame} object.
#' @param channel_specs A named list specifying channel-specific margin
#'   settings, or \code{NULL} to use the ranges stored in the FCS parameters
#'   (PeacoQC defaults) on all signal channels. Each element is a numeric
#'   vector of length 2, \code{c(lower_bound, upper_bound)}, named by channel
#'   name. When given, only the named channels are checked, and the ranges
#'   are passed to \code{PeacoQC::RemoveMargins(channel_specifications =)}.
#'   PeacoQC keeps events with \code{max(min(lower, 0), min(x)) < x <
#'   min(upper, max(x))}. So the upper bound is exact, but a positive lower
#'   bound acts as 0, and the observed minimum is always removed.
#'
#' @return A \code{flowFrame} with margin events removed. The attributes
#'   \code{events_in} and \code{events_removed} hold the event counts
#'   before filtering and the number removed.
#'
#' @details
#' Margin events are cells that saturate detector limits (floor or ceiling).
#' These can skew downstream analysis and should be removed early.
#' The function must be called before any data transformation because
#' RemoveMargins relies on raw signal intensities.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   mat <- matrix(abs(rnorm(500, 50000, 15000)), ncol = 5,
#'                 dimnames = list(NULL, c("FSC-A", "SSC-A", "CD3", "CD4", "CD8")))
#'   ff <- flowCore::flowFrame(mat)
#'   ff_clean <- sw_filter_margins(ff)
#' }
#' }
#'
#' @export
sw_filter_margins <- function(ff, channel_specs = NULL) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop(
      "Package 'PeacoQC' is required for margin removal. ",
      "Install it from Bioconductor.",
      call. = FALSE
    )
  }

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  if (!methods::is(ff, "flowFrame")) {
    stop("'ff' must be a flowFrame object.", call. = FALSE)
  }

  n_before <- nrow(ff)

  if (is.null(channel_specs)) {
    # Use signal channels only (exclude Time, Original_ID, etc.)
    signal_idx <- which(sw_channel_is_signal(ff))
    if (length(signal_idx) == 0) {
      signal_idx <- seq_len(ncol(ff))
    }
    ff_clean <- PeacoQC::RemoveMargins(ff, signal_idx)
  } else {
    if (!is.list(channel_specs) || is.null(names(channel_specs)) ||
        !all(vapply(channel_specs, function(x) is.numeric(x) && length(x) == 2,
                    logical(1)))) {
      stop("'channel_specs' must be a named list of numeric vectors of ",
           "length 2 (lower, upper).", call. = FALSE)
    }
    # Pass the ranges, not only the names, so user limits take effect.
    ff_clean <- PeacoQC::RemoveMargins(
      ff,
      channels = names(channel_specs),
      channel_specifications = lapply(channel_specs, as.numeric)
    )
  }

  n_after <- nrow(ff_clean)
  n_removed <- n_before - n_after
  pct <- round(100 * n_removed / n_before, 1)
  message(
    "RemoveMargins: removed ",
    n_removed,
    " / ",
    n_before,
    " events (",
    pct,
    "%)"
  )

  attr(ff_clean, "events_in") <- as.integer(n_before)
  attr(ff_clean, "events_removed") <- as.integer(n_removed)
  ff_clean
}
