#' @title Spectral Unmixing Wrappers
#'
#' @description
#' Functions for spectral unmixing of flow cytometry data using AutoSpectral,
#' and for loading pre-unmixed FCS files (e.g., from SpectroFlo).
#'
#' The AutoSpectral workflow is streamlined into five steps plus an orchestrator:
#' \enumerate{
#'   \item \code{\link{sw_unmix_setup}} — initialize parameters and
#'     control file
#'   \item \code{\link{sw_unmix_prepare}} — gate, clean, and extract
#'     fluorophore spectra
#'   \item \code{\link{sw_unmix_extract_af}} — per-cell autofluorescence
#'     extraction from unstained samples
#'   \item \code{\link{sw_unmix_extract_variants}} — fluorophore emission
#'     variability mapping
#'   \item \code{\link{sw_unmix_run}} — unmix fully-stained sample FCS files
#'   \item \code{\link{sw_unmix_pipeline}} — one-call end-to-end orchestrator
#' }
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

# ---- Helper: supported cytometers ----
.valid_cytometers <- function() {
  c("aurora", "auroraNL", "id7000", "s8", "a8",
    "a5se", "opteon", "mosaic", "xenith")
}

# ---- Helper: override asp output directories ----
.override_asp_dirs <- function(asp, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  asp$figure.gate.dir <- file.path(output_dir, "figure_gate")
  asp$figure.af.dir <- file.path(output_dir, "figure_autofluorescence")
  asp$figure.clean.control.dir <- file.path(output_dir,
                                             "figure_clean_controls")
  asp$figure.spectral.ribbon.dir <- file.path(output_dir,
                                               "figure_spectral_ribbon")
  asp$figure.spectra.dir <- file.path(output_dir, "figure_spectra")
  asp$figure.similarity.heatmap.dir <- file.path(output_dir,
                                                   "figure_similarity_heatmap")
  asp$figure.scatter.dir.base <- file.path(output_dir, "figure_scatter")
  asp$table.spectra.dir <- file.path(output_dir, "table_spectra")
  asp$unmixed.fcs.dir <- file.path(output_dir, "unmixed_fcs")
  asp$variant.dir <- file.path(output_dir, "figure_spectral_variants")

  asp
}

# ============================================================================
# Function 1: sw_unmix_setup
# ============================================================================

#' Initialize AutoSpectral Parameters and Control File
#'
#' Sets up the AutoSpectral parameter list for a given cytometer, and either
#' creates a draft control file CSV from the FCS filenames in
#' \code{control_dir} or validates an existing one.
#'
#' @param control_dir Character path to the directory containing single-stain
#'   control FCS files.
#' @param cytometer Character string specifying the cytometer type. One of
#'   \code{"aurora"} (default), \code{"auroraNL"}, \code{"id7000"},
#'   \code{"s8"}, \code{"a8"}, \code{"a5se"}, \code{"opteon"},
#'   \code{"mosaic"}, or \code{"xenith"}.
#' @param control_file Character path to an existing control file CSV, or
#'   \code{NULL} (default) to auto-generate one from the FCS file names.
#' @param output_dir Character path for output files (figures, tables).
#'   Default: \code{"SpectraWeaveR_unmix"}.
#' @param figures Logical; whether to generate QC plots throughout the
#'   workflow. Default: \code{TRUE}.
#'
#' @return A named list (\code{sw_setup}) with components:
#'   \describe{
#'     \item{\code{asp}}{The AutoSpectral parameter list, with output
#'       directories routed under \code{output_dir}.}
#'     \item{\code{control_file}}{The validated path to the control CSV.}
#'     \item{\code{control_dir}}{The control directory path.}
#'   }
#'
#' @details
#' This is Step 1 of the SpectraWeaveR unmixing workflow. If
#' \code{control_file = NULL}, a draft CSV is created and a message instructs
#' the user to review and edit it before proceeding to
#' \code{\link{sw_unmix_prepare}}.
#'
#' @examples
#' \dontrun{
#' setup <- sw_unmix_setup("path/to/controls", cytometer = "aurora")
#' }
#'
#' @export
sw_unmix_setup <- function(control_dir,
                                  cytometer = "aurora",
                                  control_file = NULL,
                                  output_dir = "SpectraWeaveR_unmix",
                                  figures = TRUE) {
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
    stop("'cytometer' must be one of: ",
         paste(valid, collapse = ", "), call. = FALSE)
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a single directory path.", call. = FALSE)
  }

  if (!is.logical(figures) || length(figures) != 1) {
    stop("'figures' must be TRUE or FALSE.", call. = FALSE)
  }

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
  asp <- AutoSpectral::get.autospectral.param(
    cytometer = cytometer,
    figures = figures
  )

  # Override output directories to route under output_dir
  asp <- .override_asp_dirs(asp, output_dir)

  # Handle control file
  if (is.null(control_file)) {
    message("Creating draft control file from FCS files in: ", control_dir)
    # Temporarily change to output_dir so the CSV is written there;
    # on.exit restores wd even if create.control.file() errors
    old_wd <- setwd(output_dir)
    on.exit(setwd(old_wd), add = TRUE)
    AutoSpectral::create.control.file(control_dir, asp)
    setwd(old_wd)

    # Find the generated file
    csv_files <- list.files(output_dir, pattern = "^fcs_control_file.*\\.csv$",
                            full.names = TRUE)
    if (length(csv_files) == 0) {
      stop("Failed to create control file.", call. = FALSE)
    }
    control_file <- csv_files[length(csv_files)]  # latest

    message("\n*** IMPORTANT ***")
    message("A draft control file has been created at: ", control_file)
    message("Please review and edit it before proceeding to ",
            "sw_unmix_prepare().")
    message("See AutoSpectral documentation for control file details:")
    message("  https://drcytometer.github.io/AutoSpectral/articles/",
            "02_Control_File_example.html")
  } else {
    # Validate and check existing control file
    message("Validating control file: ", control_file)
    AutoSpectral::check.control.file(control_dir, control_file, asp)
    message("Control file validation passed.")
  }

  structure(
    list(
      asp = asp,
      control_file = control_file,
      control_dir = control_dir
    ),
    class = "sw_setup"
  )
}

# ============================================================================
# Function 2: sw_unmix_prepare
# ============================================================================

#' Prepare Controls: Gate, Clean, and Extract Fluorophore Spectra
#'
#' Loads single-stain control FCS files, applies automated gating, cleans
#' the controls (autofluorescence removal, universal negatives,
#' downsampling), and extracts fluorophore spectral signatures via robust
#' linear modeling.
#'
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param gating_system Character; the gating algorithm to use.
#'   \code{"landmarks"} (default) or \code{"density"}.
#' @param gate_list Optional named list of pre-defined gates (created via
#'   \code{AutoSpectral::define.gate.landmarks()} or
#'   \code{AutoSpectral::define.gate.density()}). Default: \code{NULL}
#'   (auto-generate gates).
#' @param clean Logical; whether to run control cleaning. Default:
#'   \code{TRUE}.
#' @param af_remove Logical; whether to remove autofluorescence intrusions
#'   from cell-based controls during cleaning. Default: \code{TRUE}.
#' @param parallel Logical; whether to enable parallel processing. Default:
#'   \code{FALSE}.
#' @param threads Integer or \code{NULL}; number of threads for parallel
#'   processing. Default: \code{NULL} (auto-detect).
#'
#' @return A named list (\code{sw_controls}) with components:
#'   \describe{
#'     \item{\code{flow_control}}{The full \code{flow.control} data
#'       structure from AutoSpectral.}
#'     \item{\code{spectra}}{Fluorophore spectral signature matrix
#'       (fluorophores in rows, detectors in columns), L-infinity
#'       normalized.}
#'     \item{\code{setup}}{Pass-through of the \code{sw_setup} object.}
#'   }
#'
#' @details
#' This is Step 2 of the SpectraWeaveR unmixing workflow. It wraps
#' \code{AutoSpectral::define.flow.control()},
#' \code{AutoSpectral::clean.controls()}, and
#' \code{AutoSpectral::get.fluorophore.spectra()} into a single call.
#'
#' @export
sw_unmix_prepare <- function(setup,
                                gating_system = c("landmarks", "density"),
                                gate_list = NULL,
                                clean = TRUE,
                                af_remove = TRUE,
                                parallel = FALSE,
                                threads = NULL) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop("'setup' must be an sw_setup object from sw_unmix_setup().",
         call. = FALSE)
  }

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
  if (!is.null(threads) && (!is.numeric(threads) || length(threads) != 1 ||
                            threads < 1)) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }

  .check_autospectral()

  asp <- setup$asp
  control_dir <- setup$control_dir
  control_file <- setup$control_file

  # Step 2a: Load controls, apply gating, organize into flow.control
  message("\n=== Loading and gating controls ===")
  flow_control <- AutoSpectral::define.flow.control(
    control.dir = control_dir,
    control.def.file = control_file,
    asp = asp,
    gate = TRUE,
    gating.system = gating_system,
    gate.list = gate_list,
    parallel = parallel,
    verbose = TRUE,
    threads = threads
  )

  # Step 2b: Clean controls (optional)
  if (clean) {
    message("\n=== Cleaning controls ===")
    flow_control <- AutoSpectral::clean.controls(
      flow.control = flow_control,
      asp = asp,
      af.remove = af_remove,
      universal.negative = TRUE,
      downsample = TRUE,
      parallel = parallel,
      verbose = TRUE,
      threads = threads
    )
  }

  # Step 2c: Extract fluorophore spectra
  message("\n=== Extracting fluorophore spectra ===")
  spectra <- AutoSpectral::get.fluorophore.spectra(
    flow.control = flow_control,
    asp = asp,
    use.clean.expr = clean,
    figures = asp$figures
  )

  message("Extracted spectra for ", nrow(spectra), " fluorophores across ",
          ncol(spectra), " detectors")

  structure(
    list(
      flow_control = flow_control,
      spectra = spectra,
      setup = setup
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
#' using SOM clustering. For multi-tissue experiments, provide a named list
#' of unstained files to get tissue-matched AF spectra.
#'
#' @param unstained_fcs Character path to an unstained FCS file, or a named
#'   list of paths for tissue-specific AF (e.g.,
#'   \code{list(spleen = "spleen_unstained.fcs",
#'              lung = "lung_unstained.fcs")}).
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param refine Logical; whether to perform a second-pass refinement
#'   targeting high-error cells. Default: \code{TRUE}. Recommended for
#'   tissue samples with complex autofluorescence.
#' @param som_dim Integer; SOM grid dimension. Default: \code{10} (up to
#'   100 AF clusters).
#' @param parallel Logical; whether to enable parallel processing. Default:
#'   \code{TRUE}.
#' @param threads Integer or \code{NULL}; number of threads. Default:
#'   \code{NULL} (auto-detect).
#'
#' @return If \code{unstained_fcs} is a single path: a matrix of AF spectra
#'   (clusters in rows, detectors in columns). If a named list: a named list
#'   of such matrices.
#'
#' @details
#' This is Step 3 of the SpectraWeaveR unmixing workflow. Per-cell AF
#' extraction allows AutoSpectral to model cell-to-cell variation in
#' autofluorescence, producing better unmixing with less spread.
#'
#' For different tissue types (e.g., spleen, lung, liver), provide separate
#' unstained samples and match the resulting AF spectra to the corresponding
#' stained samples during unmixing.
#'
#' Installation of \code{AutoSpectralRcpp} is strongly recommended for
#' faster processing.
#'
#' @export
sw_unmix_extract_af <- function(unstained_fcs,
                                  setup,
                                  spectra,
                                  refine = TRUE,
                                  som_dim = 10,
                                  parallel = TRUE,
                                  threads = NULL) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop("'setup' must be an sw_setup object from sw_unmix_setup().",
         call. = FALSE)
  }

  if (!is.matrix(spectra)) {
    stop("'spectra' must be a matrix of fluorophore spectral signatures.",
         call. = FALSE)
  }

  if (!is.logical(refine) || length(refine) != 1) {
    stop("'refine' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(som_dim) || length(som_dim) != 1 || som_dim < 2) {
    stop("'som_dim' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("'parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(threads) && (!is.numeric(threads) || length(threads) != 1 ||
                            threads < 1)) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }

  # Validate unstained_fcs type early (before AutoSpectral check)
  if (!(is.character(unstained_fcs) && length(unstained_fcs) == 1) &&
      !is.list(unstained_fcs)) {
    stop("'unstained_fcs' must be a single file path or a named list of ",
         "file paths.", call. = FALSE)
  }

  .check_autospectral()

  # Warn if AutoSpectralRcpp is not available
  if (!requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    message("Note: Install AutoSpectralRcpp for faster AF extraction: ",
            "remotes::install_github('DrCytometer/AutoSpectralRcpp')")
  }

  asp <- setup$asp

  # Handle single path vs named list
  if (is.character(unstained_fcs) && length(unstained_fcs) == 1) {
    # Single unstained file
    if (!file.exists(unstained_fcs)) {
      stop("Unstained FCS file does not exist: ", unstained_fcs,
           call. = FALSE)
    }
    message("Extracting AF spectra from: ", basename(unstained_fcs))
    af_spectra <- AutoSpectral::get.af.spectra(
      unstained.sample = unstained_fcs,
      asp = asp,
      spectra = spectra,
      som.dim = as.integer(som_dim),
      figures = asp$figures,
      refine = refine,
      parallel = parallel,
      threads = threads
    )
    message("Extracted ", nrow(af_spectra), " AF spectral signatures")
    return(af_spectra)

  } else if (is.list(unstained_fcs)) {
    # Named list of unstained files (multi-tissue)
    if (is.null(names(unstained_fcs)) ||
        any(names(unstained_fcs) == "")) {
      stop("'unstained_fcs' list must have non-empty names ",
           "(e.g., tissue types).", call. = FALSE)
    }

    af_list <- list()
    for (tissue in names(unstained_fcs)) {
      fcs_path <- unstained_fcs[[tissue]]
      if (!is.character(fcs_path) || length(fcs_path) != 1) {
        stop("Each element of 'unstained_fcs' must be a single file path. ",
             "Problem with: ", tissue, call. = FALSE)
      }
      if (!file.exists(fcs_path)) {
        stop("Unstained FCS file does not exist for '", tissue, "': ",
             fcs_path, call. = FALSE)
      }

      message("Extracting AF spectra for '", tissue, "' from: ",
              basename(fcs_path))
      af_list[[tissue]] <- AutoSpectral::get.af.spectra(
        unstained.sample = fcs_path,
        asp = asp,
        spectra = spectra,
        som.dim = as.integer(som_dim),
        figures = asp$figures,
        title = paste("AF", tissue),
        refine = refine,
        parallel = parallel,
        threads = threads
      )
      message("  -> ", nrow(af_list[[tissue]]), " AF signatures for '",
              tissue, "'")
    }
    return(af_list)

  } else {
    stop("'unstained_fcs' must be a single file path or a named list of ",
         "file paths.", call. = FALSE)
  }
}

# ============================================================================
# Function 4: sw_unmix_extract_variants
# ============================================================================

#' Extract Spectral Variants for Per-Cell Fluorophore Optimization
#'
#' Maps fluorophore emission variability from single-stain controls using
#' SOM clustering. The resulting spectral variants are used during unmixing
#' to optimize fluorophore signatures on a per-cell basis, reducing spillover
#' spread.
#'
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param af_spectra A matrix of AF spectra matching the control cell type
#'   (from \code{\link{sw_unmix_extract_af}}).
#' @param refine Logical; whether to perform a second pass on high-error
#'   cells. Default: \code{TRUE}.
#' @param som_dim Integer; SOM grid dimension. Default: \code{10} (up to
#'   100 variants per fluorophore).
#' @param n_cells Integer; number of brightest positive cells to use per
#'   fluorophore. Default: \code{2000}.
#' @param parallel Logical; enable parallel processing. Default:
#'   \code{FALSE}.
#' @param threads Integer or \code{NULL}; number of threads. Default:
#'   \code{NULL} (auto-detect).
#'
#' @return A list (class \code{sw_variants}) containing:
#'   \describe{
#'     \item{\code{thresholds}}{Per-channel positivity thresholds.}
#'     \item{\code{variants}}{Named list of variant spectral matrices.}
#'     \item{\code{delta.list}}{Pre-computed spectral deltas.}
#'     \item{\code{delta.norms}}{Pre-computed delta norms.}
#'   }
#'
#' @details
#' This is Step 4 of the SpectraWeaveR unmixing workflow. It is optional
#' but recommended: per-cell fluorophore optimization can significantly
#' reduce spillover spread, especially for tandem dyes.
#'
#' Results are saved as an RDS file in the output directory for reuse.
#'
#' @export
sw_unmix_extract_variants <- function(setup,
                                         spectra,
                                         af_spectra,
                                         refine = TRUE,
                                         som_dim = 10,
                                         n_cells = 2000,
                                         parallel = FALSE,
                                         threads = NULL) {
  # Validate inputs before checking dependency
  if (!inherits(setup, "sw_setup")) {
    stop("'setup' must be an sw_setup object from sw_unmix_setup().",
         call. = FALSE)
  }
  if (!is.matrix(spectra)) {
    stop("'spectra' must be a matrix of fluorophore spectral signatures.",
         call. = FALSE)
  }
  if (!is.matrix(af_spectra) || nrow(af_spectra) < 2) {
    stop("'af_spectra' must be a matrix with at least 2 rows of AF spectra.",
         call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1) {
    stop("'refine' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(som_dim) || length(som_dim) != 1 || som_dim < 2) {
    stop("'som_dim' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(n_cells) || length(n_cells) != 1 || n_cells < 1) {
    stop("'n_cells' must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("'parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(threads) && (!is.numeric(threads) || length(threads) != 1 ||
                            threads < 1)) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }

  .check_autospectral()

  asp <- setup$asp
  control_dir <- setup$control_dir
  control_file <- setup$control_file

  message("\n=== Extracting spectral variants ===")
  variants <- AutoSpectral::get.spectral.variants(
    control.dir = control_dir,
    control.def.file = control_file,
    asp = asp,
    spectra = spectra,
    af.spectra = af_spectra,
    n.cells = as.integer(n_cells),
    som.dim = as.integer(som_dim),
    figures = asp$figures,
    parallel = parallel,
    verbose = TRUE,
    threads = threads,
    refine = refine
  )

  n_fluors <- length(variants$variants)
  message("Extracted spectral variants for ", n_fluors, " fluorophores")

  structure(variants, class = "sw_variants")
}

# ============================================================================
# Function 5: sw_unmix_run (rewritten)
# ============================================================================

#' Unmix Spectral Flow Cytometry Data
#'
#' Performs spectral unmixing on fully-stained sample FCS files using
#' AutoSpectral. Supports OLS, WLS, Poisson, and full AutoSpectral unmixing
#' with per-cell autofluorescence extraction and fluorophore optimization.
#'
#' @param input Character path to a single FCS file, a character vector of
#'   FCS file paths, or a directory containing FCS files.
#' @param spectra A fluorophore spectral matrix (from
#'   \code{\link{sw_unmix_prepare}}).
#' @param setup An \code{sw_setup} object from
#'   \code{\link{sw_unmix_setup}}.
#' @param flow_control The \code{flow.control} list from
#'   \code{\link{sw_unmix_prepare}}.
#' @param method Character; unmixing algorithm. One of
#'   \code{"AutoSpectral"} (default), \code{"OLS"}, \code{"WLS"},
#'   \code{"Poisson"}, or \code{"Automatic"}.
#' @param af_spectra A matrix of AF spectra (from
#'   \code{\link{sw_unmix_extract_af}}). Required for
#'   \code{method = "AutoSpectral"}.
#' @param spectra_variants A spectral variants list (from
#'   \code{\link{sw_unmix_extract_variants}}). Optional; enables per-cell
#'   fluorophore optimization.
#' @param speed Character; precision-speed trade-off for per-cell fluorophore
#'   optimization. \code{"fast"} (1 variant/cell, default), \code{"medium"}
#'   (3), or \code{"slow"} (10).
#' @param output_dir Character path for output unmixed FCS files. Default:
#'   \code{NULL} (uses setup output directory).
#' @param file_suffix Character; optional suffix appended to output file
#'   names. Default: \code{NULL}.
#' @param parallel Logical; enable parallel processing. Default:
#'   \code{TRUE}.
#' @param threads Integer or \code{NULL}; number of threads. Default:
#'   \code{NULL} (auto-detect).
#' @param chunk_size Integer; number of events per processing chunk for
#'   memory management. Default: \code{2e6}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{output_dir}}{Path to the directory containing unmixed
#'       FCS files.}
#'     \item{\code{output_files}}{Character vector of paths to the unmixed
#'       FCS files written to disk.}
#'     \item{\code{spectra}}{The spectra matrix used for unmixing.}
#'     \item{\code{method}}{The unmixing method used.}
#'   }
#'
#' @details
#' This is Step 5 of the SpectraWeaveR unmixing workflow. The function
#' auto-detects whether \code{input} is a single file, multiple files, or
#' a directory.
#'
#' When \code{method = "AutoSpectral"}, \code{af_spectra} is required.
#' Providing \code{spectra_variants} additionally enables per-cell
#' fluorophore optimization, which can reduce spillover spread.
#'
#' Installation of \code{AutoSpectralRcpp} is strongly recommended for
#' faster processing.
#'
#' @export
sw_unmix_run <- function(input,
                     spectra,
                     setup,
                     flow_control,
                     method = c("AutoSpectral", "OLS", "WLS",
                                "Poisson", "Automatic"),
                     af_spectra = NULL,
                     spectra_variants = NULL,
                     speed = c("fast", "medium", "slow"),
                     output_dir = NULL,
                     file_suffix = NULL,
                     parallel = TRUE,
                     threads = NULL,
                     chunk_size = 2e6) {
  # Validate inputs before checking dependency
  if (!is.character(input) || length(input) == 0) {
    stop("'input' must be a non-empty character vector (file path(s) ",
         "or directory).", call. = FALSE)
  }
  if (!is.matrix(spectra)) {
    stop("'spectra' must be a matrix of fluorophore spectral signatures.",
         call. = FALSE)
  }
  if (!inherits(setup, "sw_setup")) {
    stop("'setup' must be an sw_setup object from sw_unmix_setup().",
         call. = FALSE)
  }
  if (!is.list(flow_control)) {
    stop("'flow_control' must be the flow.control list from ",
         "sw_unmix_prepare().", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("'parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(threads) && (!is.numeric(threads) || length(threads) != 1 ||
                            threads < 1)) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 1) {
    stop("'chunk_size' must be a positive number.", call. = FALSE)
  }
  if (!is.null(file_suffix) &&
      (!is.character(file_suffix) || length(file_suffix) != 1)) {
    stop("'file_suffix' must be NULL or a single character string.",
         call. = FALSE)
  }

  method <- match.arg(method)
  speed <- match.arg(speed)

  # Validate method-specific requirements
  if (method == "AutoSpectral" && is.null(af_spectra)) {
    stop("'af_spectra' is required for method = 'AutoSpectral'. ",
         "Use sw_unmix_extract_af() to obtain it, or choose ",
         "method = 'OLS' or 'WLS'.", call. = FALSE)
  }

  .check_autospectral()

  if (method == "AutoSpectral" && is.null(spectra_variants)) {
    message("Note: Providing 'spectra_variants' via ",
            "sw_unmix_extract_variants() may improve unmixing quality.")
  }

  # Warn if AutoSpectralRcpp is not available
  if (method == "AutoSpectral" &&
      !requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    message("Note: Install AutoSpectralRcpp for faster unmixing: ",
            "remotes::install_github('DrCytometer/AutoSpectralRcpp')")
  }

  asp <- setup$asp

  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- asp$unmixed.fcs.dir
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Determine input type: directory vs file(s)
  is_dir <- length(input) == 1 && dir.exists(input)

  message("\n=== Unmixing with method: ", method, " ===")

  if (is_dir) {
    # Unmix all files in directory
    fcs_files <- list.files(input, pattern = "\\.fcs$",
                            full.names = TRUE, ignore.case = TRUE)
    if (length(fcs_files) == 0) {
      stop("No FCS files found in directory: ", input, call. = FALSE)
    }
    message("Found ", length(fcs_files), " FCS files in: ", input)

    AutoSpectral::unmix.folder(
      fcs.dir = input,
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
      threads = threads,
      verbose = TRUE,
      chunk.size = as.integer(chunk_size)
    )

  } else {
    # Unmix individual file(s)
    # Validate files exist
    missing <- input[!file.exists(input)]
    if (length(missing) > 0) {
      stop("FCS file(s) not found: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }

    for (fcs_file in input) {
      message("Unmixing: ", basename(fcs_file))
      AutoSpectral::unmix.fcs(
        fcs.file = fcs_file,
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
        threads = threads,
        verbose = TRUE,
        chunk.size = as.integer(chunk_size)
      )
    }
  }

  # Collect output files
  output_files <- list.files(output_dir, pattern = "\\.fcs$",
                              full.names = TRUE, ignore.case = TRUE)

  message("\nUnmixing complete. ", length(output_files),
          " file(s) written to: ", output_dir)

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
#' One-call orchestrator that runs the complete AutoSpectral unmixing
#' workflow: setup, control preparation, AF extraction, spectral variant
#' mapping, and unmixing.
#'
#' @param control_dir Character path to single-stain control FCS files.
#' @param sample_input Character path to fully-stained FCS file(s) or a
#'   directory.
#' @param unstained_fcs Character path to unstained FCS file(s), or a named
#'   list for tissue-specific AF (see \code{\link{sw_unmix_extract_af}}).
#'   Required when \code{method = "AutoSpectral"}.
#' @param cytometer Character; cytometer type. Default: \code{"aurora"}.
#' @param control_file Character path to an existing control CSV, or
#'   \code{NULL} to auto-generate. Default: \code{NULL}.
#' @param method Character; unmixing method. Default:
#'   \code{"AutoSpectral"}.
#' @param speed Character; variant testing speed. Default: \code{"fast"}.
#' @param refine Logical; whether to refine AF and variant extraction.
#'   Default: \code{TRUE}.
#' @param output_dir Character path for all output files. Default:
#'   \code{"SpectraWeaveR_unmix"}.
#' @param parallel Logical; enable parallel processing. Default:
#'   \code{FALSE}.
#' @param threads Integer or \code{NULL}; number of threads. Default:
#'   \code{NULL}.
#' @param gating_system Character; gating algorithm. Default:
#'   \code{"landmarks"}.
#' @param som_dim Integer; SOM grid dimension for AF and variant extraction.
#'   Default: \code{10}.
#' @param chunk_size Integer; events per unmixing chunk. Default:
#'   \code{2e6}.
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
#' in sequence. Steps 3--4 are skipped when \code{method} is \code{"OLS"}
#' or \code{"WLS"} (no per-cell optimization needed).
#'
#' \strong{Multi-tissue experiments}: When \code{unstained_fcs} is a named
#' list of tissue-specific unstained files, the pipeline extracts AF spectra
#' for each tissue but uses only the first tissue's AF spectra for spectral
#' variant extraction and for the final unmixing call. This is because
#' spectral variants are derived from single-stain controls (which come from
#' one tissue type), and the pipeline cannot automatically match sample
#' files to tissue types. For tissue-specific unmixing with matched AF
#' spectra, call \code{\link{sw_unmix_run}} separately for each tissue group,
#' passing the corresponding AF spectra from the named list returned by
#' \code{\link{sw_unmix_extract_af}}.
#'
#' @export
sw_unmix_pipeline <- function(control_dir,
                              sample_input,
                              unstained_fcs = NULL,
                              cytometer = "aurora",
                              control_file = NULL,
                              method = c("AutoSpectral", "OLS", "WLS",
                                         "Poisson", "Automatic"),
                              speed = c("fast", "medium", "slow"),
                              refine = TRUE,
                              output_dir = "SpectraWeaveR_unmix",
                              parallel = FALSE,
                              threads = NULL,
                              gating_system = c("landmarks", "density"),
                              som_dim = 10,
                              chunk_size = 2e6) {
  method <- match.arg(method)
  speed <- match.arg(speed)
  gating_system <- match.arg(gating_system)

  # Validate inputs
  if (!is.character(control_dir) || length(control_dir) != 1) {
    stop("'control_dir' must be a single directory path.", call. = FALSE)
  }
  if (!is.character(sample_input) || length(sample_input) == 0) {
    stop("'sample_input' must be a non-empty character vector (file path(s) ",
         "or directory).", call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1) {
    stop("'refine' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("'parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(threads) && (!is.numeric(threads) || length(threads) != 1 ||
                            threads < 1)) {
    stop("'threads' must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.numeric(som_dim) || length(som_dim) != 1 || som_dim < 2) {
    stop("'som_dim' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 1) {
    stop("'chunk_size' must be a positive number.", call. = FALSE)
  }

  # Check that unstained is provided for AutoSpectral
  needs_percell <- method %in% c("AutoSpectral", "Automatic")
  if (needs_percell && is.null(unstained_fcs)) {
    stop("'unstained_fcs' is required for method = '", method, "'. ",
         "Provide path(s) to unstained FCS file(s).", call. = FALSE)
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
    figures = TRUE
  )
  results$setup <- setup

  # --- Step 2: Prepare controls ---
  message("\n", strrep("=", 60))
  message("Step 2/5: Prepare Controls")
  message(strrep("=", 60))
  controls <- sw_unmix_prepare(
    setup = setup,
    gating_system = gating_system,
    parallel = parallel,
    threads = threads
  )
  results$flow_control <- controls$flow_control
  results$spectra <- controls$spectra

  # --- Steps 3-4: Per-cell optimization (conditional) ---
  af_spectra <- NULL
  spectra_variants <- NULL

  if (needs_percell) {
    # Step 3: AF spectra extraction
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

    # Step 4: Spectral variants
    message("\n", strrep("=", 60))
    message("Step 4/5: Extract Spectral Variants")
    message(strrep("=", 60))

    # For multi-tissue AF, use the first tissue's AF for variant
    # extraction. Variants come from single-stain controls which share
    # the same tissue source, so the first entry is representative.
    af_for_variants <- if (is.list(af_spectra) && !is.matrix(af_spectra)) {
      af_spectra[[1]]
    } else {
      af_spectra
    }

    spectra_variants <- sw_unmix_extract_variants(
      setup = setup,
      spectra = controls$spectra,
      af_spectra = af_for_variants,
      refine = refine,
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

  # For multi-tissue AF with AutoSpectral method, use first tissue's AF.
  # For tissue-specific unmixing, call sw_unmix_run() separately per tissue.
  unmix_af <- if (is.list(af_spectra) && !is.matrix(af_spectra)) {
    message("Note: Multiple tissue AF spectra provided. Using '",
            names(af_spectra)[1], "' for unmixing. ",
            "For tissue-specific AF, call sw_unmix_run() separately per tissue.")
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

  fcs_files <- list.files(fcs_dir, pattern = pattern,
                          full.names = TRUE, ignore.case = TRUE)

  if (length(fcs_files) == 0) {
    stop("No FCS files matching pattern '", pattern, "' found in: ", fcs_dir,
         call. = FALSE)
  }

  message("Loading ", length(fcs_files), " unmixed FCS files from: ", fcs_dir)

  fs <- sw_io_read_fcs(fcs_files, ...)

  message("Loaded flowSet with ", length(fs), " samples, ",
          ncol(fs[[1]]), " channels each")

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
#'   settings, or \code{NULL} to use PeacoQC defaults. Each element should
#'   be a numeric vector of length 2: \code{c(lower_bound, upper_bound)}.
#'
#' @return A \code{flowFrame} with margin events removed.
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
    stop("Package 'PeacoQC' is required for margin removal. ",
         "Install it from Bioconductor.", call. = FALSE)
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
    # channel_specs should list channel indices or names
    ff_clean <- PeacoQC::RemoveMargins(ff, names(channel_specs))
  }

  n_after <- nrow(ff_clean)
  n_removed <- n_before - n_after
  pct <- round(100 * n_removed / n_before, 1)
  message("RemoveMargins: removed ", n_removed, " / ", n_before,
          " events (", pct, "%)")

  ff_clean
}
