#' @title Automated Gating with openCyto
#'
#' @description
#' Functions for building gating templates and applying automated gating
#' using the openCyto framework.
#'
#' @name gate
#' @keywords internal
NULL

#' Build a Gating Template CSV
#'
#' Writes a default CSV gating template file suitable for a standard
#' lymphocyte panel. The template follows the openCyto CSV format and can be
#' edited before use.
#'
#' @param output_file Character path where the CSV template will be written.
#' @param template_type Character string specifying the template type.
#'   Currently supported: \code{"lymphocyte"} (default).
#'
#' @return The file path (invisibly).
#'
#' @details
#' The generated template includes gating strategies for:
#' \itemize{
#'   \item Singlets (FSC-A vs FSC-H)
#'   \item Live cells (viability dye negative)
#'   \item Lymphocytes (FSC-A vs SSC-A)
#' }
#'
#' Users should review and customize the template for their specific panel
#' before applying it with \code{\link{sw_gate}}.
#'
#' @export
sw_build_gating_template <- function(output_file, template_type = "lymphocyte") {
  if (!is.character(output_file) || length(output_file) != 1) {
    stop("'output_file' must be a single file path.", call. = FALSE)
  }

  valid_types <- c("lymphocyte")
  if (!template_type %in% valid_types) {
    stop("Unknown template_type '", template_type, "'. Must be one of: ",
         paste(valid_types, collapse = ", "), call. = FALSE)
  }

  if (template_type == "lymphocyte") {
    template <- data.frame(
      alias = c("nonDebris", "singlets", "lymphocytes"),
      pop = c("+", "+", "+"),
      parent = c("root", "nonDebris", "singlets"),
      dims = c("FSC-A", "FSC-A,FSC-H", "FSC-A,SSC-A"),
      gating_method = c("mindensity", "singletGate",
                         "flowClust"),
      gating_args = c("min = 50000",
                       "wider_gate = TRUE",
                       "K = 2, target = 1"),
      collapseDataForGating = c(TRUE, TRUE, TRUE),
      groupBy = c(NA, NA, NA),
      preprocessing_method = c(NA, NA, NA),
      preprocessing_args = c(NA, NA, NA),
      stringsAsFactors = FALSE
    )
  }

  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir) && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE)
  }

  utils::write.csv(template, output_file, row.names = FALSE)
  message("Gating template written to: ", output_file)

  invisible(output_file)
}

#' Apply Automated Gating
#'
#' Full openCyto gating workflow: loads FCS files into a GatingSet, applies
#' optional transformation, and runs hierarchical gating using a CSV template.
#'
#' @param fcs_files Character vector of FCS file paths.
#' @param gating_template Character path to a CSV gating template file
#'   (see \code{\link{sw_build_gating_template}}).
#' @param transform_channels Character vector of channel names to transform.
#'   If \code{NULL}, no transformation is applied.
#' @param transform_func Transformation function to apply (default:
#'   \code{flowCore::arcsinhTransform} with the specified cofactor).
#' @param cofactor Numeric; cofactor for the default arcsinh transformation
#'   (default: 6000, appropriate for spectral flow cytometry). Ignored if
#'   \code{transform_func} is provided.
#' @param ... Additional arguments.
#'
#' @return A \code{GatingSet} object with gates applied.
#'
#' @export
sw_gate <- function(fcs_files, gating_template, transform_channels = NULL,
                    transform_func = NULL, cofactor = 6000, ...) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }
  if (!requireNamespace("flowWorkspace", quietly = TRUE)) {
    stop("Package 'flowWorkspace' is required.", call. = FALSE)
  }
  if (!requireNamespace("openCyto", quietly = TRUE)) {
    stop("Package 'openCyto' is required.", call. = FALSE)
  }

  if (!is.character(fcs_files) || length(fcs_files) == 0) {
    stop("'fcs_files' must be a non-empty character vector.", call. = FALSE)
  }

  if (!file.exists(gating_template)) {
    stop("Gating template file not found: ", gating_template, call. = FALSE)
  }

  # Load FCS data
  message("Loading FCS files...")
  fs <- sw_read_fcs(fcs_files)

  # Convert to GatingSet
  message("Creating GatingSet...")
  gs <- flowWorkspace::GatingSet(fs)

  # Apply transformation if specified
  if (!is.null(transform_channels)) {
    message("Applying transformation to ", length(transform_channels),
            " channels...")

    if (is.null(transform_func)) {
      # Default: arcsinh for spectral flow
      trans_list <- flowCore::transformList(
        transform_channels,
        flowCore::arcsinhTransform(a = 0, b = 1 / cofactor, c = 0)
      )
    } else {
      trans_list <- flowCore::transformList(transform_channels, transform_func)
    }

    gs <- flowWorkspace::transform(gs, trans_list)
  }

  # Load and apply gating template
  message("Loading gating template: ", gating_template)
  gt <- openCyto::gatingTemplate(gating_template)

  message("Applying gates...")
  openCyto::gt_gating(gt, gs)

  message("Gating complete. Nodes: ",
          paste(flowWorkspace::gs_get_pop_paths(gs), collapse = ", "))

  gs
}

#' Extract Gated Populations
#'
#' Extracts \code{flowFrame} objects from a specific gate node of a
#' \code{GatingSet}.
#'
#' @param gs A \code{GatingSet} object.
#' @param node Character string specifying the gate node path
#'   (e.g., \code{"/singlets/lymphocytes"}).
#'
#' @return A named list of \code{flowFrame} objects, one per sample in the
#'   GatingSet.
#'
#' @export
sw_extract_gated <- function(gs, node) {
  if (!requireNamespace("flowWorkspace", quietly = TRUE)) {
    stop("Package 'flowWorkspace' is required.", call. = FALSE)
  }

  if (!methods::is(gs, "GatingSet")) {
    stop("'gs' must be a GatingSet object.", call. = FALSE)
  }

  if (!is.character(node) || length(node) != 1) {
    stop("'node' must be a single character string.", call. = FALSE)
  }

  # Extract data from the specified gate node
  sample_names <- flowWorkspace::sampleNames(gs)
  ff_list <- stats::setNames(
    lapply(sample_names, function(sn) {
      flowWorkspace::gh_pop_get_data(gs[[sn]], node)
    }),
    sample_names
  )

  message("Extracted ", length(ff_list), " samples from node '", node, "'")
  ff_list
}
