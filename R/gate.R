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
#'   Supported: \code{"lymphocyte"} (default), \code{"myeloid"},
#'   \code{"nk"}, \code{"treg"}, \code{"full_pbmc"}.
#'
#' @return The file path (invisibly).
#'
#' @details
#' Templates use marker names in the \code{dims} column (e.g., \code{"CD3"},
#' \code{"CD4"}) except for scatter channels (\code{FSC-A}, \code{SSC-A}).
#' Ensure marker names are set on your \code{flowFrame}/\code{flowSet}
#' via \code{\link{sw_set_marker_names}} before gating.
#'
#' Available templates:
#' \describe{
#'   \item{\code{lymphocyte}}{nonDebris -> singlets -> lymphocytes
#'     (FSC/SSC-based, 3 gates)}
#'   \item{\code{myeloid}}{nonDebris -> singlets -> CD3neg -> HLA-DR+
#'     (4 gates, targets monocytes/DCs)}
#'   \item{\code{nk}}{nonDebris -> singlets -> CD3neg -> CD56+
#'     (4 gates, targets NK cells)}
#'   \item{\code{treg}}{nonDebris -> singlets -> CD3+ -> CD4+ ->
#'     CD25hiCD127lo (5 gates, targets regulatory T cells)}
#'   \item{\code{full_pbmc}}{nonDebris -> singlets -> branching into
#'     T cells, B cells, NK cells, monocytes (8 gates)}
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

  valid_types <- c("lymphocyte", "myeloid", "nk", "treg", "full_pbmc")
  if (!template_type %in% valid_types) {
    stop("Unknown template_type '", template_type, "'. Must be one of: ",
         paste(valid_types, collapse = ", "), call. = FALSE)
  }

  template <- .build_template(template_type)

  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir) && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE)
  }

  utils::write.csv(template, output_file, row.names = FALSE)
  message("Gating template written to: ", output_file)

  invisible(output_file)
}


# ---------------------------------------------------------------------------
# Internal: template builders
# ---------------------------------------------------------------------------

.make_gate_row <- function(alias, pop, parent, dims, gating_method,
                           gating_args = NA_character_,
                           collapseDataForGating = TRUE,
                           groupBy = NA_character_,
                           preprocessing_method = NA_character_,
                           preprocessing_args = NA_character_) {
  data.frame(
    alias = alias, pop = pop, parent = parent, dims = dims,
    gating_method = gating_method, gating_args = gating_args,
    collapseDataForGating = collapseDataForGating,
    groupBy = groupBy,
    preprocessing_method = preprocessing_method,
    preprocessing_args = preprocessing_args,
    stringsAsFactors = FALSE
  )
}

.build_template <- function(type) {
  # Common initial gates (debris + singlets)
  debris <- .make_gate_row("nonDebris", "+", "root", "FSC-A",
                           "mindensity", "min = 50000")
  singlets <- .make_gate_row("singlets", "+", "nonDebris", "FSC-A,FSC-H",
                             "singletGate", "wider_gate = TRUE")

  if (type == "lymphocyte") {
    lymph <- .make_gate_row("lymphocytes", "+", "singlets", "FSC-A,SSC-A",
                            "flowClust", "K = 2, target = 1")
    return(do.call(rbind, list(debris, singlets, lymph)))
  }

  if (type == "myeloid") {
    cd3neg <- .make_gate_row("CD3neg", "-", "singlets", "CD3",
                             "mindensity")
    hladr <- .make_gate_row("HLADRpos", "+", "CD3neg", "HLA-DR",
                            "mindensity")
    return(do.call(rbind, list(debris, singlets, cd3neg, hladr)))
  }

  if (type == "nk") {
    cd3neg <- .make_gate_row("CD3neg", "-", "singlets", "CD3",
                             "mindensity")
    cd56pos <- .make_gate_row("CD56pos", "+", "CD3neg", "CD56",
                              "mindensity")
    return(do.call(rbind, list(debris, singlets, cd3neg, cd56pos)))
  }

  if (type == "treg") {
    cd3pos <- .make_gate_row("CD3pos", "+", "singlets", "CD3",
                             "mindensity")
    cd4pos <- .make_gate_row("CD4pos", "+", "CD3pos", "CD4",
                             "mindensity")
    treg <- .make_gate_row("Treg", "+", "CD4pos", "CD25,CD127",
                           "flowClust", "K = 2, target = 1")
    return(do.call(rbind, list(debris, singlets, cd3pos, cd4pos, treg)))
  }

  if (type == "full_pbmc") {
    # T cells branch
    cd3pos <- .make_gate_row("CD3pos", "+", "singlets", "CD3",
                             "mindensity")
    cd4pos <- .make_gate_row("CD4pos", "+", "CD3pos", "CD4",
                             "mindensity")
    cd8pos <- .make_gate_row("CD8pos", "+", "CD3pos", "CD8",
                             "mindensity")
    # B cells
    cd3neg <- .make_gate_row("CD3neg", "-", "singlets", "CD3",
                             "mindensity")
    cd19pos <- .make_gate_row("CD19pos", "+", "CD3neg", "CD19",
                              "mindensity")
    # NK cells
    cd56pos <- .make_gate_row("CD56pos", "+", "CD3neg", "CD56",
                              "mindensity")
    # Monocytes
    cd14pos <- .make_gate_row("CD14pos", "+", "CD3neg", "CD14",
                              "mindensity")
    return(do.call(rbind, list(debris, singlets,
                               cd3pos, cd4pos, cd8pos,
                               cd3neg, cd19pos, cd56pos, cd14pos)))
  }
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
