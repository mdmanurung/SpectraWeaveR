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
#' @param debris_min Numeric; required, no default. Lower FSC-A cutoff for
#'   the \code{nonDebris} gate, written as the \code{min} argument of the
#'   openCyto \code{mindensity} method. Events below it are never taken as
#'   the valley, so it sets the debris boundary. The right value depends on
#'   the instrument and the FSC-A scale.
#' @param singlet_wider_gate Logical; required, no default. Written as
#'   \code{wider_gate} of the openCyto \code{singletGate} method. \code{TRUE}
#'   widens the singlet gate and keeps more events near the edge.
#' @param flowclust_k Integer; required for templates that use flowClust
#'   (\code{"lymphocyte"}, \code{"treg"}), ignored otherwise. Number of
#'   clusters \code{K} fitted by \code{flowClust}.
#' @param flowclust_target Integer; required for templates that use
#'   flowClust (\code{"lymphocyte"}, \code{"treg"}), ignored otherwise.
#'   Index of the cluster kept as the population (\code{target} in
#'   \code{flowClust}). Clusters are not ordered by size, so this is a
#'   biological choice.
#' @param pool_samples Logical; required, no default. Written to the
#'   \code{collapseDataForGating} column of every row, with \code{groupBy}
#'   left empty. \code{TRUE}: openCyto pools the events of all samples in
#'   the GatingSet, fits one gate, and applies that gate to every sample.
#'   \code{FALSE}: each sample gets its own gate fitted on its own events.
#'
#' @return The file path (invisibly).
#'
#' @details
#' Templates use marker names in the \code{dims} column (e.g., \code{"CD3"},
#' \code{"CD4"}) except for scatter channels (\code{FSC-A}, \code{SSC-A}).
#' Ensure marker names are set on your \code{flowFrame}/\code{flowSet}
#' via \code{\link{sw_channel_set_markers}} before gating.
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
#' before applying it with \code{\link{sw_gate_run}}.
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' sw_gate_template(tmp, template_type = "lymphocyte", debris_min = 50000,
#'                  singlet_wider_gate = TRUE, flowclust_k = 2,
#'                  flowclust_target = 1, pool_samples = TRUE)
#' head(read.csv(tmp))
#' unlink(tmp)
#'
#' @export
sw_gate_template <- function(output_file, template_type = "lymphocyte",
                             debris_min, singlet_wider_gate,
                             flowclust_k, flowclust_target, pool_samples) {
  if (!is.character(output_file) || length(output_file) != 1) {
    stop("'output_file' must be a single file path.", call. = FALSE)
  }
  if (missing(debris_min) || missing(singlet_wider_gate) ||
      missing(pool_samples)) {
    stop("'debris_min', 'singlet_wider_gate' and 'pool_samples' are ",
         "required (no defaults).", call. = FALSE)
  }
  if (!is.numeric(debris_min) || length(debris_min) != 1 ||
      !is.finite(debris_min)) {
    stop("'debris_min' must be a single finite number.", call. = FALSE)
  }
  if (!is.logical(singlet_wider_gate) || length(singlet_wider_gate) != 1 ||
      is.na(singlet_wider_gate)) {
    stop("'singlet_wider_gate' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(pool_samples) || length(pool_samples) != 1 ||
      is.na(pool_samples)) {
    stop("'pool_samples' must be TRUE or FALSE.", call. = FALSE)
  }

  valid_types <- c("lymphocyte", "myeloid", "nk", "treg", "full_pbmc")
  if (!template_type %in% valid_types) {
    stop("Unknown template_type '", template_type, "'. Must be one of: ",
         paste(valid_types, collapse = ", "), call. = FALSE)
  }

  if (template_type %in% c("lymphocyte", "treg")) {
    if (missing(flowclust_k) || missing(flowclust_target)) {
      stop("'flowclust_k' and 'flowclust_target' are required for the '",
           template_type, "' template (no defaults).", call. = FALSE)
    }
    for (nm in c("flowclust_k", "flowclust_target")) {
      v <- get(nm)
      if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 1 ||
          v != round(v)) {
        stop("'", nm, "' must be a single positive integer.", call. = FALSE)
      }
    }
  } else {
    flowclust_k <- NA_integer_
    flowclust_target <- NA_integer_
  }

  template <- .build_template(
    template_type,
    debris_min = debris_min,
    singlet_wider_gate = singlet_wider_gate,
    flowclust_k = flowclust_k,
    flowclust_target = flowclust_target,
    pool_samples = pool_samples
  )

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

.make_gate_row_raw <- function(alias, pop, parent, dims, gating_method,
                           gating_args = NA_character_,
                           collapseDataForGating,
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

.build_template <- function(type, debris_min, singlet_wider_gate,
                            flowclust_k, flowclust_target, pool_samples) {
  # Every row takes the pooling choice; bind it to .make_gate_row.
  .make_gate_row <- function(...) {
    .make_gate_row_raw(..., collapseDataForGating = pool_samples)
  }
  fc_args <- paste0("K = ", flowclust_k, ", target = ", flowclust_target)
  # Common initial gates (debris + singlets)
  debris <- .make_gate_row("nonDebris", "+", "root", "FSC-A",
                           "mindensity", paste0("min = ", debris_min))
  singlets <- .make_gate_row("singlets", "+", "nonDebris", "FSC-A,FSC-H",
                             "singletGate",
                             paste0("wider_gate = ",
                                    toupper(as.character(singlet_wider_gate))))

  if (type == "lymphocyte") {
    lymph <- .make_gate_row("lymphocytes", "+", "singlets", "FSC-A,SSC-A",
                            "flowClust", fc_args)
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
                           "flowClust", fc_args)
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
#'   (see \code{\link{sw_gate_template}}).
#' @param transform_channels Character vector of channel names to transform.
#'   If \code{NULL}, no transformation is applied.
#' @param transform_func Transformation function to apply (default:
#'   \code{flowCore::arcsinhTransform} with the specified cofactor).
#' @param cofactor Numeric; no default. Cofactor of the arcsinh
#'   transformation (events are divided by this value before \code{asinh}).
#'   Required whenever \code{transform_channels} is given and
#'   \code{transform_func} is \code{NULL}; unused otherwise. For spectral
#'   flow data a value around 6000 is common, but the choice changes every
#'   downstream gate, so pass it explicitly.
#'
#' @details
#' Per-sample or pooled gating is set in the template, not here. The
#' \code{collapseDataForGating} and \code{groupBy} columns of the CSV are
#' read by \code{openCyto::gatingTemplate()}; see the \code{pool_samples}
#' argument of \code{\link{sw_gate_template}}.
#'
#' @return A \code{GatingSet} object with gates applied.
#'
#' @examples
#' \dontrun{
#' # Build template and apply
#' tmpl <- tempfile(fileext = ".csv")
#' sw_gate_template(tmpl, "lymphocyte", debris_min = 50000,
#'                  singlet_wider_gate = TRUE, flowclust_k = 2,
#'                  flowclust_target = 1, pool_samples = TRUE)
#' gs <- sw_gate_run(c("sample1.fcs", "sample2.fcs"), gating_template = tmpl,
#'                    transform_channels = c("FSC-A", "FSC-H", "SSC-A"),
#'                    cofactor = 6000)
#' }
#'
#' @export
sw_gate_run <- function(fcs_files, gating_template, transform_channels = NULL,
                    transform_func = NULL, cofactor) {
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

  if (!is.null(transform_channels) && is.null(transform_func) &&
      missing(cofactor)) {
    stop("'cofactor' is required (no default) when 'transform_channels' ",
         "is given without 'transform_func'.", call. = FALSE)
  }

  # openCyto calls flowCore functions such as rectangleGate() and exprs() unqualified,
  # so flowCore has to be attached while gating.
  if (!"package:flowCore" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("flowCore"))
  }

  # Load FCS data
  message("Loading FCS files...")
  fs <- sw_io_read_fcs(fcs_files)
  # A single file comes back as a flowFrame; GatingSet() needs a flowSet.
  if (methods::is(fs, "flowFrame")) {
    fs <- flowCore::flowSet(fs)
  }

  # Convert to GatingSet
  message("Creating GatingSet...")
  gs <- flowWorkspace::GatingSet(fs)

  # Apply transformation if specified
  if (!is.null(transform_channels)) {
    message("Applying transformation to ", length(transform_channels),
            " channels...")

    # flowWorkspace::transform() on a GatingSet needs a transformerList (it
    # rejects flowCore::transformList), so the default arcsinh is a
    # flow_trans() with its inverse; the inverse lets gates be drawn back on
    # the raw scale.
    if (is.null(transform_func)) {
      # Default: arcsinh for spectral flow, asinh(x / cofactor)
      transform_func <- flowWorkspace::flow_trans(
        "arcsinh",
        function(x) asinh(x / cofactor),
        function(x) sinh(x) * cofactor
      )
    }
    if (!inherits(transform_func, "transform")) {
      stop("'transform_func' must be a transformer, e.g. ",
           "flowWorkspace::flow_trans().", call. = FALSE)
    }
    trans_list <- flowWorkspace::transformerList(transform_channels,
                                                 transform_func)

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
#' @examples
#' \dontrun{
#' ff_list <- sw_gate_extract(gs, node = "/singlets/lymphocytes")
#' length(ff_list)
#' }
#'
#' @export
sw_gate_extract <- function(gs, node) {
  if (!requireNamespace("flowWorkspace", quietly = TRUE)) {
    stop("Package 'flowWorkspace' is required.", call. = FALSE)
  }

  if (!is.character(node) || length(node) != 1) {
    stop("'node' must be a single character string.", call. = FALSE)
  }

  if (!methods::is(gs, "GatingSet")) {
    stop("'gs' must be a GatingSet object.", call. = FALSE)
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
