# R/gate_wrap.R
# Wrappers that move sw_gate objects in and out of external gating tools:
# flowCore gates, GatingML 2.0 files, FlowJo workspaces (CytoML, AGPL-3, wrapped
# only), flowGate (MIT, manual gating in Shiny) and flowMagic (Apache-2.0,
# template-based automated gating). Nothing is copied. External packages are
# Suggests; tests replace the lookup helper .sw_gate_wrap_fn().

# Fetch an exported function from a Suggests package, or stop with install help.
.sw_gate_wrap_fn <- function(name, pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    hint <- switch(pkg,
      flowMagic = "remotes::install_github('semontante/flowMagic')",
      paste0("BiocManager::install('", pkg, "')")
    )
    stop("Package '", pkg, "' is required. Install it with ", hint, ".", call. = FALSE)
  }
  getExportedValue(pkg, name)
}

# ---- sw_gate <-> flowCore ---------------------------------------------------

#' Convert an sw_gate to a flowCore gate
#'
#' A polygon gate becomes a `flowCore::polygonGate`. A range gate becomes a
#' one-dimensional `flowCore::rectangleGate`. The parent is not stored in a
#' flowCore filter; keep it in the `sw_gate`.
#'
#' @param gate An `sw_gate`.
#' @return A `polygonGate` or `rectangleGate`, with `filterId` set to the gate id.
#' @seealso [sw_gate_from_flowcore()]
#' @export
sw_gate_to_flowcore <- function(gate) {
  if (!inherits(gate, "sw_gate")) {
    stop("`gate` must be an sw_gate.", call. = FALSE)
  }
  if (gate$type == "range") {
    bounds <- stats::setNames(list(c(gate$vertices["min", 1], gate$vertices["max", 1])), gate$channels)
    return(do.call(flowCore::rectangleGate, c(list(filterId = gate$id), bounds)))
  }
  v <- gate$vertices
  colnames(v) <- gate$channels
  rownames(v) <- NULL
  flowCore::polygonGate(.gate = v, filterId = gate$id)
}

#' Convert a flowCore gate to an sw_gate
#'
#' Accepts a `polygonGate` or a `rectangleGate`. A one-dimensional rectangle
#' becomes a range gate. A two-dimensional rectangle becomes a polygon with four
#' vertices (open sides are not allowed there). Other filters are refused.
#'
#' @param filter A `flowCore::polygonGate` or `flowCore::rectangleGate`.
#' @param id Character, the new gate id. Required; the filter's own id is not used.
#' @param parent Character parent gate id, or `NA`. Required; flowCore filters
#'   do not hold a parent.
#' @return An `sw_gate`.
#' @seealso [sw_gate_to_flowcore()]
#' @export
sw_gate_from_flowcore <- function(filter, id, parent) {
  if (methods::is(filter, "polygonGate")) {
    v <- filter@boundaries
    return(sw_gate_polygon(id, colnames(v), v, parent))
  }
  if (methods::is(filter, "rectangleGate")) {
    chans <- as.character(flowCore::parameters(filter))
    lo <- unname(filter@min)
    hi <- unname(filter@max)
    if (length(chans) == 1) {
      return(sw_gate_range(id, chans, lo, hi, parent))
    }
    if (length(chans) == 2) {
      if (any(!is.finite(c(lo, hi)))) {
        stop("A 2D rectangle with an open side cannot become a polygon.", call. = FALSE)
      }
      v <- cbind(c(lo[1], hi[1], hi[1], lo[1]), c(lo[2], lo[2], hi[2], hi[2]))
      return(sw_gate_polygon(id, chans, v, parent))
    }
    stop("Rectangle gates on more than 2 channels are not supported.", call. = FALSE)
  }
  stop("`filter` must be a flowCore polygonGate or rectangleGate.", call. = FALSE)
}

# ---- GatingML 2.0 (own writer and reader, xml2) -----------------------------

.sw_gml_ns <- c(
  gating = "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
  dt = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
)
.sw_gml_ext <- "urn:spectraweaver:gatingml-ext"

.sw_xml_esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}
.sw_num <- function(x) format(x, digits = 15, scientific = FALSE, trim = TRUE)

.sw_gml_dim <- function(channel, lo = NULL, hi = NULL) {
  attrs <- ""
  if (!is.null(lo) && is.finite(lo)) attrs <- paste0(attrs, " gating:min=\"", .sw_num(lo), "\"")
  if (!is.null(hi) && is.finite(hi)) attrs <- paste0(attrs, " gating:max=\"", .sw_num(hi), "\"")
  paste0(
    "  <gating:dimension gating:compensation-ref=\"uncompensated\"", attrs, ">",
    "<dt:fcs-dimension dt:name=\"", .sw_xml_esc(channel), "\"/></gating:dimension>"
  )
}

#' Write gates as a GatingML 2.0 file
#'
#' Writes `sw_gate` objects as `RectangleGate` (range) and `PolygonGate`
#' (polygon) elements. The parent id goes in an extension attribute
#' `sw:parent`, because GatingML 2.0 has no parent field; other tools ignore it.
#' The writer is the package's own (xml2 is needed only to read). Which
#' installed function could write GatingML instead: none. `flowUtils`, which
#' had `write.gatingML`, is not on Bioconductor 3.20, and CytoML only exports
#' Cytobank and FlowJo writers that need a GatingSet.
#'
#' @param gates A list of `sw_gate` objects, or one `sw_gate`.
#' @param path Output file path.
#' @return `path`, invisibly.
#' @seealso [sw_gate_read_gatingml()]
#' @export
sw_gate_write_gatingml <- function(gates, path) {
  if (inherits(gates, "sw_gate")) gates <- list(gates)
  ok <- is.list(gates) && length(gates) > 0 && all(vapply(gates, inherits, logical(1), "sw_gate"))
  if (!ok) stop("`gates` must be a non-empty list of sw_gate objects.", call. = FALSE)
  ids <- vapply(gates, function(g) g$id, character(1))
  if (anyDuplicated(ids)) stop("Gate ids must be unique.", call. = FALSE)
  one <- function(g) {
    par <- if (is.na(g$parent)) "" else paste0(" sw:parent=\"", .sw_xml_esc(g$parent), "\"")
    if (g$type == "range") {
      return(paste0(
        " <gating:RectangleGate gating:id=\"", .sw_xml_esc(g$id), "\"", par, ">\n",
        .sw_gml_dim(g$channels, g$vertices["min", 1], g$vertices["max", 1]), "\n",
        " </gating:RectangleGate>"
      ))
    }
    dims <- vapply(g$channels, .sw_gml_dim, character(1), USE.NAMES = FALSE)
    verts <- vapply(seq_len(nrow(g$vertices)), function(i) {
      paste0(
        "  <gating:vertex><gating:coordinate dt:value=\"", .sw_num(g$vertices[i, 1]),
        "\"/><gating:coordinate dt:value=\"", .sw_num(g$vertices[i, 2]), "\"/></gating:vertex>"
      )
    }, character(1))
    paste0(
      " <gating:PolygonGate gating:id=\"", .sw_xml_esc(g$id), "\"", par, ">\n",
      paste(dims, collapse = "\n"), "\n", paste(verts, collapse = "\n"), "\n",
      " </gating:PolygonGate>"
    )
  }
  txt <- paste0(
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
    "<gating:Gating-ML xmlns:gating=\"", .sw_gml_ns[["gating"]], "\" xmlns:dt=\"",
    .sw_gml_ns[["dt"]], "\" xmlns:sw=\"", .sw_gml_ext, "\">\n",
    paste(vapply(gates, one, character(1)), collapse = "\n"), "\n</gating:Gating-ML>\n"
  )
  writeLines(txt, path, useBytes = TRUE)
  invisible(path)
}

#' Read gates from a GatingML 2.0 file
#'
#' Reads `RectangleGate` (one dimension, or two dimensions with finite bounds)
#' and `PolygonGate` elements. Other gate types are skipped and named in the
#' `skipped` attribute. The parent comes from the `sw:parent` attribute that
#' [sw_gate_write_gatingml()] writes; files from other tools give `NA`.
#'
#' @param path GatingML file path.
#' @return A named list of `sw_gate` objects, with attribute `skipped` (ids of
#'   elements that could not be read).
#' @seealso [sw_gate_write_gatingml()]
#' @export
sw_gate_read_gatingml <- function(path) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required to read GatingML.", call. = FALSE)
  }
  ns <- c(.sw_gml_ns, sw = .sw_gml_ext)
  doc <- xml2::read_xml(path)
  nodes <- xml2::xml_find_all(doc, "/*/*", ns)
  gates <- list()
  skipped <- character(0)
  for (nd in nodes) {
    kind <- xml2::xml_name(nd)
    gid <- xml2::xml_attr(nd, "gating:id", ns)
    if (is.na(gid)) gid <- xml2::xml_attr(nd, "id")
    par <- xml2::xml_attr(nd, "sw:parent", ns)
    dims <- xml2::xml_find_all(nd, "gating:dimension", ns)
    chans <- vapply(dims, function(d) {
      xml2::xml_attr(xml2::xml_find_first(d, "dt:fcs-dimension", ns), "dt:name", ns)
    }, character(1))
    num <- function(d, a, dflt) {
      v <- xml2::xml_attr(d, paste0("gating:", a), ns)
      if (is.na(v)) dflt else as.numeric(v)
    }
    g <- tryCatch(
      {
        if (kind == "RectangleGate" && length(dims) == 1) {
          sw_gate_range(gid, chans, num(dims[[1]], "min", -Inf), num(dims[[1]], "max", Inf), par)
        } else if (kind == "RectangleGate" && length(dims) == 2) {
          lo <- vapply(dims, num, numeric(1), "min", -Inf)
          hi <- vapply(dims, num, numeric(1), "max", Inf)
          if (any(!is.finite(c(lo, hi)))) stop("open side")
          sw_gate_polygon(gid, chans, cbind(c(lo[1], hi[1], hi[1], lo[1]), c(lo[2], lo[2], hi[2], hi[2])), par)
        } else if (kind == "PolygonGate" && length(dims) == 2) {
          co <- xml2::xml_attr(xml2::xml_find_all(nd, "gating:vertex/gating:coordinate", ns), "dt:value", ns)
          sw_gate_polygon(gid, chans, matrix(as.numeric(co), ncol = 2, byrow = TRUE), par)
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )
    if (is.null(g)) skipped <- c(skipped, gid) else gates[[gid]] <- g
  }
  attr(gates, "skipped") <- skipped
  gates
}

# ---- GatingSet readers -------------------------------------------------------

# Turn the gates of one GatingHierarchy-like sample into a named list of
# sw_gate. Supports polygon and rectangle filters; others are recorded.
.sw_gs_gates <- function(gs, sample, nodes = NULL) {
  gp <- function(name) .sw_gate_wrap_fn(name, "flowWorkspace")
  full <- gp("gs_get_pop_paths")(gs, path = "full")
  auto <- gp("gs_get_pop_paths")(gs, path = "auto")
  keep <- full != "root"
  full <- full[keep]
  auto <- auto[keep]
  to_auto <- stats::setNames(auto, full)
  gates <- list()
  skipped <- character(0)
  for (i in seq_along(full)) {
    if (!is.null(nodes) && !(auto[i] %in% nodes)) next
    flt <- gp("gs_pop_get_gate")(gs, full[i])
    flt <- if (is.list(flt)) flt[[sample]] else flt
    par_full <- dirname(full[i])
    par <- if (par_full %in% c("/", ".", "root")) NA_character_ else unname(to_auto[par_full])
    g <- tryCatch(sw_gate_from_flowcore(flt, auto[i], par), error = function(e) NULL)
    if (is.null(g)) skipped <- c(skipped, auto[i]) else gates[[auto[i]]] <- g
  }
  attr(gates, "skipped") <- skipped
  gates
}

#' Read the gates of a FlowJo workspace
#'
#' Wraps CytoML (AGPL-3, not copied): `open_flowjo_xml()` and
#' `flowjo_to_gatingset()` build a GatingSet without loading events
#' (`execute = FALSE`), then flowWorkspace getters return each sample's gates.
#' Only polygon and rectangle gates are converted; others are named in the
#' `skipped` attribute of each sample's list.
#'
#' Gate coordinates are in the workspace's own scale. FlowJo stores them on the
#' transformed (for example biexponential) scale, so apply them only to data on
#' the same scale.
#'
#' @param wsp_path Path to a FlowJo `.wsp` file.
#' @param group Character or number, the workspace group to import (the
#'   `name` argument of `flowjo_to_gatingset()`). Required: groups can hold
#'   different gates.
#' @param transform `TRUE` or `FALSE`, passed to `flowjo_to_gatingset()`. It
#'   decides whether the workspace transformations are applied to the gates.
#'   Required.
#' @param ... Further arguments for `flowjo_to_gatingset()`.
#'
#' @return A named list (one entry per sample) of named lists of `sw_gate`.
#' @export
sw_gate_read_flowjo <- function(wsp_path, group, transform, ...) {
  if (!file.exists(wsp_path)) stop("`wsp_path` does not exist: ", wsp_path, call. = FALSE)
  force(group)
  force(transform)
  ws <- .sw_gate_wrap_fn("open_flowjo_xml", "CytoML")(wsp_path)
  gs <- .sw_gate_wrap_fn("flowjo_to_gatingset", "CytoML")(
    ws,
    name = group, execute = FALSE, transform = transform, ...
  )
  samples <- .sw_gate_wrap_fn("sampleNames", "flowWorkspace")(gs)
  stats::setNames(lapply(samples, function(s) .sw_gs_gates(gs, s)), samples)
}

#' Draw gates by hand with flowGate and return them as sw_gate objects
#'
#' Wraps `flowGate::gs_gate_interactive()`, which opens a Shiny app on a
#' GatingSet. The gate is added to the GatingSet (in place) under `filter_id`;
#' this function then reads it back per sample. The app needs a browser, so it
#' runs only in an interactive session. Gate coordinates are on the scale the
#' GatingSet holds. For a transformed axis use `flowGate::gs_gate_transform_interactive()`
#' yourself and read the gates with [sw_gate_from_gatingset()].
#'
#' @param gs A `GatingSet`, or a `flowSet`/`cytoset` (wrapped with
#'   `flowWorkspace::GatingSet()`).
#' @param filter_id Character, name of the gate to draw. Required.
#' @param dims List of two channel names (x, y), or one for a range gate.
#'   Required.
#' @param sample Number or name of the sample shown while drawing. Required.
#' @param subset Character, the parent population the gate is drawn in
#'   (`"root"` for all events). Required.
#' @param regate `TRUE` to redraw an existing gate. Required.
#' @param overlay_gates Character vector of existing gates to show, or `NULL`.
#'
#' @return A named list (per sample) of `sw_gate`, and the GatingSet in
#'   attribute `gs`.
#' @export
sw_gate_draw <- function(gs, filter_id, dims, sample, subset, regate, overlay_gates = NULL) {
  force(filter_id)
  force(dims)
  force(sample)
  force(subset)
  force(regate)
  if (!interactive() && !isTRUE(getOption("sw.gate_draw_allow_noninteractive"))) {
    stop("sw_gate_draw() opens a Shiny app and needs an interactive session.", call. = FALSE)
  }
  if (!inherits(gs, "GatingSet")) {
    gs <- .sw_gate_wrap_fn("GatingSet", "flowWorkspace")(gs)
  }
  draw <- .sw_gate_wrap_fn("gs_gate_interactive", "flowGate")
  draw(
    gs,
    filterId = filter_id, sample = sample, dims = as.list(dims), subset = subset,
    regate = regate, overlayGates = overlay_gates
  )
  samples <- .sw_gate_wrap_fn("sampleNames", "flowWorkspace")(gs)
  out <- stats::setNames(lapply(samples, function(s) .sw_gs_gates(gs, s, nodes = filter_id)), samples)
  attr(out, "gs") <- gs
  out
}

#' Read the gates of a GatingSet as sw_gate objects
#'
#' @param gs A `GatingSet`.
#' @param nodes Character vector of population names to keep, or `NULL` for all.
#' @return A named list (per sample) of named lists of `sw_gate`.
#' @export
sw_gate_from_gatingset <- function(gs, nodes) {
  force(nodes)
  samples <- .sw_gate_wrap_fn("sampleNames", "flowWorkspace")(gs)
  stats::setNames(lapply(samples, function(s) .sw_gs_gates(gs, s, nodes = nodes)), samples)
}

# ---- flowMagic ---------------------------------------------------------------

#' Train a flowMagic template model from hand-drawn template gates
#'
#' Wraps `flowMagic::magic_template_train()`. Unless `use_existing_templates`
#' is `TRUE`, flowMagic opens an interactive window to draw each template gate
#' on each sample in `samples_id`. The result is the flowMagic list
#' (`ref_model_info`, `ref_train`, `list_out_final`, `all_list_out_obj`). Save it
#' with `saveRDS()`: pass it back as `input_template_obj` with
#' `use_existing_templates = TRUE` to retrain without redrawing.
#'
#' Every argument below changes the model, so none has a default. No pretrained
#' model is needed for this template workflow. flowMagic's separate
#' "generalized" models (ALL_GP_models.tar.gz from FRDR) are not used here.
#'
#' @param gs A `GatingSet`.
#' @param parent_node Character, population the template gates are drawn in.
#' @param samples_id Character, samples to draw templates on.
#' @param channel_x,channel_y Character, the two channels.
#' @param n_gates Whole number of template gates per sample.
#' @param train_model Character, model family (flowMagic: `"rf"`, `"knn"`,
#'   `"nb"`, `"nnet"`, `"dt"`).
#' @param n_points_per_plot Whole number of events per sample, or `NULL`.
#' @param prop_down Proportion of events kept (0 to 1), or `NULL`. Give exactly
#'   one of `n_points_per_plot` and `prop_down`.
#' @param seed Whole number, seed set before downsampling.
#' @param input_template_obj A saved result of this function, or `NULL`.
#' @param use_existing_templates `TRUE` to reuse `input_template_obj` without
#'   drawing.
#' @param ... Further arguments for `magic_template_train()`.
#'
#' @return The flowMagic training list.
#' @seealso [sw_gate_flowmagic()]
#' @export
sw_gate_flowmagic_train <- function(gs, parent_node, samples_id, channel_x, channel_y,
                                    n_gates, train_model, n_points_per_plot, prop_down,
                                    seed, input_template_obj, use_existing_templates, ...) {
  force(parent_node)
  force(samples_id)
  force(n_gates)
  force(train_model)
  force(seed)
  force(input_template_obj)
  force(use_existing_templates)
  if (is.null(n_points_per_plot) == is.null(prop_down)) {
    stop("Give exactly one of `n_points_per_plot` and `prop_down`; set the other to NULL.", call. = FALSE)
  }
  train <- .sw_gate_wrap_fn("magic_template_train", "flowMagic")
  train(
    gs_input = gs, parent_node = parent_node, samples_id = samples_id,
    channel_x = channel_x, channel_y = channel_y, n_gates = n_gates,
    input_template_obj = input_template_obj,
    use_existing_templates = use_existing_templates, train_model = train_model,
    prop_down = prop_down, n_points_per_plot = n_points_per_plot, seed = seed, ...
  )
}

#' Predict gates on new samples with a flowMagic template model
#'
#' Wraps `flowMagic::magic_template_predict()`. flowMagic labels events, then
#' draws a concave-hull polygon around each requested label for each sample.
#' This function returns those polygons as `sw_gate` objects.
#'
#' flowMagic adds the predicted gates to `gs` in place (by reference). Pass a
#' copy (`flowWorkspace::gs_clone()`) if `gs` must stay unchanged. Event labels
#' for one sample come from [sw_gate_apply()] on the returned polygon.
#'
#' @param gs A `GatingSet` holding the samples to gate.
#' @param parent_node Character, population events are taken from and gates
#'   are added under (`"root"` for all events).
#' @param channel_x,channel_y Character, the two channels. Must match training.
#' @param ref_train,ref_model_info The `ref_train` and `ref_model_info`
#'   elements from [sw_gate_flowmagic_train()].
#' @param label_to_node Named character vector, predicted label to new gate
#'   name, for example `c("1" = "Live")`. Required: it picks which labels
#'   become gates.
#' @param concavity_val Number, flowMagic's concavity for the hull. Lower
#'   values follow the events more tightly. Required: it changes the polygon.
#' @param n_cores Whole number of cores. Does not change results.
#' @param ... Further arguments for `flowMagic::export_raw_gs_plots()`.
#'
#' @return A named list (per sample) of named lists (per gate name) of
#'   `sw_gate`, with the modified GatingSet in attribute `gs`.
#' @seealso [sw_gate_flowmagic_train()]
#' @export
sw_gate_flowmagic <- function(gs, parent_node, channel_x, channel_y, ref_train, ref_model_info,
                              label_to_node, concavity_val, n_cores = 1, ...) {
  force(parent_node)
  force(label_to_node)
  force(concavity_val)
  if (is.null(names(label_to_node))) {
    stop("`label_to_node` must be a named character vector (label = gate name).", call. = FALSE)
  }
  predict_fn <- .sw_gate_wrap_fn("magic_template_predict", "flowMagic")
  res <- predict_fn(
    gs_input = gs, parent_node = parent_node, channel_x = channel_x, channel_y = channel_y,
    ref_train = ref_train, ref_model_info = ref_model_info, label_to_node = label_to_node,
    n_cores = n_cores, concavity_val = concavity_val, add_to_gs = TRUE,
    return_gs_only = FALSE, ...
  )
  parent <- if (parent_node == "root") NA_character_ else basename(parent_node)
  by_node <- res$list_poly_gates
  samples <- unique(unlist(lapply(by_node, names)))
  out <- stats::setNames(lapply(samples, function(s) {
    got <- lapply(names(by_node), function(node) {
      flt <- by_node[[node]][[s]]
      if (is.null(flt)) NULL else sw_gate_from_flowcore(flt, node, parent)
    })
    names(got) <- names(by_node)
    got[!vapply(got, is.null, logical(1))]
  }), samples)
  attr(out, "gs") <- res$gs
  out
}
