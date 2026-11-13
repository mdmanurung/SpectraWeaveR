# R/report_import.R
# Purpose: write the import-step results (I1 to I6 and I9) into one new run
#   folder: tables, figures with their source data, provenance and session info.
# Run: sw_report_import(result, out_dir) with a named list of function outputs.
# Inputs:
#   result  named list; any subset of bound, header_audit, integrity, signatures,
#           compensation, batch_date, dye_map, channel_summary
# Outputs (under out_dir):
#   tables/*.csv         every table of every input, nothing filtered
#   figures/*.png        one per figure panel
#   source_data/*.csv    the data frame passed to ggplot, one per figure
#   provenance.json      package, R and input provenance, inputs present and absent
#   sessionInfo.txt      sessionInfo() written last
# Reused: .sw_stop_missing() (R/corrector.R), .sw_pkg_info() (R/unmix_report.R),
#   .sw_fig_theme() and .sw_fig_wrap() (R/unmix_figures.R). sw_unmix_save_figures()
#   is tied to panels F1 to F8, so the PNG and source CSV writing is repeated here.

# Settings ----------------------------------------------------------------------
# Figure sizes are layout only; none of them changes a number in a table.

# Two-column journal width, as in the unmixing figures.
.sw_ri_width_mm <- 180
# One file row needs about 2.6 mm for a 7 pt label to stay readable.
.sw_ri_row_mm <- 2.6
# Axis titles, legend and caption.
.sw_ri_base_mm <- 45
.sw_ri_min_height_mm <- 70
.sw_ri_dpi <- 300
# Height of one facet row in the channel summary figure.
.sw_ri_facet_mm <- 50

# Names a result may hold, and the class each one must carry.
.sw_ri_inputs <- c(
  "bound",
  "header_audit",
  "integrity",
  "signatures",
  "compensation",
  "batch_date",
  "dye_map",
  "channel_summary"
)
.sw_ri_classes <- c(
  header_audit = "sw_fcs_header_audit",
  integrity = "sw_fcs_integrity_audit",
  signatures = "sw_fcs_header_signatures",
  compensation = "sw_fcs_compensation_state",
  batch_date = "sw_fcs_batch_date_check",
  dye_map = "sw_fcs_marker_dye_map",
  channel_summary = "sw_fcs_signal_audit"
)
.sw_ri_bound_parts <- c(
  "bound",
  "manifest_no_file",
  "file_no_manifest",
  "multi_dataset",
  "not_fcs"
)

# Okabe-Ito colours keyed by label; an unknown label stops the figure.
.sw_ri_pal_header <- c(
  match = "grey85",
  differs = "#D55E00",
  missing_in_file = "#0072B2",
  absent_in_reference = "#E69F00"
)
.sw_ri_pal_status <- c(
  ok = "grey70",
  assessed = "grey70",
  match = "grey70",
  differs = "#D55E00",
  not_assessed = "#CC79A7"
)
.sw_ri_pal_marker <- c(
  missing = "#0072B2",
  extra = "#E69F00",
  no_deviation = "grey85",
  not_assessed = "#CC79A7"
)
.sw_ri_pal_state <- c(
  "TRUE" = "#009E73",
  "FALSE" = "#D55E00",
  "no value" = "grey80",
  "not_assessed" = "#CC79A7"
)


#' Write the import-step results into a new report folder
#'
#' @description
#' Takes the outputs of the import-step functions and writes them as tables,
#' figures with their source data, a provenance file and the session info into
#' a new folder. The report holds numbers only. It draws no cutoff and gives no
#' pass or fail.
#'
#' @details
#' **What is written.** Under `out_dir`: `tables/*.csv` (every table of every
#' input, including unmatched and not-assessed rows), `figures/*.png`,
#' `source_data/<panel>.csv` (exactly the data frame passed to ggplot, one per
#' figure), `provenance.json` and `sessionInfo.txt`.
#'
#' **Nothing is skipped silently.** An element of `result` that you leave out
#' is listed in `absent` of the return value and in `provenance.json`. A panel
#' that cannot be drawn gets the text `table fallback: <reason>` in `figures`
#' instead of a path. A file that could not be assessed stays in the plotted
#' data with an empty value, so its axis label shows but it has no point.
#'
#' **Figures.** Tile maps (header audit, signatures, batch by date, marker by
#' batch, dye map), dot panels (`$TOT` difference, channel summary) and a file
#' by state grid (compensation). Figures are proposals: inspect each one before
#' you draw a conclusion from it.
#'
#' **No overwrite.** `out_dir` must not exist. There is no overwrite argument.
#'
#' @param result A named list holding any subset of: `bound`
#'   ([sw_manifest_bind_files()]), `header_audit` ([sw_fcs_header_audit()]),
#'   `integrity` ([sw_fcs_integrity_audit()]), `signatures`
#'   ([sw_fcs_header_signatures()]), `compensation`
#'   ([sw_fcs_compensation_state()]), `batch_date`
#'   ([sw_fcs_batch_date_check()]), `dye_map` ([sw_fcs_marker_dye_map()]) and
#'   `channel_summary` ([sw_fcs_signal_audit()] run with `summary` set).
#'   Required.
#' @param out_dir Path of a folder that does not exist yet. Required.
#'
#' @return A list with `tables` (named character vector of CSV paths),
#'   `figures` (named list; a PNG path, or `table fallback: <reason>`),
#'   `source_data` (named character vector of CSV paths, one per drawn figure),
#'   `dir` (`out_dir`) and `absent` (names of inputs that were not given).
#'
#' @seealso [sw_fcs_header_audit()], [sw_fcs_integrity_audit()]
#'
#' @examples
#' \dontrun{
#' res <- list(
#'   dye_map = sw_fcs_marker_dye_map(paths, batch = batch, separator = " "),
#'   integrity = sw_fcs_integrity_audit(
#'     paths,
#'     tot_tolerance = 0,
#'     time_keyword = "Time",
#'     header_only = FALSE
#'   )
#' )
#' sw_report_import(res, out_dir = "runs/2026-11-13_import")
#' }
#'
#' @export
sw_report_import <- function(result, out_dir) {
  .sw_stop_missing(c(result = missing(result), out_dir = missing(out_dir)))
  .sw_ri_check_result(result)
  if (
    !is.character(out_dir) ||
      length(out_dir) != 1L ||
      is.na(out_dir) ||
      !nzchar(out_dir)
  ) {
    stop("`out_dir` must be one folder path.", call. = FALSE)
  }
  if (file.exists(out_dir)) {
    stop(
      sprintf(
        "`out_dir` exists: %s. The report never overwrites; choose a new folder.",
        out_dir
      ),
      call. = FALSE
    )
  }
  for (pkg in c("ggplot2", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' is required for sw_report_import().", pkg),
        call. = FALSE
      )
    }
  }

  # Build every table and figure object in memory first. Rendering happens in the
  # write block below, so that block removes the folder again if it fails.
  present <- intersect(.sw_ri_inputs, names(result))
  absent <- setdiff(.sw_ri_inputs, present)
  tables <- unlist(
    lapply(present, \(nm) .sw_ri_tables(nm, result[[nm]])),
    recursive = FALSE
  )
  panels <- unlist(
    lapply(present, \(nm) .sw_ri_panels(nm, result[[nm]])),
    recursive = FALSE
  )

  # Write run folder ------------------------------------------------------------
  made <- suppressWarnings(
    dir.create(file.path(out_dir, "tables"), recursive = TRUE)
  )
  if (!made) {
    stop(
      sprintf("Cannot create `out_dir`: %s.", out_dir),
      call. = FALSE
    )
  }
  finished <- FALSE
  on.exit(if (!finished) unlink(out_dir, recursive = TRUE), add = TRUE)
  for (sub in c("figures", "source_data")) {
    if (!dir.create(file.path(out_dir, sub))) {
      stop(sprintf("Cannot create %s in `out_dir`.", sub), call. = FALSE)
    }
  }
  table_paths <- vapply(
    names(tables),
    \(nm) {
      path <- file.path(out_dir, "tables", paste0(nm, ".csv"))
      utils::write.csv(.sw_ri_flat(tables[[nm]]), path, row.names = FALSE)
      path
    },
    character(1)
  )
  is_plot <- vapply(panels, is.list, logical(1))
  figure_paths <- lapply(names(panels), \(nm) {
    if (!is_plot[[nm]]) {
      return(panels[[nm]])
    }
    .sw_ri_save(panels[[nm]], nm, out_dir)
  })
  names(figure_paths) <- names(panels)
  source_paths <- vapply(
    names(panels)[is_plot],
    \(nm) {
      path <- file.path(out_dir, "source_data", paste0(nm, ".csv"))
      utils::write.csv(
        .sw_ri_flat(panels[[nm]]$frame),
        path,
        row.names = FALSE
      )
      path
    },
    character(1)
  )

  prov <- list(
    package_version = .sw_pkg_info("SpectraWeaveR")$version,
    r_version = R.version.string,
    platform = R.version$platform,
    datetime = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    seed = NA,
    packages = do.call(
      dplyr::bind_rows,
      lapply(
        c("SpectraWeaveR", "flowCore", "ggplot2", "jsonlite"),
        .sw_pkg_info
      )
    ),
    inputs = I(present),
    absent = I(absent),
    figures = lapply(figure_paths, \(f) if (file.exists(f)) basename(f) else f),
    input_provenance = stats::setNames(
      lapply(present, \(nm) .sw_ri_input_prov(nm, result[[nm]])),
      present
    )
  )
  writeLines(
    jsonlite::toJSON(
      prov,
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null",
      na = "null",
      digits = NA,
      force = TRUE
    ),
    file.path(out_dir, "provenance.json")
  )
  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
  finished <- TRUE

  list(
    tables = table_paths,
    figures = figure_paths,
    source_data = source_paths,
    dir = out_dir,
    absent = absent
  )
}


# Validation --------------------------------------------------------------------

.sw_ri_check_result <- function(result) {
  if (is.list(result) && length(result) == 0L) {
    stop(
      "`result` has no recognised element. Accepted names: ",
      paste(.sw_ri_inputs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  keys <- names(result)
  if (!is.list(result) || is.null(keys) || anyNA(keys) || any(keys == "")) {
    stop("`result` must be a named list of import-step outputs.", call. = FALSE)
  }
  unknown <- setdiff(keys, .sw_ri_inputs)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown element(s) in `result`: %s. Accepted: %s.",
        paste(unknown, collapse = ", "),
        paste(.sw_ri_inputs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(keys)) {
    stop("Names in `result` must be unique.", call. = FALSE)
  }
  for (nm in setdiff(keys, "bound")) {
    if (!inherits(result[[nm]], .sw_ri_classes[[nm]])) {
      stop(
        sprintf(
          "`result$%s` must be the output of %s().",
          nm,
          .sw_ri_classes[[nm]]
        ),
        call. = FALSE
      )
    }
  }
  if (
    "bound" %in%
      keys &&
      !(is.list(result$bound) &&
        all(.sw_ri_bound_parts %in% names(result$bound)))
  ) {
    stop(
      "`result$bound` must be the output of sw_manifest_bind_files().",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# Tables ------------------------------------------------------------------------

# Every table of one input, named for its CSV file. NULL parts are not written.
.sw_ri_tables <- function(name, x) {
  parts <- switch(
    name,
    bound = list(
      bound_files = x$bound,
      bound_manifest_no_file = x$manifest_no_file,
      bound_file_no_manifest = x$file_no_manifest,
      bound_multi_dataset = x$multi_dataset,
      bound_not_fcs = x$not_fcs
    ),
    header_audit = list(header_audit = x),
    integrity = list(integrity = x),
    signatures = list(
      signatures = x$files,
      signatures_counts = x$counts,
      signatures_keywords = x$keywords,
      signatures_panel_match = x$panel_match
    ),
    compensation = list(compensation_state = x),
    batch_date = list(
      batch_date = x$files,
      batch_date_crosstab = x$crosstab,
      batch_date_groups = x$groups,
      batch_date_unmatched = x$unmatched
    ),
    dye_map = list(
      dye_map = x$map,
      dye_map_marker_batch = x$marker_batch,
      dye_map_unsplit = x$unsplit
    ),
    channel_summary = list(channel_summary = x)
  )
  Filter(Negate(is.null), parts)
}

# A list column cannot go into a CSV; join its elements with "; ".
.sw_ri_flat <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  attr(df, "provenance") <- NULL
  is_list <- vapply(df, is.list, logical(1))
  for (col in names(df)[is_list]) {
    df[[col]] <- vapply(
      df[[col]],
      \(v) paste(format(v), collapse = "; "),
      character(1)
    )
  }
  df
}

# The provenance block an input carries, or NULL when it has none.
.sw_ri_input_prov <- function(name, x) {
  prov <- if (name %in% c("bound", "signatures", "batch_date")) {
    x$provenance
  } else {
    attr(x, "provenance")
  }
  if (is.null(prov)) {
    return(list(note = "this input carries no provenance"))
  }
  prov
}


# Figures -----------------------------------------------------------------------

# Panels of one input: a named list; each element is list(frame, plot, width,
# height) or a "table fallback: <reason>" string.
.sw_ri_panels <- function(name, x) {
  switch(
    name,
    bound = list(
      bound_files = "table fallback: names and counts, nothing to plot"
    ),
    header_audit = list(header_audit_tiles = .sw_ri_fig_header(x)),
    integrity = list(integrity_tot = .sw_ri_fig_integrity(x)),
    signatures = list(
      signatures_files = .sw_ri_fig_sig_files(x),
      signatures_markers = .sw_ri_fig_sig_markers(x)
    ),
    compensation = list(compensation_grid = .sw_ri_fig_comp(x)),
    batch_date = list(batch_date_tiles = .sw_ri_fig_date(x)),
    dye_map = list(dye_map_tiles = .sw_ri_fig_dye(x)),
    channel_summary = list(channel_summary_dots = .sw_ri_fig_channel(x))
  )
}

# Colour scale from a named palette; a label without a colour stops the figure.
.sw_ri_fill <- function(values, palette, title, aesthetic = "fill") {
  unknown <- setdiff(unique(values), names(palette))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "No colour for %s: %s.",
        title,
        paste(shQuote(unknown), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (aesthetic == "fill") {
    ggplot2::scale_fill_manual(values = palette, name = title)
  } else {
    ggplot2::scale_colour_manual(values = palette, name = title)
  }
}

.sw_ri_height <- function(n_rows) {
  max(.sw_ri_min_height_mm, .sw_ri_base_mm + n_rows * .sw_ri_row_mm)
}
.sw_ri_width <- function(n_cols) {
  max(.sw_ri_width_mm, 30 + n_cols * .sw_ri_row_mm)
}

# File x header item tile map; fill = status against the reference file.
.sw_ri_fig_header <- function(x) {
  frame <- tibble::tibble(
    file = x$file,
    item = trimws(paste(
      x$kind,
      x$key,
      ifelse(is.na(x$channel), "", paste0("[", x$channel, "]"))
    )),
    status = x$status
  )
  if (nrow(frame) == 0L) {
    return("table fallback: the header audit has no rows")
  }
  caption <- .sw_fig_wrap(sprintf(
    "Reference file: %s. Tiles are compared as text after trimming white space.",
    attr(x, "reference")
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = item, y = file, fill = status)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    .sw_ri_fill(frame$status, .sw_ri_pal_header, "Status vs reference") +
    ggplot2::scale_x_discrete(limits = unique(frame$item)) +
    ggplot2::scale_y_discrete(limits = rev(unique(frame$file))) +
    ggplot2::labs(x = "Header item", y = "File", caption = caption) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 4
      )
    )
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(length(unique(frame$item))),
    height = .sw_ri_height(length(unique(frame$file)))
  )
}

# One dot per file for $TOT minus rows read, ordered by value.
.sw_ri_fig_integrity <- function(x) {
  rows <- x[x$item == "tot_vs_rows", ]
  rows <- rows[order(rows$value, rows$file, na.last = TRUE), ]
  frame <- tibble::tibble(
    file = rows$file,
    value = rows$value,
    status = rows$status,
    reason = rows$reason
  )
  if (all(is.na(frame$value))) {
    return(paste0(
      "table fallback: no $TOT value was assessed ",
      "(header_only is TRUE or no file could be read)"
    ))
  }
  caption <- .sw_fig_wrap(paste(
    "Files without a dot were not assessed; see tables/integrity.csv.",
    "header_only =",
    attr(x, "header_only")
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = file, y = value, colour = status)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
    ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
    .sw_ri_fill(frame$status, .sw_ri_pal_status, "Status", "colour") +
    ggplot2::scale_x_discrete(limits = frame$file) +
    ggplot2::labs(
      x = "File, ordered by value",
      y = "$TOT minus rows read (events)",
      caption = caption
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 4
      )
    )
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(nrow(frame)),
    height = 90
  )
}

# File x signature tile map; files in the same signature sit in the same column.
.sw_ri_fig_sig_files <- function(x) {
  f <- x$files
  frame <- tibble::tibble(
    file = f$file,
    signature = ifelse(is.na(f$signature), "not_assessed", f$signature),
    n_channels = f$n_channels,
    status = f$status
  )
  sig_levels <- c(
    sort(setdiff(frame$signature, "not_assessed")),
    intersect("not_assessed", frame$signature)
  )
  frame <- frame[order(match(frame$signature, sig_levels), frame$file), ]
  caption <- .sw_fig_wrap(sprintf(
    "%d signature(s); each definition is in tables/signatures_counts.csv.",
    nrow(x$counts)
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = signature, y = file, fill = status)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    .sw_ri_fill(frame$status, .sw_ri_pal_status, "Status") +
    ggplot2::scale_x_discrete(limits = sig_levels) +
    ggplot2::scale_y_discrete(limits = rev(frame$file)) +
    ggplot2::labs(x = "Header signature", y = "File", caption = caption) +
    .sw_fig_theme()
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width_mm,
    height = .sw_ri_height(nrow(frame))
  )
}

# File x marker tile map from the panel match: markers missing or extra per panel.
.sw_ri_fig_sig_markers <- function(x) {
  pm <- x$panel_match
  if (is.null(pm)) {
    return("table fallback: panel was NULL, so no panel match was run")
  }
  split_names <- \(v) {
    if (is.na(v) || !nzchar(v)) {
      character()
    } else {
      strsplit(v, "; ", fixed = TRUE)[[1]]
    }
  }
  rows <- lapply(seq_len(nrow(pm)), \(i) {
    if (pm$status[i] == "not_assessed") {
      markers <- "(not assessed)"
      state <- "not_assessed"
      return(tibble::tibble(
        file = pm$file[i],
        panel = pm$panel[i],
        marker = markers,
        state = state
      ))
    }
    miss <- split_names(pm$missing[i])
    extra <- split_names(pm$extra[i])
    if (length(miss) + length(extra) == 0L) {
      return(tibble::tibble(
        file = pm$file[i],
        panel = pm$panel[i],
        marker = "(no difference)",
        state = "no_deviation"
      ))
    }
    tibble::tibble(
      file = pm$file[i],
      panel = pm$panel[i],
      marker = c(miss, extra),
      state = c(rep("missing", length(miss)), rep("extra", length(extra)))
    )
  })
  frame <- dplyr::bind_rows(rows)
  caption <- .sw_fig_wrap(paste(
    "Only markers that differ from the panel are drawn.",
    "Counts of matched markers are in tables/signatures_panel_match.csv."
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = marker, y = file, fill = state)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(panel),
      scales = "free_x",
      space = "free_x"
    ) +
    .sw_ri_fill(frame$state, .sw_ri_pal_marker, "Marker vs panel") +
    ggplot2::scale_y_discrete(limits = rev(unique(pm$file))) +
    ggplot2::labs(x = "Marker", y = "File", caption = caption) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 4
      )
    )
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(length(unique(paste(frame$panel, frame$marker)))),
    height = .sw_ri_height(length(unique(pm$file)))
  )
}

# File x state grid: one point per file and item.
.sw_ri_fig_comp <- function(x) {
  items <- c("matrix_present", "matrix_identity", "apply_parsed")
  frame <- dplyr::bind_rows(lapply(items, \(it) {
    val <- x[[it]]
    tibble::tibble(
      file = x$file,
      item = it,
      state = ifelse(
        x$header_status != "read",
        "not_assessed",
        ifelse(is.na(val), "no value", as.character(val))
      )
    )
  }))
  caption <- .sw_fig_wrap(paste(
    "no value: the item is empty for a readable file (for example the matrix",
    "is absent, so it cannot be identity). not_assessed: the header was not read."
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = item, y = file, colour = state)
  ) +
    ggplot2::geom_point(shape = 15, size = 2) +
    .sw_ri_fill(frame$state, .sw_ri_pal_state, "State", "colour") +
    ggplot2::scale_x_discrete(limits = items) +
    ggplot2::scale_y_discrete(limits = rev(unique(frame$file))) +
    ggplot2::labs(x = "Compensation item", y = "File", caption = caption) +
    .sw_fig_theme()
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width_mm,
    height = .sw_ri_height(length(unique(frame$file)))
  )
}

# Declared batch x acquisition date tile map; fill = number of files.
.sw_ri_fig_date <- function(x) {
  ct <- x$crosstab
  if (nrow(ct) == 0L) {
    return("table fallback: the batch by date table has no rows")
  }
  frame <- tibble::tibble(
    batch = ct$batch,
    date = ifelse(is.na(ct$date), "no date", format(ct$date)),
    n_files = ct$n_files
  )
  dates <- c(
    sort(setdiff(unique(frame$date), "no date")),
    intersect("no date", frame$date)
  )
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = date, y = batch, fill = n_files)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = n_files), size = 2.2) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#0072B2",
      limits = c(0, NA),
      name = "Files"
    ) +
    ggplot2::scale_x_discrete(limits = dates) +
    ggplot2::scale_y_discrete(limits = rev(unique(frame$batch))) +
    ggplot2::labs(
      x = "Acquisition date",
      y = "Declared batch",
      caption = .sw_fig_wrap(
        "Zero counts are drawn; a file without a date is counted under no date."
      )
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(length(dates)),
    height = .sw_ri_height(length(unique(frame$batch)))
  )
}

# Marker x batch tile map; text = dyes, fill = distinct dyes for the marker.
.sw_ri_fig_dye <- function(x) {
  mb <- x$marker_batch
  if (nrow(mb) == 0L) {
    return("table fallback: no channel name was split into marker and dye")
  }
  frame <- tibble::tibble(
    marker = mb$marker,
    batch = mb$batch,
    dyes = mb$dyes,
    n_dyes_marker = mb$n_dyes_marker
  )
  caption <- .sw_fig_wrap(sprintf(
    "%d channel name(s) did not split and are listed in tables/dye_map_unsplit.csv.",
    nrow(x$unsplit)
  ))
  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = batch, y = marker, fill = factor(n_dyes_marker))
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = dyes), size = 2) +
    ggplot2::scale_fill_viridis_d(
      end = 0.85,
      name = "Distinct dyes\nfor the marker"
    ) +
    ggplot2::scale_y_discrete(limits = rev(sort(unique(frame$marker)))) +
    ggplot2::labs(x = "Batch", y = "Marker", caption = caption) +
    .sw_fig_theme()
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(length(unique(frame$batch))),
    height = .sw_ri_height(length(unique(frame$marker)))
  )
}

# Per-file location and IQR, one panel per channel with a free y axis.
.sw_ri_fig_channel <- function(x) {
  d <- x[startsWith(x$check, "summary_") & !is.na(x$channel), ]
  loc_check <- intersect(c("summary_median", "summary_mean"), unique(d$check))
  if (length(loc_check) != 1L) {
    return(paste0(
      "table fallback: no summary rows; ",
      "sw_fcs_signal_audit() was run with summary = NULL"
    ))
  }
  stat <- sub("summary_", "", loc_check, fixed = TRUE)
  key <- unique(d[c("file", "channel")])
  pick <- \(chk, col = "value") {
    s <- d[d$check == chk, ]
    s[[col]][match(paste(key$file, key$channel), paste(s$file, s$channel))]
  }
  frame <- tibble::tibble(
    file = key$file,
    channel = key$channel,
    location = pick(loc_check),
    q25 = pick("summary_q0.25"),
    q75 = pick("summary_q0.75"),
    n_events = pick("summary_n"),
    detail = pick(loc_check, "detail")
  )
  files <- unique(x$file)
  caption <- .sw_fig_wrap(paste(
    "Bars span the 25th to 75th percentile and need probs 0.25 and 0.75 in summary.",
    "Files are in the order they were passed. A file with no dot has too few finite events;",
    "see tables/channel_summary.csv."
  ))
  plot <- ggplot2::ggplot(frame, ggplot2::aes(x = file, y = location)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = q25, ymax = q75),
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 1, na.rm = TRUE) +
    ggplot2::facet_wrap(
      ggplot2::vars(channel),
      ncol = 2,
      scales = "free_y",
      axes = "all",
      axis.labels = "margins"
    ) +
    ggplot2::scale_x_discrete(limits = files) +
    ggplot2::labs(
      x = "File, in input order",
      y = sprintf("Per-file %s (units of the data given to the audit)", stat),
      caption = caption
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 4
      )
    )
  list(
    frame = frame,
    plot = plot,
    width = .sw_ri_width(length(files) / 2),
    height = .sw_ri_base_mm +
      ceiling(length(unique(frame$channel)) / 2) * .sw_ri_facet_mm
  )
}

# Save one panel as PNG; its source CSV is written by the caller.
.sw_ri_save <- function(panel, name, out_dir) {
  path <- file.path(out_dir, "figures", paste0(name, ".png"))
  ggplot2::ggsave(
    path,
    plot = panel$plot,
    width = panel$width,
    height = panel$height,
    units = "mm",
    dpi = .sw_ri_dpi,
    bg = "white",
    limitsize = FALSE
  )
  path
}

# Bare column names used inside ggplot2::aes().
utils::globalVariables(c(
  "item",
  "file",
  "status",
  "value",
  "signature",
  "marker",
  "state",
  "date",
  "batch",
  "n_files",
  "dyes",
  "n_dyes_marker",
  "location",
  "q25",
  "q75",
  "label"
))
