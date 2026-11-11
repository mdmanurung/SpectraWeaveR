# R/faust_wrap.R
# Wrapper for the GPL (>= 3) package faust (RGLab/FAUST, Greene et al. 2021).
# Nothing is copied: each function calls the installed package and reads the
# files it writes. Tests replace the lookup helper .sw_faust_fn().
#
# Use: annotate the same events before and after batch correction, then ask
# whether correction kept the biology (label agreement, per-phenotype shares and
# marker medians, phenotype shares by batch). This is not a fault detector.

# Fetch an exported function from faust, or stop with install help.
.sw_faust_fn <- function(name, pkg = "faust") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '",
      pkg,
      "' is required. Install it with ",
      "remotes::install_github('RGLab/scamp') and ",
      "remotes::install_github('RGLab/FAUST').",
      call. = FALSE
    )
  }
  getExportedValue(pkg, name)
}

# FAUST names its unclassified bin with zero tokens, "0_0_0_0_0"; the length
# follows the marker count, so match the pattern and not the literal.
.sw_faust_is_unclassified <- function(x) {
  grepl("^0(_0)*$", x)
}

# A named list of numeric matrices with column names.
.sw_faust_check_mats <- function(data, arg = "data") {
  ok <- is.list(data) &&
    length(data) > 0 &&
    !is.null(names(data)) &&
    !anyNA(names(data)) &&
    !anyDuplicated(names(data)) &&
    all(vapply(
      data,
      \(m) is.matrix(m) && is.numeric(m) && !is.null(colnames(m)),
      logical(1)
    ))
  if (!ok) {
    stop(
      "`",
      arg,
      "` must be a uniquely named list of numeric matrices with column names.",
      call. = FALSE
    )
  }
}

#' Build a GatingSet from per-file event matrices, ready for FAUST
#'
#' @description
#' Subsample each file, wrap the events in `flowFrame`s whose `desc` field is
#' the marker name, and return a `GatingSet` that [sw_faust_run()] can pass to
#' FAUST. The wrapper never transforms the data.
#'
#' @details
#' **No defaults.** `n_events` and `subsample_seed` change which cells FAUST
#' sees, so you choose them.
#'
#' **Scale.** Pass data already on the analysis scale (arcsinh, cofactor 6000
#' for spectral flow). FAUST reads the values as given.
#'
#' **Row index.** A file with `n_events` rows or fewer keeps all its rows in
#' their order. A larger file keeps `n_events` rows drawn with
#' `subsample_seed`, sorted into original order. The kept row numbers are the
#' attribute `kept_rows`. Two data sets with the same row counts and the same
#' `subsample_seed` keep the same rows, which is what a before and after
#' comparison needs.
#'
#' @param data Uniquely named list of numeric matrices, one per file, events in
#'   rows and channels in columns.
#' @param markers Character vector of the channels to keep. Names are column
#'   names of the matrices and values are the marker names given to FAUST (its
#'   `desc` field). Unnamed values are used as both.
#' @param sample_table Data frame with columns `sample`, `batch`, `donor` and
#'   `role`. `sample` must equal `names(data)`, in the same order.
#' @param n_events Whole number, events kept per file.
#' @param subsample_seed Number, seed for the row draw.
#'
#' @return A `GatingSet` with attributes `kept_rows` (named list of row
#'   indices), `markers` and `sample_table`.
#'
#' @seealso [sw_faust_run()]
#'
#' @examples
#' \dontrun{
#' gs <- sw_faust_gatingset(
#'   data = events, markers = c(`BUV395-A` = "CD45RA", `BV480-A` = "CD3"),
#'   sample_table = tbl, n_events = 5000, subsample_seed = 1
#' )
#' }
#'
#' @export
sw_faust_gatingset <- function(
  data,
  markers,
  sample_table,
  n_events,
  subsample_seed
) {
  .sw_stop_missing(c(
    data = missing(data),
    markers = missing(markers),
    sample_table = missing(sample_table),
    n_events = missing(n_events),
    subsample_seed = missing(subsample_seed)
  ))
  .sw_faust_check_mats(data)
  if (!is.character(markers) || length(markers) == 0 || anyNA(markers)) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (is.null(names(markers))) {
    names(markers) <- markers
  }
  if (anyDuplicated(markers) || anyDuplicated(names(markers))) {
    stop("`markers` must not repeat a channel or a marker name.", call. = FALSE)
  }
  needed <- c("sample", "batch", "donor", "role")
  if (!is.data.frame(sample_table) || !all(needed %in% names(sample_table))) {
    stop(
      "`sample_table` needs columns ",
      paste(needed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  .sw_check_count(n_events, "n_events")
  if (!is.numeric(subsample_seed) || length(subsample_seed) != 1) {
    stop("`subsample_seed` must be one number.", call. = FALSE)
  }
  stopifnot(
    "names(data) equals sample_table$sample" = identical(
      names(data),
      as.character(sample_table$sample)
    )
  )
  absent <- setdiff(names(markers), unique(unlist(lapply(data, colnames))))
  if (length(absent) > 0) {
    stop(
      "`markers` names channels missing from `data`: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Row draw, sorted so the original order is kept.
  kept_rows <- .sw_with_seed(
    subsample_seed,
    lapply(data, function(m) {
      if (nrow(m) <= n_events) {
        seq_len(nrow(m))
      } else {
        sort(sample.int(nrow(m), n_events))
      }
    })
  )

  frames <- lapply(names(data), function(s) {
    m <- data[[s]][kept_rows[[s]], names(markers), drop = FALSE]
    ff <- flowCore::flowFrame(m)
    pd <- flowCore::pData(flowCore::parameters(ff))
    pd$desc <- unname(markers[pd$name])
    flowCore::pData(flowCore::parameters(ff)) <- pd
    ff
  })
  names(frames) <- names(data)
  gs <- flowWorkspace::GatingSet(flowCore::flowSet(frames))

  meta <- data.frame(
    batch = as.character(sample_table$batch),
    donor = as.character(sample_table$donor),
    role = as.character(sample_table$role),
    row.names = names(data)
  )
  meta$name <- rownames(meta)
  flowWorkspace::pData(gs) <- meta[, c("name", "batch", "donor", "role")]

  attr(gs, "kept_rows") <- kept_rows
  attr(gs, "markers") <- markers
  attr(gs, "sample_table") <- sample_table
  gs
}

#' Run FAUST on a GatingSet
#'
#' @description
#' Wrap `faust::faust()`. FAUST writes its results to disk under
#' `project_path`; this function returns that path as a classed object that
#' [sw_faust_labels()] reads.
#'
#' @details
#' **No defaults.** Every argument below changes the phenotypes FAUST finds, so
#' none has one. FAUST's own values are conventions: depth score threshold
#' 0.01, selection quantile 0.5, name occurrence 0, density subsample defaults
#' 1e6, 1e6 and 1.
#'
#' **Run folder.** The function stops if `project_path/faustData` exists. FAUST
#' check-points its intermediates there, so a re-run with new settings needs a
#' new folder.
#'
#' **Fixed inside.** `plottingDevice = "png"`, `debugFlag = FALSE` and
#' `archDescriptionList = list(targetArch = "singleCPU")`.
#'
#' **Experimental unit.** Per-cell labels are read back by file row order,
#' which was checked for `experimental_unit = ""` only. [sw_faust_labels()]
#' stops for any other value.
#'
#' @param gs `GatingSet` from [sw_faust_gatingset()].
#' @param project_path Folder for FAUST output. It must not hold `faustData`.
#' @param starting_pop `"root"` or a gate node name.
#' @param markers Marker names (the `desc` field) to give FAUST. The same
#'   vector for the before and after runs.
#' @param depth_score_threshold Number, FAUST depth score threshold.
#' @param selection_quantile Number in (0, 1), FAUST selection quantile.
#' @param experimental_unit `""` (one unit per sample) or a `pData` column.
#' @param imputation_hierarchy `""` or a `pData` column.
#' @param channel_bounds `""` (empirical bounds) or a 2 x markers matrix.
#' @param name_occurrence_num Whole number >= 0, FAUST `nameOccuranceNum`.
#' @param supervised_list `NA` or a list of boundary preferences per marker.
#' @param seed Number, FAUST `seedValue`.
#' @param threads Whole number, FAUST `threadNum`. Whether the thread count
#'   changes results is unverified.
#' @param density_subsample_threshold,density_subsample_size,density_subsample_iterations
#'   Numbers passed to the `densitySubSample*` arguments.
#' @param annotation_forest_depth Whole number >= 1, FAUST
#'   `annotationForestDepth` (its default is 3). It sets how deep the
#'   annotation forest grows, so which phenotypes FAUST can find.
#' @param annotations_approved `TRUE` runs on to phenotype discovery. `FALSE`
#'   stops after the thresholds so you can review `plotData/scoreLines.*`.
#' @param draw_histograms `TRUE` or `FALSE`, passed to
#'   `drawAnnotationHistograms`. Changes plots only.
#'
#' @return `project_path`, with class `sw_faust_run`. Attribute `status` is
#'   `"complete"` or `"thresholds_only"`; `params` echoes every argument plus
#'   `faust_version`; `kept_rows` and `markers` come from `gs`.
#'
#' @seealso [sw_faust_gatingset()], [sw_faust_labels()]
#'
#' @examples
#' \dontrun{
#' run <- sw_faust_run(
#'   gs, project_path = "faust_before", starting_pop = "root", markers = mk,
#'   depth_score_threshold = 0.01, selection_quantile = 0.5,
#'   experimental_unit = "", imputation_hierarchy = "", channel_bounds = "",
#'   name_occurrence_num = 0, supervised_list = NA, seed = 1, threads = 8,
#'   density_subsample_threshold = 1e6, density_subsample_size = 1e6,
#'   density_subsample_iterations = 1, annotation_forest_depth = 3,
#'   annotations_approved = TRUE,
#'   draw_histograms = FALSE
#' )
#' }
#'
#' @export
sw_faust_run <- function(
  gs,
  project_path,
  starting_pop,
  markers,
  depth_score_threshold,
  selection_quantile,
  experimental_unit,
  imputation_hierarchy,
  channel_bounds,
  name_occurrence_num,
  supervised_list,
  seed,
  threads,
  density_subsample_threshold,
  density_subsample_size,
  density_subsample_iterations,
  annotation_forest_depth,
  annotations_approved,
  draw_histograms
) {
  .sw_stop_missing(c(
    gs = missing(gs),
    project_path = missing(project_path),
    starting_pop = missing(starting_pop),
    markers = missing(markers),
    depth_score_threshold = missing(depth_score_threshold),
    selection_quantile = missing(selection_quantile),
    experimental_unit = missing(experimental_unit),
    imputation_hierarchy = missing(imputation_hierarchy),
    channel_bounds = missing(channel_bounds),
    name_occurrence_num = missing(name_occurrence_num),
    supervised_list = missing(supervised_list),
    seed = missing(seed),
    threads = missing(threads),
    density_subsample_threshold = missing(density_subsample_threshold),
    density_subsample_size = missing(density_subsample_size),
    density_subsample_iterations = missing(density_subsample_iterations),
    annotation_forest_depth = missing(annotation_forest_depth),
    annotations_approved = missing(annotations_approved),
    draw_histograms = missing(draw_histograms)
  ))
  if (!inherits(gs, "GatingSet")) {
    stop("`gs` must be a GatingSet from sw_faust_gatingset().", call. = FALSE)
  }
  one_string <- function(x) is.character(x) && length(x) == 1 && !is.na(x)
  one_number <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x)
  one_flag <- function(x) is.logical(x) && length(x) == 1 && !is.na(x)
  if (!one_string(project_path)) {
    stop("`project_path` must be one path.", call. = FALSE)
  }
  if (file.exists(file.path(project_path, "faustData"))) {
    stop(
      "`project_path` already holds 'faustData'. Use a new folder; FAUST ",
      "check-points there and would reuse old intermediates.",
      call. = FALSE
    )
  }
  if (!one_string(starting_pop)) {
    stop("`starting_pop` must be one string.", call. = FALSE)
  }
  absent <- setdiff(markers, flowWorkspace::markernames(gs))
  if (!is.character(markers) || length(markers) == 0 || length(absent) > 0) {
    stop(
      "`markers` must match the marker names (desc) of `gs`; not found: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  for (nm in c("depth_score_threshold", "seed")) {
    if (!one_number(get(nm))) {
      stop("`", nm, "` must be one number.", call. = FALSE)
    }
  }
  if (
    !one_number(selection_quantile) ||
      selection_quantile <= 0 ||
      selection_quantile >= 1
  ) {
    stop("`selection_quantile` must be one number in (0, 1).", call. = FALSE)
  }
  for (nm in c("experimental_unit", "imputation_hierarchy")) {
    v <- get(nm)
    if (
      !one_string(v) || (v != "" && !v %in% names(flowWorkspace::pData(gs)))
    ) {
      stop(
        "`",
        nm,
        "` must be \"\" or a column of pData(gs).",
        call. = FALSE
      )
    }
  }
  if (
    !(identical(channel_bounds, "") ||
      (is.matrix(channel_bounds) && ncol(channel_bounds) == length(markers)))
  ) {
    stop(
      "`channel_bounds` must be \"\" or a 2 x markers matrix.",
      call. = FALSE
    )
  }
  if (
    !one_number(name_occurrence_num) ||
      name_occurrence_num < 0 ||
      name_occurrence_num != round(name_occurrence_num)
  ) {
    stop("`name_occurrence_num` must be one whole number >= 0.", call. = FALSE)
  }
  .sw_check_count(threads, "threads")
  .sw_check_count(annotation_forest_depth, "annotation_forest_depth")
  for (nm in c(
    "density_subsample_threshold",
    "density_subsample_size",
    "density_subsample_iterations"
  )) {
    if (!one_number(get(nm)) || get(nm) < 1) {
      stop("`", nm, "` must be one number >= 1.", call. = FALSE)
    }
  }
  if (!one_flag(annotations_approved) || !one_flag(draw_histograms)) {
    stop(
      "`annotations_approved` and `draw_histograms` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  dir.create(project_path, recursive = TRUE, showWarnings = FALSE)
  .sw_faust_fn("faust")(
    gatingSet = gs,
    startingCellPop = starting_pop,
    activeChannels = markers,
    channelBounds = channel_bounds,
    experimentalUnit = experimental_unit,
    imputationHierarchy = imputation_hierarchy,
    projectPath = project_path,
    depthScoreThreshold = depth_score_threshold,
    selectionQuantile = selection_quantile,
    nameOccuranceNum = name_occurrence_num,
    supervisedList = supervised_list,
    debugFlag = FALSE,
    threadNum = threads,
    seedValue = seed,
    drawAnnotationHistograms = draw_histograms,
    annotationsApproved = annotations_approved,
    densitySubSampleThreshold = density_subsample_threshold,
    densitySubSampleSize = density_subsample_size,
    densitySubSampleIterations = density_subsample_iterations,
    annotationForestDepth = annotation_forest_depth,
    archDescriptionList = list(targetArch = "singleCPU"),
    plottingDevice = "png"
  )

  complete <- file.exists(file.path(
    project_path,
    "faustData",
    "faustCountMatrix.rds"
  ))
  out <- project_path
  attr(out, "status") <- if (complete) "complete" else "thresholds_only"
  attr(out, "params") <- list(
    starting_pop = starting_pop,
    markers = markers,
    depth_score_threshold = depth_score_threshold,
    selection_quantile = selection_quantile,
    experimental_unit = experimental_unit,
    imputation_hierarchy = imputation_hierarchy,
    channel_bounds = channel_bounds,
    name_occurrence_num = name_occurrence_num,
    supervised_list = supervised_list,
    seed = seed,
    threads = threads,
    density_subsample_threshold = density_subsample_threshold,
    density_subsample_size = density_subsample_size,
    density_subsample_iterations = density_subsample_iterations,
    annotations_approved = annotations_approved,
    draw_histograms = draw_histograms,
    annotation_forest_depth = annotation_forest_depth,
    faust_version = tryCatch(
      as.character(utils::packageVersion("faust")),
      error = \(e) NA_character_
    )
  )
  attr(out, "kept_rows") <- attr(gs, "kept_rows")
  attr(out, "markers") <- attr(gs, "markers")
  class(out) <- c("sw_faust_run", "character")
  out
}

#' Read FAUST phenotype labels
#'
#' @description
#' Read the per-cell phenotype of every kept event and the per-file phenotype
#' shares from a finished [sw_faust_run()].
#'
#' @details
#' Per-cell labels come from `faustData/sampleData/<sample>/faustAnnotation.csv`,
#' one row per kept event in input order (checked against the saved expression
#' matrix on a four-file test, `experimental_unit = ""`). The one run-time
#' check: for each file the sorted label counts must equal the non-zero
#' entries of `faustCountMatrix.rds`, else the function stops.
#'
#' Phenotype strings are FAUST's per-cell tokens (`CD3~2~2~`, level 2 of 2).
#' The unclassified bin is kept as its own phenotype (all-zero tokens).
#'
#' @param run A `sw_faust_run` with status `"complete"`.
#' @param data The per-file matrices given to [sw_faust_gatingset()]. They
#'   only fix the row count of each file.
#'
#' @return A list of class `sw_faust_labels` with `cells` (tibble `sample`,
#'   `cell` as row of the original matrix, `phenotype`), `proportions` (tibble
#'   `sample`, `phenotype`, `n`, `proportion`) and `kept_rows`.
#'
#' @seealso [sw_faust_compare()]
#'
#' @examples
#' \dontrun{
#' labels <- sw_faust_labels(run, events)
#' }
#'
#' @export
sw_faust_labels <- function(run, data) {
  .sw_stop_missing(c(run = missing(run), data = missing(data)))
  if (!inherits(run, "sw_faust_run")) {
    stop("`run` must come from sw_faust_run().", call. = FALSE)
  }
  if (!identical(attr(run, "status"), "complete")) {
    stop(
      "FAUST stopped after the thresholds (status 'thresholds_only'). Review ",
      "'plotData/scoreLines.*' in the project folder, then re-run in a new ",
      "folder with `annotations_approved = TRUE`.",
      call. = FALSE
    )
  }
  if (!identical(attr(run, "params")$experimental_unit, "")) {
    stop(
      "Per-cell row order was verified only for `experimental_unit = \"\"`.",
      call. = FALSE
    )
  }
  .sw_faust_check_mats(data)
  kept <- attr(run, "kept_rows")
  if (!identical(names(kept), names(data))) {
    stop("`data` must hold the same files as the run.", call. = FALSE)
  }
  stopifnot(
    "kept rows lie within data" = all(mapply(
      \(k, m) max(k) <= nrow(m),
      kept,
      data
    ))
  )

  root <- unclass(run)
  counts <- readRDS(file.path(root, "faustData", "faustCountMatrix.rds"))
  cells <- lapply(names(kept), function(s) {
    ann <- utils::read.csv(
      file.path(root, "faustData", "sampleData", s, "faustAnnotation.csv"),
      header = FALSE,
      colClasses = "character"
    )[[1]]
    if (length(ann) != length(kept[[s]])) {
      stop(
        "A file has ",
        length(ann),
        " labels for ",
        length(kept[[s]]),
        " kept events.",
        call. = FALSE
      )
    }
    # Per-file label counts must equal the FAUST count matrix row.
    row <- counts[s, ]
    stopifnot(
      "label counts equal faustCountMatrix.rds" = identical(
        sort(as.numeric(table(ann))),
        sort(as.numeric(row[row > 0]))
      )
    )
    tibble::tibble(sample = s, cell = kept[[s]], phenotype = ann)
  })
  cells <- dplyr::bind_rows(cells)
  proportions <- dplyr::bind_rows(
    lapply(split(cells$phenotype, cells$sample), function(p) {
      tab <- table(p)
      tibble::tibble(
        phenotype = names(tab),
        n = as.integer(tab),
        proportion = as.integer(tab) / length(p)
      )
    }),
    .id = "sample"
  )

  structure(
    list(cells = cells, proportions = proportions, kept_rows = kept),
    class = "sw_faust_labels"
  )
}

# Rewrite each phenotype with its marker~level~levels~ tokens sorted by
# marker, so two runs that list markers in a different order match. The
# unclassified bin is kept as is. Stops if two phenotypes of one run collapse
# to the same sorted string.
.sw_faust_sort_tokens <- function(labels) {
  phenos <- unique(labels$cells$phenotype)
  sorted <- vapply(
    phenos,
    \(p) {
      if (.sw_faust_is_unclassified(p)) {
        return(p)
      }
      tok <- regmatches(p, gregexpr("[^~]+~[0-9]+~[0-9]+~", p))[[1]]
      if (!identical(paste(tok, collapse = ""), p)) {
        stop(
          sprintf(
            "Phenotype `%s` is not a run of marker~level~levels~ tokens.",
            p
          ),
          call. = FALSE
        )
      }
      marker <- sub("~.*$", "", tok)
      paste(tok[order(marker, method = "radix")], collapse = "")
    },
    character(1)
  )
  if (anyDuplicated(sorted)) {
    stop(
      "Two phenotypes of one run share the same marker states after sorting.",
      call. = FALSE
    )
  }
  labels$cells$phenotype <- unname(sorted[labels$cells$phenotype])
  labels$proportions$phenotype <- unname(sorted[labels$proportions$phenotype])
  labels
}

# Median of each marker over the cells of one phenotype, per file or pooled.
.sw_faust_medians <- function(labels, data, markers, scope) {
  cells <- labels$cells
  key <- if (scope == "file") {
    paste(cells$sample, cells$phenotype, sep = "\r")
  } else {
    cells$phenotype
  }
  out <- lapply(split(cells, key), function(d) {
    m <- do.call(
      rbind,
      lapply(split(d, d$sample), \(x) {
        data[[x$sample[[1]]]][x$cell, markers, drop = FALSE]
      })
    )
    tibble::tibble(
      subject = if (scope == "file") d$sample[[1]] else "pooled",
      phenotype = d$phenotype[[1]],
      item = markers,
      value = apply(m, 2, stats::median),
      n_cells = nrow(m)
    )
  })
  dplyr::bind_rows(out)
}


#' Compare FAUST labels before and after batch correction
#'
#' @description
#' Return, on the same events, the label agreement (adjusted Rand index),
#' per-phenotype shares, per-phenotype marker medians, and the spread of each
#' phenotype's share across batches. It reports numbers and gives no pass or
#' fail.
#'
#' @details
#' **No defaults.** `match_by`, `ari_scope` and `median_scope` change the
#' numbers, so you choose them.
#'
#' **Same events.** `before` and `after` must come from the same kept rows,
#' else the function stops. `data_after` must keep the row order of
#' `data_before`; that holds for the corrected tibble only if the correction
#' keeps row order, which the caller must confirm.
#'
#' **Phenotypes on one side only** get `NA` on the other side; they are never
#' dropped. The unclassified bin stays as a row.
#'
#' **Batch and donor.** Donor is nested in batch, so a phenotype that varies by
#' batch among donor files is confounded with donor. `prop_by_batch_*` is
#' split by `role` so reference files can be read apart from donor files.
#'
#' @param before,after `sw_faust_labels` objects.
#' @param data_before,data_after Named lists of per-file matrices, the same
#'   files as the labels.
#' @param batch,role Named character vectors (names = files) giving the batch
#'   and the role (for example reference or donor) of each file.
#' @param markers Character vector of columns of the matrices for the medians.
#' @param match_by `"exact"`: phenotype strings must be equal.
#'   `"marker_state"`: the `marker~level~levels~` tokens of each phenotype are
#'   sorted by marker before matching, so the same marker states match
#'   whatever order each FAUST run lists them in. Phenotypes whose marker sets
#'   differ still get `NA`. The output then holds the sorted strings.
#' @param ari_scope `"file"` or `"pooled"`.
#' @param median_scope `"file"` (per phenotype per file) or `"pooled"` (per
#'   phenotype over all files).
#'
#' @return A tibble of class `sw_faust_compare` with columns `module`
#'   (`"faust"`), `subject`, `check`, `item`, `phenotype`, `value`, `text`.
#'   Checks: `ari`; `prop_before`, `prop_after`, `prop_diff` (item =
#'   phenotype); `median_before`, `median_after`, `median_diff` (item =
#'   marker); `prop_by_batch_before`, `prop_by_batch_after` (subject = role,
#'   value = range over batches of the per-batch median share).
#'
#' @seealso [sw_faust_labels()], [sw_diagnose()]
#'
#' @examples
#' \dontrun{
#' sw_faust_compare(
#'   lab_before, lab_after, ev_before, ev_after,
#'   batch = batch, role = role, markers = chans,
#'   match_by = "exact", ari_scope = "file", median_scope = "file"
#' )
#' }
#'
#' @export
sw_faust_compare <- function(
  before,
  after,
  data_before,
  data_after,
  batch,
  role,
  markers,
  match_by,
  ari_scope,
  median_scope
) {
  .sw_stop_missing(c(
    before = missing(before),
    after = missing(after),
    data_before = missing(data_before),
    data_after = missing(data_after),
    batch = missing(batch),
    role = missing(role),
    markers = missing(markers),
    match_by = missing(match_by),
    ari_scope = missing(ari_scope),
    median_scope = missing(median_scope)
  ))
  if (
    !inherits(before, "sw_faust_labels") || !inherits(after, "sw_faust_labels")
  ) {
    stop(
      "`before` and `after` must come from sw_faust_labels().",
      call. = FALSE
    )
  }
  if (length(match_by) != 1 || !match_by %in% c("exact", "marker_state")) {
    stop("`match_by` must be \"exact\" or \"marker_state\".", call. = FALSE)
  }
  if (match_by == "marker_state") {
    before <- .sw_faust_sort_tokens(before)
    after <- .sw_faust_sort_tokens(after)
  }
  if (!ari_scope %in% c("file", "pooled") || length(ari_scope) != 1) {
    stop("`ari_scope` must be \"file\" or \"pooled\".", call. = FALSE)
  }
  if (!median_scope %in% c("file", "pooled") || length(median_scope) != 1) {
    stop("`median_scope` must be \"file\" or \"pooled\".", call. = FALSE)
  }
  stopifnot(
    "before and after carry the same kept rows" = identical(
      before$kept_rows,
      after$kept_rows
    )
  )
  .sw_faust_check_mats(data_before, "data_before")
  .sw_faust_check_mats(data_after, "data_after")
  files <- names(before$kept_rows)
  if (
    !identical(names(data_before), files) ||
      !identical(names(data_after), files) ||
      !all(files %in% names(batch)) ||
      !all(files %in% names(role))
  ) {
    stop(
      "`data_before`, `data_after`, `batch` and `role` must cover the labelled files.",
      call. = FALSE
    )
  }
  if (
    !all(vapply(
      files,
      \(s) nrow(data_before[[s]]) == nrow(data_after[[s]]),
      logical(1)
    ))
  ) {
    stop(
      "`data_after` must have the row count of `data_before`.",
      call. = FALSE
    )
  }
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a character vector.", call. = FALSE)
  }

  # Label agreement on the same cells.
  stopifnot(
    "before and after label the same cells" = identical(
      before$cells[c("sample", "cell")],
      after$cells[c("sample", "cell")]
    )
  )
  ari_by <- function(idx) {
    .sw_ari(before$cells$phenotype[idx], after$cells$phenotype[idx])
  }
  ari_rows <- if (ari_scope == "file") {
    by_file <- split(seq_len(nrow(before$cells)), before$cells$sample)
    tibble::tibble(
      subject = names(by_file),
      check = "ari",
      item = NA_character_,
      phenotype = NA_character_,
      value = unname(vapply(by_file, ari_by, numeric(1))),
      text = paste0("cells=", lengths(by_file))
    )
  } else {
    tibble::tibble(
      subject = "pooled",
      check = "ari",
      item = NA_character_,
      phenotype = NA_character_,
      value = ari_by(seq_len(nrow(before$cells))),
      text = paste0("cells=", nrow(before$cells))
    )
  }

  # Share matrix, files by phenotypes. A phenotype that never occurs on a side
  # is NA there; one that occurs elsewhere but not in a file is 0.
  phenos <- union(
    unique(before$proportions$phenotype),
    unique(after$proportions$phenotype)
  )
  share <- function(lab) {
    m <- matrix(
      NA_real_,
      length(files),
      length(phenos),
      dimnames = list(files, phenos)
    )
    m[, unique(lab$proportions$phenotype)] <- 0
    p <- lab$proportions
    m[cbind(p$sample, p$phenotype)] <- p$proportion
    m
  }
  share_b <- share(before)
  share_a <- share(after)
  as_rows <- function(m, check) {
    tibble::tibble(
      subject = rep(rownames(m), times = ncol(m)),
      check = check,
      item = rep(colnames(m), each = nrow(m)),
      phenotype = rep(colnames(m), each = nrow(m)),
      value = as.vector(m),
      text = NA_character_
    )
  }
  prop_rows <- dplyr::bind_rows(
    as_rows(share_b, "prop_before"),
    as_rows(share_a, "prop_after"),
    as_rows(share_a - share_b, "prop_diff")
  )

  # Marker medians per phenotype, matched on file, phenotype and marker.
  med_b <- .sw_faust_medians(before, data_before, markers, median_scope)
  med_a <- .sw_faust_medians(after, data_after, markers, median_scope)
  keys <- unique(rbind(
    med_b[c("subject", "phenotype", "item")],
    med_a[c("subject", "phenotype", "item")]
  ))
  key_of <- function(d) paste(d$subject, d$phenotype, d$item, sep = "\r")
  pick <- function(d, col) d[[col]][match(key_of(keys), key_of(d))]
  med_rows <- dplyr::bind_rows(
    lapply(
      list(
        median_before = pick(med_b, "value"),
        median_after = pick(med_a, "value"),
        median_diff = pick(med_a, "value") - pick(med_b, "value")
      ),
      \(v) tibble::tibble(keys, value = v)
    ),
    .id = "check"
  )
  n_b <- pick(med_b, "n_cells")
  n_a <- pick(med_a, "n_cells")
  med_rows$text <- rep(
    paste0("cells_before=", n_b, " cells_after=", n_a),
    times = 3
  )

  # Range over batches of the per-batch median share, split by role.
  by_batch <- function(m, check) {
    groups <- unique(data.frame(
      role = unname(role[files]),
      batch = unname(batch[files])
    ))
    rows <- lapply(unique(groups$role), function(r) {
      bs <- groups$batch[groups$role == r]
      med <- vapply(
        bs,
        \(b) {
          f <- files[role[files] == r & batch[files] == b]
          apply(m[f, , drop = FALSE], 2, \(x) {
            if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
          })
        },
        numeric(ncol(m))
      )
      med <- matrix(med, nrow = ncol(m), dimnames = list(colnames(m), bs))
      tibble::tibble(
        subject = r,
        check = check,
        item = colnames(m),
        phenotype = colnames(m),
        value = apply(med, 1, \(x) {
          if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
        }),
        text = paste0("batches=", length(bs))
      )
    })
    dplyr::bind_rows(rows)
  }
  batch_rows <- dplyr::bind_rows(
    by_batch(share_b, "prop_by_batch_before"),
    by_batch(share_a, "prop_by_batch_after")
  )

  out <- dplyr::bind_rows(ari_rows, prop_rows, med_rows, batch_rows)
  out$module <- "faust"
  out <- out[c(
    "module",
    "subject",
    "check",
    "item",
    "phenotype",
    "value",
    "text"
  )]
  class(out) <- c("sw_faust_compare", class(out))
  out
}
