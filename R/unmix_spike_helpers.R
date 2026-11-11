# Helpers around the spike-in design: gated rows to raw rows, the job table
# for unmixing, and the table of unmixed files for scoring. See unmix_spike.R.

# Population indices of one sample: a GatingSet, or a plain list for tests and
# for gates kept outside flowWorkspace (list(sample = list(node = rows))).
.sw_spike_gate_rows <- function(gs, sample, node) {
  if (methods::is(gs, "GatingSet")) {
    return(which(flowWorkspace::gh_pop_get_indices(gs[[sample]], node)))
  }
  idx <- gs[[sample]][[node]]
  if (is.null(idx)) {
    stop(
      sprintf("No rows for node '%s' in sample '%s'.", node, sample),
      call. = FALSE
    )
  }
  if (is.logical(idx)) which(idx) else as.integer(idx)
}

#' Map gated populations to raw-file rows for a spike design
#'
#' Gates are usually drawn on cleaned files (margins or QC removed), while
#' spiking works on the raw files that are unmixed. This maps each gated
#' (cleaned) row to its raw row through the key columns. Set operations such as
#' "live minus CD3+" stay with the caller, on the returned integer vectors.
#'
#' @param gs A `flowWorkspace::GatingSet` built on `cleaned_files`, or a named
#'   list `list(sample = list(node = rows))` of logical or integer indices on
#'   the cleaned rows.
#' @param nodes Named character: population name = gate path in `gs`.
#' @param cleaned_files Named (by sample) FCS paths or matrices the gates index.
#' @param raw_files Named (by sample) FCS paths or matrices to be spiked.
#' @param key_cols Columns that identify a row in both; must be unique per file.
#'
#' @return A named list `populations[[sample]][[pop]]` of integer raw rows, the
#'   `populations` argument of [sw_unmix_spike_build()].
#'
#' @export
sw_unmix_spike_rows <- function(gs, nodes, cleaned_files, raw_files, key_cols) {
  .sw_stop_missing(c(
    gs = missing(gs),
    nodes = missing(nodes),
    cleaned_files = missing(cleaned_files),
    raw_files = missing(raw_files),
    key_cols = missing(key_cols)
  ))
  if (!is.character(nodes) || is.null(names(nodes))) {
    stop(
      "'nodes' must be a named character vector (population = gate).",
      call. = FALSE
    )
  }
  cleaned_files <- as.list(cleaned_files)
  raw_files <- as.list(raw_files)
  if (
    is.null(names(raw_files)) ||
      !setequal(names(raw_files), names(cleaned_files))
  ) {
    stop(
      "'cleaned_files' and 'raw_files' must be named by the same samples.",
      call. = FALSE
    )
  }
  out <- list()
  for (s in names(raw_files)) {
    raw_keys <- .sw_spike_check_unique_keys(
      .sw_spike_read(raw_files[[s]]),
      key_cols,
      sprintf("Raw file of '%s'", s)
    )
    clean_keys <- .sw_spike_check_unique_keys(
      .sw_spike_read(cleaned_files[[s]]),
      key_cols,
      sprintf("Cleaned file of '%s'", s)
    )
    raw_row <- match(clean_keys, raw_keys)
    if (anyNA(raw_row)) {
      stop(
        sprintf(
          "Sample '%s': %d cleaned row(s) have no raw row with the same key.",
          s,
          sum(is.na(raw_row))
        ),
        call. = FALSE
      )
    }
    out[[s]] <- lapply(nodes, function(node) {
      raw_row[.sw_spike_gate_rows(gs, s, node)]
    })
  }
  out
}

# Folder holding the files of one set: written by build, or the raw folder.
.sw_spike_set_dir <- function(design, set) {
  paths <- design$files$path[design$files$set == set]
  if (length(paths) == 0L && set == "base") {
    paths <- unlist(design$raw)
  }
  dirs <- unique(dirname(paths))
  if (length(dirs) != 1L || anyNA(paths)) {
    stop(
      sprintf("Set '%s' has no single folder of FCS files on disk.", set),
      call. = FALSE
    )
  }
  dirs
}

#' Job table for unmixing a spike design
#'
#' One row per setup and set (the design's sets plus `"base"`): the folder of
#' raw input files and the folder the unmixed files should go to,
#' `<setup root>/<set>`. Map it to your scheduler; the package runs no
#' unmixing.
#'
#' @param design An `sw_spike_design` written to disk.
#' @param setups Named character: setup name = root folder for its unmixed
#'   output.
#'
#' @return A tibble with `setup`, `set`, `raw_dir` and `out_dir`.
#'
#' @export
sw_unmix_spike_jobs <- function(design, setups) {
  .sw_stop_missing(c(design = missing(design), setups = missing(setups)))
  if (!inherits(design, "sw_spike_design")) {
    stop("'design' must be an sw_spike_design.", call. = FALSE)
  }
  if (!is.character(setups) || is.null(names(setups))) {
    stop(
      "'setups' must be a named character vector (setup = root folder).",
      call. = FALSE
    )
  }
  sets <- c(design$sets$set, "base")
  grid <- expand.grid(
    setup = names(setups),
    set = sets,
    stringsAsFactors = FALSE
  )
  tibble::tibble(
    setup = grid$setup,
    set = grid$set,
    raw_dir = vapply(
      grid$set,
      function(s) .sw_spike_set_dir(design, s),
      character(1)
    ),
    out_dir = file.path(unname(setups[grid$setup]), grid$set)
  )
}

#' Table of unmixed files for scoring a spike design
#'
#' Builds the `unmixed` argument of [sw_unmix_spike_score()] for the layout
#' `<setup root>/<set>/<raw name><file_suffix>.fcs`.
#'
#' @param design An `sw_spike_design` written to disk.
#' @param setups Named character: setup name = root folder of its unmixed
#'   output.
#' @param file_suffix Text AutoSpectral appends to the raw file name.
#'
#' @return A tibble with `setup`, `sample`, `set`, `path` and `exists`.
#'
#' @export
sw_unmix_spike_pairs <- function(
  design,
  setups,
  file_suffix = " AutoSpectral"
) {
  .sw_stop_missing(c(design = missing(design), setups = missing(setups)))
  if (!inherits(design, "sw_spike_design")) {
    stop("'design' must be an sw_spike_design.", call. = FALSE)
  }
  if (!all(vapply(design$raw, is.character, logical(1)))) {
    stop("'design' must be built from FCS files.", call. = FALSE)
  }
  if (!is.character(setups) || is.null(names(setups))) {
    stop(
      "'setups' must be a named character vector (setup = root folder).",
      call. = FALSE
    )
  }
  sets <- c(design$sets$set, "base")
  grid <- expand.grid(
    setup = names(setups),
    sample = names(design$raw),
    set = sets,
    stringsAsFactors = FALSE
  )
  stem <- tools::file_path_sans_ext(basename(unlist(design$raw)[grid$sample]))
  path <- file.path(
    unname(setups[grid$setup]),
    grid$set,
    paste0(stem, file_suffix, ".fcs")
  )
  tibble::tibble(
    setup = grid$setup,
    sample = grid$sample,
    set = grid$set,
    path = path,
    exists = file.exists(path)
  )
}
