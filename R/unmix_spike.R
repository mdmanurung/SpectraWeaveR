#' @title Spike-in recovery for unmixing
#'
#' @description
#' Add a known amount of one dye's spectrum to raw detector values of chosen
#' cells, let the user unmix the spiked and unspiked files with any setup, and
#' measure how much of the spike each setup gives back.
#'
#' * [sw_unmix_spike_build()] writes spiked copies of raw files, one set per
#'   marker, and a key of which cell got which spike.
#' * [sw_unmix_spike_cospike()] spikes a second dye into the same cells of one
#'   level, as benchmark 05 script 88d did.
#' * [sw_unmix_spike_score()] reads the unmixed files and returns recovery per
#'   cell, summaries with quantiles, and an independence check.
#'
#' @details
#' The package never runs or submits unmixing. Unmix every set of the design,
#' and the unspiked `base` set, with each setup you want to compare, then pass
#' the file table to [sw_unmix_spike_score()]. [sw_unmix_spike_jobs()] and
#' [sw_unmix_spike_pairs()] build those tables for the folder layout
#' `<setup root>/<set>/`.
#'
#' Recovery of a cell is `(spiked unmix - base unmix) / k` in the dye's
#' unmixed channel, with `k = cofactor * sinh(level)`. This difference measures
#' the spike only if a cell unmixes independently of the other cells in its
#' file. Score checks that assumption on every unspiked cell, and marks each
#' file `independent` or `coupled`. Recovery values are kept for coupled files;
#' read them with the status.
#'
#' The spike is noise-free: no photon noise is added with the signal, so
#' recovery is an optimistic test for noise-weighted methods.
#'
#' Method source: benchmark 05 of this repository (scripts 62, 64, 88, 88d,
#' 88e, 88h). Every analysis choice is a required argument (ROADMAP item B11,
#' register rows "Spike levels" to "Co-spike level").
#'
#' @name unmix_spike
#' @keywords internal
NULL

# ---- Constants ---------------------------------------------------------------

# Channel name patterns that are not fluorescence (scatter and time).
.sw_spike_nonfluor_pattern <- "^(FSC|SSC|Time)"

# Tolerance for "each spiked spectrum row has max 1".
.sw_spike_max_tol <- 1e-6

# ---- Internal helpers --------------------------------------------------------

# Read one raw input: an FCS path or a matrix with column names.
.sw_spike_read <- function(x) {
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("Every raw matrix needs column names.", call. = FALSE)
    }
    return(x)
  }
  if (!is.character(x) || length(x) != 1L) {
    stop("Each raw input must be an FCS path or a matrix.", call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required to read FCS files.", call. = FALSE)
  }
  flowCore::exprs(flowCore::read.FCS(
    x,
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
}

# One string per row from the alignment key columns.
.sw_spike_keys <- function(x, key_cols) {
  do.call(paste, c(as.data.frame(x[, key_cols, drop = FALSE]), sep = "|"))
}

# Stop when the key columns do not identify every row of a file.
.sw_spike_check_unique_keys <- function(x, key_cols, label) {
  keys <- .sw_spike_keys(x, key_cols)
  n_dup <- sum(duplicated(keys))
  if (n_dup > 0L) {
    stop(
      sprintf(
        "%s: %d row(s) share a key on %s; the key must be unique per file.",
        label,
        n_dup,
        paste(key_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(keys)
}

# Stop unless each named spectrum row peaks at 1, so k is the added peak.
.sw_spike_check_spectrum_rows <- function(spectra, rows) {
  peaks <- apply(abs(spectra[rows, , drop = FALSE]), 1, max)
  bad <- names(peaks)[abs(peaks - 1) > .sw_spike_max_tol]
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Spike spectrum rows must have max 1; not so for: %s.",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Run expr with a seed, leaving the caller's random stream as it was.
.sw_spike_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  set.seed(seed)
  expr
}

# Levels per population: round(fraction * n) cells get a level above 0, spread
# evenly over the non-zero levels; the rest get 0. Drawn per population (R4).
.sw_spike_draw <- function(populations, levels, spike_fraction, seed) {
  nonzero <- sort(levels[levels > 0])
  .sw_spike_with_seed(seed, {
    rows <- list()
    for (s in sort(names(populations))) {
      for (p in names(populations[[s]])) {
        r <- populations[[s]][[p]]
        n_spike <- round(spike_fraction * length(r))
        lv <- c(rep_len(nonzero, n_spike), rep(0, length(r) - n_spike))
        rows[[length(rows) + 1L]] <- tibble::tibble(
          sample = s,
          raw_row = as.integer(r),
          pop = p,
          level = sample(lv)
        )
      }
    }
    dplyr::bind_rows(rows)
  })
}

# Raise $PnR and the parameter range where the spiked maximum exceeds it.
.sw_spike_raise_range <- function(ff) {
  x <- flowCore::exprs(ff)
  pd <- flowCore::parameters(ff)
  pdata <- methods::slot(pd, "data")
  kw <- flowCore::keyword(ff)
  for (i in seq_len(nrow(pdata))) {
    top <- max(x[, pdata$name[i]], na.rm = TRUE)
    if (is.finite(top) && top >= pdata$range[i]) {
      new_range <- ceiling(top) + 1
      pdata$range[i] <- new_range
      pdata$maxRange[i] <- new_range - 1
      kw[[paste0("$P", i, "R")]] <- as.character(new_range)
    }
  }
  methods::slot(pd, "data") <- pdata
  flowCore::parameters(ff) <- pd
  flowCore::keyword(ff) <- kw
  ff
}

# Write one spiked matrix as FCS, reusing the raw file's header when there is one.
.sw_spike_write_fcs <- function(x, raw_input, path, keywords) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required to write FCS files.", call. = FALSE)
  }
  if (is.character(raw_input)) {
    ff <- flowCore::read.FCS(
      raw_input,
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    flowCore::exprs(ff) <- x
  } else {
    ff <- flowCore::flowFrame(x)
  }
  ff <- .sw_spike_raise_range(ff)
  kw <- flowCore::keyword(ff)
  kw[names(keywords)] <- keywords
  flowCore::keyword(ff) <- kw
  flowCore::write.FCS(ff, path)
  invisible(path)
}

# File name of one sample's spiked copy.
.sw_spike_file_name <- function(raw_input, sample) {
  if (is.character(raw_input)) basename(raw_input) else paste0(sample, ".fcs")
}

# Write or keep the spiked matrices of one set; return the file table rows.
.sw_spike_emit_set <- function(spiked, raw, set, output_dir, keywords) {
  if (is.null(output_dir)) {
    return(tibble::tibble(
      set = set,
      sample = names(spiked),
      path = NA_character_
    ))
  }
  set_dir <- file.path(output_dir, set)
  dir.create(set_dir, recursive = TRUE)
  paths <- vapply(
    names(spiked),
    function(s) {
      path <- file.path(set_dir, .sw_spike_file_name(raw[[s]], s))
      .sw_spike_write_fcs(spiked[[s]], raw[[s]], path, keywords)
      path
    },
    character(1)
  )
  tibble::tibble(set = set, sample = names(spiked), path = unname(paths))
}

# Spiked share per population and per file, from one key.
.sw_spike_fraction <- function(key, n_events) {
  per_pop <- key |>
    dplyr::summarise(
      n = dplyr::n(),
      n_spiked = sum(.data$k > 0),
      .by = c("set", "sample", "pop")
    ) |>
    dplyr::mutate(fraction = .data$n_spiked / .data$n)
  per_file <- per_pop |>
    dplyr::summarise(
      file_spiked = sum(.data$n_spiked),
      .by = c("set", "sample")
    ) |>
    dplyr::mutate(
      file_events = unname(n_events[.data$sample]),
      file_fraction = .data$file_spiked / .data$file_events
    )
  dplyr::left_join(
    per_pop,
    per_file,
    by = c("set", "sample"),
    relationship = "many-to-one"
  )
}

# Package and R versions for the provenance block.
.sw_spike_provenance <- function(...) {
  c(
    list(
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      package_version = as.character(utils::packageVersion("SpectraWeaveR")),
      r_version = R.version.string
    ),
    list(...)
  )
}

# ---- Build -------------------------------------------------------------------

#' Build spiked raw files for a spike-in recovery test
#'
#' Adds `k = cofactor * sinh(level)` times a dye's spectrum to the raw detector
#' values of chosen cells. Each marker in `dyes` gets its own set of files; all
#' sets spike the same cells at the same levels, so the sets differ only in the
#' dye.
#'
#' @param raw Named raw inputs, one per sample: a character vector of FCS paths,
#'   or a list of matrices (events by channels, with column names). Names are
#'   sample ids.
#' @param populations Named list over samples; each element a named list of
#'   integer raw row indices per population, for example
#'   `list(A3 = list(CD4T = c(...), CD8T = c(...)))`. Populations must not share
#'   rows. [sw_unmix_spike_rows()] builds this from a GatingSet.
#' @param spectra Matrix of spectra, fluorophores in rows and detectors in
#'   columns, with names. Rows used for spiking must have max 1.
#' @param dyes Named character: marker name = row name of `spectra`. One set is
#'   built per marker.
#' @param levels Numeric spike levels in arcsinh units; must contain 0.
#' @param cofactor Positive arcsinh cofactor; sets `k`.
#' @param spike_fraction Share of each population's cells that get a level above
#'   0 (`round(spike_fraction * n)` cells, spread evenly over the non-zero
#'   levels). `NULL` when `assignment` is given.
#' @param seed Seed for the level draw. `NULL` when `assignment` is given.
#' @param key_cols Columns that identify a row in raw and unmixed files (for
#'   example time and scatter). They must not be spectrum detectors and must be
#'   unique per file.
#' @param output_dir Folder for the spiked FCS files; must not exist. `NULL`
#'   keeps the spiked matrices in memory (matrix input only).
#' @param assignment Optional data frame with columns `sample`, `raw_row`, `pop`
#'   and `level` that replaces the draw (for example a saved spike key).
#' @param write_base If `TRUE` and `output_dir` is given, also write an unspiked
#'   copy of each raw file as set `"base"`, so base and spiked files pass through
#'   the same writer.
#'
#' @return An `sw_spike_design`: a list with `design_id`, `type` (`"single"`),
#'   `key` (one row per set and population cell: `set`, `sample`, `raw_row`,
#'   `pop`, `marker`, `dye`, `level`, `k`), `sets`, `files` (`set`, `sample`,
#'   `path`), `data` (in-memory spiked matrices or `NULL`), `fraction` (realised
#'   spiked share per population and per file), the inputs and `provenance`.
#'   With `output_dir`, `spike_key.csv`, `fraction.csv` and `design.rds` are
#'   written there.
#'
#' @examples
#' set.seed(1)
#' det <- paste0("D", 1:4)
#' spectra <- rbind(A = c(1, 0.4, 0.1, 0), B = c(0, 0.2, 1, 0.5))
#' colnames(spectra) <- det
#' raw <- lapply(c(s1 = 1, s2 = 2), function(i) {
#'   cbind(Time = seq_len(50) + i / 10, matrix(abs(rnorm(200, 100, 10)), 50,
#'     dimnames = list(NULL, det)))
#' })
#' design <- sw_unmix_spike_build(
#'   raw = raw,
#'   populations = list(s1 = list(T = 1:40), s2 = list(T = 1:40)),
#'   spectra = spectra,
#'   dyes = c(markerA = "A"),
#'   levels = c(0, 1, 2),
#'   cofactor = 6000,
#'   spike_fraction = 0.5,
#'   seed = 42,
#'   key_cols = "Time",
#'   output_dir = NULL
#' )
#' design
#'
#' @export
sw_unmix_spike_build <- function(
  raw,
  populations,
  spectra,
  dyes,
  levels,
  cofactor,
  spike_fraction,
  seed,
  key_cols,
  output_dir,
  assignment = NULL,
  write_base = TRUE
) {
  .sw_stop_missing(c(
    raw = missing(raw),
    populations = missing(populations),
    spectra = missing(spectra),
    dyes = missing(dyes),
    levels = missing(levels),
    cofactor = missing(cofactor),
    spike_fraction = missing(spike_fraction),
    seed = missing(seed),
    key_cols = missing(key_cols),
    output_dir = missing(output_dir)
  ))
  if (is.character(raw)) {
    raw <- as.list(raw)
  }
  if (!is.list(raw) || is.null(names(raw)) || anyDuplicated(names(raw))) {
    stop("'raw' must be named by sample, with unique names.", call. = FALSE)
  }
  file_input <- all(vapply(raw, is.character, logical(1)))
  if (!file_input && !all(vapply(raw, is.matrix, logical(1)))) {
    stop("'raw' must be all FCS paths or all matrices.", call. = FALSE)
  }
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L) {
      stop("'output_dir' must be one path or NULL.", call. = FALSE)
    }
    if (dir.exists(output_dir)) {
      stop(
        sprintf("'output_dir' already exists: %s", output_dir),
        call. = FALSE
      )
    }
  } else if (file_input) {
    stop("'output_dir' may be NULL only for matrix input.", call. = FALSE)
  }

  # Spectrum, dyes, levels and cofactor -------------------------------------
  if (
    !is.matrix(spectra) ||
      is.null(rownames(spectra)) ||
      is.null(colnames(spectra))
  ) {
    stop("'spectra' must be a matrix with row and column names.", call. = FALSE)
  }
  if (
    !is.character(dyes) || is.null(names(dyes)) || anyDuplicated(names(dyes))
  ) {
    stop(
      "'dyes' must be a named character vector (marker = spectrum row).",
      call. = FALSE
    )
  }
  unknown <- setdiff(dyes, rownames(spectra))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Dye(s) not in spectra rows: %s.",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  .sw_spike_check_spectrum_rows(spectra, unname(dyes))
  if (
    !is.numeric(levels) || anyNA(levels) || any(levels < 0) || !0 %in% levels
  ) {
    stop("'levels' must be non-negative numbers and contain 0.", call. = FALSE)
  }
  levels <- sort(unique(levels))
  if (length(levels) < 2L) {
    stop("'levels' needs at least one level above 0.", call. = FALSE)
  }
  if (!is.numeric(cofactor) || length(cofactor) != 1L || !(cofactor > 0)) {
    stop("'cofactor' must be one positive number.", call. = FALSE)
  }
  if (!is.character(key_cols) || length(key_cols) == 0L) {
    stop("'key_cols' must name at least one column.", call. = FALSE)
  }
  if (length(intersect(key_cols, colnames(spectra))) > 0L) {
    stop("'key_cols' must not include spectrum detectors.", call. = FALSE)
  }

  # Populations -------------------------------------------------------------
  if (!is.list(populations) || !setequal(names(populations), names(raw))) {
    stop(
      "'populations' must be named by the same samples as 'raw'.",
      call. = FALSE
    )
  }
  for (s in names(populations)) {
    pops <- populations[[s]]
    if (!is.list(pops) || is.null(names(pops)) || any(!nzchar(names(pops)))) {
      stop(
        sprintf("populations[['%s']] must be a named list.", s),
        call. = FALSE
      )
    }
    all_rows <- unlist(pops, use.names = FALSE)
    if (anyDuplicated(all_rows)) {
      stop(
        sprintf(
          "Populations of sample '%s' share %d row(s); they must be disjoint.",
          s,
          sum(duplicated(all_rows))
        ),
        call. = FALSE
      )
    }
  }

  # Assignment: drawn from spike_fraction and seed, or given ---------------
  if (is.null(assignment)) {
    if (
      !is.numeric(spike_fraction) ||
        length(spike_fraction) != 1L ||
        !(spike_fraction > 0 && spike_fraction <= 1)
    ) {
      stop("'spike_fraction' must be one number in (0, 1].", call. = FALSE)
    }
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop(
        "'seed' must be one number when 'assignment' is NULL.",
        call. = FALSE
      )
    }
    assign_tbl <- .sw_spike_draw(populations, levels, spike_fraction, seed)
  } else {
    if (!is.null(spike_fraction) || !is.null(seed)) {
      stop(
        "With 'assignment', set 'spike_fraction' and 'seed' to NULL.",
        call. = FALSE
      )
    }
    need <- c("sample", "raw_row", "pop", "level")
    if (!is.data.frame(assignment) || !all(need %in% names(assignment))) {
      stop(
        "'assignment' needs columns sample, raw_row, pop, level.",
        call. = FALSE
      )
    }
    assign_tbl <- tibble::as_tibble(assignment[need])
    assign_tbl$raw_row <- as.integer(assign_tbl$raw_row)
    if (!all(assign_tbl$level %in% levels)) {
      stop("'assignment' has levels outside 'levels'.", call. = FALSE)
    }
    pop_tbl <- dplyr::bind_rows(lapply(names(populations), function(s) {
      dplyr::bind_rows(lapply(names(populations[[s]]), function(p) {
        tibble::tibble(
          sample = s,
          pop = p,
          raw_row = as.integer(populations[[s]][[p]])
        )
      }))
    }))
    if (
      nrow(dplyr::anti_join(
        assign_tbl,
        pop_tbl,
        by = c("sample", "pop", "raw_row")
      )) >
        0L ||
        nrow(dplyr::anti_join(
          pop_tbl,
          assign_tbl,
          by = c("sample", "pop", "raw_row")
        )) >
          0L
    ) {
      stop("'assignment' rows must match 'populations' exactly.", call. = FALSE)
    }
  }

  # Read raw, check keys and rows -------------------------------------------
  raw_x <- lapply(raw, .sw_spike_read)
  for (s in names(raw_x)) {
    x <- raw_x[[s]]
    missing_cols <- setdiff(c(key_cols, colnames(spectra)), colnames(x))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "Sample '%s' lacks column(s): %s.",
          s,
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    .sw_spike_check_unique_keys(x, key_cols, sprintf("Sample '%s'", s))
    if (
      any(unlist(populations[[s]]) > nrow(x)) ||
        any(unlist(populations[[s]]) < 1)
    ) {
      stop(
        sprintf("Population rows out of range in sample '%s'.", s),
        call. = FALSE
      )
    }
  }
  n_events <- vapply(raw_x, nrow, integer(1))

  design_id <- if (is.null(output_dir)) "in-memory" else basename(output_dir)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # One set per marker: same cells, same levels, that marker's dye ----------
  detectors <- colnames(spectra)
  keys <- list()
  files <- list()
  data <- if (is.null(output_dir)) list() else NULL
  for (m in names(dyes)) {
    key_m <- assign_tbl |>
      dplyr::mutate(
        set = m,
        marker = m,
        dye = unname(dyes[[m]]),
        k = cofactor * sinh(.data$level)
      )
    spiked <- lapply(stats::setNames(names(raw_x), names(raw_x)), function(s) {
      x <- raw_x[[s]]
      ks <- key_m[key_m$sample == s & key_m$k > 0, ]
      x[ks$raw_row, detectors] <- x[ks$raw_row, detectors] +
        outer(ks$k, spectra[dyes[[m]], detectors])
      x
    })
    keywords <- list(
      `$SW_SPIKE_DESIGN` = design_id,
      `$SW_SPIKE_DYE` = unname(dyes[[m]]),
      `$SW_SPIKE_LEVELS` = paste(levels, collapse = ",")
    )
    files[[m]] <- .sw_spike_emit_set(spiked, raw, m, output_dir, keywords)
    if (is.null(output_dir)) {
      data[[m]] <- spiked
    }
    keys[[m]] <- key_m
  }
  if (!is.null(output_dir) && isTRUE(write_base)) {
    base_kw <- list(
      `$SW_SPIKE_DESIGN` = design_id,
      `$SW_SPIKE_DYE` = "none",
      `$SW_SPIKE_LEVELS` = "0"
    )
    files[["base"]] <- .sw_spike_emit_set(
      raw_x,
      raw,
      "base",
      output_dir,
      base_kw
    )
  }
  key <- dplyr::bind_rows(keys) |>
    dplyr::select(
      "set",
      "sample",
      "raw_row",
      "pop",
      "marker",
      "dye",
      "level",
      "k"
    )
  fraction <- .sw_spike_fraction(key, n_events)

  design <- structure(
    list(
      design_id = design_id,
      type = "single",
      key = key,
      sets = tibble::tibble(
        set = names(dyes),
        marker = names(dyes),
        dye = unname(dyes)
      ),
      files = dplyr::bind_rows(files),
      data = data,
      fraction = fraction,
      raw = raw,
      n_events = n_events,
      populations = populations,
      spectra = spectra,
      dyes = dyes,
      levels = levels,
      cofactor = cofactor,
      spike_fraction = spike_fraction,
      seed = seed,
      key_cols = key_cols,
      output_dir = output_dir,
      provenance = .sw_spike_provenance(assignment_given = !is.null(assignment))
    ),
    class = "sw_spike_design"
  )
  if (!is.null(output_dir)) {
    utils::write.csv(
      key,
      file.path(output_dir, "spike_key.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      fraction,
      file.path(output_dir, "fraction.csv"),
      row.names = FALSE
    )
    saveRDS(design, file.path(output_dir, "design.rds"))
  }
  design
}

#' @export
print.sw_spike_design <- function(x, ...) {
  cat(
    "<sw_spike_design>",
    x$design_id,
    sprintf("(%s)", x$type),
    "\n  sets:",
    paste(x$sets$set, collapse = ", "),
    "\n  samples:",
    length(x$n_events),
    "| levels:",
    paste(x$levels, collapse = ", "),
    "| cofactor:",
    x$cofactor,
    "\n  spiked share per file:",
    sprintf(
      "%.3f to %.3f",
      min(x$fraction$file_fraction),
      max(x$fraction$file_fraction)
    ),
    "\n"
  )
  invisible(x)
}

# ---- Co-spike ----------------------------------------------------------------

#' Co-spike a second dye into the cells of one level
#'
#' Rebuilds from the raw files of a single-spike design. Only the cells the
#' design gave `level1` are spiked, with the first dye at `level1` and a second
#' dye at each of `levels2`; every other cell stays unspiked. This is the
#' construction of benchmark 05 script 88d, so the spiked share per file is
#' lower than in the single-spike design. `levels2 = 0` gives a first-dye-only
#' set at that lower share, which checks the effect of the spiked share.
#'
#' @param design An `sw_spike_design` of type `"single"` from
#'   [sw_unmix_spike_build()].
#' @param marker Marker of `design` whose dye is the first dye.
#' @param dye2 Named character of length 1: second marker name = row of the
#'   design's `spectra`. The row must have max 1.
#' @param levels2 Numeric levels of the second dye, in arcsinh units; may
#'   contain 0.
#' @param level1 The design level whose cells are co-spiked.
#' @param output_dir Folder for the co-spike files; must not exist. `NULL`
#'   only when the design holds matrices.
#'
#' @return An `sw_spike_design` of type `"cospike"`. Sets are named
#'   `cospike_<marker2>_<level2>`; `key` adds `marker2`, `dye2`, `level2` and
#'   `k2` (0 for unspiked cells).
#'
#' @export
sw_unmix_spike_cospike <- function(
  design,
  marker,
  dye2,
  levels2,
  level1,
  output_dir
) {
  .sw_stop_missing(c(
    design = missing(design),
    marker = missing(marker),
    dye2 = missing(dye2),
    levels2 = missing(levels2),
    level1 = missing(level1),
    output_dir = missing(output_dir)
  ))
  if (!inherits(design, "sw_spike_design") || design$type != "single") {
    stop("'design' must be a single-spike sw_spike_design.", call. = FALSE)
  }
  if (
    !is.character(marker) ||
      length(marker) != 1L ||
      !marker %in% names(design$dyes)
  ) {
    stop("'marker' must be one marker of the design.", call. = FALSE)
  }
  if (
    !is.character(dye2) ||
      length(dye2) != 1L ||
      is.null(names(dye2)) ||
      !dye2 %in% rownames(design$spectra)
  ) {
    stop(
      "'dye2' must be one named spectrum row (marker2 = row).",
      call. = FALSE
    )
  }
  .sw_spike_check_spectrum_rows(design$spectra, unname(dye2))
  if (!is.numeric(levels2) || anyNA(levels2) || any(levels2 < 0)) {
    stop("'levels2' must be non-negative numbers.", call. = FALSE)
  }
  if (
    !is.numeric(level1) ||
      length(level1) != 1L ||
      !(level1 %in% design$levels) ||
      level1 == 0
  ) {
    stop("'level1' must be one non-zero level of the design.", call. = FALSE)
  }
  file_input <- all(vapply(design$raw, is.character, logical(1)))
  if (is.null(output_dir) && file_input) {
    stop("'output_dir' may be NULL only for matrix input.", call. = FALSE)
  }
  if (!is.null(output_dir)) {
    if (dir.exists(output_dir)) {
      stop(
        sprintf("'output_dir' already exists: %s", output_dir),
        call. = FALSE
      )
    }
    dir.create(output_dir, recursive = TRUE)
  }

  # Cells of level1 get both dyes; the rest of the population stays unspiked --
  dye1 <- design$dyes[[marker]]
  base_key <- design$key[design$key$set == marker, ]
  detectors <- colnames(design$spectra)
  raw_x <- lapply(design$raw, .sw_spike_read)
  k1 <- design$cofactor * sinh(level1)
  marker2 <- names(dye2)
  design_id <- if (is.null(output_dir)) "in-memory" else basename(output_dir)
  keys <- list()
  files <- list()
  data <- if (is.null(output_dir)) list() else NULL
  for (l2 in sort(unique(levels2))) {
    set_name <- sprintf("cospike_%s_%s", marker2, format(l2))
    key_s <- base_key |>
      dplyr::mutate(
        set = set_name,
        level = dplyr::if_else(.data$level == level1, level1, 0),
        k = dplyr::if_else(.data$level > 0, k1, 0),
        marker2 = marker2,
        dye2 = unname(dye2),
        level2 = dplyr::if_else(.data$level > 0, l2, 0),
        k2 = design$cofactor * sinh(.data$level2)
      )
    spiked <- lapply(stats::setNames(names(raw_x), names(raw_x)), function(s) {
      x <- raw_x[[s]]
      ks <- key_s[key_s$sample == s & key_s$k > 0, ]
      x[ks$raw_row, detectors] <- x[ks$raw_row, detectors] +
        outer(ks$k, design$spectra[dye1, detectors]) +
        outer(ks$k2, design$spectra[dye2, detectors])
      x
    })
    keywords <- list(
      `$SW_SPIKE_DESIGN` = design_id,
      `$SW_SPIKE_DYE` = paste(dye1, unname(dye2), sep = "+"),
      `$SW_SPIKE_LEVELS` = sprintf("%s,%s", format(level1), format(l2))
    )
    files[[set_name]] <- .sw_spike_emit_set(
      spiked,
      design$raw,
      set_name,
      output_dir,
      keywords
    )
    if (is.null(output_dir)) {
      data[[set_name]] <- spiked
    }
    keys[[set_name]] <- key_s
  }
  key <- dplyr::bind_rows(keys)
  sets <- tibble::tibble(
    set = names(keys),
    marker = marker,
    dye = unname(dye1),
    marker2 = marker2,
    dye2 = unname(dye2),
    level2 = sort(unique(levels2))
  )
  fraction <- .sw_spike_fraction(key, design$n_events)
  out <- design
  out$design_id <- design_id
  out$type <- "cospike"
  out$key <- key
  out$sets <- sets
  out$files <- dplyr::bind_rows(files)
  out$data <- data
  out$fraction <- fraction
  out$level1 <- level1
  out$output_dir <- output_dir
  out$provenance <- .sw_spike_provenance(from_design = design$design_id)
  if (!is.null(output_dir)) {
    utils::write.csv(
      key,
      file.path(output_dir, "spike_key.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      fraction,
      file.path(output_dir, "fraction.csv"),
      row.names = FALSE
    )
    saveRDS(out, file.path(output_dir, "design.rds"))
  }
  out
}

# ---- Score -------------------------------------------------------------------

# Unmixed matrix of one row of the unmixed table: a path or a matrix list column.
.sw_spike_unmixed_x <- function(row) {
  if ("data" %in% names(row) && !is.null(row$data[[1]])) {
    return(row$data[[1]])
  }
  .sw_spike_read(row$path[[1]])
}

# Fluorescence channels shared by two unmixed matrices, minus the key columns.
.sw_spike_fluor_channels <- function(a, b, key_cols) {
  shared <- intersect(colnames(a), colnames(b))
  shared <- setdiff(shared, key_cols)
  shared[!grepl(.sw_spike_nonfluor_pattern, shared)]
}

# Shift of unspiked rows per channel in arcsinh units, one row per channel.
.sw_spike_baseline_rows <- function(
  x_spk,
  x_base,
  rows,
  channels,
  cofactor,
  tol
) {
  if (length(rows) == 0L) {
    return(tibble::tibble(
      channel = channels,
      n = 0L,
      max_abs_diff = NA_real_,
      share_moved = NA_real_
    ))
  }
  d <- abs(
    asinh(x_spk[rows, channels, drop = FALSE] / cofactor) -
      asinh(x_base[rows, channels, drop = FALSE] / cofactor)
  )
  tibble::tibble(
    channel = channels,
    n = length(rows),
    max_abs_diff = unname(apply(d, 2, max)),
    share_moved = unname(colMeans(d > tol))
  )
}

# Cosine between two named numeric vectors on their shared names.
.sw_spike_cosine <- function(a, b) {
  shared <- intersect(names(a), names(b))
  if (length(shared) == 0L) {
    return(NA_real_)
  }
  a <- a[shared]
  b <- b[shared]
  sum(a * b) / sqrt(sum(a^2) * sum(b^2))
}

#' Score spike-in recovery from unmixed files
#'
#' For each setup, sample and set, compares the unmixed spiked file with the
#' unmixed base file cell by cell. Spiked cells give recovery; unspiked cells
#' give the independence check.
#'
#' @param design An `sw_spike_design` from [sw_unmix_spike_build()] or
#'   [sw_unmix_spike_cospike()].
#' @param unmixed Data frame with columns `setup`, `sample`, `set` and either
#'   `path` (unmixed FCS) or a list column `data` (unmixed matrices). Each setup
#'   needs every design set plus `"base"` for every sample. Unmixed rows must be
#'   in raw row order, as AutoSpectral writes them; the key columns are checked.
#' @param loss_cut A cell counts as lost when its recovery is below this value.
#' @param baseline_tol Shift in arcsinh units above which an unspiked cell
#'   counts as moved. A file with any moved cell in any fluorescence channel is
#'   `coupled`.
#' @param results Optional named list (setup = `sw_unmix_result`). Gives the
#'   cosine between each spike dye's row and the setup's spectrum, and a note
#'   when a setup records `n.af.passes` above 1.
#' @param channel_suffix Suffix that turns a dye name into its unmixed channel
#'   name (AutoSpectral writes `<dye>-A`).
#'
#' @return An `sw_spike_score`: a list with
#' \describe{
#'   \item{cells}{one row per spiked cell: `setup`, `set`, `sample`, `pop`,
#'     `raw_row`, `marker`, `level`, `k`, `base_value`, `recovered`, `recovery`
#'     (and `level2`, `k2`, `recovered_2`, `recovery_2` for co-spikes).}
#'   \item{summary}{per `setup`, `set`, `marker`, `pop`, `sample`, `level`:
#'     `n`, `share_lost`, `median_recovery`, `q05`, `q25`, `q50`, `q75`, `q95`,
#'     `file_status` (and `median_recovery_2` for co-spikes).}
#'   \item{by_setup}{donor medians and ranges of `share_lost` and
#'     `median_recovery` per `setup`, `set`, `marker`, `pop`, `level`, with the
#'     number of donors and of coupled files.}
#'   \item{baseline_check}{per `setup`, `sample`, `set`, `group` (`"level0"`:
#'     unspiked cells of the populations; `"untouched"`: cells outside them)
#'     and `channel`: `n`, `max_abs_diff`, `share_moved`.}
#'   \item{file_status}{per `setup`, `sample`, `set`: `independent` or
#'     `coupled`.}
#'   \item{spike_cosine}{per `setup` and spike dye, when `results` is given.}
#'   \item{notes}{character vector.}
#' }
#'
#' @export
sw_unmix_spike_score <- function(
  design,
  unmixed,
  loss_cut,
  baseline_tol,
  results = NULL,
  channel_suffix = "-A"
) {
  .sw_stop_missing(c(
    design = missing(design),
    unmixed = missing(unmixed),
    loss_cut = missing(loss_cut),
    baseline_tol = missing(baseline_tol)
  ))
  if (!inherits(design, "sw_spike_design")) {
    stop("'design' must be an sw_spike_design.", call. = FALSE)
  }
  if (!is.numeric(loss_cut) || length(loss_cut) != 1L || is.na(loss_cut)) {
    stop("'loss_cut' must be one number.", call. = FALSE)
  }
  if (
    !is.numeric(baseline_tol) ||
      length(baseline_tol) != 1L ||
      !(baseline_tol >= 0)
  ) {
    stop("'baseline_tol' must be one non-negative number.", call. = FALSE)
  }
  if (
    !is.data.frame(unmixed) ||
      !all(c("setup", "sample", "set") %in% names(unmixed)) ||
      !any(c("path", "data") %in% names(unmixed))
  ) {
    stop(
      "'unmixed' needs columns setup, sample, set and path or data.",
      call. = FALSE
    )
  }
  unmixed <- tibble::as_tibble(unmixed)

  # Every setup needs every set and the base for every sample, once ----------
  samples <- names(design$n_events)
  wanted <- expand.grid(
    setup = unique(unmixed$setup),
    sample = samples,
    set = c(design$sets$set, "base"),
    stringsAsFactors = FALSE
  )
  have <- paste(unmixed$setup, unmixed$sample, unmixed$set)
  if (anyDuplicated(have)) {
    stop("'unmixed' has duplicated setup, sample and set rows.", call. = FALSE)
  }
  absent <- setdiff(paste(wanted$setup, wanted$sample, wanted$set), have)
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "'unmixed' lacks %d file(s), e.g. %s.",
        length(absent),
        paste(utils::head(absent, 3), collapse = "; ")
      ),
      call. = FALSE
    )
  }

  key_cols <- design$key_cols
  cof <- design$cofactor
  is_cospike <- design$type == "cospike"
  cells <- list()
  baseline <- list()
  for (setup in unique(unmixed$setup)) {
    for (s in samples) {
      base_row <- unmixed[
        unmixed$setup == setup & unmixed$sample == s & unmixed$set == "base",
      ]
      x_base <- .sw_spike_unmixed_x(base_row)
      if (nrow(x_base) != design$n_events[[s]]) {
        stop(
          sprintf(
            "Setup '%s', sample '%s': base unmix has %d rows, raw has %d.",
            setup,
            s,
            nrow(x_base),
            design$n_events[[s]]
          ),
          call. = FALSE
        )
      }
      pop_rows <- unique(design$key$raw_row[design$key$sample == s])
      untouched <- setdiff(seq_len(nrow(x_base)), pop_rows)
      for (set in design$sets$set) {
        spk_row <- unmixed[
          unmixed$setup == setup & unmixed$sample == s & unmixed$set == set,
        ]
        x_spk <- .sw_spike_unmixed_x(spk_row)
        if (
          !isTRUE(all.equal(
            x_spk[, key_cols, drop = FALSE],
            x_base[, key_cols, drop = FALSE],
            check.attributes = FALSE
          ))
        ) {
          stop(
            sprintf(
              "Setup '%s', sample '%s', set '%s': key columns differ from the base unmix.",
              setup,
              s,
              set
            ),
            call. = FALSE
          )
        }
        key_s <- design$key[design$key$sample == s & design$key$set == set, ]
        dye <- design$sets$dye[design$sets$set == set]
        channel <- paste0(dye, channel_suffix)
        if (!channel %in% colnames(x_spk)) {
          stop(
            sprintf("Unmixed channel '%s' not found.", channel),
            call. = FALSE
          )
        }

        # Recovery of the spiked cells ----------------------------------------
        sp <- key_s[key_s$k > 0, ]
        out <- sp |>
          dplyr::mutate(
            setup = .env$setup,
            base_value = x_base[sp$raw_row, channel],
            recovered = x_spk[sp$raw_row, channel] -
              x_base[sp$raw_row, channel],
            recovery = .data$recovered / .data$k
          )
        if (is_cospike) {
          channel2 <- paste0(
            design$sets$dye2[design$sets$set == set],
            channel_suffix
          )
          if (!channel2 %in% colnames(x_spk)) {
            stop(
              sprintf("Unmixed channel '%s' not found.", channel2),
              call. = FALSE
            )
          }
          out <- out |>
            dplyr::mutate(
              recovered_2 = x_spk[sp$raw_row, channel2] -
                x_base[sp$raw_row, channel2],
              recovery_2 = dplyr::if_else(
                .data$k2 > 0,
                .data$recovered_2 / .data$k2,
                NA_real_
              )
            )
        }
        cells[[length(cells) + 1L]] <- out

        # Independence: unspiked cells of the populations and outside them ----
        channels <- .sw_spike_fluor_channels(x_spk, x_base, key_cols)
        unspiked <- key_s$k == 0
        if (is_cospike) {
          unspiked <- unspiked & key_s$k2 == 0
        }
        level0 <- key_s$raw_row[unspiked]
        for (group in c("level0", "untouched")) {
          rows <- if (group == "level0") level0 else untouched
          baseline[[length(baseline) + 1L]] <- .sw_spike_baseline_rows(
            x_spk,
            x_base,
            rows,
            channels,
            cof,
            baseline_tol
          ) |>
            dplyr::mutate(
              setup = .env$setup,
              sample = .env$s,
              set = .env$set,
              group = .env$group,
              .before = 1
            )
        }
      }
    }
  }
  cells <- dplyr::bind_rows(cells) |>
    dplyr::relocate("setup", "set", "sample", "pop", "raw_row")
  baseline_check <- dplyr::bind_rows(baseline)
  file_status <- baseline_check |>
    dplyr::summarise(
      max_shift = suppressWarnings(max(.data$max_abs_diff, na.rm = TRUE)),
      file_status = dplyr::if_else(
        any(.data$share_moved > 0, na.rm = TRUE),
        "coupled",
        "independent"
      ),
      .by = c("setup", "sample", "set")
    )

  # Summaries: per donor, then across donors --------------------------------
  summary_tbl <- cells |>
    dplyr::summarise(
      n = dplyr::n(),
      share_lost = mean(.data$recovery < loss_cut),
      median_recovery = stats::median(.data$recovery),
      q05 = unname(stats::quantile(.data$recovery, 0.05)),
      q25 = unname(stats::quantile(.data$recovery, 0.25)),
      q50 = unname(stats::quantile(.data$recovery, 0.50)),
      q75 = unname(stats::quantile(.data$recovery, 0.75)),
      q95 = unname(stats::quantile(.data$recovery, 0.95)),
      median_recovery_2 = if (is_cospike) {
        stats::median(.data$recovery_2, na.rm = TRUE)
      } else {
        NA_real_
      },
      .by = c("setup", "set", "marker", "pop", "sample", "level")
    ) |>
    dplyr::left_join(
      file_status[c("setup", "sample", "set", "file_status")],
      by = c("setup", "sample", "set"),
      relationship = "many-to-one"
    )
  if (!is_cospike) {
    summary_tbl$median_recovery_2 <- NULL
  }
  by_setup <- summary_tbl |>
    dplyr::summarise(
      n_donors = dplyr::n(),
      n_coupled = sum(.data$file_status == "coupled"),
      share_lost_median = stats::median(.data$share_lost),
      share_lost_min = min(.data$share_lost),
      share_lost_max = max(.data$share_lost),
      recovery_median = stats::median(.data$median_recovery),
      recovery_min = min(.data$median_recovery),
      recovery_max = max(.data$median_recovery),
      .by = c("setup", "set", "marker", "pop", "level")
    )

  # Optional: spike-to-setup cosine and pass notes --------------------------
  spike_cosine <- NULL
  notes <- character()
  if (!is.null(results)) {
    if (!is.list(results) || is.null(names(results))) {
      stop("'results' must be a named list of sw_unmix_result.", call. = FALSE)
    }
    spike_dyes <- unique(c(design$sets$dye, if (is_cospike) design$sets$dye2))
    spike_cosine <- dplyr::bind_rows(lapply(names(results), function(setup) {
      res <- results[[setup]]
      dplyr::bind_rows(lapply(spike_dyes, function(d) {
        cos <- if (d %in% rownames(res$spectra)) {
          .sw_spike_cosine(design$spectra[d, ], res$spectra[d, ])
        } else {
          NA_real_
        }
        tibble::tibble(setup = setup, dye = d, cosine = cos)
      }))
    }))
    for (setup in names(results)) {
      af_passes <- results[[setup]]$provenance$args[["n.af.passes"]]
      if (!is.null(af_passes) && isTRUE(af_passes > 1)) {
        notes <- c(
          notes,
          sprintf(
            "Setup '%s' records n.af.passes = %s; extra AF passes select cells per file (refine.af.quantile), so unmixing may not be per cell. Read recovery with file_status.",
            setup,
            af_passes
          )
        )
      }
    }
  }
  n_coupled <- sum(file_status$file_status == "coupled")
  if (n_coupled > 0L) {
    notes <- c(
      notes,
      sprintf(
        "%d of %d unmixed file(s) are coupled: unspiked cells moved by more than %s arcsinh.",
        n_coupled,
        nrow(file_status),
        format(baseline_tol)
      )
    )
  }

  structure(
    list(
      cells = cells,
      summary = summary_tbl,
      by_setup = by_setup,
      baseline_check = baseline_check,
      file_status = file_status,
      spike_cosine = spike_cosine,
      notes = notes,
      loss_cut = loss_cut,
      baseline_tol = baseline_tol,
      design_id = design$design_id,
      type = design$type,
      provenance = .sw_spike_provenance(design = design$design_id)
    ),
    class = "sw_spike_score"
  )
}

#' @export
print.sw_spike_score <- function(x, ...) {
  status <- table(x$file_status$file_status)
  cat(
    "<sw_spike_score>",
    x$design_id,
    sprintf("(%s)", x$type),
    "\n  files:",
    paste(names(status), status, sep = " = ", collapse = ", "),
    "\n  setups:",
    paste(unique(x$summary$setup), collapse = ", "),
    "| loss cut:",
    x$loss_cut,
    "| baseline tol:",
    x$baseline_tol,
    "arcsinh\n"
  )
  for (n in x$notes) {
    cat("  note:", n, "\n")
  }
  invisible(x)
}

# dplyr data-mask pronoun used throughout this file.
utils::globalVariables(c(".data", ".env"))
