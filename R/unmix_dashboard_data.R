# Data model for the unmixing dashboard (sw_unmix_dashboard, sw_unmix_app).
# Reads the files a finished run folder holds and returns plain lists and
# data frames. Nothing here draws. Event-level tables are binned here, so the
# HTML file stays small.

# Status values and the run-folder files read (see sw_unmix_report and
# sw_unmix_save_figures for the writers).
.sw_dash_report_files <- c("checks.csv", "spectra.csv", "provenance.csv")
.sw_dash_panel_files <- c(
  spectra = "spectra-overlay.csv",
  similarity = "similarity-matrix.csv",
  distributions = "unmixed-distributions.csv",
  biplot = "biplot.csv"
)

# Argument checks shared by sw_unmix_dashboard() and sw_unmix_app(): the run
# folders and the two display settings that have no default.
.sw_dash_check_settings <- function(runs, sigma, bins) {
  if (!is.character(runs) || length(runs) < 1L || anyNA(runs)) {
    stop("`runs` must be a character vector of run folders.", call. = FALSE)
  }
  if (
    !is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0
  ) {
    stop("`pseudo_log_sigma` must be one positive number.", call. = FALSE)
  }
  if (
    !is.numeric(bins) ||
      length(bins) != 1L ||
      !is.finite(bins) ||
      bins < 10 ||
      bins != round(bins)
  ) {
    stop("`bins` must be a whole number of at least 10.", call. = FALSE)
  }
  invisible(TRUE)
}

# Run labels: the names of `runs`, or the folder names. Must be unique.
.sw_dash_run_labels <- function(runs) {
  labels <- if (is.null(names(runs))) basename(runs) else names(runs)
  if (anyDuplicated(labels) || any(!nzchar(labels))) {
    stop("Run labels must be unique and not empty.", call. = FALSE)
  }
  labels
}

# Folder of figure panel source data in a run: the last (by name) of
# `figures*/panels/source_data` that holds the spectra-overlay table. A run
# folder with `figures` and `figures_v2` gives `figures_v2`.
.sw_dash_find_figures <- function(run_dir, figures_dir = NULL) {
  if (!is.null(figures_dir)) {
    if (
      !is.character(figures_dir) ||
        length(figures_dir) != 1L ||
        is.na(figures_dir)
    ) {
      stop("`figures_dir` must be one folder name.", call. = FALSE)
    }
    cand <- figures_dir
  } else {
    cand <- sort(list.dirs(run_dir, recursive = FALSE, full.names = FALSE))
    cand <- cand[grepl("^figures", cand)]
  }
  ok <- cand[file.exists(file.path(
    run_dir,
    cand,
    "panels",
    "source_data",
    .sw_dash_panel_files[["spectra"]]
  ))]
  if (length(ok) == 0L) {
    return(NA_character_)
  }
  ok[length(ok)]
}

# Parse the limit text written by sw_unmix_diagnose, for example
# "warn if > 4; fail if > 6", into c(warn = 4, fail = 6).
.sw_dash_parse_limit <- function(text) {
  text <- unique(stats::na.omit(as.character(text)))
  if (length(text) != 1L) {
    return(c(warn = NA_real_, fail = NA_real_))
  }
  pick <- function(word) {
    m <- regmatches(
      text,
      regexec(paste0(word, " if > ([-0-9.eE+]+)"), text)
    )[[1L]]
    if (length(m) == 2L) as.numeric(m[[2L]]) else NA_real_
  }
  c(warn = pick("warn"), fail = pick("fail"))
}

# Read the named columns of a figure source table. `what` names each column
# of the file in order with its type ("character", "numeric") or NULL to skip
# it. scan() with fixed types reads a 1.2 million row table in about half the
# time read.csv() needs and gives the same vectors.
.sw_dash_read_columns <- function(path, what) {
  header <- strsplit(readLines(path, n = 1L, warn = FALSE), ",", fixed = TRUE)[[
    1L
  ]]
  header <- gsub("\"", "", header, fixed = TRUE)
  if (length(header) != length(what)) {
    stop(
      "Unexpected columns in ",
      path,
      ": ",
      paste(header, collapse = ", "),
      call. = FALSE
    )
  }
  types <- lapply(what, function(w) if (is.null(w)) NULL else vector(w, 0L))
  cols <- scan(
    path,
    what = types,
    sep = ",",
    skip = 1L,
    quote = "\"",
    na.strings = "NA",
    quiet = TRUE
  )
  keep <- !vapply(what, is.null, logical(1))
  names(cols)[keep] <- header[keep]
  cols[keep]
}

# One 1D histogram per method and channel on a common edge set, in the
# pseudo-log space.
.sw_dash_bin_distributions <- function(path, sigma, bins) {
  d <- .sw_dash_read_columns(
    path,
    list("character", NULL, NULL, "character", "numeric")
  )
  trans <- scales::transform_pseudo_log(sigma = sigma)
  tv_all <- trans$transform(d$value)
  channels <- unique(d$channel)
  methods <- unique(d$method)
  method_id <- match(d$method, methods)
  by_channel <- split(seq_along(tv_all), factor(d$channel, levels = channels))
  out <- lapply(channels, function(ch) {
    rows <- by_channel[[ch]]
    tv <- tv_all[rows]
    mid <- method_id[rows]
    r <- range(tv, na.rm = TRUE)
    edges <- seq(r[[1L]], r[[2L]], length.out = bins + 1L)
    idx <- findInterval(tv, edges, rightmost.closed = TRUE)
    counts <- lapply(seq_along(methods), function(j) {
      i <- idx[mid == j]
      tabulate(i[i >= 1L & i <= bins], nbins = bins)
    })
    names(counts) <- methods
    list(
      channel = ch,
      edges = edges,
      counts = counts,
      n = tabulate(mid, nbins = length(methods))
    )
  })
  names(out) <- channels
  list(channels = channels, methods = methods, panels = out)
}

# Sparse 2D histogram per method, as triples (x bin, y bin, count), 0-based.
.sw_dash_bin_biplot <- function(path, sigma, bins) {
  d <- .sw_dash_read_columns(
    path,
    list("character", NULL, NULL, "numeric", "numeric")
  )
  trans <- scales::transform_pseudo_log(sigma = sigma)
  tx <- trans$transform(d$x)
  ty <- trans$transform(d$y)
  rx <- range(tx, na.rm = TRUE)
  ry <- range(ty, na.rm = TRUE)
  ex <- seq(rx[[1L]], rx[[2L]], length.out = bins + 1L)
  ey <- seq(ry[[1L]], ry[[2L]], length.out = bins + 1L)
  ix_all <- findInterval(tx, ex, rightmost.closed = TRUE)
  iy_all <- findInterval(ty, ey, rightmost.closed = TRUE)
  methods <- unique(d$method)
  method_id <- match(d$method, methods)
  cells <- lapply(seq_along(methods), function(j) {
    keep <- method_id == j
    ix <- ix_all[keep]
    iy <- iy_all[keep]
    ok <- ix >= 1L & ix <= bins & iy >= 1L & iy <= bins
    key <- (iy[ok] - 1L) * bins + ix[ok]
    n <- tabulate(key, nbins = bins * bins)
    nz <- which(n > 0L)
    list(
      x = (nz - 1L) %% bins,
      y = (nz - 1L) %/% bins,
      n = n[nz],
      total = sum(keep)
    )
  })
  names(cells) <- methods
  list(methods = methods, ex = ex, ey = ey, cells = cells)
}

# Read one run folder. `sigma` and `bins` are the display settings for the
# event-level panels; they have no default (see sw_unmix_dashboard).
.sw_dash_read_run <- function(
  run_dir,
  sigma,
  bins,
  figures_dir = NULL,
  name = NULL,
  bin_events = TRUE
) {
  if (!is.character(run_dir) || length(run_dir) != 1L || !dir.exists(run_dir)) {
    stop("`run_dir` must be an existing run folder: ", run_dir, call. = FALSE)
  }
  rep_dir <- file.path(run_dir, "report")
  missing_files <- .sw_dash_report_files[
    !file.exists(file.path(rep_dir, .sw_dash_report_files))
  ]
  if (length(missing_files) > 0L) {
    stop(
      "Run folder lacks report/",
      paste(missing_files, collapse = ", report/"),
      ": ",
      run_dir,
      call. = FALSE
    )
  }
  rd <- \(f) utils::read.csv(file.path(rep_dir, f), stringsAsFactors = FALSE)
  checks <- rd("checks.csv")
  spectra_long <- rd("spectra.csv")
  provenance <- rd("provenance.csv")
  run_id <- provenance$value[
    provenance$section == "run" & provenance$field == "run_id"
  ]
  run_id <- if (length(run_id)) run_id[[1L]] else NA_character_

  # SIF rows and the run's saved limits. A run made before the max_sif check
  # existed has none: SIF is then computed from the saved spectra with the same
  # function sw_unmix_diagnose uses, and no verdict is attached to it.
  sif_rows <- checks[checks$metric == "max_sif", , drop = FALSE]
  limits <- .sw_dash_parse_limit(sif_rows$limit)
  sif_source <- "run"
  if (nrow(sif_rows) == 0L && nrow(spectra_long) > 0L) {
    sif_rows <- .sw_diag_sif(.sw_dash_spectra_matrix(list(
      spectra_long = spectra_long
    )))
    sif_rows$status <- ifelse(is.na(sif_rows$value), "not_run", "not_judged")
    sif_source <- "computed"
  }
  partner <- regmatches(
    sif_rows$reason,
    regexec("off-diagonal entry: (.*) \\(([-0-9.eE+]+)\\)\\.$", sif_rows$reason)
  )
  sif <- data.frame(
    fluorophore = sif_rows$item,
    value = sif_rows$value,
    status = sif_rows$status,
    partner = vapply(
      partner,
      \(m) if (length(m) == 3L) m[[2L]] else NA_character_,
      ""
    ),
    partner_value = vapply(
      partner,
      \(m) if (length(m) == 3L) as.numeric(m[[3L]]) else NA_real_,
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  fig <- .sw_dash_find_figures(run_dir, figures_dir)
  src <- if (is.na(fig)) {
    NA_character_
  } else {
    file.path(run_dir, fig, "panels", "source_data")
  }
  panel_path <- \(key) {
    if (is.na(src)) {
      return(NA_character_)
    }
    p <- file.path(src, .sw_dash_panel_files[[key]])
    if (file.exists(p)) p else NA_character_
  }

  spectra_overlay <- NULL
  if (!is.na(p <- panel_path("spectra"))) {
    spectra_overlay <- utils::read.csv(
      p,
      colClasses = c(
        "character",
        "character",
        "character",
        "numeric",
        "numeric",
        "numeric",
        "NULL"
      ),
      stringsAsFactors = FALSE
    )
  }
  similarity <- NULL
  if (!is.na(p <- panel_path("similarity"))) {
    similarity <- utils::read.csv(p, stringsAsFactors = FALSE)
  }
  distributions <- NULL
  if (bin_events && !is.na(p <- panel_path("distributions"))) {
    distributions <- .sw_dash_bin_distributions(p, sigma, bins)
  }
  biplot <- NULL
  if (bin_events && !is.na(p <- panel_path("biplot"))) {
    biplot <- .sw_dash_bin_biplot(p, sigma, bins)
  }

  compare_path <- c(
    if (!is.na(fig)) file.path(run_dir, fig, "spectra_compare.csv"),
    file.path(run_dir, "spectra_compare.csv")
  )
  compare_path <- compare_path[file.exists(compare_path)]
  compare <- if (length(compare_path)) {
    utils::read.csv(compare_path[[1L]], stringsAsFactors = FALSE)
  } else {
    NULL
  }

  list(
    name = if (is.null(name)) basename(normalizePath(run_dir)) else name,
    run_dir = run_dir,
    run_id = run_id,
    figures_dir = fig,
    provenance = provenance,
    checks = checks,
    status_counts = as.list(table(checks$status)),
    limits = limits,
    sif = sif,
    sif_source = sif_source,
    spectra_long = spectra_long,
    spectra_overlay = spectra_overlay,
    similarity = similarity,
    distributions = distributions,
    biplot = biplot,
    compare_spectracle = compare,
    spike = .sw_dash_read_spike(run_dir),
    events = .sw_dash_read_events(run_dir)
  )
}

# Event-metric and held-out tables that sw_unmix_report() writes. Pairs are
# sorted by |slope| (largest first) so the page shows the worst first. NULL
# when none of the four CSVs exists.
.sw_dash_event_files <- c(
  channels = "event_channels.csv",
  pairs = "event_pairs.csv",
  heldout = "heldout_scores.csv",
  repro = "heldout_repro.csv"
)
.sw_dash_read_events <- function(run_dir) {
  paths <- file.path(run_dir, "report", .sw_dash_event_files)
  names(paths) <- names(.sw_dash_event_files)
  if (!any(file.exists(paths))) {
    return(NULL)
  }
  tables <- lapply(paths, \(p) {
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
  })
  if (!is.null(tables$pairs)) {
    tables$pairs <- tables$pairs[order(-abs(tables$pairs$slope), na.last = TRUE), ]
  }
  tables
}

# Spike-in and marker-screen tables that sw_unmix_report() writes, and the F7
# and F8 PNGs that sw_unmix_save_figures() writes (file names starting with
# "spike-summary" or "spike-ecdf" in any `figures*/panels` folder). PNGs are
# embedded as data URIs so the page stays one file. NULL when none of the
# three CSVs exists.
.sw_dash_spike_files <- c(
  summary = "spike_summary.csv",
  baseline = "baseline_check.csv",
  screen = "marker_screen.csv"
)
.sw_dash_read_spike <- function(run_dir) {
  paths <- file.path(run_dir, "report", .sw_dash_spike_files)
  names(paths) <- names(.sw_dash_spike_files)
  if (!any(file.exists(paths))) {
    return(NULL)
  }
  tables <- lapply(paths, \(p) {
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
  })
  fig_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  fig_dirs <- file.path(
    fig_dirs[grepl("^figures", basename(fig_dirs))],
    "panels"
  )
  pngs <- sort(unlist(lapply(fig_dirs[dir.exists(fig_dirs)], \(d) {
    list.files(d, pattern = "^spike-(summary|ecdf).*\\.png$", full.names = TRUE)
  })))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    pngs <- character()
  }
  figures <- lapply(pngs, \(p) {
    bytes <- readBin(p, "raw", file.size(p))
    list(
      name = basename(p),
      src = paste0("data:image/png;base64,", jsonlite::base64_enc(bytes))
    )
  })
  c(tables, list(figures = figures))
}

# Spectra matrix (fluorophores x detectors) from the long table of a run.
.sw_dash_spectra_matrix <- function(run) {
  long <- run$spectra_long
  fl <- unique(long$fluorophore)
  dt <- unique(long$detector)
  idx <- match(
    paste(rep(fl, times = length(dt)), rep(dt, each = length(fl))),
    paste(long$fluorophore, long$detector)
  )
  matrix(long$value[idx], nrow = length(fl), dimnames = list(fl, dt))
}

# Angle (degrees) between the spectra of two runs, per shared fluorophore.
.sw_dash_compare_runs <- function(run_a, run_b) {
  a <- .sw_dash_spectra_matrix(run_a)
  b <- .sw_dash_spectra_matrix(run_b)
  shared <- intersect(rownames(a), rownames(b))
  det <- intersect(colnames(a), colnames(b))
  if (length(shared) == 0L || length(det) == 0L) {
    return(NULL)
  }
  sw_unmix_compare_spectra(
    a[shared, det, drop = FALSE],
    b[shared, det, drop = FALSE]
  )
}


# ---- AutoSpectral's own outputs (the `unmix/` folder of a run) ---------------------

# Path of `to` relative to the folder `from_dir`, with "/" separators. Both must
# exist. The dashboard links images this way so the pair of folders can be moved
# together.
.sw_dash_relpath <- function(from_dir, to) {
  f <- strsplit(normalizePath(from_dir, mustWork = TRUE), "/", fixed = TRUE)[[
    1L
  ]]
  t <- strsplit(normalizePath(to, mustWork = TRUE), "/", fixed = TRUE)[[1L]]
  k <- 0L
  while (k < min(length(f), length(t)) && identical(f[[k + 1L]], t[[k + 1L]])) {
    k <- k + 1L
  }
  paste(c(rep("..", length(f) - k), t[-seq_len(k)]), collapse = "/")
}

# A CSV with row names in column 1 as a numeric matrix; NULL when unreadable or
# larger than `max_dim` on a side (a display limit for the page, not a result).
.sw_dash_read_matrix_csv <- function(path, max_dim = 150L) {
  m <- tryCatch(
    as.matrix(utils::read.csv(path, row.names = 1L, check.names = FALSE)),
    error = function(e) NULL
  )
  if (
    is.null(m) || !is.numeric(m) || any(dim(m) > max_dim) || any(dim(m) < 1L)
  ) {
    return(NULL)
  }
  m
}

.sw_dash_matrix_entry <- function(label, m) {
  list(
    label = label,
    rows = rownames(m),
    cols = colnames(m),
    values = .sw_dash_round(as.vector(t(m)), 4L)
  )
}

# What AutoSpectral wrote under <run>/unmix: the matrices it saved as CSV (hotspot,
# similarity, spillover spread), its autofluorescence spectra, the spectra of each
# stage (automated, original, cleaned), and the list of its image and PDF files
# with paths relative to `out_dir` (the folder of the dashboard file). Images are
# linked, not embedded. Returns NULL when the run has no `unmix` folder.
.sw_dash_read_autospectral <- function(run_dir, out_dir) {
  as_dir <- file.path(run_dir, "unmix")
  if (!dir.exists(as_dir)) {
    return(NULL)
  }
  csv <- \(sub, pattern) {
    d <- file.path(as_dir, sub)
    if (!dir.exists(d)) {
      return(character())
    }
    list.files(d, pattern = pattern, full.names = TRUE)
  }
  # Matrices.
  mat_files <- csv("figure_similarity_heatmap", "matrix\\.csv$")
  matrices <- Filter(
    Negate(is.null),
    lapply(mat_files, function(p) {
      m <- .sw_dash_read_matrix_csv(p)
      if (is.null(m)) {
        return(NULL)
      }
      label <- gsub("_", " ", sub("\\.csv$", "", basename(p)))
      # AutoSpectral writes "Automatedhotspot_matrix.csv": put a space after the stage.
      label <- sub(
        "^(Automated|Clean|Normalized)(?=[a-z])",
        "\\1 ",
        label,
        perl = TRUE
      )
      .sw_dash_matrix_entry(label, m)
    })
  )
  # Autofluorescence spectra (one set per unstained file).
  af_files <- csv("table_spectra", "Autofluorescence spectra.*\\.csv$")
  af <- Filter(
    Negate(is.null),
    lapply(af_files, function(p) {
      m <- .sw_dash_read_matrix_csv(p)
      if (is.null(m)) {
        return(NULL)
      }
      list(
        label = sub(" Autofluorescence spectra.*$", "", basename(p)),
        detectors = colnames(m),
        names = rownames(m),
        values = .sw_dash_round(as.vector(t(m)), 4L)
      )
    })
  )
  # Spectra of each stage.
  sp_files <- csv("table_spectra", "autospectral_spectra.*\\.csv$")
  versions <- Filter(
    Negate(is.null),
    lapply(sp_files, function(p) {
      m <- .sw_dash_read_matrix_csv(p)
      if (is.null(m)) {
        return(NULL)
      }
      label <- gsub("_", " ", sub("_autospectral_spectra.*$", "", basename(p)))
      list(
        label = label,
        detectors = colnames(m),
        fluorophores = rownames(m),
        values = .sw_dash_round(as.vector(t(m)), 4L)
      )
    })
  )
  # Image and PDF files, grouped by AutoSpectral's output folder.
  folders <- sort(list.dirs(as_dir, recursive = FALSE, full.names = FALSE))
  folders <- folders[grepl("^figure_", folders)]
  files <- do.call(
    rbind,
    lapply(folders, function(fd) {
      f <- list.files(
        file.path(as_dir, fd),
        pattern = "\\.(jpe?g|png|pdf)$",
        ignore.case = TRUE,
        full.names = TRUE
      )
      if (length(f) == 0L) {
        return(NULL)
      }
      data.frame(
        folder = fd,
        name = basename(f),
        href = vapply(
          f,
          \(x) utils::URLencode(.sw_dash_relpath(out_dir, x)),
          ""
        ),
        pdf = grepl("\\.pdf$", f, ignore.case = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  list(
    matrices = matrices,
    af = af,
    versions = versions,
    gallery = if (is.null(files)) NULL else as.list(files)
  )
}
