#' Write an interactive unmixing dashboard as one HTML file
#'
#' Reads one or more finished run folders (the output of the benchmark runner
#' or of [sw_unmix_report()] plus [sw_unmix_save_figures()]) and writes a single
#' self-contained HTML file. The file has inline CSS, JavaScript and SVG. It
#' loads nothing from the web, so it works offline and can be mailed.
#'
#' Sections: overview, checks (filter, search, sort), spreading inflation
#' factor (SIF) with the run's cutoffs, spectra per fluorophore, a clickable
#' similarity matrix, unmixed distributions, a biplot, AutoSpectral's own
#' outputs, and a comparison of two or more runs.
#'
#' **Spike-in and marker screen.** A `Spike & screen` tab appears when a run's
#' `report/` folder holds `spike_summary.csv`, `baseline_check.csv` or
#' `marker_screen.csv` (written by [sw_unmix_report()] from diagnostics made
#' with `spike` or `screen`). It shows those tables and the F7 and F8 PNGs
#' found in any `figures*/panels` folder (file names starting with
#' `spike-summary` or `spike-ecdf`), embedded in the page. Runs without these
#' files show no such tab.
#'
#' **Event metrics.** An `Event metrics` tab appears when a run's `report/`
#' folder holds `event_channels.csv`, `event_pairs.csv`, `heldout_scores.csv`
#' or `heldout_repro.csv` (written by [sw_unmix_report()] from diagnostics made
#' with `event_metrics` or `heldout`). Pairs are listed by |slope|, largest
#' first.
#'
#' **AutoSpectral outputs.** The `AutoSpectral` tab reads the folder
#' `<run>/unmix` that AutoSpectral wrote. Its CSV tables become interactive
#' views: the hotspot, similarity and spillover-spread matrices, the
#' autofluorescence spectra, and the spectra of each stage (automated, original,
#' cleaned). Its JPG and PDF figures (scatter, gates, variants, cleaned controls,
#' QC reports) are listed in a gallery as links relative to the dashboard file.
#' They are not embedded, so they show only while the dashboard file sits next
#' to the run folder. The per-fluorophore variant spectra (`Spectral_variants.rds`)
#' are not read.
#'
#' **Verdicts come from the run.** The status of every check is the one saved in
#' the run's `report/checks.csv`. The SIF tab has two sliders that recolour the
#' bars with cutoffs of the viewer's choice. They are display only, nothing is
#' saved, and a banner says so. The sliders start at the run's saved cutoffs.
#'
#' **Event-level panels are binned.** The distribution and biplot tables hold
#' tens of thousands of events. They are binned here on the pseudo-log scale
#' `scales::transform_pseudo_log(sigma = pseudo_log_sigma)` with `bins` bins per
#' axis, and the counts are embedded. Both are choices that change what the
#' panels show, so neither has a default.
#'
#' @param runs Run folders. A character vector; names become the run labels.
#'   Each folder needs `report/checks.csv`, `report/spectra.csv` and
#'   `report/provenance.csv`. Figure source tables are read from the last
#'   `figures*/panels/source_data` folder that holds them; a panel without its
#'   table is left out.
#' @param path Output file. Must not exist.
#' @param pseudo_log_sigma Positive number, the linear-zone width of the
#'   pseudo-log scale used to bin events. Required, no default.
#' @param bins Whole number of at least 10, the bins per axis for the
#'   distribution and biplot panels. Required, no default.
#' @param title Page title.
#' @param figures_dir Optional folder name inside each run (for example
#'   `"figures_v2"`) to read the figure tables from, instead of the last
#'   `figures*` folder.
#' @param axis_labels Optional character vector of length 2: names of the
#'   biplot axes (the source table does not store them). `NULL` shows
#'   "x" and "y".
#'
#' @return `path`, invisibly.
#'
#' @seealso [sw_unmix_app()] for a live version, [sw_unmix_report()] for the
#'   static report.
#'
#' @examples
#' \dontrun{
#' sw_unmix_dashboard(
#'   c(cells = "runs/cells-set-v2", beads = "runs/beads-set-v2"),
#'   path = "dashboard.html",
#'   pseudo_log_sigma = 1000,
#'   bins = 200
#' )
#' }
#'
#' @export
sw_unmix_dashboard <- function(
  runs,
  path,
  pseudo_log_sigma,
  bins,
  title = "Unmixing dashboard",
  figures_dir = NULL,
  axis_labels = NULL
) {
  .sw_fig_require(missing(pseudo_log_sigma), "pseudo_log_sigma")
  .sw_fig_require(missing(bins), "bins")
  .sw_dash_check_args(runs, path, pseudo_log_sigma, bins, title, axis_labels)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for the dashboard.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for the dashboard.", call. = FALSE)
  }

  labels <- .sw_dash_run_labels(runs)
  run_data <- Map(
    function(dir, nm) {
      .sw_dash_read_run(
        dir,
        pseudo_log_sigma,
        as.integer(bins),
        figures_dir,
        nm
      )
    },
    unname(runs),
    labels
  )
  run_data <- unname(run_data)
  out_dir <- normalizePath(dirname(path), mustWork = TRUE)
  payload <- .sw_dash_payload(
    run_data,
    title,
    pseudo_log_sigma,
    bins,
    axis_labels,
    out_dir
  )
  html <- .sw_dash_render(payload)
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}

.sw_dash_check_args <- function(runs, path, sigma, bins, title, axis_labels) {
  .sw_dash_check_settings(runs, sigma, bins)
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("`path` must be one file path.", call. = FALSE)
  }
  if (!dir.exists(dirname(path))) {
    stop("The folder of `path` does not exist: ", dirname(path), call. = FALSE)
  }
  if (file.exists(path)) {
    stop("`path` exists; refusing to overwrite: ", path, call. = FALSE)
  }
  if (!is.character(title) || length(title) != 1L || is.na(title)) {
    stop("`title` must be one string.", call. = FALSE)
  }
  if (
    !is.null(axis_labels) &&
      (!is.character(axis_labels) || length(axis_labels) != 2L)
  ) {
    stop("`axis_labels` must be NULL or two strings.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---- Payload (what the page reads as JSON) ------------------------------------

.sw_dash_round <- function(x, digits = 5L) signif(x, digits)

.sw_dash_spectra_payload <- function(run) {
  ov <- run$spectra_overlay
  if (is.null(ov)) {
    long <- run$spectra_long
    fl <- unique(long$fluorophore)
    dt <- unique(long$detector)
    mat <- .sw_dash_spectra_matrix(run)
    peak <- apply(abs(mat), 1L, max)
    items <- lapply(seq_along(fl), function(i) {
      list(ours = .sw_dash_round(mat[i, ] / peak[[i]]), angle = list(ours = 0))
    })
    names(items) <- fl
    return(list(
      detectors = dt,
      fluorophores = fl,
      methods = "ours",
      items = items
    ))
  }
  dt <- unique(ov$detector)
  fl <- unique(ov$fluorophore)
  methods <- unique(ov$method)
  items <- lapply(fl, function(f) {
    sub <- ov[ov$fluorophore == f, , drop = FALSE]
    vals <- lapply(methods, function(m) {
      s <- sub[sub$method == m, , drop = FALSE]
      if (nrow(s) == 0L) {
        return(NULL)
      }
      .sw_dash_round(s$value_norm[match(dt, s$detector)])
    })
    names(vals) <- methods
    ang <- lapply(methods, function(m) {
      a <- sub$angle_deg[sub$method == m]
      if (length(a)) .sw_dash_round(a[[1L]]) else NULL
    })
    names(ang) <- methods
    c(Filter(Negate(is.null), vals), list(angle = Filter(Negate(is.null), ang)))
  })
  names(items) <- fl
  list(detectors = dt, fluorophores = fl, methods = methods, items = items)
}

.sw_dash_similarity_payload <- function(run) {
  s <- run$similarity
  if (is.null(s)) {
    return(NULL)
  }
  fl <- unique(s$fluorophore_x)
  m <- matrix(
    NA_real_,
    length(fl),
    length(fl),
    dimnames = list(fl, fl)
  )
  m[cbind(match(s$fluorophore_y, fl), match(s$fluorophore_x, fl))] <- s$cosine
  list(fluorophores = fl, cosine = .sw_dash_round(as.vector(t(m)), 4L))
}

.sw_dash_dist_payload <- function(run) {
  d <- run$distributions
  if (is.null(d)) {
    return(NULL)
  }
  panels <- lapply(d$panels, function(p) {
    list(
      channel = p$channel,
      edges = .sw_dash_round(p$edges, 6L),
      counts = p$counts,
      n = as.list(stats::setNames(as.integer(p$n), d$methods))
    )
  })
  list(channels = d$channels, methods = d$methods, panels = unname(panels))
}

# Non-empty bins travel as delta-coded keys: key = y * bins + x, and `k` holds
# the differences between successive keys (the first from -1). Keys are sorted,
# so most differences are small integers. This is lossless; the page decodes
# it back to (x, y) and the same counts.
.sw_dash_biplot_payload <- function(run, axis_labels, bins) {
  b <- run$biplot
  if (is.null(b)) {
    return(NULL)
  }
  list(
    methods = b$methods,
    ex = .sw_dash_round(b$ex, 6L),
    ey = .sw_dash_round(b$ey, 6L),
    cells = lapply(b$cells, function(cc) {
      key <- cc$y * bins + cc$x
      list(k = diff(c(-1L, key)), n = cc$n, total = cc$total)
    }),
    labels = if (is.null(axis_labels)) c("x", "y") else axis_labels
  )
}

.sw_dash_run_payload <- function(run, axis_labels, bins, out_dir) {
  list(
    name = run$name,
    run_id = run$run_id,
    figures_dir = if (is.na(run$figures_dir)) NULL else run$figures_dir,
    status_counts = run$status_counts,
    limits = list(
      warn = if (is.na(run$limits[["warn"]])) {
        NULL
      } else {
        unname(run$limits[["warn"]])
      },
      fail = if (is.na(run$limits[["fail"]])) {
        NULL
      } else {
        unname(run$limits[["fail"]])
      }
    ),
    provenance = as.list(run$provenance),
    checks = as.list(run$checks),
    sif = as.list(run$sif),
    sif_source = run$sif_source,
    spectra = .sw_dash_spectra_payload(run),
    similarity = .sw_dash_similarity_payload(run),
    distributions = .sw_dash_dist_payload(run),
    biplot = .sw_dash_biplot_payload(run, axis_labels, bins),
    spectracle = if (is.null(run$compare_spectracle)) {
      NULL
    } else {
      as.list(run$compare_spectracle)
    },
    autospectral = .sw_dash_read_autospectral(run$run_dir, out_dir),
    spike = .sw_dash_spike_payload(run$spike),
    events = if (is.null(run$events)) NULL else lapply(run$events, \(x) if (is.null(x)) NULL else as.list(x))
  )
}

.sw_dash_spike_payload <- function(spike) {
  if (is.null(spike)) {
    return(NULL)
  }
  tbl <- \(x) if (is.null(x)) NULL else as.list(x)
  list(
    summary = tbl(spike$summary),
    baseline = tbl(spike$baseline),
    screen = tbl(spike$screen),
    figures = spike$figures
  )
}

.sw_dash_payload <- function(
  run_data,
  title,
  sigma,
  bins,
  axis_labels,
  out_dir
) {
  methods <- unique(unlist(lapply(run_data, function(r) {
    c(
      r$spectra_overlay$method,
      r$distributions$methods,
      r$biplot$methods
    )
  })))
  compare <- NULL
  if (length(run_data) >= 2L) {
    compare <- lapply(seq_along(run_data)[-1L], function(i) {
      cmp <- .sw_dash_compare_runs(run_data[[1L]], run_data[[i]])
      list(
        a = run_data[[1L]]$name,
        b = run_data[[i]]$name,
        angles = if (is.null(cmp)) NULL else as.list(cmp)
      )
    })
  }
  list(
    title = title,
    settings = list(pseudo_log_sigma = sigma, bins = as.integer(bins)),
    palette = as.list(.sw_fig_palette(methods)),
    package_version = as.character(utils::packageVersion("SpectraWeaveR")),
    runs = lapply(
      run_data,
      .sw_dash_run_payload,
      axis_labels = axis_labels,
      bins = as.integer(bins),
      out_dir = out_dir
    ),
    compare = compare
  )
}

# ---- HTML ----------------------------------------------------------------------

.sw_dash_render <- function(payload) {
  tpl_path <- system.file(
    "dashboard",
    "dashboard.html",
    package = "SpectraWeaveR"
  )
  if (!nzchar(tpl_path)) {
    # Source-tree fallback for devtools::load_all().
    tpl_path <- file.path("inst", "dashboard", "dashboard.html")
  }
  if (!file.exists(tpl_path)) {
    stop("Dashboard template not found.", call. = FALSE)
  }
  tpl <- paste(
    readLines(tpl_path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    na = "null",
    null = "null",
    digits = NA,
    dataframe = "columns"
  )
  # A script block ends at "</", so neutralise it inside the data.
  json <- gsub("</", "<\\/", as.character(json), fixed = TRUE)
  title <- .sw_html_escape(payload$title)
  tpl <- sub("__TITLE__", title, tpl, fixed = TRUE)
  # fixed = TRUE with a function-free replacement: no backslash handling.
  pos <- regexpr("__DATA__", tpl, fixed = TRUE)
  paste0(
    substr(tpl, 1L, pos - 1L),
    json,
    substr(tpl, pos + attr(pos, "match.length"), nchar(tpl))
  )
}
