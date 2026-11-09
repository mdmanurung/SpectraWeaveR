# Tests for R/unmix_figures.R (F1 to F4 builders, sampling, saving).

# ---- Mock data ---------------------------------------------------------------

mock_spectra <- function() {
  det <- paste0("D", 1:6)
  m <- rbind(
    FITC = c(0.1, 1, 0.5, 0.2, 0.1, 0),
    BB515 = c(0.1, 0.9, 0.6, 0.2, 0.1, 0),
    PE = c(0, 0.1, 0.3, 1, 0.6, 0.2)
  )
  colnames(m) <- det
  m
}

mock_events_list <- function() {
  set.seed(42)
  ch <- c("A", "B", "C")
  list(
    ours = matrix(rnorm(300, sd = 50), ncol = 3, dimnames = list(NULL, ch)),
    spectroflo = matrix(
      rnorm(150, sd = 60),
      ncol = 3,
      dimnames = list(NULL, ch)
    )
  )
}


# ---- Palette and required arguments ------------------------------------------

test_that("unknown method labels error in every builder", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  expect_error(
    sw_unmix_figure_spectra(list(ours = s, mystery = s)),
    "Unknown method label.*mystery"
  )
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 20, seed = 1)
  ev$method[ev$method == "spectroflo"] <- "SpectroFlo"
  expect_error(
    sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = 10),
    "Unknown method label.*SpectroFlo"
  )
  expect_error(
    sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 10),
    "Unknown method label"
  )
  expect_error(
    sw_unmix_sample_events(list(foo = mock_events_list()$ours), 5, 1),
    "Unknown method label.*foo"
  )
})

test_that("required arguments error by name", {
  skip_if_not_installed("ggplot2")
  x <- mock_events_list()
  expect_error(sw_unmix_sample_events(x, seed = 1), "`n_events` is required")
  expect_error(sw_unmix_sample_events(x, n_events = 5), "`seed` is required")
  ev <- sw_unmix_sample_events(x, n_events = 20, seed = 1)
  expect_error(
    sw_unmix_figure_distributions(ev, "A"),
    "`pseudo_log_sigma` is required"
  )
  expect_error(
    sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10),
    "`bins` is required"
  )
  expect_error(
    sw_unmix_figure_biplot(ev, "A", "B", bins = 10),
    "`pseudo_log_sigma` is required"
  )
  expect_error(sw_unmix_figure_spectra(), "`spectra_list` is required")
  expect_error(sw_unmix_figure_similarity(), "`spectra` is required")
  expect_error(sw_unmix_save_figures(list()), "`dir` is required")
})


# ---- F1 ----------------------------------------------------------------------

test_that("F1 aligns spectra by name, not position", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  shuffled <- s[c("PE", "FITC", "BB515"), rev(colnames(s))]
  f1 <- sw_unmix_figure_spectra(list(ours = s, spectroflo = shuffled))
  expect_equal(f1$source$angle_deg, rep(0, nrow(f1$source)))
  a <- f1$source[f1$source$method == "ours", ]
  b <- f1$source[f1$source$method == "spectroflo", ]
  expect_identical(a$fluorophore, b$fluorophore)
  expect_identical(a$detector, b$detector)
  expect_equal(a$value_raw, b$value_raw)
})

test_that("F1 errors on mismatched fluorophore or detector names", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  other <- s
  rownames(other)[1] <- "APC"
  expect_error(
    sw_unmix_figure_spectra(list(ours = s, spectroflo = other)),
    "differ from 'ours'.*APC"
  )
  other <- s[, 1:5]
  expect_error(
    sw_unmix_figure_spectra(list(ours = s, spectroflo = other)),
    "differ from 'ours'.*D6"
  )
})

test_that("F1 peak-normalises each row by max absolute value", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  s["PE", ] <- s["PE", ] * 7
  s["FITC", "D6"] <- -2
  f1 <- sw_unmix_figure_spectra(list(ours = s))
  src <- f1$source
  for (f in rownames(s)) {
    rows <- src[src$fluorophore == f, ]
    expect_equal(rows$value_norm, rows$value_raw / max(abs(rows$value_raw)))
  }
  expect_equal(
    max(abs(src$value_norm[src$fluorophore == "FITC"])),
    1
  )
  expect_equal(
    src$value_norm[src$fluorophore == "FITC" & src$detector == "D6"],
    -1
  )
})

test_that("F1 angles: identical = 0, orthogonal = 90, shown in the strip", {
  skip_if_not_installed("ggplot2")
  a <- rbind(X = c(1, 0, 0), Y = c(0, 1, 1))
  colnames(a) <- c("D1", "D2", "D3")
  b <- rbind(X = c(0, 1, 0), Y = c(0, 2, 2))
  colnames(b) <- colnames(a)
  f1 <- sw_unmix_figure_spectra(list(ours = a, spectracle = b))
  ang <- unique(f1$source[
    f1$source$method == "spectracle",
    c("fluorophore", "angle_deg")
  ])
  expect_equal(ang$angle_deg[ang$fluorophore == "X"], 90, tolerance = 1e-10)
  expect_equal(ang$angle_deg[ang$fluorophore == "Y"], 0, tolerance = 1e-10)
  expect_true(any(grepl(
    "spectracle 90.0 deg",
    f1$source$panel_label,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "spectracle 0.0 deg",
    f1$source$panel_label,
    fixed = TRUE
  )))
  # Angles match sw_unmix_compare_spectra(), the spectra_compare.csv source.
  cmp <- sw_unmix_compare_spectra(a, b)
  expect_equal(
    ang$angle_deg[match(cmp$fluorophore, ang$fluorophore)],
    cmp$angle_deg
  )
})

test_that("F1 plot has the expected layers and plots the source rows", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  f1 <- sw_unmix_figure_spectra(list(ours = s, spectroflo = s * 1.1))
  expect_s3_class(f1$plot, "ggplot")
  geoms <- vapply(f1$plot$layers, \(l) class(l$geom)[1], character(1))
  expect_identical(unname(geoms), c("GeomLine", "GeomPoint"))
  expect_equal(nrow(f1$plot$data), 2 * length(s))
  expect_equal(nrow(f1$source), nrow(f1$plot$data))
  expect_equal(f1$plot$data$value_norm, f1$source$value_norm)
  expect_identical(levels(f1$plot$data$detector), colnames(s))
  built <- ggplot2::ggplot_build(f1$plot)
  expect_equal(nrow(ggplot2::layer_data(f1$plot, 2)), nrow(f1$source))
  expect_false(is.null(built))
})


# ---- F2 ----------------------------------------------------------------------

test_that("F2 similarity has unit diagonal, is symmetric, and keeps order", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  f2 <- sw_unmix_figure_similarity(s)
  expect_equal(unname(diag(f2$matrix)), rep(1, nrow(s)))
  expect_equal(f2$matrix, t(f2$matrix))
  expect_identical(rownames(f2$matrix), rownames(s))
  expect_equal(nrow(f2$source), nrow(s)^2)
  hand <- sum(s["FITC", ] * s["PE", ]) /
    sqrt(sum(s["FITC", ]^2) * sum(s["PE", ]^2))
  expect_equal(
    f2$source$cosine[
      f2$source$fluorophore_x == "FITC" & f2$source$fluorophore_y == "PE"
    ],
    hand
  )
  expect_identical(levels(f2$plot$data$fluorophore_x), rownames(s))
})

test_that("F2 plot is a tile layer over the full matrix with 0 to 1 limits", {
  skip_if_not_installed("ggplot2")
  f2 <- sw_unmix_figure_similarity(mock_spectra())
  expect_s3_class(f2$plot, "ggplot")
  expect_identical(class(f2$plot$layers[[1]]$geom)[1], "GeomTile")
  expect_equal(nrow(ggplot2::layer_data(f2$plot)), nrow(f2$source))
  expect_equal(f2$plot$data$cosine, f2$source$cosine)
  expect_match(f2$plot$labels$caption, "FITC vs BB515|BB515 vs FITC")
})

test_that("F2 keeps a full-colour diagonal when self cosine rounds above 1", {
  skip_if_not_installed("ggplot2")
  # Seed 2 gives a self cosine of 1 + 2.2e-16 in floating point.
  set.seed(2)
  s <- matrix(runif(60, 0.1, 3), nrow = 6)
  dimnames(s) <- list(paste0("F", 1:6), paste0("D", 1:10))
  expect_gt(max(.sw_row_cosine(s)), 1)
  f2 <- sw_unmix_figure_similarity(s)
  expect_false(any(f2$plot$data$cosine > 1))
  expect_false(anyNA(ggplot2::layer_data(f2$plot)$fill))
  expect_equal(unname(diag(f2$matrix)), rep(1, 6))
})

test_that("F2 clamps negative cosines to 0, says so, and errors on zero rows", {
  skip_if_not_installed("ggplot2")
  s <- rbind(A = c(1, 0), B = c(-1, 0.1))
  colnames(s) <- c("D1", "D2")
  f2 <- sw_unmix_figure_similarity(s)
  expect_equal(
    f2$source$cosine[
      f2$source$fluorophore_x == "B" & f2$source$fluorophore_y == "A"
    ],
    0
  )
  expect_lt(min(f2$source$cosine_raw), 0)
  expect_match(f2$plot$labels$caption, "drawn as 0")
  expect_false(anyNA(ggplot2::layer_data(f2$plot)$fill))
  s <- rbind(A = c(1, 0), B = c(0, 0))
  colnames(s) <- c("D1", "D2")
  expect_error(sw_unmix_figure_similarity(s), "all-zero row")
})


# ---- Sampling ----------------------------------------------------------------

test_that("sampling is reproducible and restores the RNG state", {
  x <- mock_events_list()
  set.seed(99)
  before <- .Random.seed
  a <- sw_unmix_sample_events(x, n_events = 30, seed = 1)
  expect_identical(.Random.seed, before)
  b <- sw_unmix_sample_events(x, n_events = 30, seed = 1)
  expect_identical(a, b)
  c2 <- sw_unmix_sample_events(x, n_events = 30, seed = 2)
  expect_false(identical(a$event_id, c2$event_id))
})

test_that("sampling removes .Random.seed when none existed", {
  x <- mock_events_list()
  old <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  rm(".Random.seed", envir = globalenv())
  sw_unmix_sample_events(x, n_events = 10, seed = 1)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
})

test_that("sample_log records available and sampled counts", {
  x <- mock_events_list()
  ev <- sw_unmix_sample_events(x, n_events = 70, seed = 1)
  log <- attr(ev, "sample_log")
  expect_identical(log$file, c("ours", "spectroflo"))
  expect_identical(log$n_available, c(100L, 50L))
  expect_identical(log$n_sampled, c(70L, 50L))
  expect_equal(as.vector(table(ev$method)[log$method]), log$n_sampled)
  expect_identical(names(ev), c("file", "method", "event_id", "A", "B", "C"))
  # Values come from the sampled rows of the source matrix.
  rows <- ev[ev$method == "ours", ]
  expect_equal(
    unname(as.matrix(rows[, c("A", "B", "C")])),
    unname(x$ours[rows$event_id, ])
  )
  expect_false(anyDuplicated(rows$event_id) > 0)
})

test_that("sampling aligns channels by name and errors on a mismatch", {
  x <- mock_events_list()
  x$spectroflo <- x$spectroflo[, c("C", "A", "B")]
  ev <- sw_unmix_sample_events(x, n_events = 200, seed = 1)
  rows <- ev[ev$method == "spectroflo", ]
  expect_equal(rows$A, unname(x$spectroflo[rows$event_id, "A"]))
  colnames(x$spectroflo)[1] <- "Z"
  expect_error(sw_unmix_sample_events(x, 10, 1), "Channels of 'spectroflo'")
})

test_that("method labels can be given separately from file names", {
  x <- mock_events_list()
  names(x) <- c("lung_a.fcs", "lung_b.fcs")
  ev <- sw_unmix_sample_events(x, 10, 1, method = c("ours", "spectroflo"))
  expect_identical(unique(ev$file), c("lung_a.fcs", "lung_b.fcs"))
  expect_identical(unique(ev$method), c("ours", "spectroflo"))
})


# ---- F3 and F4 ---------------------------------------------------------------

test_that("F3 source matches the plotted data and caption counts", {
  skip_if_not_installed("ggplot2")
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 80, seed = 1)
  ev$B[1] <- NA
  f3 <- sw_unmix_figure_distributions(ev, c("A", "B"), pseudo_log_sigma = 10)
  expect_s3_class(f3$plot, "ggplot")
  expect_identical(class(f3$plot$layers[[1]]$geom)[1], "GeomDensity")
  expect_equal(nrow(f3$plot$data), nrow(f3$source))
  expect_equal(f3$plot$data$value, f3$source$value)
  expect_equal(nrow(f3$source), 2 * nrow(ev) - 1)
  n_ours_b <- f3$counts$n_events[
    f3$counts$method == "ours" & f3$counts$channel == "B"
  ]
  expect_equal(n_ours_b, 79)
  expect_match(
    f3$plot$labels$caption,
    "ours 79 to 80; spectroflo 50",
    fixed = TRUE
  )
  expect_match(f3$plot$labels$x, "sigma = 10", fixed = TRUE)
  expect_error(
    sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = -1),
    "positive"
  )
})

test_that("F3 and F4 refuse two files under one method", {
  skip_if_not_installed("ggplot2")
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 20, seed = 1)
  ev$file[1] <- "other.fcs"
  expect_error(
    sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = 10),
    "one file"
  )
  expect_error(
    sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 10),
    "one file"
  )
})

test_that("F4 plots finite pairs per method with shared axes", {
  skip_if_not_installed("ggplot2")
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 80, seed = 1)
  ev$A[2] <- Inf
  f4 <- sw_unmix_figure_biplot(ev, "A", "C", pseudo_log_sigma = 10, bins = 20)
  expect_s3_class(f4$plot, "ggplot")
  expect_identical(class(f4$plot$layers[[1]]$stat)[1], "StatBin2d")
  expect_equal(nrow(f4$source), nrow(ev) - 1)
  expect_equal(f4$plot$data$x, f4$source$x)
  expect_identical(f4$counts$n_events, c(79L, 50L))
  expect_match(
    gsub("\\s+", " ", f4$plot$labels$caption),
    "ours 79; spectroflo 50. 20 bins",
    fixed = TRUE
  )
  # Both axes break at 0 or signed powers of ten.
  built <- ggplot2::ggplot_build(f4$plot)
  trans <- scales::transform_pseudo_log(sigma = 10)
  for (sc in list(
    built$layout$panel_scales_x[[1]],
    built$layout$panel_scales_y[[1]]
  )) {
    b <- trans$inverse(stats::na.omit(sc$get_breaks()))
    b <- b[b != 0]
    expect_equal(log10(abs(b)), round(log10(abs(b))), tolerance = 1e-8)
  }
  expect_identical(f4$plot$facet$params$free, list(x = FALSE, y = FALSE))
  expect_equal(sum(ggplot2::layer_data(f4$plot)$count), nrow(f4$source))
})


test_that("pseudo-log breaks are at least 8% of the axis apart (sigma 1000)", {
  skip_if_not_installed("ggplot2")
  trans <- scales::transform_pseudo_log(sigma = 1000)
  limits <- c(-3e5, 4e6)
  b <- .sw_fig_symlog_breaks(trans)(limits)
  pos <- trans$transform(b)
  width <- diff(trans$transform(limits))
  expect_true(all(diff(sort(pos)) >= 0.08 * width))
  # +-10^3 sit almost on 0 at this sigma and are dropped.
  expect_true(0 %in% b)
  expect_false(any(c(-1e3, 1e3) %in% b))
  expect_true(all(c(-1e5, 1e4, 1e5, 1e6) %in% b))
  # The same rule holds on the built F3 axis.
  set.seed(4)
  x <- list(
    ours = cbind(A = c(rnorm(200, 0, 300), rnorm(100, 2e5, 5e4), -2e5, 3e6))
  )
  ev <- sw_unmix_sample_events(x, n_events = 1000, seed = 1)
  f3 <- sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = 1000)
  sc <- ggplot2::ggplot_build(f3$plot)$layout$panel_scales_x[[1]]
  pos <- sort(stats::na.omit(sc$get_breaks()))
  expect_true(all(diff(pos) >= 0.08 * diff(sc$get_limits())))
})

test_that("F3 legend keys draw lines, F4 keeps a colourbar", {
  skip_if_not_installed("ggplot2")
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 20, seed = 1)
  f3 <- sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = 10)
  expect_identical(
    environment(f3$plot$layers[[1]]$geom$draw_key)$f,
    ggplot2::draw_key_path
  )
  f4 <- sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 10)
  fill_scale <- f4$plot$scales$get_scales("fill")
  expect_identical(fill_scale$guide, "colourbar")
})


# ---- Saving and side effects -------------------------------------------------

test_that("builders draw nothing and write nothing", {
  skip_if_not_installed("ggplot2")
  dev_before <- grDevices::dev.list()
  files_before <- list.files(tempdir(), recursive = TRUE)
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 20, seed = 1)
  s <- mock_spectra()
  sw_unmix_figure_spectra(list(ours = s))
  sw_unmix_figure_similarity(s)
  sw_unmix_figure_distributions(ev, "A", pseudo_log_sigma = 10)
  sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 10)
  expect_identical(grDevices::dev.list(), dev_before)
  expect_identical(list.files(tempdir(), recursive = TRUE), files_before)
})

test_that("save writes non-empty files and refuses to overwrite", {
  skip_if_not_installed("ggplot2")
  s <- mock_spectra()
  ev <- sw_unmix_sample_events(mock_events_list(), n_events = 40, seed = 1)
  figs <- list(
    F1 = sw_unmix_figure_spectra(list(ours = s, spectroflo = s * 0.9)),
    F2 = sw_unmix_figure_similarity(s),
    F3 = sw_unmix_figure_distributions(ev, c("A", "B"), pseudo_log_sigma = 10),
    F4 = sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 10)
  )
  dir <- tempfile("figs")
  out <- sw_unmix_save_figures(figs, dir)
  expect_identical(out$panel, c("F1", "F2", "F3", "F4"))
  expect_identical(
    basename(out$png),
    c(
      "spectra-overlay.png",
      "similarity-matrix.png",
      "unmixed-distributions.png",
      "biplot.png"
    )
  )
  paths <- c(out$png, out$source_csv, out$pdf[!is.na(out$pdf)])
  expect_true(all(file.exists(paths)))
  expect_true(all(file.size(paths) > 0))
  csv <- utils::read.csv(out$source_csv[2])
  expect_equal(nrow(csv), out$n_rows[2])
  expect_equal(csv$cosine, figs$F2$source$cosine)
  mtimes <- file.mtime(paths)
  expect_error(sw_unmix_save_figures(figs["F2"], dir), "Refusing to overwrite")
  expect_identical(file.mtime(paths), mtimes)
  unlink(dir, recursive = TRUE)
})

test_that("save uses per-panel sizes and accepts named overrides", {
  skip_if_not_installed("ggplot2")
  f2 <- sw_unmix_figure_similarity(mock_spectra())
  dir <- tempfile("figs")
  out <- sw_unmix_save_figures(list(F2 = f2), dir)
  expect_identical(c(out$width_mm, out$height_mm), c(120, 100))
  dir2 <- tempfile("figs")
  out2 <- sw_unmix_save_figures(
    list(F2 = f2),
    dir2,
    sizes = list(F2 = c(60, 50))
  )
  expect_identical(c(out2$width_mm, out2$height_mm), c(60, 50))
  # Width ratio of the two PNGs follows the sizes (same 300 dpi).
  dims <- function(p) {
    con <- file(p, "rb")
    on.exit(close(con))
    raw <- readBin(con, "raw", 24)
    c(
      sum(as.integer(raw[17:20]) * 256^(3:0)),
      sum(as.integer(raw[21:24]) * 256^(3:0))
    )
  }
  expect_lt(abs(dims(out$png)[1] - 120 / 25.4 * 300), 2)
  expect_lt(abs(dims(out2$png)[1] - 60 / 25.4 * 300), 2)
  expect_error(
    sw_unmix_save_figures(
      list(F2 = f2),
      tempfile(),
      sizes = list(F9 = c(1, 1))
    ),
    "`sizes` must be"
  )
  unlink(c(dir, dir2), recursive = TRUE)
})

test_that("save checks panel names against the results", {
  skip_if_not_installed("ggplot2")
  f2 <- sw_unmix_figure_similarity(mock_spectra())
  expect_error(
    sw_unmix_save_figures(list(F5 = f2), tempfile()),
    "Unknown panel"
  )
  expect_error(
    sw_unmix_save_figures(list(F1 = f2), tempfile()),
    "do not match their `panel`"
  )
})
