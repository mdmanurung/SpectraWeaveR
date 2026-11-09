# tests/testthat/test-unmix_dashboard.R
# Dashboard data model, HTML writer and Shiny app, on small mock run folders.

# A run folder with the files the runner and sw_unmix_report write. SIF values
# are fixed so every status is known: A 7 (fail), B 5 (warn), C 2 (pass).
make_run <- function(dir, sif = c(A = 7, B = 5, C = 2), limit = "warn if > 4; fail if > 6") {
  dir.create(file.path(dir, "report"), recursive = TRUE)
  src <- file.path(dir, "figures", "panels", "source_data")
  dir.create(src, recursive = TRUE)
  fl <- names(sif)
  status <- ifelse(sif > 6, "fail", ifelse(sif > 4, "warn", "pass"))
  checks <- rbind(
    data.frame(
      check_id = "controls", scope = "run", item = "control_file",
      metric = "n_control_errors", value = 0, limit = NA, status = "not_judged",
      reason = "Counts of issues </script> <b>bold</b>."
    ),
    data.frame(
      check_id = "spectra_sif", scope = "fluorophore", item = fl,
      metric = "max_sif", value = unname(sif), limit = limit, status = status,
      reason = paste0("Diagonal of the hotspot matrix. Largest off-diagonal entry: ", rev(fl), " (1.5).")
    )
  )
  utils::write.csv(checks, file.path(dir, "report", "checks.csv"), row.names = FALSE)
  det <- paste0("D", 1:4)
  long <- expand.grid(fluorophore = fl, detector = det, stringsAsFactors = FALSE)
  long$value <- seq_len(nrow(long)) / 10
  utils::write.csv(long, file.path(dir, "report", "spectra.csv"), row.names = FALSE)
  prov <- data.frame(section = "run", field = "run_id", value = paste0("id-", basename(dir)))
  utils::write.csv(prov, file.path(dir, "report", "provenance.csv"), row.names = FALSE)

  ov <- expand.grid(
    method = c("ours", "spectracle"), fluorophore = fl, detector = det,
    stringsAsFactors = FALSE
  )
  ov$value_raw <- 0.5
  ov$value_norm <- 0.5
  ov$angle_deg <- ifelse(ov$method == "ours", 0, 1.5)
  ov$panel_label <- "label"
  utils::write.csv(ov, file.path(src, "spectra-overlay.csv"), row.names = FALSE)
  sim <- expand.grid(fluorophore_x = fl, fluorophore_y = fl, stringsAsFactors = FALSE)
  sim$cosine <- ifelse(sim$fluorophore_x == sim$fluorophore_y, 1, 0.3)
  sim$cosine_raw <- sim$cosine
  utils::write.csv(sim, file.path(src, "similarity-matrix.csv"), row.names = FALSE)

  set.seed(1)
  dist <- do.call(rbind, lapply(c("ols", "spectroflo"), function(m) {
    do.call(rbind, lapply(c("X-A", "Y-A"), function(ch) {
      data.frame(
        method = m, file = m, event_id = 1:200, channel = ch,
        value = c(rnorm(100, 0, 500), rnorm(100, 50000, 8000))
      )
    }))
  }))
  utils::write.csv(dist, file.path(src, "unmixed-distributions.csv"), row.names = FALSE)
  bp <- do.call(rbind, lapply(c("ols", "spectroflo"), function(m) {
    data.frame(method = m, file = m, event_id = 1:200, x = rnorm(200, 1000, 3000), y = rnorm(200, 2000, 3000))
  }))
  utils::write.csv(bp, file.path(src, "biplot.csv"), row.names = FALSE)
  dir
}

# Pull the embedded JSON payload back out of the HTML.
payload_of <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  m <- regmatches(html, regexec("const D = (\\{.*?\\});\n\\(function", html))[[1L]]
  jsonlite::fromJSON(m[[2L]], simplifyVector = FALSE)
}

test_that("a run folder is read into tables with the saved limits and SIF partners", {
  skip_if_not_installed("scales")
  run <- .sw_dash_read_run(make_run(file.path(tempdir(), "r1")), sigma = 1000, bins = 20)
  expect_equal(unname(run$limits), c(4, 6))
  expect_identical(run$sif$fluorophore, c("A", "B", "C"))
  expect_identical(run$sif$partner, c("C", "B", "A"))
  expect_equal(run$sif$partner_value, c(1.5, 1.5, 1.5))
  expect_equal(run$status_counts$fail, 1L)
  expect_equal(run$status_counts$not_judged, 1L)
  expect_identical(run$run_id, "id-r1")
})

test_that("a limit without text gives NA cutoffs", {
  expect_true(all(is.na(.sw_dash_parse_limit(NA_character_))))
  expect_equal(unname(.sw_dash_parse_limit("warn if > 4; fail if > 6")), c(4, 6))
  expect_equal(unname(.sw_dash_parse_limit("fail if > 6")[["fail"]]), 6)
})

test_that("binned counts add up to the events read", {
  skip_if_not_installed("scales")
  run <- .sw_dash_read_run(make_run(file.path(tempdir(), "r2")), sigma = 1000, bins = 20)
  d <- run$distributions
  expect_setequal(d$channels, c("X-A", "Y-A"))
  for (ch in d$channels) {
    for (m in d$methods) {
      expect_equal(sum(d$panels[[ch]]$counts[[m]]), 200L)
    }
  }
  b <- run$biplot
  for (m in b$methods) {
    expect_equal(sum(b$cells[[m]]$n), 200L)
    expect_equal(b$cells[[m]]$total, 200L)
    expect_true(all(b$cells[[m]]$x >= 0 & b$cells[[m]]$x < 20))
  }
})

# The binning before the scan()-based reader, kept here as the reference. The
# new functions must give identical counts, edges and totals.
old_bin_distributions <- function(path, sigma, bins) {
  d <- utils::read.csv(path, colClasses = c("character", "NULL", "NULL", "character", "numeric"), stringsAsFactors = FALSE)
  trans <- scales::transform_pseudo_log(sigma = sigma)
  d$t <- trans$transform(d$value)
  channels <- unique(d$channel)
  methods <- unique(d$method)
  out <- lapply(channels, function(ch) {
    sub <- d[d$channel == ch, , drop = FALSE]
    r <- range(sub$t, na.rm = TRUE)
    edges <- seq(r[[1L]], r[[2L]], length.out = bins + 1L)
    counts <- lapply(methods, function(m) {
      tv <- sub$t[sub$method == m]
      idx <- findInterval(tv, edges, rightmost.closed = TRUE)
      tabulate(idx[idx >= 1L & idx <= bins], nbins = bins)
    })
    names(counts) <- methods
    list(channel = ch, edges = edges, counts = counts, n = lengths(lapply(methods, \(m) which(sub$method == m))))
  })
  names(out) <- channels
  list(channels = channels, methods = methods, panels = out)
}
old_bin_biplot <- function(path, sigma, bins) {
  d <- utils::read.csv(path, colClasses = c("character", "NULL", "NULL", "numeric", "numeric"), stringsAsFactors = FALSE)
  trans <- scales::transform_pseudo_log(sigma = sigma)
  tx <- trans$transform(d$x)
  ty <- trans$transform(d$y)
  rx <- range(tx, na.rm = TRUE)
  ry <- range(ty, na.rm = TRUE)
  ex <- seq(rx[[1L]], rx[[2L]], length.out = bins + 1L)
  ey <- seq(ry[[1L]], ry[[2L]], length.out = bins + 1L)
  methods <- unique(d$method)
  cells <- lapply(methods, function(m) {
    keep <- d$method == m
    ix <- findInterval(tx[keep], ex, rightmost.closed = TRUE)
    iy <- findInterval(ty[keep], ey, rightmost.closed = TRUE)
    ok <- ix >= 1L & ix <= bins & iy >= 1L & iy <= bins
    key <- (iy[ok] - 1L) * bins + ix[ok]
    n <- tabulate(key, nbins = bins * bins)
    nz <- which(n > 0L)
    list(x = (nz - 1L) %% bins, y = (nz - 1L) %/% bins, n = n[nz], total = sum(keep))
  })
  names(cells) <- methods
  list(methods = methods, ex = ex, ey = ey, cells = cells)
}

test_that("the scan()-based binning gives the same result as the read.csv() version", {
  skip_if_not_installed("scales")
  dir <- make_run(file.path(tempdir(), "r_bin"))
  src <- file.path(dir, "figures", "panels", "source_data")
  # Add a third method, a third channel, a quoted name with a comma, and NA values.
  dist <- utils::read.csv(file.path(src, "unmixed-distributions.csv"), stringsAsFactors = FALSE)
  extra <- dist[dist$method == "ols", ]
  extra$method <- "autospectral_fast"
  extra$channel <- ifelse(extra$channel == "X-A", "Z, comma-A", "X-A")
  extra$value[1:5] <- NA
  dist <- rbind(dist, extra)
  utils::write.csv(dist, file.path(src, "unmixed-distributions.csv"), row.names = FALSE)
  bp <- utils::read.csv(file.path(src, "biplot.csv"), stringsAsFactors = FALSE)
  bp$x[3] <- NA
  bp <- rbind(bp, transform(bp[bp$method == "ols", ], method = "autospectral_fast", x = -x))
  utils::write.csv(bp, file.path(src, "biplot.csv"), row.names = FALSE)
  for (bins in c(10L, 37L)) {
    expect_identical(
      .sw_dash_bin_distributions(file.path(src, "unmixed-distributions.csv"), 1000, bins),
      old_bin_distributions(file.path(src, "unmixed-distributions.csv"), 1000, bins)
    )
    expect_identical(
      .sw_dash_bin_biplot(file.path(src, "biplot.csv"), 1000, bins),
      old_bin_biplot(file.path(src, "biplot.csv"), 1000, bins)
    )
  }
  expect_error(
    .sw_dash_read_columns(file.path(src, "biplot.csv"), list("character", NULL, "numeric")),
    "Unexpected columns"
  )
})

test_that("the biplot payload codes bins losslessly as key differences", {
  skip_if_not_installed("scales")
  run <- .sw_dash_read_run(make_run(file.path(tempdir(), "r_key")), sigma = 1000, bins = 20)
  p <- .sw_dash_biplot_payload(run, NULL, 20L)
  for (m in run$biplot$methods) {
    key <- cumsum(p$cells[[m]]$k) - 1L
    expect_identical(key %% 20L, run$biplot$cells[[m]]$x)
    expect_identical(key %/% 20L, run$biplot$cells[[m]]$y)
    expect_identical(p$cells[[m]]$n, run$biplot$cells[[m]]$n)
    expect_true(all(p$cells[[m]]$k >= 1L))
  }
  expect_identical(p$labels, c("x", "y"))
})

test_that("bin_events = FALSE skips the event tables", {
  run <- .sw_dash_read_run(make_run(file.path(tempdir(), "r3")), sigma = 1000, bins = 20, bin_events = FALSE)
  expect_null(run$distributions)
  expect_null(run$biplot)
})

test_that("the dashboard needs sigma and bins and refuses to overwrite", {
  skip_if_not_installed("jsonlite")
  dir <- make_run(file.path(tempdir(), "r4"))
  out <- tempfile(fileext = ".html")
  expect_error(sw_unmix_dashboard(dir, out, bins = 20), "pseudo_log_sigma")
  expect_error(sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000), "bins")
  expect_error(sw_unmix_dashboard(dir, out, pseudo_log_sigma = -1, bins = 20), "positive")
  expect_error(sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000, bins = 5), "at least 10")
  sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000, bins = 20)
  expect_true(file.exists(out))
  expect_error(sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000, bins = 20), "refusing to overwrite")
})

test_that("the HTML is self-contained and its data match the run", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("scales")
  dir <- make_run(file.path(tempdir(), "r5"))
  out <- tempfile(fileext = ".html")
  sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000, bins = 20, title = "T & <x>")
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # No resource loads: the only URL is the SVG namespace.
  urls <- regmatches(html, gregexpr("https?://[^\"'\\s<)]+", html, perl = TRUE))[[1L]]
  expect_identical(unique(urls), "http://www.w3.org/2000/svg")
  expect_false(grepl("<link|src=|@import", html))
  # Text with a script end tag in the data cannot close the script block.
  expect_equal(lengths(regmatches(html, gregexpr("</script>", html, fixed = TRUE))), 1L)
  expect_match(html, "<title>T &amp; &lt;x&gt;</title>", fixed = TRUE)

  p <- payload_of(out)
  expect_equal(p$settings$bins, 20)
  expect_equal(p$settings$pseudo_log_sigma, 1000)
  run <- p$runs[[1L]]
  expect_equal(run$limits$warn, 4)
  expect_equal(run$limits$fail, 6)
  expect_equal(unlist(run$sif$value), c(7, 5, 2))
  expect_identical(unlist(run$sif$status), c("fail", "warn", "pass"))
  expect_equal(run$status_counts$fail, 1)
  expect_equal(length(unlist(run$checks$status)), 4L)
  expect_true(all(c("ours", "spectracle", "ols", "spectroflo") %in% names(p$palette)))
  expect_null(p$compare)
  # The page decodes `k`; the old x and y arrays are not sent.
  expect_null(run$biplot$cells$ols$x)
  expect_equal(sum(unlist(run$biplot$cells$ols$n)), 200)
  expect_equal(length(run$biplot$cells$ols$k), length(run$biplot$cells$ols$n))
  # The what-if banner and the verdict source stay in the page.
  expect_match(html, "What-if cutoffs are display only", fixed = TRUE)
  expect_match(html, "Verdicts come from the run|saved verdicts in Checks")
})

test_that("two runs give a comparison and keep their labels", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("scales")
  a <- make_run(file.path(tempdir(), "ra"))
  b <- make_run(file.path(tempdir(), "rb"), sif = c(A = 3, B = 9, C = 2))
  out <- tempfile(fileext = ".html")
  sw_unmix_dashboard(c(cells = a, beads = b), out, pseudo_log_sigma = 1000, bins = 20)
  p <- payload_of(out)
  expect_identical(vapply(p$runs, function(r) r$name, ""), c("cells", "beads"))
  expect_length(p$compare, 1L)
  expect_identical(p$compare[[1L]]$a, "cells")
  expect_identical(p$compare[[1L]]$b, "beads")
  expect_equal(unlist(p$compare[[1L]]$angles$angle_deg), rep(0, 3))
  expect_error(
    sw_unmix_dashboard(c(a, a), tempfile(fileext = ".html"), pseudo_log_sigma = 1000, bins = 20),
    "unique"
  )
})

test_that("a run without figure tables still gives a dashboard", {
  skip_if_not_installed("jsonlite")
  dir <- make_run(file.path(tempdir(), "r6"))
  unlink(file.path(dir, "figures"), recursive = TRUE)
  out <- tempfile(fileext = ".html")
  sw_unmix_dashboard(dir, out, pseudo_log_sigma = 1000, bins = 20)
  run <- payload_of(out)$runs[[1L]]
  expect_null(run$distributions)
  expect_null(run$biplot)
  expect_null(run$similarity)
  expect_identical(run$spectra$methods, "ours")
})

test_that("an unknown method label is an error", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("scales")
  dir <- make_run(file.path(tempdir(), "r7"))
  p <- file.path(dir, "figures", "panels", "source_data", "biplot.csv")
  bp <- utils::read.csv(p)
  bp$method <- "mystery"
  utils::write.csv(bp, p, row.names = FALSE)
  expect_error(
    sw_unmix_dashboard(dir, tempfile(fileext = ".html"), pseudo_log_sigma = 1000, bins = 20),
    "Unknown method"
  )
})

test_that("a run folder without a report is refused", {
  d <- file.path(tempdir(), "empty_run")
  dir.create(d)
  expect_error(.sw_dash_read_run(d, 1000, 20), "report/")
  expect_error(.sw_dash_read_run(file.path(tempdir(), "nope"), 1000, 20), "existing run folder")
})

# ---- Shiny app ----------------------------------------------------------------

test_that("what-if status uses the same rule as the diagnostics", {
  # Above the cutoff is warn or fail; equal to it passes.
  expect_identical(
    .sw_app_status(c(7, 6, 5, 4, 2, NA), warn = 4, fail = 6),
    c("fail", "warn", "warn", "pass", "pass", "not_run")
  )
})

test_that("the app needs its event settings and builds from run folders", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("ggplot2")
  dir <- make_run(file.path(tempdir(), "r8"))
  expect_error(sw_unmix_app(dir, seed = 1, pseudo_log_sigma = 1000, bins = 20), "n_events")
  expect_error(sw_unmix_app(dir, n_events = 10, pseudo_log_sigma = 1000, bins = 20), "seed")
  expect_error(sw_unmix_app(dir, n_events = 10, seed = 1, bins = 20), "pseudo_log_sigma")
  expect_error(sw_unmix_app(dir, n_events = 10, seed = 1, pseudo_log_sigma = 1000), "bins")
  app <- sw_unmix_app(dir, n_events = 10, seed = 1, pseudo_log_sigma = 1000, bins = 20)
  expect_s3_class(app, "shiny.appobj")
})

test_that("the app server filters checks and recolours SIF without changing saved verdicts", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("ggplot2")
  a <- make_run(file.path(tempdir(), "r9a"))
  b <- make_run(file.path(tempdir(), "r9b"), sif = c(A = 3, B = 9, C = 2), limit = NA_character_)
  app <- sw_unmix_app(
    c(first = a, second = b),
    n_events = 10, seed = 1, pseudo_log_sigma = 1000, bins = 20
  )
  shiny::testServer(app, {
    session$setInputs(run = "first", statuses = c("pass", "warn", "fail", "not_judged", "not_run"), search = "")
    session$setInputs(sif_warn = 4, sif_fail = 6)
    d <- sif_data()
    expect_identical(d$fluorophore, c("A", "B", "C"))
    expect_identical(d$whatif, c("fail", "warn", "pass"))
    expect_identical(d$status, c("fail", "warn", "pass"))
    # Moving the what-if sliders changes only the what-if column.
    session$setInputs(sif_warn = 3, sif_fail = 8)
    d <- sif_data()
    expect_identical(d$whatif, c("warn", "warn", "pass"))
    expect_identical(d$status, c("fail", "warn", "pass"))
    # The second run saved no limits, so its saved status stays as written.
    session$setInputs(run = "second")
    expect_true(all(is.na(run()$limits)))
    expect_identical(run()$run_id, "id-r9b")
  })
})

test_that("common channels keep area channels present in every file", {
  m1 <- matrix(1, 1, 5, dimnames = list(NULL, c("FSC-A", "BV421-A", "AF-A", "BV421-H", "Time")))
  m2 <- matrix(1, 1, 3, dimnames = list(NULL, c("FSC-A", "BV421-A", "PE-A")))
  expect_identical(.sw_app_common_channels(list(m1, m2)), "BV421-A")
})

test_that("a run without a saved SIF check gets SIF computed from its spectra, unjudged", {
  skip_if_not_installed("scales")
  dir <- make_run(file.path(tempdir(), "r10"))
  ck <- utils::read.csv(file.path(dir, "report", "checks.csv"))
  utils::write.csv(ck[ck$metric != "max_sif", ], file.path(dir, "report", "checks.csv"), row.names = FALSE)
  run <- .sw_dash_read_run(dir, sigma = 1000, bins = 20, bin_events = FALSE)
  expect_identical(run$sif_source, "computed")
  expect_identical(run$sif$status, rep("not_judged", 3))
  expect_true(all(is.na(run$limits)))
  direct <- .sw_diag_sif(.sw_dash_spectra_matrix(run))
  expect_equal(run$sif$value, direct$value)
  # A run that saved the check keeps its own values and verdicts.
  run2 <- .sw_dash_read_run(make_run(file.path(tempdir(), "r11")), sigma = 1000, bins = 20, bin_events = FALSE)
  expect_identical(run2$sif_source, "run")
})

# ---- AutoSpectral's own outputs -------------------------------------------------

# An unmix/ folder shaped like AutoSpectral's: matrices and spectra as CSV, images
# and a PDF in figure_* folders, an rds that must not be read.
add_unmix_folder <- function(run_dir) {
  u <- file.path(run_dir, "unmix")
  for (d in c("figure_similarity_heatmap", "table_spectra", "figure_scatter", "figure_spectra", "figure_spectral_variants")) {
    dir.create(file.path(u, d), recursive = TRUE)
  }
  fl <- c("A", "B", "C")
  hot <- matrix(c(1.5, 0.4, 0.1, 0.4, 2.0, 0.3, 0.1, 0.3, 1.2), 3, dimnames = list(fl, fl))
  utils::write.csv(hot, file.path(u, "figure_similarity_heatmap", "Automatedhotspot_matrix.csv"))
  spill <- hot
  diag(spill) <- NA
  utils::write.csv(spill, file.path(u, "figure_similarity_heatmap", "autospectral_spillover_spread_matrix.csv"))
  big <- matrix(1, 200, 2, dimnames = list(paste0("r", 1:200), c("x", "y")))
  utils::write.csv(big, file.path(u, "figure_similarity_heatmap", "too_big_matrix.csv"))
  det <- paste0("D", 1:4)
  af <- matrix(runif(8), 2, dimnames = list(c("AF1", "AF2"), det))
  utils::write.csv(af, file.path(u, "table_spectra", "Unstained (Cells).fcs Autofluorescence spectra 2026.csv"))
  sp <- matrix(runif(12), 3, dimnames = list(fl, det))
  utils::write.csv(sp, file.path(u, "table_spectra", "Automated_autospectral_spectra_2026.csv"))
  utils::write.csv(sp[1, , drop = FALSE], file.path(u, "table_spectra", "Clean_autospectral_spectra_2026.csv"))
  file.create(file.path(u, "figure_scatter", "A universal (neg).jpg"))
  file.create(file.path(u, "figure_spectra", "qc report.pdf"))
  file.create(file.path(u, "figure_spectral_variants", "Spectral_variants.rds"))
  run_dir
}

test_that("a relative path walks up and down the folder tree", {
  base <- file.path(tempdir(), "rel")
  dir.create(file.path(base, "a", "b"), recursive = TRUE)
  dir.create(file.path(base, "c"), recursive = TRUE)
  file.create(file.path(base, "c", "x y.jpg"))
  expect_identical(
    .sw_dash_relpath(file.path(base, "a", "b"), file.path(base, "c", "x y.jpg")),
    "../../c/x y.jpg"
  )
  expect_identical(.sw_dash_relpath(base, file.path(base, "c", "x y.jpg")), "c/x y.jpg")
})

test_that("AutoSpectral tables and files are read, labelled and linked relative to the page", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("scales")
  run_dir <- add_unmix_folder(make_run(file.path(tempdir(), "r12")))
  out_dir <- file.path(tempdir(), "pages")
  dir.create(out_dir)
  out <- file.path(out_dir, "d.html")
  sw_unmix_dashboard(run_dir, out, pseudo_log_sigma = 1000, bins = 20)
  a <- payload_of(out)$runs[[1L]]$autospectral
  labels <- vapply(a$matrices, function(m) m$label, "")
  expect_true("Automated hotspot matrix" %in% labels)
  expect_true("autospectral spillover spread matrix" %in% labels)
  # A matrix larger than the display limit is left out, not cut.
  expect_false(any(grepl("too big", labels)))
  hot <- a$matrices[[which(labels == "Automated hotspot matrix")]]
  expect_equal(unlist(hot$values), c(1.5, 0.4, 0.1, 0.4, 2, 0.3, 0.1, 0.3, 1.2))
  expect_identical(unlist(hot$rows), c("A", "B", "C"))
  spill <- a$matrices[[which(labels == "autospectral spillover spread matrix")]]
  expect_null(spill$values[[1L]])  # the NA diagonal is null, not 0
  expect_identical(a$af[[1L]]$names, list("AF1", "AF2"))
  expect_identical(a$af[[1L]]$label, "Unstained (Cells).fcs")
  stages <- vapply(a$versions, function(v) v$label, "")
  expect_setequal(stages, c("Automated", "Clean"))
  # Gallery: images and PDF only, grouped by folder, links resolve from the page.
  g <- a$gallery
  expect_setequal(unlist(g$name), c("A universal (neg).jpg", "qc report.pdf"))
  expect_false(any(grepl("rds", unlist(g$name))))
  for (h in unlist(g$href)) {
    expect_true(file.exists(file.path(out_dir, utils::URLdecode(h))))
  }
  expect_true(all(grepl("%20", unlist(g$href))))
})

test_that("a run without an unmix folder has no AutoSpectral section", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("scales")
  out <- tempfile(fileext = ".html")
  sw_unmix_dashboard(make_run(file.path(tempdir(), "r13")), out, pseudo_log_sigma = 1000, bins = 20)
  expect_null(payload_of(out)$runs[[1L]]$autospectral)
  expect_error(
    sw_unmix_dashboard(make_run(file.path(tempdir(), "r14")), file.path(tempdir(), "nodir", "d.html"), pseudo_log_sigma = 1000, bins = 20),
    "folder of `path`"
  )
})
