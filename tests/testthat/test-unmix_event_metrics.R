# Unit tests for R/unmix_event_metrics.R: edges, event metrics, held-out
# wrapper and the diagnose/report hooks. AutoSpectral is mocked.

sim_unstained <- function(n = 4000, ch = c("BV421", "PE")) {
  set.seed(1)
  matrix(rnorm(n * length(ch), 0, 300), ncol = length(ch), dimnames = list(NULL, ch))
}

em_args <- function(...) {
  utils::modifyList(
    list(
      cofactor = 6000, negative_cut = -1.5, spread_probs = c(0.025, 0.975),
      min_pos = 100, hyper_k = 3, spread = NULL, kappa = NULL,
      pairs = NULL, n_bins = NULL, min_cells = NULL
    ),
    list(...)
  )
}

run_em <- function(data, edges, ...) {
  do.call(sw_unmix_event_metrics, c(list(data = data, edges = edges), em_args(...)))
}

# ---- sw_unmix_edges -------------------------------------------------------

test_that("edges from an unstained matrix are its quantiles", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  expect_equal(e$channel, c("BV421", "PE"))
  expect_equal(e$edge[[1]], unname(quantile(un[, 1], 0.995)))
  expect_equal(e$floor[[2]], unname(quantile(un[, 2], 0.005)))
})

test_that("edges from variants use AutoSpectral's stored thresholds", {
  v <- list(thresholds = c(BV421 = 900, PE = 1000), neg.thresholds = c(PE = -950, BV421 = -880))
  e <- sw_unmix_edges(variants = v)
  expect_equal(e$floor, c(-880, -950))
  expect_equal(e$source, rep("variants", 2))
})

test_that("edges need exactly one source and named probs", {
  un <- sim_unstained()
  expect_error(sw_unmix_edges(), "exactly one")
  expect_error(sw_unmix_edges(unstained = un, probs = c(0.005, 0.995)), "floor")
  expect_error(sw_unmix_edges(variants = list(thresholds = c(a = 1)), probs = c(floor = 0.1, edge = 0.9)), "only with")
  expect_error(sw_unmix_edges(variants = list(thresholds = c(a = 1), neg.thresholds = c(a = NA_real_))), "NA")
})

# ---- sw_unmix_event_metrics ---------------------------------------------------

test_that("every knob is required and missing ones are listed", {
  e <- sw_unmix_edges(unstained = sim_unstained(), probs = c(floor = 0.005, edge = 0.995))
  expect_error(sw_unmix_event_metrics(list(a = sim_unstained()), e), "missing: cofactor, negative_cut")
})

test_that("negatives, spread and stain index follow their definitions", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  set.seed(2)
  s1 <- rbind(un, cbind(rnorm(500, 4e4, 8e3), rnorm(500, 0, 300)))
  em <- run_em(list(s1 = s1), e)
  expect_s3_class(em, "sw_event_metrics")
  r <- em$channels[em$channels$channel == "BV421", ]
  a <- asinh(s1[, 1] / 6000)
  neg <- a[s1[, 1] < e$edge[[1]]]
  pos <- a[s1[, 1] >= e$edge[[1]]]
  expect_equal(r$pct_below_cut, 100 * mean(a < -1.5))
  expect_equal(r$pct_below_floor, 100 * mean(s1[, 1] < e$floor[[1]]))
  sq <- unname(diff(quantile(neg, c(0.025, 0.975)))) / (qnorm(0.975) - qnorm(0.025))
  expect_equal(r$spread_q, sq)
  expect_equal(r$pileup, mad(neg) / sq)
  expect_equal(r$si_mad, (median(pos) - median(neg)) / (2 * mad(neg)))
  expect_equal(r$n_pos, length(pos))
  # PE has too few positives for a stain index
  expect_true(is.na(em$channels$si_mad[em$channels$channel == "PE"]))
})

test_that("a heavy negative tail raises pct_below_cut and hyper_share", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  set.seed(3)
  bad <- un
  bad[1:400, 2] <- -2e4
  em <- run_em(list(good = un, bad = bad), e)
  pe <- em$channels[em$channels$channel == "PE", ]
  expect_gt(pe$pct_below_cut[pe$file == "bad"], pe$pct_below_cut[pe$file == "good"])
  expect_gt(pe$hyper_share[pe$file == "bad"], pe$hyper_share[pe$file == "good"])
})

test_that("the line check finds a leaning population and a flat one", {
  un <- sim_unstained(ch = c("A", "B", "C"))
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  set.seed(4)
  xa <- rnorm(3000, 3e4, 1.5e4)
  lean <- cbind(A = xa, B = -0.02 * xa + rnorm(3000, 0, 100), C = rnorm(3000, 0, 300))
  em <- run_em(list(s = rbind(un, lean)), e, pairs = "all", n_bins = 10, min_cells = 100)
  ab <- em$pairs[em$pairs$channel_i == "A" & em$pairs$channel_j == "B", ]
  ac <- em$pairs[em$pairs$channel_i == "A" & em$pairs$channel_j == "C", ]
  expect_lt(ab$slope, -0.05)
  expect_lt(ab$rho, -0.5)
  expect_lt(abs(ac$slope), 0.02)
  expect_equal(nrow(em$pairs), 6)
  # B has no positives, so B -> A is below min_cells and NA
  expect_true(is.na(em$pairs$slope[em$pairs$channel_i == "B" & em$pairs$channel_j == "A"]))
})

test_that("spread-aware edges call AutoSpectral's thresholds", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  spread <- matrix(0.1, 2, 2, dimnames = list(e$channel, e$channel))
  sides <- character()
  local_mocked_bindings(.sw_em_spread_edges = function(x, thresholds, spread, kappa, side) {
    sides <<- c(sides, side)
    matrix(thresholds * 10, nrow(x), ncol(x), byrow = TRUE, dimnames = dimnames(x))
  })
  em <- run_em(list(s = un), e, spread = spread, kappa = 2)
  expect_equal(sides, c("upper", "lower"))
  expect_equal(em$settings$edge_type, "spread-aware")
  expect_true(all(em$channels$n_pos == 0))
  expect_error(run_em(list(s = un), e, spread = spread), "together")
})

test_that("unmatched channels are reported, not dropped silently", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  only_pe <- un[, "PE", drop = FALSE]
  colnames(only_pe) <- "PE-A"
  em <- run_em(list(s = only_pe), e)
  expect_equal(em$missing$channel, "BV421")
  expect_equal(em$channels$column, "PE-A")
  expect_error(run_em(list(s = matrix(1, 2, 1, dimnames = list(NULL, "X"))), e), "no column")
})

test_that("pairs arguments are checked", {
  e <- sw_unmix_edges(unstained = sim_unstained(), probs = c(floor = 0.005, edge = 0.995))
  d <- list(s = sim_unstained())
  expect_error(run_em(d, e, n_bins = 10), "must be NULL")
  expect_error(run_em(d, e, pairs = "some", n_bins = 10, min_cells = 100), "pairs")
  expect_error(run_em(d, e, pairs = "all", n_bins = 10, min_cells = 5), "at least")
  expect_error(run_em(d, e, pairs = data.frame(i = "BV421", j = "X"), n_bins = 10, min_cells = 100), "not in")
})

# ---- Diagnose and report ------------------------------------------------------

test_that("diagnose adds event rows and judges them only with a limit", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  em <- run_em(list(s = un), e, pairs = "all", n_bins = 5, min_cells = 10)
  sp <- matrix(c(1, 0.2, 0.1, 1), 2, dimnames = list(c("BV421", "PE"), c("D1", "D2")))
  res <- sw_unmix_result(spectra = sp, run_id = "t")
  dx <- sw_unmix_diagnose(res, event_metrics = em)
  ev <- dx$checks[dx$checks$check_id == "event_channel", ]
  expect_true(all(ev$status == "not_judged"))
  expect_setequal(unique(ev$metric), c("pct_below_cut", "pct_below_floor", "neg_median", "neg_rsd", "spread_q", "pileup", "hyper_share", "si_mad", "si_q"))
  expect_true(any(dx$checks$check_id == "event_pair"))
  dx2 <- sw_unmix_diagnose(res, event_metrics = em, limits = list(max_pct_below_cut = c(warn = -1, fail = 100)))
  cut_rows <- dx2$checks[dx2$checks$metric == "pct_below_cut", ]
  expect_true(all(cut_rows$status %in% c("warn", "pass")))
  expect_error(sw_unmix_diagnose(res, event_metrics = list()), "sw_event_metrics")
})

test_that("report writes event and held-out files", {
  un <- sim_unstained()
  e <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
  em <- run_em(list(s = un), e, pairs = "all", n_bins = 5, min_cells = 10)
  ho <- structure(
    list(
      scores = tibble::tibble(control = "PE.fcs", type = "cells", target = "PE-A", channel = "BV421-A", n = 10L, spread_q = 1, bias = 0.5),
      repro = tibble::tibble(object = "spectra", identical = TRUE, max_abs_diff = 0),
      split = tibble::tibble(), settings = list(train_fraction = 0.5, seed = 42, pos_quantile = 0.995)
    ),
    class = "sw_heldout"
  )
  sp <- matrix(c(1, 0.2, 0.1, 1), 2, dimnames = list(c("BV421", "PE"), c("D1", "D2")))
  res <- sw_unmix_result(spectra = sp, run_id = "t")
  dx <- sw_unmix_diagnose(res, event_metrics = em, heldout = ho)
  expect_true(all(c("heldout_spread_q", "heldout_bias") %in% dx$checks$metric))
  path <- file.path(tempdir(), paste0("em_report_", Sys.getpid()))
  on.exit(unlink(path, recursive = TRUE))
  out <- sw_unmix_report(res, dx, path = path)
  expect_true(all(file.exists(file.path(path, c("event_channels.csv", "event_pairs.csv", "heldout_scores.csv", "heldout_repro.csv")))))
  html <- readLines(file.path(path, "report.html"))
  expect_true(any(grepl("Event metrics", html)))
  expect_true(any(grepl("Held-out control error", html)))
})

# ---- sw_unmix_heldout ------------------------------------------------------------

test_that("heldout lists every missing argument", {
  expect_error(sw_unmix_heldout(control_dir = "x"), "missing: control_file")
})

test_that("heldout splits, extracts, unmixes and scores with mocked AutoSpectral", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("AutoSpectral")
  base <- file.path(tempdir(), paste0("ho_", Sys.getpid()))
  on.exit(unlink(base, recursive = TRUE))
  cdir <- file.path(base, "ctrl")
  dir.create(cdir, recursive = TRUE)
  set.seed(5)
  mk <- \(m) flowCore::flowFrame(m)
  cols <- c("FSC-A", "SSC-A", "D1", "D2")
  neg <- cbind(rnorm(400, 1e6, 1e5), rnorm(400, 5e5, 5e4), rnorm(400, 0, 50), rnorm(400, 0, 50))
  pe <- cbind(rnorm(400, 1e6, 1e5), rnorm(400, 5e5, 5e4), rnorm(400, 0, 50), rnorm(400, 5e3, 500))
  colnames(neg) <- colnames(pe) <- cols
  flowCore::write.FCS(mk(neg), file.path(cdir, "Neg.fcs"))
  flowCore::write.FCS(mk(pe), file.path(cdir, "PE.fcs"))
  flowCore::write.FCS(mk(neg), file.path(base, "unstained.fcs"))
  cf <- file.path(cdir, "fcs_control_file.csv")
  utils::write.csv(
    data.frame(filename = c("Neg.fcs", "PE.fcs"), fluorophore = c("Negative", "PE"), control.type = "cells", universal.negative = "Neg.fcs"),
    cf, row.names = FALSE
  )
  calls <- character()
  local_mocked_bindings(.sw_ho_call = function(fn, args) {
    calls <<- c(calls, fn)
    switch(fn,
      get.autospectral.param = list(figure.dir = "fig"),
      get.spectra.automated = matrix(c(1, 0, 0, 1), 2, dimnames = list(c("BV421", "PE"), c("D1", "D2"))),
      get.af.spectra = matrix(c(0.5, 0.5), 1, dimnames = list("AF1", c("D1", "D2"))),
      get.spectral.variants = list(thresholds = c(BV421 = 1, PE = 1)),
      reload.flow.control = list(),
      unmix.folder = {
        for (f in list.files(args$fcs.dir, pattern = "\\.fcs$", full.names = TRUE)) {
          m <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
          um <- cbind(m[, 1:2], `BV421-A` = m[, "D1"], `PE-A` = m[, "D2"])
          flowCore::write.FCS(mk(um), file.path(args$output.dir, paste0(tools::file_path_sans_ext(basename(f)), " AutoSpectral.fcs")))
        }
        invisible(NULL)
      }
    )
  })
  ho <- sw_unmix_heldout(
    control_dir = cdir, control_file = cf, unstained_file = file.path(base, "unstained.fcs"),
    out_dir = file.path(base, "out"), cytometer = "aurora", seed = 42, train_fraction = 0.5,
    pos_quantile = 0.995, spread_probs = c(0.025, 0.975),
    cell_gate = list(`FSC-A` = c(5e5, 2e6), `SSC-A` = c(1e5, 1e6)),
    extract_args = list(), af_args = list(refine = TRUE), variant_args = list(),
    unmix_args = list(method = "AutoSpectral"), repro = TRUE
  )
  expect_s3_class(ho, "sw_heldout")
  expect_equal(sort(unique(ho$scores$control)), c("Neg.fcs", "PE.fcs", "unstained__unstained.fcs"))
  pe_row <- ho$scores[ho$scores$control == "PE.fcs", ]
  expect_equal(pe_row$target, "PE-A")
  expect_equal(pe_row$channel, "BV421-A")
  expect_gt(pe_row$n, 100)
  expect_true(all(ho$repro$identical[1:2]))
  expect_equal(sum(calls == "get.spectra.automated"), 3)
  expect_equal(ho$split$train + ho$split$test, ho$split$events)
  expect_error(
    sw_unmix_heldout(
      control_dir = cdir, control_file = cf, unstained_file = file.path(base, "unstained.fcs"),
      out_dir = file.path(base, "out"), cytometer = "aurora", seed = 42, train_fraction = 0.5,
      pos_quantile = 0.995, spread_probs = c(0.025, 0.975), cell_gate = NULL,
      extract_args = list(), af_args = list(), variant_args = list(), unmix_args = list(), repro = FALSE
    ),
    "already exists"
  )
})
