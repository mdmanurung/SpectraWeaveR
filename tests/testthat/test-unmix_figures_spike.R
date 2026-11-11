# Tests for R/unmix_figures_spike.R (F6 hexbin grid, F7 and F8 spike panels).

skip_if_not_installed("ggplot2")
skip_if_not_installed("hexbin")

# ---- Mock data ---------------------------------------------------------------

# Two setups, two donors, linear values; TNF is high in setup b.
mock_grid_events <- function(spike = FALSE) {
  set.seed(11)
  n <- 300
  ev <- tibble::tibble(
    setup = rep(c("a", "b"), each = 2 * n),
    sample = rep(rep(c("d1", "d2"), each = n), 2),
    GLUT1 = stats::rnorm(4 * n, 8000, 6000),
    TNF = c(stats::rnorm(2 * n, 0, 3000), stats::rnorm(2 * n, 12000, 3000))
  )
  if (spike) {
    ev$lost <- rep(c(TRUE, FALSE, FALSE, FALSE), length.out = nrow(ev))
  }
  ev
}

grid_fill <- c(a = "#0072B2", b = "#D55E00")

# A spike score with chosen setups: cells and per-donor summaries at one level.
mock_spike_score <- function(
  setups = c("as188", "pooled"),
  donors = c("d1", "d2", "d3")
) {
  set.seed(5)
  cells <- expand.grid(
    i = seq_len(200),
    sample = donors,
    setup = setups,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      set = "s1",
      pop = "CD8T",
      raw_row = .data$i,
      marker = "TNF",
      level = 2,
      k = 6000 * sinh(2),
      base_value = 0,
      recovered = 0,
      recovery = stats::rnorm(dplyr::n(), 0.9, 0.3)
    ) |>
    dplyr::select(-"i")
  summary <- cells |>
    dplyr::summarise(
      n = dplyr::n(),
      share_lost = mean(.data$recovery < 0.5),
      median_recovery = stats::median(.data$recovery),
      q05 = unname(stats::quantile(.data$recovery, 0.05)),
      q25 = unname(stats::quantile(.data$recovery, 0.25)),
      q50 = unname(stats::quantile(.data$recovery, 0.50)),
      q75 = unname(stats::quantile(.data$recovery, 0.75)),
      q95 = unname(stats::quantile(.data$recovery, 0.95)),
      .by = c("setup", "set", "marker", "pop", "sample", "level")
    )
  structure(
    list(cells = cells, summary = summary, loss_cut = 0.5),
    class = "sw_spike_score"
  )
}

# ---- F6 ----------------------------------------------------------------------

test_that("F6 scatter mode returns a ggplot with the share above the cut", {
  ev <- mock_grid_events()
  f6 <- sw_unmix_figure_scatter_grid(
    ev,
    "GLUT1",
    "TNF",
    mode = "scatter",
    cut = 1,
    cofactor = 6000,
    setup_fill = grid_fill,
    stem = "scatter-grid_TNF"
  )
  expect_s3_class(f6$plot, "ggplot")
  expect_equal(f6$panel, "F6")
  expect_equal(f6$stem, "scatter-grid_TNF")
  hand <- 100 *
    mean(asinh(ev$TNF[ev$setup == "b" & ev$sample == "d2"] / 6000) > 1)
  got <- f6$source$pct[f6$source$setup == "b" & f6$source$sample == "d2"]
  expect_equal(got, hand)
  expect_no_error(ggplot2::ggplot_build(f6$plot))
})

test_that("F6 spike mode needs lost and level, and counts lost cells", {
  ev <- mock_grid_events(spike = TRUE)
  args <- list(
    events = ev,
    x = "GLUT1",
    y = "TNF",
    mode = "spike",
    cofactor = 6000,
    setup_fill = grid_fill,
    stem = "spike-grid_TNF"
  )
  expect_error(
    do.call(sw_unmix_figure_scatter_grid, args),
    "`level` is required"
  )
  f6 <- do.call(sw_unmix_figure_scatter_grid, c(args, level = 2))
  expect_equal(f6$source$pct, rep(25, 4))
  args$events$lost <- NULL
  expect_error(
    do.call(sw_unmix_figure_scatter_grid, c(args, level = 2)),
    "lacks column\\(s\\): lost"
  )
})

test_that("F6 errors on a setup without a colour and on a missing argument", {
  ev <- mock_grid_events()
  expect_error(
    sw_unmix_figure_scatter_grid(
      ev,
      "GLUT1",
      "TNF",
      mode = "scatter",
      cut = 1,
      cofactor = 6000,
      setup_fill = grid_fill["a"],
      stem = "x"
    ),
    "without a colour in `setup_fill`: b"
  )
  expect_error(
    sw_unmix_figure_scatter_grid(
      ev,
      "GLUT1",
      "TNF",
      mode = "scatter",
      cofactor = 6000,
      setup_fill = grid_fill,
      stem = "x"
    ),
    "`cut` is required"
  )
})

# ---- F7 and F8 ---------------------------------------------------------------

test_that("F7 plots loss and recovery per donor for every setup", {
  score <- mock_spike_score()
  f7 <- sw_unmix_figure_spike_summary(
    score,
    readout_level = 2,
    stem = "spike-summary_TNF"
  )
  expect_s3_class(f7$plot, "ggplot")
  expect_setequal(unique(f7$source$setup), c("as188", "pooled"))
  lost <- f7$source[grepl("losing", f7$source$readout), ]
  hand <- 100 *
    mean(
      score$cells$recovery[
        score$cells$setup == "pooled" & score$cells$sample == "d1"
      ] <
        0.5
    )
  expect_equal(lost$value[lost$setup == "pooled" & lost$sample == "d1"], hand)
  expect_error(
    sw_unmix_figure_spike_summary(score, readout_level = 9, stem = "x"),
    "No spiked cells at level 9"
  )
  expect_error(
    sw_unmix_figure_spike_summary(score, stem = "x"),
    "`readout_level` is required"
  )
})

test_that("F8 draws one curve per setup and errors when a setup lacks a linetype", {
  score <- mock_spike_score()
  f8 <- sw_unmix_figure_spike_ecdf(
    score,
    readout_level = 2,
    stem = "spike-ecdf_TNF"
  )
  expect_s3_class(f8$plot, "ggplot")
  built <- ggplot2::ggplot_build(f8$plot)
  expect_equal(length(unique(built$data[[2]]$group)), 2L)
  at_cut <- f8$source[
    f8$source$setup == "pooled" & abs(f8$source$recovery - 0.5) < 1e-9,
  ]
  hand <- mean(score$cells$recovery[score$cells$setup == "pooled"] <= 0.5)
  expect_equal(at_cut$cum_share, hand)
  too_many <- mock_spike_score(setups = paste0("s", 1:7))
  expect_error(
    sw_unmix_figure_spike_ecdf(too_many, readout_level = 2, stem = "x"),
    "Setup\\(s\\) without a linetype: s7"
  )
})

# ---- Saving ------------------------------------------------------------------

test_that("F6 to F8 pass the Nature check and save under their stems", {
  skip_if_not(capabilities("cairo"), "cairo device not available")
  ev <- mock_grid_events()
  score <- mock_spike_score()
  figs <- list(
    F6 = sw_unmix_figure_scatter_grid(
      ev,
      "GLUT1",
      "TNF",
      mode = "scatter",
      cut = 1,
      cofactor = 6000,
      setup_fill = grid_fill,
      stem = "scatter-grid_TNF"
    ),
    F7 = sw_unmix_figure_spike_summary(
      score,
      readout_level = 2,
      stem = "spike-summary_TNF"
    ),
    F8 = sw_unmix_figure_spike_ecdf(
      score,
      readout_level = 2,
      stem = "spike-ecdf_TNF"
    )
  )
  out <- withr::local_tempdir()
  expect_no_error(sw_unmix_save_figures(figs, out, check = "nature"))
  files <- list.files(out, recursive = TRUE)
  for (stem in c("scatter-grid_TNF", "spike-summary_TNF", "spike-ecdf_TNF")) {
    expect_true(any(grepl(paste0(stem, "\\.png$"), files)), info = stem)
  }
})
