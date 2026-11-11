# Tests for R/unmix_figures_spike.R (F6 hexbin grid, F7 and F8 spike panels).

skip_if_not_installed("ggplot2")
skip_if_not_installed("hexbin")

# ---- Mock data ---------------------------------------------------------------

# Two setups, two donors, linear values; TNF is high in setup b.
mock_grid_events <- function() {
  set.seed(11)
  n <- 300
  tibble::tibble(
    setup = rep(c("a", "b"), each = 2 * n),
    sample = rep(rep(c("d1", "d2"), each = n), 2),
    GLUT1 = stats::rnorm(4 * n, 8000, 6000),
    TNF = c(stats::rnorm(2 * n, 0, 3000), stats::rnorm(2 * n, 12000, 3000))
  )
}

grid_fill <- c(a = "#0072B2", b = "#D55E00")

# A spike score: TNF negatives (base unmix) spiked at level 2. Each cell either
# keeps most of the spike or loses most of it; the lost share differs by setup
# and donor, so the per-panel % lost differs.
mock_spike_score <- function(
  setups = c("as188", "pooled"),
  donors = c("d1", "d2", "d3")
) {
  set.seed(5)
  k <- 6000 * sinh(2)
  p_lost <- outer(
    seq(0.05, 0.4, length.out = length(setups)),
    seq(0, 0.1, length.out = length(donors)),
    `+`
  )
  dimnames(p_lost) <- list(setups, donors)
  cells <- expand.grid(
    i = seq_len(400),
    sample = donors,
    setup = setups,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      lost = stats::runif(dplyr::n()) <
        p_lost[cbind(.data$setup, .data$sample)],
      share = dplyr::if_else(
        .data$lost,
        stats::rnorm(dplyr::n(), 0.15, 0.1),
        stats::rnorm(dplyr::n(), 0.95, 0.1)
      ),
      set = "s1",
      pop = "CD8T",
      raw_row = .data$i,
      marker = "TNF",
      level = 2,
      k = k,
      base_value = stats::rnorm(dplyr::n(), 0, 3000),
      recovered = .data$share * k,
      recovery = .data$recovered / .data$k
    ) |>
    dplyr::select(-c("i", "lost", "share"))
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

test_that("F6 spike mode plots one marker before vs after and counts lost cells", {
  score <- mock_spike_score()
  args <- list(
    events = score,
    x = "TNF",
    mode = "spike",
    cofactor = 6000,
    setup_fill = c(as188 = "#0072B2", pooled = "#D55E00"),
    stem = "spike-grid_TNF"
  )
  expect_error(
    do.call(sw_unmix_figure_scatter_grid, args),
    "`level` is required"
  )
  f6 <- do.call(sw_unmix_figure_scatter_grid, c(args, level = 2))

  # % lost per panel equals a hand count, and differs across panels.
  cells <- score$cells
  hand <- 100 *
    mean(cells$recovery[cells$setup == "pooled" & cells$sample == "d3"] < 0.5)
  got <- f6$source$pct[f6$source$setup == "pooled" & f6$source$sample == "d3"]
  expect_equal(got, hand)
  expect_gt(length(unique(round(f6$source$pct, 1))), 3L)

  # Axes are the marker before and after: x = base unmix, y = base + recovered.
  d <- f6$plot$data
  pick <- cells$setup == "as188" & cells$sample == "d1"
  expect_equal(
    d$x[d$setup == "as188" & d$sample == "d1"],
    asinh(cells$base_value[pick] / 6000)
  )
  expect_equal(
    d$y[d$setup == "as188" & d$sample == "d1"],
    asinh((cells$base_value[pick] + cells$recovered[pick]) / 6000)
  )
  expect_equal(f6$plot$labels$x, "TNF, unspiked unmix (arcsinh, cofactor 6000)")
  expect_equal(f6$plot$labels$y, "TNF, +2 arcsinh spiked unmix")

  # Wrong inputs fail: a plain data frame, a second channel, a mismatched cofactor.
  expect_error(
    sw_unmix_figure_scatter_grid(
      mock_grid_events(),
      "TNF",
      mode = "spike",
      cofactor = 6000,
      setup_fill = grid_fill,
      stem = "x",
      level = 2
    ),
    "must be an sw_spike_score"
  )
  expect_error(
    do.call(sw_unmix_figure_scatter_grid, c(args, y = "GLUT1", level = 2)),
    "`y` is not used in spike mode"
  )
  expect_error(
    do.call(
      sw_unmix_figure_scatter_grid,
      c(utils::modifyList(args, list(cofactor = 150)), level = 2)
    ),
    "does not equal cofactor"
  )
})

test_that("F6 height follows the grid of square panels", {
  one_row <- mock_grid_events() |>
    dplyr::filter(.data$sample == "d1")
  f_one <- sw_unmix_figure_scatter_grid(
    one_row,
    "GLUT1",
    "TNF",
    mode = "scatter",
    cut = 1,
    cofactor = 6000,
    setup_fill = grid_fill,
    stem = "x"
  )
  f_two <- sw_unmix_figure_scatter_grid(
    mock_grid_events(),
    "GLUT1",
    "TNF",
    mode = "scatter",
    cut = 1,
    cofactor = 6000,
    setup_fill = grid_fill,
    stem = "x"
  )
  expect_equal(f_one$size[1], 183)
  expect_lt(f_one$size[2], f_two$size[2])
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
