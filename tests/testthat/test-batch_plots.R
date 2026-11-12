# tests/testthat/test-batch_plots.R
# Unit tests for R/batch_plots.R on a mock batch matrix.

make_matrix <- function(population = FALSE) {
  s <- tibble::tibble(
    check = "rest_emd",
    marker = rep(c("A", "B", "C"), each = 3),
    batch = rep(c("b1", "b2", "b3"), 3),
    value = c(0.1, 0.9, 0.2, 0.3, 0.31, 0.29, 0.1, 0.1, 0.5)
  )
  if (population) {
    s <- rbind(transform(s, population = "T"), transform(s, population = "B"))
  }
  sw_batch_matrix(s, "rest_emd", c(rest_emd = "high"))
}

test_that("plot returns plot and source data, markers ordered by max value", {
  skip_if_not_installed("ggplot2")
  fig <- sw_plot_batch_matrix(
    make_matrix(),
    "rest_emd",
    n_markers = NULL,
    population = NULL
  )
  expect_named(fig, c("plot", "source_data"))
  expect_s3_class(fig$plot, "ggplot")
  expect_equal(nrow(fig$source_data), 9)
  # Top row of the plot is the marker with the largest value (A, 0.9).
  expect_equal(rev(levels(fig$source_data$marker)), c("A", "C", "B"))
  # Largest dot is rank 1.
  sd <- fig$source_data
  expect_equal(
    sd$batch[sd$marker == "A" & sd$rank_size == max(sd$rank_size)],
    "b2"
  )
  expect_match(fig$plot$labels$caption, "No cutoff drawn")
  expect_no_error(ggplot2::ggplot_build(fig$plot))
})

test_that("n_markers keeps the first markers in plot order", {
  skip_if_not_installed("ggplot2")
  fig <- sw_plot_batch_matrix(
    make_matrix(),
    "rest_emd",
    n_markers = 2,
    population = NULL
  )
  expect_setequal(as.character(unique(fig$source_data$marker)), c("A", "C"))
})

test_that("population is required when the matrix has one", {
  skip_if_not_installed("ggplot2")
  m <- make_matrix(population = TRUE)
  expect_error(
    sw_plot_batch_matrix(m, "rest_emd", n_markers = NULL, population = NULL),
    "population"
  )
  fig <- sw_plot_batch_matrix(m, "rest_emd", n_markers = NULL, population = "T")
  expect_equal(nrow(fig$source_data), 9)
})

test_that("errors on missing arguments and wrong inputs", {
  skip_if_not_installed("ggplot2")
  m <- make_matrix()
  expect_error(sw_plot_batch_matrix(m), "Required argument")
  expect_error(
    sw_plot_batch_matrix(m, "rest_emd", n_markers = NULL),
    "'population'"
  )
  expect_error(
    sw_plot_batch_matrix(m, "nope", n_markers = NULL, population = NULL),
    "present"
  )
  expect_error(
    sw_plot_batch_matrix(
      as.data.frame(m),
      "rest_emd",
      n_markers = NULL,
      population = NULL
    ),
    "sw_batch_matrix"
  )
})
