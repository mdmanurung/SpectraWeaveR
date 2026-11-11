test_that("sw_gate_plot returns a ggplot for polygon and range gates", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  d <- list(
    a = cbind(A = rnorm(500), B = rnorm(500)),
    b = cbind(A = rnorm(500, 1), B = rnorm(500))
  )
  poly <- sw_gate_polygon("pop", c("A", "B"), rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1)), parent = NA)
  adapted <- sw_gate_polygon("pop", c("A", "B"), rbind(c(0, -1), c(2, -1), c(2, 1), c(0, 1)), parent = NA)
  p <- sw_gate_plot(d, list(a = list(poly), b = list(adapted)), samples = c("a", "b"),
    gate = "pop", n_events = 200, seed = 1, reference_gates = list(poly))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  # event cap: 200 per sample drawn
  expect_equal(nrow(p$layers[[1]]$data), 400)

  rg <- sw_gate_range("r", "A", -1, Inf, parent = NA)
  p2 <- sw_gate_plot(d, list(rg), samples = "a", gate = "r", n_events = 100, seed = 1,
    reference_gates = NULL)
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("sw_gate_plot restores the RNG state and checks input", {
  skip_if_not_installed("ggplot2")
  d <- list(a = cbind(A = rnorm(50), B = rnorm(50)))
  g <- sw_gate_polygon("pop", c("A", "B"), rbind(c(-1, -1), c(1, -1), c(1, 1)), parent = NA)
  set.seed(5)
  before <- .Random.seed
  sw_gate_plot(d, list(g), "a", "pop", n_events = 10, seed = 99, reference_gates = NULL)
  expect_identical(.Random.seed, before)
  expect_error(sw_gate_plot(d, list(g), "zz", "pop", 10, 1, NULL), "samples")
  expect_error(sw_gate_plot(d, list(g), "a", "nope", 10, 1, NULL), "not in")
  expect_error(sw_gate_plot(d, list(g), "a", "pop", 0, 1, NULL), "n_events")
})
