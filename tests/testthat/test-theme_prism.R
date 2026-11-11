# Tests for R/theme_prism.R (standalone prism-style theme and the Nature theme).

test_that("the prism theme is complete and needs base_size", {
  skip_if_not_installed("ggplot2")
  th <- .sw_theme_prism(base_size = 10)
  expect_true(attr(th, "complete"))
  expect_error(.sw_theme_prism(), "`base_size` is required")
})

test_that("the prism theme sets the black-and-white elements", {
  skip_if_not_installed("ggplot2")
  th <- .sw_theme_prism(base_size = 14)
  expect_identical(th$line$colour, "#000000")
  expect_identical(th$text$face, "bold")
  expect_identical(th$text$size, 14)
  expect_s3_class(th$panel.grid, "element_blank")
  expect_s3_class(th$strip.background, "element_blank")
  expect_s3_class(th$legend.title, "element_blank")
  expect_identical(th$legend.position, "right")
  expect_equal(as.numeric(th$axis.ticks.length), 14 / 2.5)
})

test_that("axis angle and border are checked", {
  skip_if_not_installed("ggplot2")
  expect_error(
    .sw_theme_prism(base_size = 10, axis_text_angle = 30),
    "must be one of"
  )
  expect_error(.sw_theme_prism(base_size = 10, border = NA), "TRUE or FALSE")
  with_border <- .sw_theme_prism(base_size = 10, border = TRUE)
  expect_s3_class(with_border$panel.border, "element_rect")
  expect_s3_class(with_border$axis.line, "element_blank")
  no_border <- .sw_theme_prism(base_size = 10, border = FALSE)
  expect_s3_class(no_border$panel.border, "element_blank")
  expect_s3_class(no_border$axis.line, "element_line")
})

test_that("the prism theme matches ggprism on the elements it sets", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggprism")
  ours <- .sw_theme_prism(
    base_size = 12,
    base_family = "sans",
    axis_text_angle = 45
  )
  theirs <- ggprism::theme_prism(
    base_size = 12,
    base_family = "sans",
    axis_text_angle = 45
  )
  set_here <- c(
    "line",
    "rect",
    "text",
    "axis.line",
    "axis.text",
    "axis.text.x",
    "axis.text.y",
    "axis.ticks",
    "axis.ticks.length",
    "axis.title",
    "axis.title.x",
    "axis.title.y",
    "legend.key.width",
    "legend.text",
    "legend.title",
    "legend.position",
    "panel.border",
    "panel.grid",
    "panel.spacing",
    "strip.background",
    "strip.text",
    "strip.text.x",
    "strip.text.y",
    "plot.title",
    "plot.caption",
    "plot.margin"
  )
  for (el in set_here) {
    expect_equal(ours[[el]], theirs[[el]], info = el, ignore_attr = TRUE)
  }
})

test_that("the Nature theme draws text between 5 and 7 pt with the legend right", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = 1:4, y = c(1, 3, 2, 4), g = c("a", "a", "b", "b"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(g)) +
    ggplot2::labs(title = "Title", caption = "Caption") +
    .sw_theme_nature()
  pt <- .sw_fig_drawn_font_sizes(p)
  expect_gte(min(pt), 5)
  expect_lte(max(pt), 7)
  th <- .sw_theme_nature()
  expect_identical(th$legend.position, "right")
  expect_identical(th$strip.clip, "off")
  expect_identical(th$legend.title$face, "bold")
})
