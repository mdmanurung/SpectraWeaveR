# tests/testthat/test-qc_app.R
# Tests for R/qc_app.R: pure helpers, app construction and server logic.

make_diag <- function() {
  d <- tibble::tibble(
    module = c(
      rep("signal", 6),
      rep("batch", 2),
      "header",
      "header",
      "anomaly",
      "anomaly"
    ),
    subject = c(
      "f1",
      "f2",
      "f3",
      "f1",
      "f2",
      "f3",
      "b1 ~ b2",
      "b1 ~ b3",
      "f1",
      "f2",
      "f1",
      "f2"
    ),
    check = c(
      rep("saturation_share", 3),
      rep("rate_cv", 3),
      rep("pair_emd", 2),
      "keyword $CYT",
      "keyword $CYT",
      "novelty_log_density",
      "novelty_log_density"
    ),
    item = c("A", "A", "A", NA, NA, NA, "A", "A", NA, NA, NA, NA),
    value = c(0.1, 0.4, 0.0, 0.5, 0.2, 0.9, 1.1, 0.2, 0, 1, -3, -30),
    text = c(rep(NA, 8), "match: X vs X", "differs: Y vs X", NA, NA)
  )
  class(d) <- c("sw_diagnosis", class(d))
  d
}

test_that("sw_qc_app needs every argument and a diagnosis", {
  expect_error(sw_qc_app(), "no default")
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("ggiraph")
  expect_error(
    sw_qc_app(data.frame(), n = 5, qc_report = NULL),
    "must come from"
  )
  expect_error(sw_qc_app(make_diag(), n = 5), "qc_report")
})

test_that("data helper sorts by value and ranks novelty from the lowest", {
  d <- make_diag()
  top <- .sw_qc_app_data(d, "signal", "saturation_share", 2)
  expect_equal(top$subject, c("f2", "f1"))
  expect_equal(top$label, c("f2 / A", "f1 / A"))
  nov <- .sw_qc_app_data(d, "anomaly", "novelty_log_density", 5)
  expect_equal(nov$subject, c("f2", "f1"))
  expect_equal(nrow(.sw_qc_app_data(d, "signal", "no_such_check", 3)), 0)
})

test_that("the plot helper returns an interactive widget", {
  skip_if_not_installed("ggiraph")
  skip_if_not_installed("ggplot2")
  top <- .sw_qc_app_data(make_diag(), "signal", "rate_cv", 3)
  w <- .sw_qc_app_plot(top, "rate_cv")
  expect_s3_class(w, "girafe")
})

test_that("the app builds and its server reacts to the check choice", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("ggiraph")
  app <- sw_qc_app(make_diag(), n = 2, qc_report = NULL)
  expect_s3_class(app, "shiny.appobj")
  shiny::testServer(app, {
    session$setInputs(check = "signal : rate_cv", n = 2)
    expect_equal(nrow(chosen()), 2)
    expect_equal(chosen()$subject, c("f3", "f1"))
    session$setInputs(check = "batch : pair_emd", n = 5)
    expect_equal(chosen()$subject, c("b1 ~ b2", "b1 ~ b3"))
  })
})
