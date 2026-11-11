# Tests for R/unmix_dashboard_quarto.R. Mock figures only.

# One section with a PNG figure and its source table beside it.
make_sections <- function(dir) {
  dir.create(file.path(dir, "panels", "source_data"), recursive = TRUE)
  png <- file.path(dir, "panels", "spike-summary_TNF_T.png")
  grDevices::png(png, width = 200, height = 150)
  graphics::plot(1:3)
  grDevices::dev.off()
  utils::write.csv(
    data.frame(setup = c("a", "b"), share_lost = c(0.0123456, 0.2)),
    file.path(dir, "panels", "source_data", "spike-summary_TNF_T.csv"),
    row.names = FALSE
  )
  list(list(
    heading = "Spike-in recovery",
    finding = "Setup b loses more of the spike.",
    figures = data.frame(path = png, caption = "Share lost [mock]."),
    notes = "AutoSpectral: a quoted note."
  ))
}

test_that("arguments are checked before anything is written", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("reactable")
  root <- withr::local_tempdir()
  secs <- make_sections(root)
  expect_error(sw_unmix_dashboard_quarto(secs, title = "x"), "path")
  expect_error(
    sw_unmix_dashboard_quarto(secs, path = root, title = "x"),
    "refusing to overwrite"
  )
  out <- file.path(root, "out")
  expect_error(
    sw_unmix_dashboard_quarto(
      secs,
      out,
      "x",
      verdicts = data.frame(label = "a", value = "1", note = "", status = "ok")
    ),
    "Unknown verdict status"
  )
  bad <- secs
  bad[[1]]$figures$path <- file.path(root, "missing.png")
  expect_error(sw_unmix_dashboard_quarto(bad, out, "x"), "not found")
  expect_false(dir.exists(out))
})

test_that("the dashboard renders one self-contained file with tabs and data", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if_not_installed("reactable")
  skip_if(is.null(quarto::quarto_path()), "Quarto CLI not found")
  root <- withr::local_tempdir()
  out_dir <- file.path(root, "dash")
  html_path <- sw_unmix_dashboard_quarto(
    make_sections(root),
    out_dir,
    title = "Mock <dashboard>",
    verdicts = data.frame(
      label = "TNF loss",
      value = "2% vs 14%",
      note = "mock",
      status = "good"
    ),
    footer = "mock footer"
  )
  expect_true(file.exists(html_path))
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  expect_match(html, "verdict good", fixed = TRUE)
  expect_match(html, "Mock &lt;dashboard&gt;", fixed = TRUE)
  expect_match(html, "Spike-in recovery", fixed = TRUE)
  expect_match(html, "AutoSpectral notes", fixed = TRUE)
  expect_match(html, "data:image/png;base64", fixed = TRUE)
  expect_match(html, "2 rows; file:", fixed = TRUE)
  expect_match(html, "mock footer", fixed = TRUE)
  log <- utils::read.csv(file.path(out_dir, "img_log.csv"))
  expect_identical(log$source, "png")
  expect_true(file.exists(file.path(out_dir, "sessionInfo.txt")))
})
