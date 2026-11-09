# R/qc_app.R
# Shiny explorer for a diagnosis (sw_diagnose()) and an optional PeacoQC report.
#
# Portions of the page layout are adapted from iidx
# (https://github.com/saeyslab/iidx, 06_Dashboard.R), Copyright 2025 David Novak,
# licensed under the Apache License, Version 2.0. Adapted parts: the page CSS for
# the card look (.container-fluid background, .card border) and the spinner
# wrapper idea (.sw_app_spinner). Changed: rewritten for bslib cards and for the
# sw_diagnosis table; nothing of iidx's differential-expression code is used. The
# licence text is in inst/licenses/iidx-Apache-2.0.txt and the notice in NOTICE.

# CSS adapted from iidx: pale page background and thin card borders.
.sw_app_css <- "
  .container-fluid {background-color: #fafeff;}
  .card {box-shadow: 0px 0px 0px; border: 1px solid #d4d4d4;}
"

# Wrap an output in a loading spinner when shinycssloaders is installed (idea
# adapted from iidx's f_spinner()).
.sw_app_spinner <- function(ui) {
  if (requireNamespace("shinycssloaders", quietly = TRUE)) {
    shinycssloaders::withSpinner(
      ui,
      type = 8,
      color = "#3579bc",
      size = 0.75,
      hide.ui = FALSE
    )
  } else {
    ui
  }
}

# Rows of one check, largest value first (novelty log-density: lowest first).
.sw_qc_app_data <- function(diagnosis, module, check, n) {
  d <- diagnosis[
    diagnosis$module == module &
      diagnosis$check == check &
      !is.na(diagnosis$value),
  ]
  if (nrow(d) == 0) {
    return(d)
  }
  d$label <- ifelse(
    is.na(d$item),
    d$subject,
    paste(d$subject, d$item, sep = " / ")
  )
  d <- d[order(if (check == "novelty_log_density") d$value else -d$value), ]
  utils::head(d, n)
}

# Interactive dot plot of one check. Hover shows the subject; click selects it.
.sw_qc_app_plot <- function(d, check) {
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = .data$value,
      y = .data$label,
      tooltip = .data$subject,
      data_id = .data$subject
    )
  ) +
    ggiraph::geom_point_interactive(size = 2.5, colour = "#3579bc") +
    ggplot2::labs(x = check, y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  ggiraph::girafe(
    ggobj = p,
    options = list(
      ggiraph::opts_selection(type = "single"),
      ggiraph::opts_hover(css = "fill:#d95f02;stroke:#d95f02;")
    )
  )
}

#' Shiny explorer for a technical diagnosis
#'
#' @description
#' Open a Shiny app over an [sw_diagnose()] result. Pick a check to see the
#' files, batches or markers with the largest value, click a point to see every
#' value of that subject, and list the header differences. If you give the path
#' of a PeacoQC report, a heatmap of the QC results is shown as well.
#'
#' The app shows numbers only. It colours nothing as good or bad, because a
#' cutoff is your choice and the package sets none.
#'
#' @param diagnosis An `sw_diagnosis` result from [sw_diagnose()].
#' @param n Whole number, subjects shown per check at start. Required; you can
#'   change it in the app.
#' @param qc_report `NULL`, or the path of a `PeacoQC_report.txt` (or the QC
#'   output folder) for [sw_qc_heatmap()]. Required; pass `NULL` to skip.
#'
#' @return A `shiny.appobj`. Print it, or pass it to `shiny::runApp()`.
#'
#' @seealso [sw_diagnose()], [sw_diagnosis_top()], [sw_qc_heatmap()]
#'
#' @examples
#' \dontrun{
#' shiny::runApp(sw_qc_app(dx, n = 15, qc_report = NULL))
#' }
#'
#' @export
sw_qc_app <- function(diagnosis, n, qc_report) {
  .sw_stop_missing(c(
    diagnosis = missing(diagnosis),
    n = missing(n),
    qc_report = missing(qc_report)
  ))
  for (pkg in c("shiny", "bslib", "ggiraph", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for the app.", call. = FALSE)
    }
  }
  if (!inherits(diagnosis, "sw_diagnosis")) {
    stop("`diagnosis` must come from sw_diagnose().", call. = FALSE)
  }
  .sw_check_count(n, "n")
  if (
    !is.null(qc_report) && (!is.character(qc_report) || length(qc_report) != 1)
  ) {
    stop("`qc_report` must be NULL or one path.", call. = FALSE)
  }

  combos <- unique(diagnosis[!is.na(diagnosis$value), c("module", "check")])
  combos$id <- paste(combos$module, combos$check, sep = " : ")
  has_header <- "header" %in% diagnosis$module

  ui <- bslib::page_sidebar(
    title = "Technical diagnosis",
    theme = bslib::bs_theme(version = 5),
    sidebar = bslib::sidebar(
      shiny::selectInput("check", "Check", choices = combos$id),
      shiny::numericInput("n", "Subjects shown", value = n, min = 1, step = 1),
      shiny::helpText("Values only. No cutoff is applied.")
    ),
    shiny::tags$style(.sw_app_css),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        bslib::card_header("Largest values of the selected check"),
        .sw_app_spinner(ggiraph::girafeOutput("plot"))
      ),
      bslib::card(
        bslib::card_header("All values of the selected subject"),
        shiny::tableOutput("subject_table")
      ),
      if (has_header) {
        bslib::card(
          bslib::card_header("Header differences from the reference"),
          shiny::tableOutput("header_table")
        )
      },
      if (!is.null(qc_report)) {
        bslib::card(
          bslib::card_header("PeacoQC results per file"),
          .sw_app_spinner(shiny::plotOutput("qc_heatmap"))
        )
      }
    )
  )

  server <- function(input, output, session) {
    chosen <- shiny::reactive({
      row <- combos[combos$id == input$check, ]
      .sw_qc_app_data(diagnosis, row$module, row$check, input$n)
    })
    output$plot <- ggiraph::renderGirafe({
      d <- chosen()
      shiny::req(nrow(d) > 0)
      .sw_qc_app_plot(d, combos$check[combos$id == input$check])
    })
    output$subject_table <- shiny::renderTable({
      sel <- input$plot_selected
      shiny::req(length(sel) == 1)
      s <- diagnosis[
        diagnosis$subject == sel & !is.na(diagnosis$value),
        c("module", "check", "item", "value")
      ]
      s[order(s$module, s$check), ]
    })
    output$header_table <- shiny::renderTable({
      h <- diagnosis[
        diagnosis$module == "header" & diagnosis$value == 1,
        c("subject", "check", "item", "text")
      ]
      h
    })
    if (!is.null(qc_report)) {
      output$qc_heatmap <- shiny::renderPlot({
        sw_qc_heatmap(qc_report)
      })
    }
  }
  shiny::shinyApp(ui, server)
}
