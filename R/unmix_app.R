#' Live Shiny explorer for unmixing runs
#'
#' A Shiny app over one or more finished run folders. It shows the saved
#' checks, the spreading inflation factor (SIF) with what-if cutoffs, the
#' spectra of a fluorophore, and unmixed events redrawn for any channel pair.
#' It uses the same run-folder reader as [sw_unmix_dashboard()].
#'
#' **Verdicts come from the run.** The Checks tab shows the saved statuses. The
#' two SIF sliders recolour the bars with cutoffs of your choice. They are
#' display only: nothing is saved and a note says so. They start at the run's
#' saved cutoffs.
#'
#' **Events.** The Events tab reads the unmixed FCS files in
#' `<run>/unmixed/<method>/` with `flowCore`, draws `n_events` per file with
#' `seed` ([sw_unmix_sample_events()]), and redraws with
#' [sw_unmix_figure_distributions()] or [sw_unmix_figure_biplot()]. The event
#' count, seed, pseudo-log sigma and bins change what the panels show, so none
#' has a default.
#'
#' @param runs Run folders, a character vector; names become run labels.
#' @param n_events Positive whole number, events drawn per FCS file. Required.
#' @param seed Integer seed for the draws. Required.
#' @param pseudo_log_sigma Positive number, pseudo-log linear-zone width.
#'   Required.
#' @param bins Whole number of at least 10, bins per axis for the biplot.
#'   Required.
#' @param figures_dir Optional folder name inside each run to read the figure
#'   tables from (for example `"figures_v2"`).
#'
#' @return A `shiny.appobj`. Print it, or pass it to `shiny::runApp()`.
#'
#' @seealso [sw_unmix_dashboard()] for a single shareable HTML file.
#'
#' @examples
#' \dontrun{
#' shiny::runApp(sw_unmix_app(
#'   c(cells = "runs/cells-set-v2", beads = "runs/beads-set-v2"),
#'   n_events = 50000, seed = 1, pseudo_log_sigma = 1000, bins = 200
#' ))
#' }
#'
#' @export
sw_unmix_app <- function(
  runs,
  n_events,
  seed,
  pseudo_log_sigma,
  bins,
  figures_dir = NULL
) {
  .sw_fig_require(missing(n_events), "n_events")
  .sw_fig_require(missing(seed), "seed")
  .sw_fig_require(missing(pseudo_log_sigma), "pseudo_log_sigma")
  .sw_fig_require(missing(bins), "bins")
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for the app.", call. = FALSE)
  }
  .sw_fig_need_ggplot2()
  .sw_dash_check_settings(runs, pseudo_log_sigma, bins)
  labels <- .sw_dash_run_labels(runs)
  run_data <- Map(
    function(dir, nm) {
      .sw_dash_read_run(dir, pseudo_log_sigma, as.integer(bins), figures_dir, nm, bin_events = FALSE)
    },
    unname(runs),
    labels
  )
  names(run_data) <- labels
  shiny::shinyApp(
    ui = .sw_app_ui(labels),
    server = .sw_app_server(run_data, n_events, seed, pseudo_log_sigma, as.integer(bins))
  )
}

# Status of a value against what-if cutoffs: "warn" above `warn`, "fail" above
# `fail`, as in sw_unmix_diagnose limits. NA values are "not_run".
.sw_app_status <- function(value, warn, fail) {
  ifelse(
    is.na(value),
    "not_run",
    ifelse(value > fail, "fail", ifelse(value > warn, "warn", "pass"))
  )
}

.sw_app_status_colours <- c(
  pass = "#1b7f3b",
  warn = "#a15c00",
  fail = "#b42318",
  not_judged = "#5b6475",
  not_run = "#8a91a0"
)

# Events of one FCS file as a matrix, unmixed values as written (no transform).
.sw_app_read_exprs <- function(path, n_lines = NULL) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required to read FCS files.", call. = FALSE)
  }
  flowCore::exprs(flowCore::read.FCS(
    path,
    transformation = FALSE,
    truncate_max_range = FALSE,
    which.lines = n_lines
  ))
}

# Area channels present in every file, without scatter and time. OLS output
# has no AF channels, so files of different methods differ in the rest.
.sw_app_common_channels <- function(mats) {
  common <- Reduce(intersect, lapply(mats, colnames))
  common[grepl("-A$", common) & !grepl("^(FSC|SSC|Time)", common)]
}

.sw_app_ui <- function(labels) {
  shiny::fluidPage(
    shiny::titlePanel("Unmixing explorer"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("run", "Run", choices = labels),
        shiny::uiOutput("run_note"),
        width = 3
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Checks",
            shiny::checkboxGroupInput(
              "statuses",
              "Status",
              choices = c("pass", "warn", "fail", "not_judged", "not_run"),
              selected = c("pass", "warn", "fail", "not_judged", "not_run"),
              inline = TRUE
            ),
            shiny::textInput("search", "Search item or reason", ""),
            shiny::tableOutput("checks_table")
          ),
          shiny::tabPanel(
            "SIF",
            shiny::p(
              shiny::strong("What-if cutoffs are display only."),
              "Nothing is saved; saved verdicts are in Checks."
            ),
            shiny::sliderInput("sif_warn", "warn above", 0, 20, 4, step = 0.1),
            shiny::sliderInput("sif_fail", "fail above", 0, 20, 6, step = 0.1),
            shiny::plotOutput("sif_plot", height = "600px"),
            shiny::tableOutput("sif_table")
          ),
          shiny::tabPanel(
            "Spectra",
            shiny::selectInput("fluor", "Fluorophore", choices = NULL),
            shiny::plotOutput("spectra_plot")
          ),
          shiny::tabPanel(
            "Events",
            shiny::selectInput("panel", "Panel", c("Distributions", "Biplot")),
            shiny::selectInput("x_channel", "Channel (x)", choices = NULL),
            shiny::selectInput("y_channel", "Channel (y, biplot)", choices = NULL),
            shiny::actionButton("draw", "Read events and draw"),
            shiny::textOutput("events_note"),
            shiny::plotOutput("events_plot", height = "600px")
          )
        )
      )
    )
  )
}

.sw_app_server <- function(run_data, n_events, seed, sigma, bins) {
  function(input, output, session) {
    run <- shiny::reactive({
      shiny::req(input$run)
      run_data[[input$run]]
    })

    output$run_note <- shiny::renderUI({
      r <- run()
      lim <- r$limits
      shiny::tagList(
        shiny::p(shiny::strong("Run id: "), r$run_id),
        shiny::p(
          shiny::strong("Saved SIF cutoffs: "),
          if (all(is.na(lim))) "none" else paste0("warn above ", lim[["warn"]], ", fail above ", lim[["fail"]])
        )
      )
    })

    # Sliders start at the run's saved cutoffs (display only).
    shiny::observeEvent(run(), {
      lim <- run()$limits
      shiny::updateSliderInput(session, "sif_warn", value = if (is.na(lim[["warn"]])) 4 else lim[["warn"]])
      shiny::updateSliderInput(session, "sif_fail", value = if (is.na(lim[["fail"]])) 6 else lim[["fail"]])
      fl <- unique(c(run()$spectra_overlay$fluorophore, run()$sif$fluorophore))
      shiny::updateSelectInput(session, "fluor", choices = fl)
    })

    output$checks_table <- shiny::renderTable({
      ck <- run()$checks
      keep <- ck$status %in% input$statuses
      if (nzchar(input$search)) {
        hay <- tolower(paste(ck$item, ck$reason))
        keep <- keep & grepl(tolower(input$search), hay, fixed = TRUE)
      }
      utils::head(ck[keep, , drop = FALSE], 500L)
    })

    sif_data <- shiny::reactive({
      d <- run()$sif
      d <- d[!is.na(d$value), , drop = FALSE]
      d$whatif <- .sw_app_status(d$value, input$sif_warn, input$sif_fail)
      d[order(-d$value), , drop = FALSE]
    })

    output$sif_plot <- shiny::renderPlot({
      d <- sif_data()
      d$fluorophore <- factor(d$fluorophore, levels = rev(d$fluorophore))
      ggplot2::ggplot(d, ggplot2::aes(x = value, y = fluorophore, fill = whatif)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept = c(input$sif_warn, input$sif_fail), linetype = "dashed") +
        ggplot2::scale_fill_manual(values = .sw_app_status_colours, drop = FALSE) +
        ggplot2::labs(
          x = "Spreading inflation factor (what-if cutoffs, display only)",
          y = NULL,
          fill = "What-if"
        ) +
        ggplot2::theme_minimal()
    })

    output$sif_table <- shiny::renderTable({
      d <- sif_data()
      data.frame(
        fluorophore = d$fluorophore,
        SIF = d$value,
        saved = d$status,
        what_if = d$whatif,
        partner = d$partner
      )
    })

    output$spectra_plot <- shiny::renderPlot({
      shiny::req(input$fluor)
      ov <- run()$spectra_overlay
      shiny::req(ov)
      d <- ov[ov$fluorophore == input$fluor, , drop = FALSE]
      d$detector <- factor(d$detector, levels = unique(ov$detector))
      cols <- .sw_fig_palette(unique(d$method))
      ggplot2::ggplot(d, ggplot2::aes(detector, value_norm, colour = method, group = method)) +
        ggplot2::geom_line() +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::labs(x = "Detector", y = "Peak-normalised intensity", title = input$fluor) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
    })

    # Events: read FCS files only when asked.
    fcs_files <- shiny::reactive({
      f <- list.files(file.path(run()$run_dir, "unmixed"), pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE)
      names(f) <- basename(dirname(f))
      f
    })
    events <- shiny::eventReactive(input$draw, {
      f <- fcs_files()
      shiny::req(length(f) > 0L)
      mats <- lapply(f, .sw_app_read_exprs)
      common <- .sw_app_common_channels(mats)
      sw_unmix_sample_events(
        lapply(mats, \(m) m[, common, drop = FALSE]),
        n_events = n_events,
        seed = seed
      )
    })
    shiny::observeEvent(run(), {
      f <- fcs_files()
      if (length(f) > 0L && requireNamespace("flowCore", quietly = TRUE)) {
        heads <- lapply(f, .sw_app_read_exprs, n_lines = 1L)
        ch <- .sw_app_common_channels(heads)
        shiny::updateSelectInput(session, "x_channel", choices = ch)
        shiny::updateSelectInput(session, "y_channel", choices = ch)
      }
    })
    output$events_note <- shiny::renderText({
      paste0(
        "Draws ", n_events, " events per file, seed ", seed,
        "; pseudo-log sigma ", sigma, "; ", bins, " bins (set when the app was created)."
      )
    })
    output$events_plot <- shiny::renderPlot({
      ev <- events()
      if (input$panel == "Distributions") {
        sw_unmix_figure_distributions(ev, input$x_channel, pseudo_log_sigma = sigma)$plot
      } else {
        sw_unmix_figure_biplot(ev, input$x_channel, input$y_channel, pseudo_log_sigma = sigma, bins = bins)$plot
      }
    })
  }
}
