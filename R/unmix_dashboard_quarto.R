# sw_unmix_dashboard_quarto(): the benchmark 05 Quarto dashboard pattern
# (96_dashboard.qmd) as a package template. Verdict cards, then one section
# per finding with Figures / Data / notes tabs. The page is a local HTML file;
# it is never published.

#' Write a Quarto HTML dashboard of figures, their data and notes
#'
#' Renders one self-contained HTML page with the template in
#' `inst/quarto_dashboard`: a header band, verdict cards, an optional intro,
#' then one section per entry of `sections`. Each section shows its finding
#' and three tabs: the figures, the plotted data behind each figure (sortable,
#' searchable tables) and notes (for example quotes from the AutoSpectral
#' documentation). The page uses its own stylesheet with light and dark mode
#' (`theme: none`) and embeds every image and script (`embed-resources`), so
#' it works offline.
#'
#' **Images.** Each figure is converted from its vector PDF to a PNG at `dpi`
#' with `pdftoppm` (poppler). The PDF is the figure path itself, or the PDF
#' with the same stem beside a PNG path. A PNG without a PDF, or any PNG when
#' `pdftoppm` is missing, is copied as it is; `img_log.csv` records which
#' source each image used.
#'
#' **Data tab.** The table of a figure is the CSV with the same stem beside
#' it, or `source_data/<stem>.csv` in its folder (the layout of
#' [sw_unmix_save_figures()]). Only the first `max_rows` rows are shown, with
#' numbers rounded to 3 significant digits; the note under the table names
#' the full file.
#'
#' The function writes a local file only. It does not publish anything.
#'
#' @param sections List of sections. Each is a list with `heading` (one
#'   string), `finding` (Markdown text, may be `NULL`), `figures` (optional
#'   data frame with `path` to a PNG or PDF and `caption`) and `notes`
#'   (optional character vector).
#' @param path Folder to create. It must not exist. Gets `dashboard.html`,
#'   `img/` with the web images, `img_log.csv` and `sessionInfo.txt`.
#' @param title Page title.
#' @param verdicts Optional data frame of verdict cards with columns `label`,
#'   `value`, `note` and `status` (`"good"`, `"warn"` or `"bad"`, which sets
#'   the card colour). The status is the caller's reading, not a computed
#'   verdict.
#' @param intro Optional Markdown shown under the cards.
#' @param subtitle Optional subtitle.
#' @param kicker Small text above the title in the header band.
#' @param footer Optional footer text.
#' @param notes_label Name of the notes tab.
#' @param dpi Resolution of the PNGs made from PDFs.
#' @param max_rows Rows shown per data table.
#'
#' @return The path of `dashboard.html`, invisibly.
#'
#' @seealso [sw_unmix_dashboard()] for the interactive single-file dashboard
#'   of run folders.
#'
#' @examples
#' \dontrun{
#' sections <- list(list(
#'   heading = "Spike-in recovery",
#'   finding = "Setup b loses more of the TNF spike than setup a.",
#'   figures = data.frame(
#'     path = "figures/panels/spike-summary_TNF_CD8T.png",
#'     caption = "Share of spiked cells that lose most of the spike."
#'   ),
#'   notes = "AutoSpectral: unmixing is per cell."
#' ))
#' sw_unmix_dashboard_quarto(sections, path = "dashboard_out", title = "Unmix dashboard")
#' }
#'
#' @export
sw_unmix_dashboard_quarto <- function(
  sections,
  path,
  title,
  verdicts = NULL,
  intro = NULL,
  subtitle = NULL,
  kicker = "SpectraWeaveR",
  footer = NULL,
  notes_label = "AutoSpectral notes",
  dpi = 300,
  max_rows = 500
) {
  .sw_stop_missing(c(
    sections = missing(sections),
    path = missing(path),
    title = missing(title)
  ))
  for (pkg in c("quarto", "reactable", "knitr", "rmarkdown")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' is required for the Quarto dashboard.", pkg),
        call. = FALSE
      )
    }
  }
  one_string <- \(x) is.character(x) && length(x) == 1L && !is.na(x)
  if (!one_string(path)) {
    stop("`path` must be one folder path.", call. = FALSE)
  }
  if (file.exists(path)) {
    stop("`path` exists; refusing to overwrite: ", path, call. = FALSE)
  }
  if (!one_string(title)) {
    stop("`title` must be one string.", call. = FALSE)
  }
  for (arg in c("subtitle", "intro", "footer")) {
    v <- get(arg)
    if (!is.null(v) && !one_string(v)) {
      stop(sprintf("`%s` must be NULL or one string.", arg), call. = FALSE)
    }
  }
  if (!one_string(kicker) || !one_string(notes_label)) {
    stop("`kicker` and `notes_label` must be one string each.", call. = FALSE)
  }
  for (arg in c("dpi", "max_rows")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1) {
      stop(sprintf("`%s` must be one positive number.", arg), call. = FALSE)
    }
  }
  verdicts <- .sw_qdash_check_verdicts(verdicts)
  sections <- .sw_qdash_check_sections(sections)

  # Web images and tables go to `path`; the render runs in a scratch folder.
  dir.create(file.path(path, "img"), recursive = TRUE)
  path <- normalizePath(path)
  work <- tempfile("sw_qdash_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  has_pdftoppm <- nzchar(Sys.which("pdftoppm"))
  img_log <- list()
  for (k in seq_along(sections)) {
    figs <- sections[[k]]$figures
    if (is.null(figs) || nrow(figs) == 0L) {
      next
    }
    tables <- vector("list", nrow(figs))
    figs$label <- basename(figs$path)
    figs$web <- NA_character_
    figs$has_data <- FALSE
    figs$data_note <- ""
    for (i in seq_len(nrow(figs))) {
      stem <- sub("\\.(png|pdf)$", "", figs$path[i], ignore.case = TRUE)
      web <- file.path(path, "img", sprintf("s%02d_f%02d.png", k, i))
      pdf <- paste0(stem, ".pdf")
      png <- paste0(stem, ".png")
      source <- if (has_pdftoppm && file.exists(pdf)) "pdf" else "png"
      if (source == "pdf") {
        status <- system2(
          "pdftoppm",
          c(
            "-r",
            dpi,
            "-png",
            "-singlefile",
            shQuote(pdf),
            shQuote(sub("\\.png$", "", web))
          )
        )
        if (!identical(as.integer(status), 0L) || !file.exists(web)) {
          stop("pdftoppm failed on ", pdf, call. = FALSE)
        }
      } else {
        if (!file.exists(png)) {
          stop(
            "No PNG for ",
            figs$path[i],
            if (!has_pdftoppm) " (pdftoppm is not installed)",
            call. = FALSE
          )
        }
        file.copy(png, web)
      }
      # Quarto embeds images by paths relative to the page.
      figs$web[i] <- file.path("img", basename(web))
      img_log[[length(img_log) + 1L]] <- data.frame(
        section = sections[[k]]$heading,
        figure = figs$path[i],
        source = source,
        web = web
      )

      csv <- c(
        paste0(stem, ".csv"),
        file.path(dirname(stem), "source_data", paste0(basename(stem), ".csv"))
      )
      csv <- csv[file.exists(csv)]
      if (length(csv) == 0L) {
        figs$data_note[i] <- "No source table was saved beside this figure."
        next
      }
      tab <- utils::read.csv(csv[[1L]], check.names = FALSE)
      n <- nrow(tab)
      tab <- utils::head(tab, max_rows)
      is_num <- vapply(tab, is.double, logical(1))
      tab[is_num] <- lapply(tab[is_num], signif, digits = 3)
      tables[[i]] <- tab
      figs$has_data[i] <- TRUE
      figs$data_note[i] <- if (n > max_rows) {
        sprintf(
          "First %d of %d rows; full table: `%s`.",
          max_rows,
          n,
          csv[[1L]]
        )
      } else {
        sprintf("%d rows; file: `%s`.", n, csv[[1L]])
      }
    }
    figs$caption <- gsub("([][])", "\\\\\\1", figs$caption)
    sections[[k]]$figures <- figs
    sections[[k]]$tables <- tables
  }
  utils::write.csv(
    do.call(rbind, img_log),
    file.path(path, "img_log.csv"),
    row.names = FALSE
  )

  # Template files, with the title and partials filled in.
  tpl_dir <- system.file("quarto_dashboard", package = "SpectraWeaveR")
  if (!nzchar(tpl_dir)) {
    tpl_dir <- file.path("inst", "quarto_dashboard")
  }
  file.copy(file.path(tpl_dir, "dashboard.css"), work)
  file.copy(file.path(path, "img"), work, recursive = TRUE)
  yaml_str <- \(x) gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", x))
  qmd <- readLines(file.path(tpl_dir, "dashboard.qmd"), warn = FALSE)
  qmd <- sub("__TITLE__", yaml_str(title), qmd, fixed = TRUE)
  if (is.null(subtitle)) {
    qmd <- qmd[!grepl("__SUBTITLE__", qmd, fixed = TRUE)]
  } else {
    qmd <- sub("__SUBTITLE__", yaml_str(subtitle), qmd, fixed = TRUE)
  }
  writeLines(qmd, file.path(work, "dashboard.qmd"))
  writeLines(
    c(
      "<div class=\"sw-header\"><div class=\"sw-header-inner\">",
      sprintf("<div class=\"sw-kicker\">%s</div>", .sw_html_escape(kicker)),
      sprintf("<div class=\"sw-title\">%s</div>", .sw_html_escape(title)),
      "</div></div>"
    ),
    file.path(work, "header.html")
  )
  writeLines(
    sprintf(
      "<footer class=\"sw-footer\">%s</footer>",
      .sw_html_escape(
        if (is.null(footer)) {
          paste("SpectraWeaveR", utils::packageVersion("SpectraWeaveR"))
        } else {
          footer
        }
      )
    ),
    file.path(work, "footer.html")
  )
  saveRDS(
    list(
      verdicts = verdicts,
      intro = intro,
      sections = sections,
      notes_label = notes_label
    ),
    file.path(work, "spec.rds")
  )

  quarto::quarto_render(
    file.path(work, "dashboard.qmd"),
    execute_params = list(spec = file.path(work, "spec.rds")),
    quiet = TRUE
  )
  out <- file.path(path, "dashboard.html")
  if (!file.copy(file.path(work, "dashboard.html"), out)) {
    stop("Quarto did not write dashboard.html.", call. = FALSE)
  }
  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(path, "sessionInfo.txt")
  )
  invisible(out)
}

# Verdict cards: label, value, note, status; text is HTML-escaped here because
# the template writes it as raw HTML.
.sw_qdash_check_verdicts <- function(verdicts) {
  if (is.null(verdicts)) {
    return(NULL)
  }
  need <- c("label", "value", "note", "status")
  if (!is.data.frame(verdicts) || !all(need %in% names(verdicts))) {
    stop(
      "`verdicts` must be a data frame with label, value, note and status.",
      call. = FALSE
    )
  }
  bad <- setdiff(verdicts$status, c("good", "warn", "bad"))
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Unknown verdict status: %s. Use good, warn or bad.",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  verdicts <- verdicts[need]
  for (col in c("label", "value", "note")) {
    verdicts[[col]] <- .sw_html_escape(as.character(verdicts[[col]]))
  }
  verdicts
}

.sw_qdash_check_sections <- function(sections) {
  if (!is.list(sections) || length(sections) == 0L) {
    stop("`sections` must be a non-empty list.", call. = FALSE)
  }
  lapply(sections, \(sec) {
    if (
      !is.list(sec) ||
        !is.character(sec$heading) ||
        length(sec$heading) != 1L
    ) {
      stop("Each section needs one string `heading`.", call. = FALSE)
    }
    if (is.null(sec$finding)) {
      sec$finding <- ""
    }
    figs <- sec$figures
    if (!is.null(figs)) {
      if (!is.data.frame(figs) || !all(c("path", "caption") %in% names(figs))) {
        stop(
          sprintf(
            "Section '%s': `figures` must be a data frame with path and caption.",
            sec$heading
          ),
          call. = FALSE
        )
      }
      figs <- data.frame(
        path = as.character(figs$path),
        caption = as.character(figs$caption)
      )
      gone <- figs$path[
        !file.exists(figs$path) &
          !file.exists(sub("\\.png$", ".pdf", figs$path))
      ]
      if (length(gone) > 0L) {
        stop(
          sprintf(
            "Section '%s': figure(s) not found: %s.",
            sec$heading,
            paste(gone, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      sec$figures <- figs
    }
    sec$notes <- as.character(sec$notes)
    sec
  })
}
