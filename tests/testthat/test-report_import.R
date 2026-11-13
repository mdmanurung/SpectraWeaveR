# tests/testthat/test-report_import.R
# Unit tests for R/report_import.R. Builds every import-step result from tiny
# synthetic FCS files, then writes reports into temp folders.

# Write one FCS file with chosen channels, markers and keywords.
make_rep_fcs <- function(dir, name, marker, kw = list()) {
  skip_if_not_installed("flowCore")
  chan <- paste0("B", seq_along(marker), "-A")
  m <- matrix(
    stats::rnorm(20 * length(chan)),
    ncol = length(chan),
    dimnames = list(NULL, chan)
  )
  ff <- flowCore::flowFrame(m)
  flowCore::pData(flowCore::parameters(ff))$desc <- marker
  desc <- flowCore::keyword(ff)
  for (k in names(kw)) {
    desc[[k]] <- kw[[k]]
  }
  flowCore::keyword(ff) <- desc
  path <- file.path(dir, name)
  flowCore::write.FCS(ff, path)
  path
}

# Four files in two batches (one with another dye), a junk file and a missing one.
# Returns every function output the report accepts, in one named list.
make_rep_result <- function() {
  skip_if_not_installed("withr")
  skip_if_not_installed("flowCore")
  d <- withr::local_tempdir(.local_envir = parent.frame())
  kw <- function(cyt, date) {
    list("$CYT" = cyt, "$DATE" = date, "$BTIM" = "09:00:00", "$CYTSN" = "SN1")
  }
  p <- c(
    a = make_rep_fcs(
      d,
      "a.fcs",
      c("CD3 BV421", "CD4 BB515"),
      kw("A", "2026-03-01")
    ),
    b = make_rep_fcs(
      d,
      "b.fcs",
      c("CD3 BV421", "CD4 BB515"),
      kw("A", "2026-03-01")
    ),
    c = make_rep_fcs(
      d,
      "c.fcs",
      c("CD3 BV421", "CD4 PE"),
      kw("B", "2026-03-05")
    ),
    d = make_rep_fcs(
      d,
      "d.fcs",
      c("CD3 BV421", "CD8-APC"),
      kw("B", "2026-03-05")
    )
  )
  writeLines("not an fcs file", junk <- file.path(d, "junk.fcs"))
  all_paths <- c(unname(p), junk, file.path(d, "gone.fcs"))
  man <- sw_sample_manifest(
    data.frame(
      id = c("a", "b", "c", "d", "ghost"),
      b = c("B1", "B1", "B2", "B2", "B2")
    ),
    sample_id = "id",
    batch = "b"
  )
  set.seed(5)
  sig <- list(
    f1 = cbind(A = stats::rnorm(50, 1), B = stats::rnorm(50)),
    f2 = cbind(A = stats::rnorm(50, 2), B = stats::rnorm(50)),
    f3 = cbind(A = stats::rnorm(5, 3), B = stats::rnorm(5))
  )
  list(
    bound = sw_manifest_bind_files(
      manifest = man,
      paths = c(unname(p), junk),
      key = "file_stem",
      case = "exact",
      recursive = FALSE,
      on_duplicate_file = "unbind",
      on_extra_file = "list",
      nextdata = "bind"
    ),
    header_audit = sw_fcs_header_audit(
      unname(p),
      keywords = "$CYT",
      reference = p[["a"]]
    ),
    integrity = sw_fcs_integrity_audit(
      all_paths,
      tot_tolerance = 0,
      time_keyword = "Time",
      header_only = FALSE
    ),
    signatures = sw_fcs_header_signatures(
      all_paths,
      keywords = "$CYT",
      panel = list(T = c("CD3 BV421", "CD4 BB515"))
    ),
    compensation = sw_fcs_compensation_state(
      all_paths,
      matrix_keywords = "$SPILLOVER",
      apply_keywords = "APPLY COMPENSATION",
      apply_true = "TRUE",
      apply_false = "FALSE",
      zero_tol = 1e-6,
      expected_state = NULL
    ),
    batch_date = sw_fcs_batch_date_check(
      stats::setNames(unname(p), c("a", "b", "c", "d")),
      manifest = man,
      date_keywords = "$DATE",
      time_keywords = "$BTIM",
      serial_keyword = "$CYTSN",
      date_format = "%Y-%m-%d",
      time_zone = "UTC",
      run_gap_hours = 4
    ),
    dye_map = sw_fcs_marker_dye_map(
      unname(p),
      batch = c(a.fcs = "B1", b.fcs = "B1", c.fcs = "B2", d.fcs = "B2"),
      separator = " "
    ),
    channel_summary = sw_fcs_signal_audit(
      sig,
      channels = c("A", "B"),
      range_max = NULL,
      range_min = NULL,
      time_channel = NULL,
      time_bin = NULL,
      area_height = NULL,
      area_height_quantile = NULL,
      reference_group = NULL,
      reference = NULL,
      summary = list(
        probs = c(0.25, 0.75),
        min_events = 10,
        stat = "median"
      )
    )
  )
}

new_out <- function() {
  file.path(withr::local_tempdir(.local_envir = parent.frame()), "run1")
}

test_that("required arguments have no default", {
  expect_error(sw_report_import(), "no default.*'result'.*'out_dir'")
  expect_error(sw_report_import(result = list()), "no default.*'out_dir'")
  expect_error(sw_report_import(out_dir = "x"), "no default.*'result'")
})

test_that("an existing out_dir stops and nothing is written", {
  res <- make_rep_result()
  out <- new_out()
  dir.create(out)
  writeLines("keep", file.path(out, "mine.txt"))
  expect_error(sw_report_import(res, out), "exists")
  expect_identical(list.files(out), "mine.txt")
})

test_that("result is checked before anything is written", {
  res <- make_rep_result()
  out <- new_out()
  expect_error(sw_report_import(1, out), "named list")
  expect_error(sw_report_import(unname(res), out), "named list")
  expect_error(sw_report_import(list(), out), "no recognised")
  expect_error(sw_report_import(list(nope = 1), out), "nope")
  expect_error(
    sw_report_import(c(res["integrity"], res["integrity"]), out),
    "unique"
  )
  # A result under the wrong name is refused, not skipped.
  expect_error(
    sw_report_import(list(integrity = res$dye_map), out),
    "integrity.*sw_fcs_integrity_audit"
  )
  expect_error(sw_report_import(res, c(out, out)), "out_dir")
  expect_false(dir.exists(out))
})

test_that("a full report writes tables, one PNG and source CSV per figure, provenance and session", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res, out)
  expect_named(rep, c("tables", "figures", "source_data", "dir", "absent"))
  expect_identical(rep$dir, out)
  expect_length(rep$absent, 0)
  expect_true(all(file.exists(unlist(rep$tables))))
  expect_true(file.exists(file.path(out, "tables", "bound_files.csv")))
  expect_true(file.exists(file.path(out, "tables", "integrity.csv")))
  expect_true(file.exists(file.path(out, "sessionInfo.txt")))

  fig_png <- Filter(\(f) file.exists(f), unlist(rep$figures))
  fallback <- Filter(
    \(f) startsWith(f, "table fallback: "),
    unlist(rep$figures)
  )
  expect_gt(length(fig_png), 0)
  expect_true("bound_files" %in% names(fallback))
  # Every figure has a PNG and a source CSV of the same panel name.
  expect_setequal(names(fig_png), names(rep$source_data))
  for (panel in names(fig_png)) {
    expect_true(file.exists(rep$source_data[[panel]]))
    expect_identical(
      rep$source_data[[panel]],
      file.path(out, "source_data", paste0(panel, ".csv"))
    )
    expect_gt(file.size(fig_png[[panel]]), 0)
    expect_gt(nrow(utils::read.csv(rep$source_data[[panel]])), 0)
  }
  # Every panel is either a figure or an explained fallback.
  expect_true(all(
    vapply(
      rep$figures,
      \(f) file.exists(f) || startsWith(f, "table fallback: "),
      logical(1)
    )
  ))

  prov <- jsonlite::read_json(file.path(out, "provenance.json"))
  expect_true(all(
    c(
      "package_version",
      "r_version",
      "platform",
      "datetime",
      "inputs",
      "absent"
    ) %in%
      names(prov)
  ))
})

test_that("source data holds one row per file, including files that could not be assessed", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res, out)
  sd <- utils::read.csv(rep$source_data[["integrity_tot"]])
  # five paths: a-d, junk and gone make six files in the audit.
  expect_equal(nrow(sd), length(unique(res$integrity$file)))
  expect_true(any(is.na(sd$value)))
  expect_true(all(c("junk.fcs", "gone.fcs") %in% sd$file))
  tab <- utils::read.csv(rep$tables[["integrity"]])
  expect_equal(nrow(tab), nrow(res$integrity))
  # Unmatched manifest rows and files stay on disk next to the bound table.
  nm <- utils::read.csv(rep$tables[["bound_manifest_no_file"]])
  expect_true("ghost" %in% nm$sample_id)
  expect_equal(
    nrow(utils::read.csv(rep$tables[["bound_not_fcs"]])),
    nrow(res$bound$not_fcs)
  )
})

test_that("absent inputs are listed in the result and in provenance.json", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res["dye_map"], out)
  expect_setequal(
    rep$absent,
    setdiff(
      c(
        "bound",
        "header_audit",
        "integrity",
        "signatures",
        "compensation",
        "batch_date",
        "dye_map",
        "channel_summary"
      ),
      "dye_map"
    )
  )
  expect_false(file.exists(file.path(out, "tables", "integrity.csv")))
  expect_true(file.exists(file.path(out, "tables", "dye_map.csv")))
  prov <- jsonlite::read_json(file.path(out, "provenance.json"))
  expect_setequal(unlist(prov$absent), rep$absent)
  expect_identical(unlist(prov$inputs), "dye_map")
})

test_that("panels with no data to draw say why instead of vanishing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  res$signatures <- sw_fcs_header_signatures(
    res$integrity$path[!duplicated(res$integrity$path)][1:4],
    keywords = "$CYT",
    panel = NULL
  )
  out <- new_out()
  rep <- sw_report_import(res["signatures"], out)
  expect_match(rep$figures[["signatures_markers"]], "^table fallback: .*panel")
  expect_true(file.exists(rep$figures[["signatures_files"]]))
  expect_false("signatures_markers" %in% names(rep$source_data))
})

test_that("channel summary figure keeps every file and channel in its source data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res["channel_summary"], out)
  sd <- utils::read.csv(rep$source_data[["channel_summary_dots"]])
  expect_setequal(sd$file, c("f1", "f2", "f3"))
  expect_setequal(sd$channel, c("A", "B"))
  expect_equal(nrow(sd), 6)
  # f3 is thin (5 events): kept, with an empty value.
  expect_true(all(is.na(sd$location[sd$file == "f3"])))
})

test_that("out_dir that is empty or cannot be created stops with a clear message", {
  res <- make_rep_result()
  expect_error(sw_report_import(res["dye_map"], ""), "out_dir")
  blocker <- withr::local_tempfile()
  writeLines("x", blocker)
  expect_no_warning(
    expect_error(
      sw_report_import(res["dye_map"], file.path(blocker, "run")),
      "Cannot create `out_dir`"
    )
  )
})

test_that("a render error leaves no half-built folder behind", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  testthat::local_mocked_bindings(
    .sw_ri_save = function(panel, name, out_dir) stop("render failed")
  )
  expect_error(sw_report_import(res["dye_map"], out), "render failed")
  expect_false(dir.exists(out))
})

test_that("provenance.json carries input provenance, packages and a session file", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res, out)
  prov <- jsonlite::read_json(file.path(out, "provenance.json"))
  expect_setequal(names(prov$input_provenance), unlist(prov$inputs))
  expect_setequal(unlist(prov$inputs), setdiff(names(res), rep$absent))
  expect_true("packages" %in% names(prov))
  expect_gt(length(prov$packages), 0)
  expect_gt(length(readLines(file.path(out, "sessionInfo.txt"))), 0)
})

test_that("an unknown status label stops the figure and names the label", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  expect_error(
    .sw_ri_fill(c("ok", "weird"), .sw_ri_pal_status, "Status"),
    "weird"
  )
  res <- make_rep_result()
  res$signatures$files$status[1] <- "weird"
  out <- new_out()
  expect_error(sw_report_import(res["signatures"], out), "weird")
  expect_false(dir.exists(out))
})

test_that("integrity from header_only = TRUE gives a table fallback and no source CSV", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  paths <- unique(res$integrity$path)
  res$integrity <- sw_fcs_integrity_audit(
    paths,
    tot_tolerance = 0,
    time_keyword = "Time",
    header_only = TRUE
  )
  out <- new_out()
  rep <- sw_report_import(res["integrity"], out)
  expect_match(rep$figures[["integrity_tot"]], "^table fallback: ")
  expect_false("integrity_tot" %in% names(rep$source_data))
  expect_false(file.exists(file.path(out, "source_data", "integrity_tot.csv")))
  expect_true(file.exists(rep$tables[["integrity"]]))
})

test_that("every table of every input is written with the row count of its source", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res, out)
  expected <- unlist(
    lapply(names(res), \(nm) .sw_ri_tables(nm, res[[nm]])),
    recursive = FALSE
  )
  expect_setequal(names(rep$tables), names(expected))
  for (nm in names(expected)) {
    expect_equal(
      nrow(utils::read.csv(rep$tables[[nm]])),
      nrow(expected[[nm]]),
      info = nm
    )
    # No verdict or roll-up column in anything written.
    expect_false(
      any(grepl(
        "pass|fail|score|verdict",
        names(utils::read.csv(rep$tables[[nm]])),
        ignore.case = TRUE
      )),
      info = nm
    )
  }
  for (nm in names(rep$source_data)) {
    expect_false(
      any(grepl(
        "pass|fail|score|verdict",
        names(utils::read.csv(rep$source_data[[nm]])),
        ignore.case = TRUE
      )),
      info = nm
    )
  }
})

test_that("files that were not assessed stay in the compensation and marker figures", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  out <- new_out()
  rep <- sw_report_import(res[c("compensation", "signatures")], out)
  comp <- utils::read.csv(rep$source_data[["compensation_grid"]])
  expect_true(all(c("junk.fcs", "gone.fcs") %in% comp$file))
  expect_true(all(
    comp$state[comp$file %in% c("junk.fcs", "gone.fcs")] == "not_assessed"
  ))
  mk <- utils::read.csv(rep$source_data[["signatures_markers"]])
  expect_true(all(c("junk.fcs", "gone.fcs") %in% mk$file))
  expect_true(all(
    mk$state[mk$file %in% c("junk.fcs", "gone.fcs")] == "not_assessed"
  ))
})

test_that("the batch date figure labels a missing date and the header figure falls back on no rows", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  res$batch_date$crosstab$date[1] <- NA
  res$header_audit <- res$header_audit[0, ]
  out <- new_out()
  rep <- sw_report_import(res[c("batch_date", "header_audit")], out)
  dt <- utils::read.csv(rep$source_data[["batch_date_tiles"]])
  expect_true("no date" %in% dt$date)
  expect_match(rep$figures[["header_audit_tiles"]], "^table fallback: ")
  expect_true(file.exists(rep$tables[["header_audit"]]))
})

test_that("the report does not touch the input FCS files", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  res <- make_rep_result()
  paths <- unique(res$integrity$path[file.exists(res$integrity$path)])
  before <- tools::md5sum(paths)
  sw_report_import(res, new_out())
  expect_identical(tools::md5sum(paths), before)
})
