# tests/testthat/test-clean_gate.R
# Unit tests for R/clean_gate.R: margins, transform, PeacoQC, gates.

# Values for the tests only. They are not defaults of sw_clean_run().
clean_args <- list(
  debris_window = c(2e5, 1.5e6),
  singlet_wider_gate = TRUE,
  pool_samples = TRUE,
  viability_channel = "LIVE DEAD NIR-A",
  viability_method = "mindensity",
  viability_cut = NULL,
  IT_limit = 0.55,
  MAD = 6,
  cofactor = 6000
)

skip_if_no_clean_stack <- function() {
  for (p in c("flowCore", "flowWorkspace", "openCyto", "PeacoQC")) {
    skip_if_not_installed(p)
  }
}

# Synthetic frame with one planted event class per cleaning step:
# debris (low FSC-A), doublets (low FSC-H for their FSC-A), dead cells (high
# viability dye) and one event saturated in BV421-A. FSC-A is unique per event
# and is the key that tags each event with its class.
make_clean_ff <- function(
  seed = 1,
  n_cells = 3000,
  n_debris = 200,
  n_doublet = 200,
  n_dead = 300
) {
  set.seed(seed)
  a <- c(
    rnorm(n_cells, 6e5, 3e4),
    rnorm(n_debris, 5e4, 1.5e4),
    rnorm(n_doublet, 6e5, 3e4),
    rnorm(n_dead, 6e5, 3e4),
    5.5e5
  )
  # flowFrames hold single precision: put FSC-A on a 0.25 grid and make every
  # value unique, so FSC-A tags each event with its class.
  a <- round(abs(a) * 4) / 4
  while (anyDuplicated(a) > 0) {
    dup <- duplicated(a)
    a[dup] <- a[dup] + 0.25
  }
  class <- c(
    rep("cell", n_cells),
    rep("debris", n_debris),
    rep("doublet", n_doublet),
    rep("dead", n_dead),
    "margin"
  )
  n <- length(a)
  h <- a * 0.8 + rnorm(n, 0, 8000)
  is_dbl <- class == "doublet"
  h[is_dbl] <- a[is_dbl] * 0.45 + rnorm(n_doublet, 0, 8000)
  nir <- rnorm(n, 1000, 600)
  nir[class == "dead"] <- rnorm(n_dead, 3e5, 5e4)
  bv <- rnorm(n, 2000, 800)
  bv[class == "margin"] <- 4e6
  mat <- cbind(
    "FSC-A" = a,
    "FSC-H" = h,
    "SSC-A" = abs(rnorm(n, 1e5, 2e4)),
    "Time" = seq_len(n) / 100,
    "BV421-A" = bv,
    "PE-A" = rnorm(n, 3000, 900),
    "LIVE DEAD NIR-A" = nir
  )
  ord <- sample(n)
  mat <- mat[ord, ]
  mat[, "Time"] <- sort(mat[, "Time"])
  ff <- flowCore::flowFrame(mat)
  pd <- flowCore::pData(flowCore::parameters(ff))
  pd$minRange <- -1e6
  pd$maxRange <- 4e6
  flowCore::pData(flowCore::parameters(ff)) <- pd
  list(ff = ff, tags = stats::setNames(class[ord], sprintf("%.2f", a[ord])))
}

# Class of each event of a frame, by its FSC-A key.
tag_events <- function(ff, tags) {
  key <- sprintf("%.2f", flowCore::exprs(ff)[, "FSC-A"])
  unname(tags[key])
}

run_clean <- function(x, ...) {
  args <- utils::modifyList(clean_args, list(...), keep.null = TRUE)
  out <- NULL
  utils::capture.output(
    suppressWarnings(suppressMessages(
      out <- do.call(sw_clean_run, c(list(x), args))
    ))
  )
  out
}

msg_of <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

# =========================================================================
# Missing required arguments
# =========================================================================

test_that("the required-argument table matches the arguments without a default", {
  f <- formals(sw_clean_run)
  no_default <- names(f)[vapply(
    f,
    function(v) identical(v, quote(expr = )),
    logical(1)
  )]
  expect_setequal(names(.sw_clean_args), setdiff(no_default, "x"))
  expect_true(all(
    vapply(
      .sw_clean_args,
      function(a) {
        all(c("controls", "how", "reference", "example") %in% names(a))
      },
      logical(1)
    )
  ))
})

test_that("with nothing supplied one message names every missing argument", {
  msg <- msg_of(sw_clean_run("sample.fcs"))
  expect_match(msg, "needs 9 required argument", fixed = TRUE)
  for (nm in names(.sw_clean_args)) {
    expect_match(msg, paste0("\n", nm, "\n"), fixed = TRUE)
    a <- .sw_clean_args[[nm]]
    expect_match(msg, paste0("How to choose: ", a$how), fixed = TRUE)
    expect_match(msg, paste0("Controls: ", a$controls), fixed = TRUE)
    expect_match(msg, paste0("Example: ", a$example), fixed = TRUE)
    expect_match(
      msg,
      paste0(
        "Reference: benchmark 05 used ",
        a$reference,
        " on 40-colour Aurora PBMC; not a recommendation."
      ),
      fixed = TRUE
    )
  }
  # One how-to line and one labelled reference per argument.
  expect_equal(lengths(regmatches(msg, gregexpr("How to choose:", msg))), 9)
  expect_equal(
    lengths(regmatches(msg, gregexpr("not a recommendation", msg))),
    9
  )
})

test_that("each argument alone is reported by name", {
  for (nm in names(.sw_clean_args)) {
    args <- clean_args[setdiff(names(clean_args), nm)]
    args["viability_cut"] <- list(NULL)
    if (nm == "viability_cut") {
      args <- args[setdiff(names(args), "viability_cut")]
    }
    msg <- msg_of(do.call(sw_clean_run, c(list("sample.fcs"), args)))
    expect_match(msg, "needs 1 required argument", fixed = TRUE)
    expect_match(msg, paste0("\n", nm, "\n"), fixed = TRUE)
  }
})

test_that("a partial call reports only the remaining arguments", {
  msg <- msg_of(sw_clean_run(
    "sample.fcs",
    debris_window = c(2e5, 1.5e6),
    pool_samples = TRUE,
    cofactor = 6000,
    IT_limit = 0.55
  ))
  expect_match(msg, "needs 5 required argument", fixed = TRUE)
  for (nm in c(
    "singlet_wider_gate",
    "viability_channel",
    "viability_method",
    "viability_cut",
    "MAD"
  )) {
    expect_match(msg, paste0("\n", nm, "\n"), fixed = TRUE)
  }
  for (nm in c("debris_window", "pool_samples", "cofactor", "IT_limit")) {
    expect_no_match(msg, paste0("\n", nm, "\n"), fixed = TRUE)
  }
  expect_equal(lengths(regmatches(msg, gregexpr("How to choose:", msg))), 5)
})

test_that("a missing argument stops before any processing step", {
  called <- character()
  spy <- function(step) {
    force(step)
    function(...) {
      called <<- c(called, step)
      stop("step ran: ", step)
    }
  }
  local_mocked_bindings(
    sw_io_read_fcs = spy("read"),
    sw_filter_margins_peacoqc = spy("margins"),
    sw_transform_apply = spy("transform"),
    sw_qc_run = spy("qc_run"),
    sw_qc_batch = spy("qc_batch"),
    sw_gate_template = spy("template"),
    sw_gate_run = spy("gate")
  )
  expect_error(
    sw_clean_run("sample.fcs", debris_window = c(2e5, 1.5e6)),
    "needs 8 required argument"
  )
  expect_error(sw_clean_run(), "needs 9 required argument")
  expect_length(called, 0)
})

test_that("a missing argument is an error, never a warning that continues", {
  expect_no_warning(try(sw_clean_run("sample.fcs"), silent = TRUE))
  expect_error(sw_clean_run("sample.fcs"), class = "error")
})

# =========================================================================
# Argument validation, before dependency checks
# =========================================================================

test_that("invalid argument types error before any dependency check", {
  checked <- FALSE
  # Records only the requireNamespace() calls made by the package itself.
  local_mocked_bindings(
    requireNamespace = function(...) {
      if (identical(topenv(parent.frame()), asNamespace("SpectraWeaveR"))) {
        checked <<- TRUE
      }
      TRUE
    },
    .package = "base"
  )
  bad <- list(
    list(list(x = 42), "'x' must be"),
    list(list(debris_window = c(5, 1)), "'debris_window' must be"),
    list(list(debris_window = "a"), "'debris_window' must be"),
    list(list(singlet_wider_gate = "yes"), "'singlet_wider_gate' must be"),
    list(list(pool_samples = NA), "'pool_samples' must be"),
    list(list(viability_channel = c("a", "b")), "'viability_channel' must be"),
    list(list(viability_method = "gmm"), "'viability_method' must be"),
    list(
      list(viability_method = "boundary", viability_cut = NULL),
      "'viability_cut' must be a single finite number"
    ),
    list(
      list(viability_method = "mindensity", viability_cut = 2),
      "'viability_cut' must be NULL"
    ),
    list(list(IT_limit = 1.5), "'IT_limit' must be"),
    list(list(MAD = -1), "'MAD' must be"),
    list(list(cofactor = 0), "'cofactor' must be")
  )
  for (b in bad) {
    args <- utils::modifyList(
      c(list(x = "sample.fcs"), clean_args),
      b[[1]],
      keep.null = TRUE
    )
    expect_error(do.call(sw_clean_run, args), b[[2]], fixed = TRUE)
  }
  expect_false(checked)
})

test_that("a missing FCS file is reported by name", {
  args <- c(list(x = "/nonexistent/sample.fcs"), clean_args)
  expect_error(do.call(sw_clean_run, args), "not found")
})

# =========================================================================
# Step order and event accounting on synthetic data
# =========================================================================

test_that("steps run in order: margins, transform, PeacoQC, gates", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  log <- character()
  seen <- list()
  orig_margins <- sw_filter_margins_peacoqc
  orig_transform <- sw_transform_apply
  orig_qc <- sw_qc_run
  orig_gate <- sw_gate_run
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(...) {
      log <<- c(log, "margins")
      out <- orig_margins(...)
      seen$margins <<- out
      out
    },
    sw_transform_apply = function(...) {
      log <<- c(log, "transform")
      orig_transform(...)
    },
    sw_qc_run = function(ff, ...) {
      log <<- c(log, "qc")
      seen$qc_in <<- ff
      orig_qc(ff, ...)
    },
    sw_gate_run = function(fcs_files, ...) {
      log <<- c(log, "gate")
      seen$gate_in <<- fcs_files[[1]]
      orig_gate(fcs_files, ...)
    }
  )
  res <- run_clean(m$ff)

  expect_equal(log, c("margins", "transform", "qc", "gate"))

  # PeacoQC and the gates see arcsinh fluorescence and raw scatter.
  e <- flowCore::exprs(seen$qc_in)
  expect_lt(max(abs(e[, "LIVE DEAD NIR-A"])), 20)
  expect_gt(max(e[, "FSC-A"]), 1e5)
  # The gates see the QC-cleaned events, not the raw ones.
  expect_equal(nrow(seen$gate_in), nrow(seen$qc_in))

  # The saturated event is gone after the margin step, before transformation.
  expect_false("margin" %in% tag_events(seen$margins, m$tags))
  expect_equal(
    res$audit$step,
    c("raw", "margins", "qc", "non_debris", "singlets", "live")
  )
  expect_equal(res$audit$n_events[2], unname(nrow(seen$margins)))
})

test_that("each planted class is removed once and the audit adds up", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  margins_out <- NULL
  orig_margins <- sw_filter_margins_peacoqc
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(...) {
      margins_out <<- orig_margins(...)
      margins_out
    }
  )
  res <- run_clean(m$ff)

  # Only cells survive all steps.
  live <- tag_events(res$events[[1]], m$tags)
  expect_true(all(live == "cell"))
  expect_gt(length(live), 0.95 * 3000)

  # Class counts after margins, then the step that must remove each class.
  after_margins <- table(factor(
    tag_events(margins_out, m$tags),
    levels = c("cell", "debris", "doublet", "dead", "margin")
  ))
  removed <- stats::setNames(res$audit$n_removed, res$audit$step)
  expect_equal(unname(after_margins[["margin"]]), 0L)
  expect_gte(removed[["non_debris"]], after_margins[["debris"]])
  expect_gte(removed[["singlets"]], after_margins[["doublet"]])
  expect_gte(removed[["live"]], after_margins[["dead"]])
  # No class is removed twice: the cells lost to gate edges stay small.
  cells_lost <- sum(removed[c("non_debris", "singlets", "live")]) -
    sum(after_margins[c("debris", "doublet", "dead")])
  expect_lt(cells_lost, 0.03 * 3000)

  # Counts add up from raw to live.
  expect_equal(res$audit$n_events[1], unname(nrow(m$ff)))
  expect_equal(
    res$audit$n_events[1] - res$audit$n_events[6],
    sum(res$audit$n_removed, na.rm = TRUE)
  )
  expect_equal(res$audit$n_events[6], unname(nrow(res$events[[1]])))
  expect_equal(res$audit_total$n_events, res$audit$n_events)
  expect_setequal(
    names(res$audit),
    c(
      "sample",
      "step",
      "n_events",
      "pct_of_initial",
      "pct_of_previous",
      "n_removed"
    )
  )
})

test_that("audit and gates are reported per sample", {
  skip_if_no_clean_stack()
  fr <- list(a = make_clean_ff(1)$ff, b = make_clean_ff(2, n_cells = 2500)$ff)
  res <- run_clean(fr)

  expect_equal(unique(res$audit$sample), c("a", "b"))
  expect_equal(nrow(res$audit), 12)
  expect_equal(names(res$events), c("a", "b"))
  summed <- tapply(res$audit$n_events, res$audit$step, sum)
  expect_equal(
    as.vector(summed[res$audit_total$step]),
    res$audit_total$n_events
  )
  expect_equal(res$qc$sample, c("a", "b"))
  expect_equal(unique(res$qc$IT_limit), 0.55)
  expect_equal(unique(res$qc$MAD), 6)
  expect_equal(nrow(res$gates), 6)
  expect_setequal(res$gates$gate, c("nonDebris", "singlets", "live"))
  expect_equal(
    res$gates$keeps[res$gates$gate == "live"],
    c("outside", "outside")
  )
  expect_equal(res$parameters$cofactor, 6000)
  expect_equal(res$template$alias, c("nonDebris", "singlets", "live"))
  expect_false(any(c("CD3neg", "lymphocytes") %in% res$template$alias))
})

test_that("a fixed boundary on the viability channel gates live cells", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  res <- run_clean(m$ff, viability_method = "boundary", viability_cut = 3)
  expect_true(all(tag_events(res$events[[1]], m$tags) == "cell"))
  expect_match(res$template$gating_args[3], "max = 3", fixed = TRUE)
})

test_that("a sample emptied by a step stops the run", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(ff, ...) ff[0, ]
  )
  expect_error(run_clean(m$ff), "No events left after step 'margins'")
})

test_that("a gating failure is reported with the settings to check", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  expect_error(
    run_clean(m$ff, debris_window = c(1e9, Inf)),
    "Gating failed (check debris_window",
    fixed = TRUE
  )
})

test_that("a viability channel absent from a sample is an error", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  expect_error(run_clean(m$ff, viability_channel = "ZOMBIE-A"), "lacks channel")
})

test_that("sw_gate_run accepts an in-memory flowSet", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  fs <- flowCore::flowSet(list(s1 = m$ff))
  tmpl <- tempfile(fileext = ".csv")
  on.exit(unlink(tmpl))
  suppressMessages(sw_gate_template(
    tmpl,
    "lymphocyte",
    debris_min = 2e5,
    singlet_wider_gate = TRUE,
    flowclust_k = 2,
    flowclust_target = 1,
    pool_samples = TRUE
  ))
  gs <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    gs <- sw_gate_run(fs, tmpl, transform_channels = NULL)
  )))
  expect_true(methods::is(gs, "GatingSet"))
})
