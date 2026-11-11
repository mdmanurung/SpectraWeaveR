# tests/testthat/test-clean_gate.R
# Unit tests for R/clean_gate.R: margins, transform, PeacoQC, gates.

# Values for the tests only. They are not defaults of sw_clean_run().
clean_args <- list(
  margin_range = c(-21761, 1e7),
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
  n_dead = 300,
  fsc_scale = 1
) {
  set.seed(seed)
  # fsc_scale moves every FSC-A peak, so samples can have different valleys.
  a <- c(
    rnorm(n_cells, 6e5 * fsc_scale, 3e4),
    rnorm(n_debris, 5e4 * fsc_scale, 1.5e4),
    rnorm(n_doublet, 6e5 * fsc_scale, 3e4),
    rnorm(n_dead, 6e5 * fsc_scale, 3e4),
    5.5e5 * fsc_scale
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

# Same as run_clean() but warnings from sw_clean_run() reach the caller.
run_clean_warn <- function(x, ...) {
  args <- utils::modifyList(clean_args, list(...), keep.null = TRUE)
  out <- NULL
  utils::capture.output(
    suppressMessages(out <- do.call(sw_clean_run, c(list(x), args)))
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
  expect_match(msg, "needs 10 required argument", fixed = TRUE)
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
  expect_equal(lengths(regmatches(msg, gregexpr("How to choose:", msg))), 10)
  expect_equal(
    lengths(regmatches(msg, gregexpr("not a recommendation", msg))),
    10
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
  expect_match(msg, "needs 6 required argument", fixed = TRUE)
  for (nm in c(
    "margin_range",
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
  expect_equal(lengths(regmatches(msg, gregexpr("How to choose:", msg))), 6)
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
    "needs 9 required argument"
  )
  expect_error(sw_clean_run(), "needs 10 required argument")
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
    list(list(margin_range = c(5, 1)), "'margin_range' must be"),
    list(list(margin_range = c(5, 5)), "'margin_range' must be"),
    list(list(margin_range = c(-1, Inf)), "'margin_range' must be"),
    list(list(margin_range = c(NA, 1)), "'margin_range' must be"),
    list(list(margin_range = 1), "'margin_range' must be"),
    list(list(margin_range = c(0, 1, 2)), "'margin_range' must be"),
    list(list(margin_range = c("a", "b")), "'margin_range' must be"),
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
      seen$qc_channels <<- list(...)$channels
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
  # PeacoQC evaluates the fluorescence channels, not the Original_ID counter.
  expect_setequal(seen$qc_channels, c("BV421-A", "PE-A", "LIVE DEAD NIR-A"))

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

test_that("margin_range reaches the margin step as the fluorescence range", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  seen <- NULL
  orig_margins <- sw_filter_margins_peacoqc
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(
      x,
      channel_specifications = NULL,
      ...
    ) {
      seen <<- channel_specifications
      orig_margins(x, channel_specifications = channel_specifications, ...)
    }
  )
  res <- run_clean(m$ff)
  expect_equal(seen, list(AllFluoChannels = clean_args$margin_range))
  expect_equal(res$parameters$margin_range, clean_args$margin_range)
})

test_that("margin_range sets the range of fluorescence channels, scatter keeps the FCS range", {
  skip_if_no_clean_stack()
  # PeacoQC::RemoveMargins() always drops the observed minimum and maximum of
  # every channel (strict comparisons). Two planted events, hi_all and lo_all,
  # hold the extremes of every channel, so no other event is an observed
  # extreme and only margin_range or the FCS range can remove the rest.
  m <- make_clean_ff()
  mat <- flowCore::exprs(m$ff)
  planted <- function(...) {
    row <- mat[1, , drop = FALSE]
    row[1, ] <- c(
      "FSC-A" = 6e5,
      "FSC-H" = 4.8e5,
      "SSC-A" = 1e5,
      "Time" = 0,
      "BV421-A" = 2000,
      "PE-A" = 3000,
      "LIVE DEAD NIR-A" = 1000
    )
    over <- c(...)
    row[1, names(over)] <- over
    row
  }
  extra <- rbind(
    # Observed extremes of every channel: removed whatever margin_range is.
    planted(
      "FSC-A" = 3e6, "FSC-H" = 3e6, "SSC-A" = 3e6,
      "BV421-A" = 5e7, "PE-A" = 5e7, "LIVE DEAD NIR-A" = 5e7
    ),
    planted(
      "FSC-A" = 10, "FSC-H" = 10, "SSC-A" = 10,
      "BV421-A" = -3e5, "PE-A" = -3e5, "LIVE DEAD NIR-A" = -3e5
    ),
    # Scatter above the FCS range of FSC-A (1e6), not the observed maximum:
    # removed whatever margin_range is.
    planted("FSC-A" = 2e6),
    # Scatter inside the FCS range: kept, also when margin_range is narrow.
    planted("FSC-A" = 9e5),
    # Fluorescence above the narrow maximum (5e5), inside the wide one (1e7).
    planted("BV421-A" = 1e6),
    planted("BV421-A" = 1e6),
    planted("BV421-A" = 1e6),
    # Fluorescence above both maxima, not the observed maximum.
    planted("PE-A" = 2e7),
    planted("PE-A" = 2e7),
    # Fluorescence below the narrow minimum (-1e4), above the wide one (-1e5).
    planted("PE-A" = -5e4),
    planted("PE-A" = -5e4),
    # Fluorescence below both minima, not the observed minimum.
    planted("BV421-A" = -2e5),
    planted("BV421-A" = -2e5),
    # Inside both ranges: kept.
    planted("BV421-A" = 5e4),
    planted("BV421-A" = 5e4)
  )
  extra[, "Time"] <- max(mat[, "Time"]) + seq_len(nrow(extra)) / 100
  ff <- flowCore::flowFrame(rbind(mat, extra))
  pd <- flowCore::pData(flowCore::parameters(ff))
  pd$minRange <- -1e6
  pd$maxRange <- 4e6
  pd$maxRange[pd$name == "FSC-A"] <- 1e6
  flowCore::pData(flowCore::parameters(ff)) <- pd

  # Run sw_clean_run() up to the margin step: its real call of
  # sw_filter_margins_peacoqc() counts the events removed, then the stub that
  # replaces the transformation stops the run.
  removed_by_margins <- function(margin_range) {
    n_kept <- NULL
    orig_margins <- sw_filter_margins_peacoqc
    local_mocked_bindings(
      sw_filter_margins_peacoqc = function(x, ...) {
        out <- orig_margins(x, ...)
        n_kept <<- nrow(out)
        out
      },
      sw_transform_apply = function(...) stop("stop after the margin step")
    )
    expect_error(
      run_clean(ff, margin_range = margin_range),
      "stop after the margin step"
    )
    unname(nrow(ff) - n_kept)
  }
  # Wide range: the two extreme events, the FSC-A event above the FCS range,
  # the two PE events at 2e7 and the two BV421 events at -2e5. The FCS range
  # (-1e6, 4e6) would remove the two PE events but keep the two BV421 events
  # at -2e5 and remove the make_clean_ff event at BV421 4e6.
  expect_equal(removed_by_margins(c(-1e5, 1e7)), 7)
  # Narrow range: also the three BV421 events at 1e6, the make_clean_ff event
  # at BV421 4e6 and the two PE events at -5e4. The 3900 scatter events near
  # 6e5 lie above the narrow maximum, so they show that scatter keeps the FCS
  # range.
  expect_equal(removed_by_margins(c(-1e4, 5e5)), 13)
})

test_that("a margin_range minimum above 0 warns and gives the same audit as 0", {
  skip_if_no_clean_stack()
  m <- make_clean_ff()
  warnings_of <- function(margin_range) {
    msgs <- character()
    res <- withCallingHandlers(
      run_clean_warn(m$ff, margin_range = margin_range),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(res = res, msgs = msgs)
  }
  positive <- warnings_of(c(500, 1e7))
  zero <- warnings_of(c(0, 1e7))

  # The warning names margin_range and the value, and says PeacoQC reads it as 0.
  hit <- grepl("'margin_range' lower limit is 500", positive$msgs, fixed = TRUE)
  expect_equal(sum(hit), 1)
  expect_match(positive$msgs[hit], "reads a lower limit above 0 as 0", fixed = TRUE)
  # A lower limit of 0 raises no such warning.
  expect_false(any(grepl("'margin_range' lower limit", zero$msgs, fixed = TRUE)))

  # The run continues and the counts equal those of margin_range[1] = 0.
  expect_equal(positive$res$audit, zero$res$audit)
  expect_equal(positive$res$parameters$margin_range, c(500, 1e7))
})

# Frame lists for the margin and PeacoQC error tests: b is the one that fails.
two_samples <- function() {
  list(a = make_clean_ff(1)$ff, b = make_clean_ff(2, n_cells = 2500)$ff)
}

test_that("a sample emptied by the margin step is named with counts and margin_range", {
  skip_if_no_clean_stack()
  fr <- two_samples()
  n_b <- nrow(fr$b)
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(x, ...) {
      if (nrow(x) == n_b) x[0, ] else x
    }
  )
  msg <- msg_of(run_clean(fr))
  expect_match(
    msg,
    "No events left after step 'margins' in sample(s): b",
    fixed = TRUE
  )
  expect_match(
    msg,
    sprintf("b: %d events before, 0 after, 100%% removed", n_b),
    fixed = TRUE
  )
  expect_match(msg, "margin_range = c(-21761, 1e+07)", fixed = TRUE)
  expect_no_match(msg, "a: ", fixed = TRUE)
})

test_that("every sample the margin step emptied is listed in one error", {
  skip_if_no_clean_stack()
  fr <- two_samples()
  local_mocked_bindings(sw_filter_margins_peacoqc = function(x, ...) x[0, ])
  msg <- msg_of(run_clean(fr))
  expect_match(msg, "in sample(s): a, b", fixed = TRUE)
  expect_match(msg, "a: ", fixed = TRUE)
  expect_match(msg, "b: ", fixed = TRUE)
})

test_that("a PeacoQC error names the sample, the event count and the original message", {
  skip_if_no_clean_stack()
  fr <- two_samples()
  n_b <- nrow(fr$b)
  n_to_qc <- NULL
  orig_qc <- sw_qc_run
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(x, ...) x,
    sw_qc_run = function(ff, ...) {
      if (nrow(ff) == n_b) {
        n_to_qc <<- nrow(ff)
        stop("need at least four unique x values", call. = FALSE)
      }
      orig_qc(ff, ...)
    }
  )
  msg <- msg_of(run_clean(fr))
  expect_match(msg, "PeacoQC failed on sample 'b'", fixed = TRUE)
  expect_match(msg, sprintf("with %d events handed to it", n_b), fixed = TRUE)
  expect_match(
    msg,
    "PeacoQC said: need at least four unique x values",
    fixed = TRUE
  )
  expect_equal(n_to_qc, n_b)
})

test_that("a sample left with a few events reaches PeacoQC; no event-count threshold stops it", {
  skip_if_no_clean_stack()
  fr <- two_samples()
  n_b <- nrow(fr$b)
  handed <- NULL
  orig_qc <- sw_qc_run
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(x, ...) {
      if (nrow(x) == n_b) x[1:5, ] else x
    },
    sw_qc_run = function(ff, ...) {
      if (nrow(ff) == 5) {
        handed <<- nrow(ff)
        stop("need at least four unique x values", call. = FALSE)
      }
      orig_qc(ff, ...)
    }
  )
  msg <- msg_of(run_clean(fr))
  expect_equal(handed, 5L)
  expect_match(msg, "PeacoQC failed on sample 'b' with 5 events", fixed = TRUE)
  expect_no_match(msg, "No events left", fixed = TRUE)
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

# Cut of the nonDebris gate of one sample, from the gates table.
debris_cut_of <- function(res, sample = 1) {
  d <- res$gates[res$gates$gate == "nonDebris", ]
  g <- d$definition[sample]
  as.numeric(sub("^FSC-A: \\[([^,]+),.*$", "\\1", g))
}

test_that("debris_window bounds the data of the density: the cut is a valley inside it", {
  skip_if_no_clean_stack()
  # Debris peak at 5e4 (SD 1.5e4) and cell peak at 6e5 (SD 3e4): the valley
  # lies between them. The window holds both peaks.
  m <- make_clean_ff()
  window <- c(1e4, 1.5e6)
  res <- suppressWarnings(run_clean(m$ff, debris_window = window))

  cut <- debris_cut_of(res)
  expect_gt(cut, window[1])
  expect_lt(cut, window[2])
  # The cut sits between the two planted peaks, in the valley.
  expect_gt(cut, 1e5)
  expect_lt(cut, 5e5)
  # mindensity was called with both limits.
  expect_match(
    res$template$gating_args[1],
    "min = 10000, max = 1500000",
    fixed = TRUE
  )
  # A valley was found, so the window floor is not the cut.
  d <- res$gates[res$gates$gate == "nonDebris", ]
  expect_false(d$cut_at_window_floor)
  expect_true(all(is.na(res$gates$cut_at_window_floor[
    res$gates$gate != "nonDebris"
  ])))
  # All planted debris is below the cut.
  expect_true(all(tag_events(res$events[[1]], m$tags) == "cell"))
  # No floor warning is raised.
  floor_warned <- FALSE
  withCallingHandlers(
    run_clean_warn(m$ff, debris_window = window),
    warning = function(w) {
      if (grepl("lowest event inside debris_window", conditionMessage(w))) {
        floor_warned <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_false(floor_warned)
})

test_that("events above the window max are kept by the debris gate", {
  skip_if_no_clean_stack()
  # max = 5.9e5 truncates the data of the density only. The gate is
  # c(cut, Inf), and the cells above 5.9e5 stay in the live events.
  m <- make_clean_ff()
  window <- c(1e4, 5.9e5)
  res <- suppressWarnings(run_clean(m$ff, debris_window = window))
  expect_match(res$template$gating_args[1], "max = 590000", fixed = TRUE)
  d <- res$gates[res$gates$gate == "nonDebris", ]
  expect_match(d$definition, ", Inf]", fixed = TRUE)
  live_fsc <- flowCore::exprs(res$events[[1]])[, "FSC-A"]
  expect_gt(sum(live_fsc > window[2]), 1000)
})

test_that("a window that leaves out the valley is reported, not silent", {
  skip_if_no_clean_stack()
  # The window starts above the debris peak and the valley, so the density
  # has one peak. openCyto then puts the cut on the lowest event in the window.
  m <- make_clean_ff()
  window <- c(3e5, 1.5e6)
  warned <- NULL
  res <- withCallingHandlers(
    run_clean_warn(m$ff, debris_window = window),
    warning = function(w) {
      if (grepl("lowest event inside debris_window", conditionMessage(w))) {
        warned <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warned, "no density valley inside the window", fixed = TRUE)
  d <- res$gates[res$gates$gate == "nonDebris", ]
  expect_true(d$cut_at_window_floor)
  x <- flowCore::exprs(m$ff)[, "FSC-A"]
  expect_equal(
    debris_cut_of(res),
    min(x[x >= window[1] & x <= window[2]]),
    tolerance = 1e-4
  )
})

# =========================================================================
# Pooling: one shared cut or one cut per sample
# =========================================================================

# Definition column of one gate, one entry per sample.
gate_definitions <- function(res, gate) {
  res$gates$definition[res$gates$gate == gate]
}

test_that("the template pools through groupBy: the number of samples when TRUE, empty when FALSE", {
  for (method in c("mindensity", "boundary")) {
    cut <- if (method == "boundary") 3 else NULL
    pooled <- .sw_clean_template(
      c(1e4, 1.5e6),
      TRUE,
      TRUE,
      "LIVE DEAD NIR-A",
      method,
      cut,
      n_samples = 4
    )
    expect_true(all(pooled$collapseDataForGating))
    expect_equal(as.integer(pooled$groupBy), rep(4L, 3))

    single <- .sw_clean_template(
      c(1e4, 1.5e6),
      TRUE,
      FALSE,
      "LIVE DEAD NIR-A",
      method,
      cut,
      n_samples = 4
    )
    expect_false(any(single$collapseDataForGating))
    expect_true(all(is.na(single$groupBy)))
  }
})

test_that("pool_samples = TRUE gives one shared cut for every gate", {
  skip_if_no_clean_stack()
  # Peaks of the three samples sit at different FSC-A, so their own valleys
  # differ; a pooled fit must give them the same cut.
  fr <- list(
    a = make_clean_ff(1, fsc_scale = 0.9)$ff,
    b = make_clean_ff(2)$ff,
    c = make_clean_ff(3, fsc_scale = 1.1)$ff
  )
  res <- suppressWarnings(run_clean(
    fr,
    debris_window = c(1e4, 1.5e6),
    pool_samples = TRUE
  ))
  expect_equal(unique(as.integer(res$template$groupBy)), 3L)
  for (g in c("nonDebris", "singlets", "live")) {
    expect_length(unique(gate_definitions(res, g)), 1)
  }
})

test_that("pool_samples = FALSE gives each sample its own cut", {
  skip_if_no_clean_stack()
  fr <- list(
    a = make_clean_ff(1, fsc_scale = 0.9)$ff,
    b = make_clean_ff(2)$ff,
    c = make_clean_ff(3, fsc_scale = 1.1)$ff
  )
  res <- suppressWarnings(run_clean(
    fr,
    debris_window = c(1e4, 1.5e6),
    pool_samples = FALSE
  ))
  expect_true(all(is.na(res$template$groupBy)))
  for (g in c("nonDebris", "singlets", "live")) {
    expect_length(unique(gate_definitions(res, g)), 3)
  }
})

test_that("the window floor is one pooled value when pooling and one value per sample otherwise", {
  skip_if_not_installed("flowCore")
  frame_of <- function(fsc) {
    flowCore::flowFrame(cbind("FSC-A" = fsc, "FSC-H" = fsc))
  }
  fr <- list(a = frame_of(c(10, 20, 30, 500)), b = frame_of(c(5, 25, 40, 90)))
  # Pooled: one value for every sample, the lowest event of all samples.
  expect_equal(.sw_clean_window_floor(fr, c(8, 90), TRUE), c(10, 10))
  expect_equal(.sw_clean_window_floor(fr, c(1, 90), TRUE), c(5, 5))
  # Not pooled: the lowest event of each sample.
  expect_equal(.sw_clean_window_floor(fr, c(1, 90), FALSE), c(a = 10, b = 5))
  # Window limits are inclusive, as in the openCyto truncation.
  expect_equal(.sw_clean_window_floor(fr, c(10, 90), FALSE), c(a = 10, b = 25))
})

# Frames handed to the gates, and the lowest FSC-A of each inside a window.
spy_gate_input <- function() {
  env <- new.env()
  orig_gate <- sw_gate_run
  list(
    env = env,
    mock = function(fcs_files, ...) {
      env$fs <- fcs_files
      orig_gate(fcs_files, ...)
    }
  )
}

floor_of <- function(fs, window, pool) {
  v <- lapply(flowCore::sampleNames(fs), function(sn) {
    x <- flowCore::exprs(fs[[sn]])[, "FSC-A"]
    x[x >= window[1] & x <= window[2]]
  })
  if (pool) rep(min(unlist(v)), length(v)) else vapply(v, min, numeric(1))
}

test_that("cut_at_window_floor compares the one pooled cut with the pooled floor", {
  skip_if_no_clean_stack()
  # The window starts above the debris peak, so each sample and their pool
  # have one peak: the cut sits on the lowest event inside the window.
  fr <- list(
    a = make_clean_ff(1)$ff,
    b = make_clean_ff(2)$ff,
    c = make_clean_ff(3)$ff
  )
  window <- c(3e5, 1.5e6)
  spy <- spy_gate_input()
  local_mocked_bindings(sw_gate_run = spy$mock)
  res <- run_clean(fr, debris_window = window, pool_samples = TRUE)

  pooled_floor <- floor_of(spy$env$fs, window, pool = TRUE)[1]
  own_floor <- floor_of(spy$env$fs, window, pool = FALSE)
  # The samples have different own floors, so the pooled floor tells the cases apart.
  expect_gt(length(unique(own_floor)), 1)

  d <- res$gates[res$gates$gate == "nonDebris", ]
  cuts <- vapply(
    seq_len(nrow(d)),
    function(i) debris_cut_of(res, i),
    numeric(1)
  )
  expect_equal(cuts, rep(pooled_floor, 3), tolerance = 1e-4)
  expect_true(all(d$cut_at_window_floor))
})

test_that("cut_at_window_floor compares each per-sample cut with that sample's floor", {
  skip_if_no_clean_stack()
  fr <- list(
    a = make_clean_ff(1)$ff,
    b = make_clean_ff(2)$ff,
    c = make_clean_ff(3)$ff
  )
  window <- c(3e5, 1.5e6)
  spy <- spy_gate_input()
  local_mocked_bindings(sw_gate_run = spy$mock)
  res <- run_clean(fr, debris_window = window, pool_samples = FALSE)

  own_floor <- floor_of(spy$env$fs, window, pool = FALSE)
  pooled_floor <- floor_of(spy$env$fs, window, pool = TRUE)[1]
  d <- res$gates[res$gates$gate == "nonDebris", ]
  cuts <- vapply(
    seq_len(nrow(d)),
    function(i) debris_cut_of(res, i),
    numeric(1)
  )
  expect_equal(cuts, own_floor, tolerance = 1e-4)
  # Against the pooled floor, samples whose own floor is higher would look
  # like valley cuts. Against their own floor all three are at the floor.
  expect_gt(sum(own_floor > pooled_floor), 0)
  expect_true(all(d$cut_at_window_floor))
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

test_that("the dependency-check probe sees the package's own calls", {
  skip_if_not_installed("flowCore")
  checked <- FALSE
  local_mocked_bindings(
    requireNamespace = function(...) {
      if (identical(topenv(parent.frame()), asNamespace("SpectraWeaveR"))) {
        checked <<- TRUE
      }
      TRUE
    },
    .package = "base"
  )
  local_mocked_bindings(
    sw_filter_margins_peacoqc = function(...) stop("margins reached")
  )
  expect_error(
    do.call(sw_clean_run, c(list(x = make_clean_ff()$ff), clean_args)),
    "margins reached"
  )
  expect_true(checked)
})
