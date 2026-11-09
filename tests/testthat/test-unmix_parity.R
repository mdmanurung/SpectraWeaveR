# tests/testthat/test-unmix_parity.R
# Parity of the AutoSpectral wrappers in R/unmix.R with AutoSpectral 1.8.8.
# Each .sw_args_*() builder returns the argument list one AutoSpectral call
# gets. These tests check every name against the target function's formals:
# parsed from an AutoSpectral source clone (SW_AUTOSPECTRAL_SRC), and from the
# installed package when there is one. `...` is never counted as a match,
# because it would accept any name.

# Source clone --------------------------------------------------------------

AS_SRC_DEFAULT <- paste0(
  "/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-SpectraWeaveR/",
  "2e2f34af-592d-4e69-8db9-a71ae68d552e/scratchpad/src/AutoSpectral"
)
as_src <- Sys.getenv("SW_AUTOSPECTRAL_SRC", AS_SRC_DEFAULT)
# Accept the repository root or its R/ folder
as_r_dir <- if (dir.exists(file.path(as_src, "R"))) {
  file.path(as_src, "R")
} else {
  as_src
}

# Index every `name <- function(...)` in the clone, parsed, not run.
# Each entry is the `function` call: [[2]] the formals, [[3]] the body.
parse_as_defs <- function(r_dir) {
  files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)
  out <- list()
  for (f in files) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    for (e in exprs) {
      is_assign <- is.call(e) &&
        (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))
      if (!is_assign || !is.name(e[[2]])) {
        next
      }
      rhs <- e[[3]]
      if (is.call(rhs) && identical(rhs[[1]], as.name("function"))) {
        out[[as.character(e[[2]])]] <- rhs
      }
    }
  }
  out
}

as_defs <- if (dir.exists(as_r_dir)) parse_as_defs(as_r_dir) else list()

skip_if_no_src <- function() {
  skip_if(
    length(as_defs) == 0,
    paste("AutoSpectral source clone not readable at", as_r_dir)
  )
}

# Formal names an AutoSpectral function accepts, `...` excluded
src_accepts <- function(fn) {
  f <- as_defs[[fn]]
  if (is.null(f)) {
    stop("No definition of ", fn, "() in ", as_r_dir, call. = FALSE)
  }
  setdiff(names(f[[2]]), "...")
}

# Evaluate a formal's default, e.g. the choices vector of `method`
src_default <- function(fn, arg) {
  eval(as_defs[[fn]][[2]][[arg]], baseenv())
}

# Builder calls with every optional argument set --------------------------

mock_asp <- list(
  figures = FALSE,
  unmixed.fcs.dir = tempdir(),
  cytometer = "Aurora"
)
mock_setup <- structure(
  list(
    asp = mock_asp,
    control_file = "f.csv",
    control_dir = ".",
    workflow = "automated"
  ),
  class = "sw_setup"
)
mock_spectra <- matrix(
  1,
  2,
  3,
  dimnames = list(c("A", "B"), c("d1", "d2", "d3"))
)
mock_af <- matrix(1, 2, 3)

unmix_full <- list(
  spectra = mock_spectra,
  asp = mock_asp,
  flow_control = list(),
  method = "AutoSpectral",
  af_spectra = mock_af,
  spectra_variants = list(),
  output_dir = tempdir(),
  file_suffix = "_x",
  speed = "fast",
  parallel = TRUE,
  threads = 2,
  chunk_size = 1e5
)

# One entry per AutoSpectral call the wrappers make
built_full <- list(
  get.autospectral.param = .sw_args_get_autospectral_param("aurora", FALSE),
  create.control.file = .sw_args_create_control_file(
    ".",
    mock_asp,
    tempdir(),
    TRUE
  ),
  check.control.file = .sw_args_check_control_file(
    ".",
    "f.csv",
    mock_asp,
    TRUE
  ),
  get.spectra.automated = .sw_args_get_spectra_automated(mock_setup),
  reload.flow.control = .sw_args_reload_flow_control(mock_setup),
  define.flow.control = .sw_args_define_flow_control(
    mock_setup,
    "landmarks",
    list(),
    FALSE,
    2
  ),
  clean.controls = .sw_args_clean_controls(list(), mock_asp, TRUE, FALSE, 2),
  get.fluorophore.spectra = .sw_args_get_fluorophore_spectra(
    list(),
    mock_asp,
    TRUE
  ),
  get.af.spectra = .sw_args_get_af_spectra(
    "u.fcs",
    mock_asp,
    mock_spectra,
    title = "AF lung",
    som_dim = 10,
    refine = TRUE,
    parallel = TRUE,
    threads = 2
  ),
  get.spectral.variants = .sw_args_get_spectral_variants(
    mock_setup,
    mock_spectra,
    n_cells = 5000,
    som_dim = 5,
    parallel = TRUE,
    threads = 2
  ),
  unmix.fcs = do.call(
    .sw_args_unmix_fcs,
    c(list(fcs_file = "s.fcs"), unmix_full)
  ),
  unmix.folder = do.call(
    .sw_args_unmix_folder,
    c(list(fcs_dir = "."), unmix_full)
  )
)

# Tests that need no AutoSpectral -----------------------------------------

test_that("builders drop unset arguments so AutoSpectral defaults apply", {
  af_min <- .sw_args_get_af_spectra("u.fcs", list(figures = TRUE), mock_spectra)
  expect_named(af_min, c("unstained.sample", "asp", "spectra"))

  var_min <- .sw_args_get_spectral_variants(
    structure(
      list(asp = list(figures = TRUE), control_file = "f", control_dir = "."),
      class = "sw_setup"
    ),
    mock_spectra
  )
  expect_named(var_min, c("control.dir", "control.def.file", "asp", "spectra"))

  run_min <- .sw_args_unmix_fcs(
    "s.fcs",
    spectra = mock_spectra,
    asp = mock_asp,
    flow_control = list(),
    method = "OLS"
  )
  expect_named(
    run_min,
    c("fcs.file", "spectra", "asp", "flow.control", "method")
  )
  expect_false("threads" %in% names(run_min))

  sa_min <- .sw_args_get_spectra_automated(
    structure(
      list(asp = list(figures = TRUE), control_file = "f", control_dir = "."),
      class = "sw_setup"
    )
  )
  expect_named(sa_min, c("control.dir", "control.def.file", "asp"))

  for (fn in names(built_full)) {
    expect_false(any(vapply(built_full[[fn]], is.null, logical(1))), info = fn)
  }
})

test_that("figures = FALSE from the setup reaches AutoSpectral, figures = TRUE is left out", {
  expect_false(built_full$get.spectra.automated$figures)
  expect_false(built_full$get.af.spectra$figures)
  user_fig <- .sw_args_get_spectra_automated(mock_setup, list(figures = TRUE))
  expect_true(user_fig$figures)
})

test_that("variants builder passes no af.spectra or refine", {
  expect_false(any(
    c("af.spectra", "refine") %in% names(built_full$get.spectral.variants)
  ))
})

test_that("dots must be named and must not clash with wrapper arguments", {
  expect_error(.sw_add_dots(list(a = 1), list(2), "f"), "must be named")
  expect_error(
    .sw_add_dots(list(a = 1), list(a = 2), "f"),
    "set by the wrapper"
  )
  expect_identical(
    .sw_add_dots(list(a = 1), list(b = 2), "f"),
    list(a = 1, b = 2)
  )
})

test_that("the call guard rejects names a function with ... would swallow", {
  # stats::median has `...`, so R itself would accept `bogus`
  expect_error(
    .sw_call_autospectral("median", list(x = 1:3, bogus = 1), ns = "stats"),
    "has no argument\\(s\\): bogus"
  )
  expect_identical(
    .sw_call_autospectral("median", list(x = c(1, 3, 5)), ns = "stats"),
    3
  )
})

# Tests against the source clone ------------------------------------------

test_that("every argument the wrappers pass exists in the 1.8.8 source formals", {
  skip_if_no_src()
  for (fn in names(built_full)) {
    unknown <- setdiff(names(built_full[[fn]]), src_accepts(fn))
    expect_identical(unknown, character(0), info = fn)
  }
})

test_that("wrapper method and speed choices match unmix.fcs and unmix.folder", {
  skip_if_no_src()
  for (fn in c("unmix.fcs", "unmix.folder")) {
    expect_identical(
      eval(formals(sw_unmix_run)$method),
      src_default(fn, "method"),
      info = fn
    )
    expect_identical(
      eval(formals(sw_unmix_run)$speed),
      src_default(fn, "speed"),
      info = fn
    )
  }
  expect_identical(
    eval(formals(sw_unmix_pipeline)$method),
    src_default("unmix.fcs", "method")
  )
})

test_that("legacy gating choices are the ones define.flow.control accepts", {
  skip_if_no_src()
  # Same set; the wrapper's order (landmarks first) differs on purpose
  expect_setequal(
    eval(formals(sw_unmix_prepare)$gating_system),
    src_default("define.flow.control", "gating.system")
  )
})

test_that("every wrapper cytometer is one get.autospectral.param accepts", {
  skip_if_no_src()
  # The valid set sits in the function body: quoted names (cytometer.direct)
  # and list names (cytometer.synonym, e.g. `a8 = list(...)`)
  body_txt <- paste(
    deparse(as_defs[["get.autospectral.param"]][[3]]),
    collapse = "\n"
  )
  for (cy in .valid_cytometers()) {
    as_string <- grepl(paste0('"', cy, '"'), body_txt, fixed = TRUE)
    as_list_name <- grepl(paste0("\\b", cy, " = list\\("), body_txt)
    expect_true(as_string || as_list_name, info = cy)
  }
})

# Tests against the installed package -------------------------------------

test_that("every argument the wrappers pass exists in installed AutoSpectral", {
  skip_if_not_installed("AutoSpectral")
  for (fn in names(built_full)) {
    accepted <- setdiff(
      names(formals(getExportedValue("AutoSpectral", fn))),
      "..."
    )
    unknown <- setdiff(names(built_full[[fn]]), accepted)
    expect_identical(unknown, character(0), info = fn)
  }
  expect_identical(
    eval(formals(sw_unmix_run)$method),
    eval(formals(AutoSpectral::unmix.fcs)$method)
  )
})

# Mocked end-to-end run ---------------------------------------------------

test_that("sw_unmix_pipeline makes the 1.8.8 call sequence with valid names", {
  skip_if_no_src()
  tmp <- withr::local_tempdir()
  ctrl_dir <- file.path(tmp, "controls")
  sample_dir <- file.path(tmp, "samples")
  dir.create(ctrl_dir)
  dir.create(sample_dir)
  ctrl_file <- file.path(tmp, "controls.csv")
  writeLines("filename", ctrl_file)
  unstained <- file.path(tmp, "unstained.fcs")
  writeLines("x", unstained)
  writeLines("x", file.path(sample_dir, "s1.fcs"))

  # Record each AutoSpectral call instead of running it
  calls <- list()
  fake_call <- function(fn, args, ns = "AutoSpectral") {
    calls[[length(calls) + 1]] <<- list(fn = fn, args = args)
    switch(
      fn,
      get.autospectral.param = list(figures = TRUE, cytometer = "Aurora"),
      check.control.file = invisible(TRUE),
      get.spectra.automated = mock_spectra,
      reload.flow.control = list(),
      get.af.spectra = mock_af,
      get.spectral.variants = list(variants = list(A = mock_spectra)),
      NULL
    )
  }
  mocks <- list(
    .check_autospectral = function() invisible(TRUE),
    .sw_call_autospectral = fake_call
  )
  wrapper_env <- environment(sw_unmix_pipeline)
  if (isNamespace(wrapper_env)) {
    # Package context (devtools::test, R CMD check)
    do.call(local_mocked_bindings, mocks)
  } else {
    # Files sourced into an environment: swap and restore by hand
    originals <- mget(names(mocks), envir = wrapper_env)
    list2env(mocks, envir = wrapper_env)
    withr::defer(list2env(originals, envir = wrapper_env))
  }

  suppressMessages(sw_unmix_pipeline(
    control_dir = ctrl_dir,
    sample_input = sample_dir,
    unstained_fcs = unstained,
    control_file = ctrl_file,
    output_dir = file.path(tmp, "out")
  ))

  fns <- vapply(calls, `[[`, character(1), "fn")
  expect_identical(
    fns,
    c(
      "get.autospectral.param",
      "check.control.file",
      "get.spectra.automated",
      "reload.flow.control",
      "get.af.spectra",
      "get.spectral.variants",
      "unmix.folder"
    )
  )
  for (cl in calls) {
    unknown <- setdiff(names(cl$args), src_accepts(cl$fn))
    expect_identical(unknown, character(0), info = cl$fn)
  }
  # Defaults left to AutoSpectral: no threads, parallel, som.dim or n.cells
  passed <- unlist(lapply(calls, \(cl) names(cl$args)))
  expect_false(any(
    c("threads", "parallel", "som.dim", "n.cells", "chunk.size") %in% passed
  ))
})
