# tests/testthat/test-targets_backend.R
# Unit tests for R/targets_backend.R — targets execution backend.
#
# These tests are skipped unless both `S7` and `targets` are installed.
# Each test uses an isolated tempdir for its store/script so reruns and
# parallel test runs do not interfere.

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# Top-level (non-closure) functions so they can be referenced from a
# `targets` callr child without dragging in the test environment.
.swt_times2 <- function(x) x * 2
.swt_plus1  <- function(x) x + 1
.swt_plus10 <- function(x) x + 10
.swt_sqrt   <- function(x) sqrt(x)
.swt_boom   <- function(x) stop("kaboom")

skip_targets_unavailable <- function() {
  skip_if_not_installed("S7")
  skip_if_not_installed("withr")
  skip_if_not_installed("targets")
}

local_targets_dir <- function(envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  withr::local_dir(dir, .local_envir = envir)
  dir
}

# -------------------------------------------------------------------------
# sw_pipeline_to_targets()
# -------------------------------------------------------------------------

test_that("sw_pipeline_to_targets writes script and definition RDS files", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times2", .swt_times2),
    sw_step("plus1",  .swt_plus1)
  ))

  info <- sw_pipeline_to_targets(pip, input = 5)

  expect_true(file.exists(info$script))
  expect_true(file.exists(info$input_rds))
  for (path in info$step_rds) expect_true(file.exists(path))

  # Script should be executable on its own via tar_make().
  targets::tar_make(script = info$script, store = "_targets",
                    callr_function = NULL, reporter = "silent")
  expect_equal(targets::tar_read_raw("plus1", store = "_targets"), 11)
})

test_that("sw_pipeline_to_targets refuses to overwrite by default", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(sw_step("id", identity)))
  sw_pipeline_to_targets(pip, input = 1)

  expect_error(
    sw_pipeline_to_targets(pip, input = 1),
    "already exists"
  )

  # overwrite = TRUE succeeds.
  expect_silent(
    sw_pipeline_to_targets(pip, input = 1, overwrite = TRUE)
  )
})

test_that("sw_pipeline_to_targets rejects empty pipelines", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("empty")
  expect_error(
    sw_pipeline_to_targets(pip, input = 1),
    "no steps"
  )
})

# -------------------------------------------------------------------------
# sw_pipeline_run_targets() — round-trip equivalence with sequential backend
# -------------------------------------------------------------------------

test_that("targets backend produces the same result as the sequential backend", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times2", .swt_times2),
    sw_step("plus1",  .swt_plus1),
    sw_step("sqrtit", .swt_sqrt)
  ))

  seq_res <- sw_pipeline_run(pip, input = 8, trace = FALSE)
  tgt_res <- sw_pipeline_run_targets(
    pip, input = 8,
    callr_function = NULL,
    trace = FALSE
  )

  expect_equal(tgt_res$result, seq_res$result)
  expect_equal(names(tgt_res$intermediates), names(seq_res$intermediates))
  expect_equal(tgt_res$intermediates, seq_res$intermediates)
  expect_equal(tgt_res$steps_completed, 3L)
})

# -------------------------------------------------------------------------
# Incremental rerun: only changed step + downstream re-execute
# -------------------------------------------------------------------------

test_that("changing one step's ARGS only invalidates that step and its downstream", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times2", .swt_times2),
    sw_step("plus1",  .swt_plus1),
    sw_step("sqrtit", .swt_sqrt)
  ))

  sw_pipeline_run_targets(pip, input = 8,
                          callr_function = NULL, trace = FALSE)

  # First run: every step should have built.
  prog1 <- targets::tar_progress(store = "_targets")
  built1 <- prog1$name[prog1$progress %in% c("built", "completed")]
  expect_true(all(c("times2", "plus1", "sqrtit") %in% built1))

  # Replace the middle step with a different function.
  pip2 <- sw_pipeline_replace(
    pip, "plus1",
    sw_step("plus1", .swt_plus10)
  )

  res2 <- sw_pipeline_run_targets(pip2, input = 8,
                                  callr_function = NULL, trace = FALSE)

  # `times2` is upstream of the change → should be skipped on the second run.
  prog2 <- targets::tar_progress(store = "_targets")
  status_for <- function(nm) prog2$progress[prog2$name == nm]
  expect_true(status_for("times2") %in% c("skipped", "skip"))
  expect_true(status_for("plus1")  %in% c("built", "completed"))
  expect_true(status_for("sqrtit") %in% c("built", "completed"))

  # And the new result reflects the replacement (sqrt(8 * 2 + 10) == sqrt(26)).
  expect_equal(res2$result, sqrt(26))
})

# -------------------------------------------------------------------------
# backend = "targets" dispatch via sw_pipeline_run()
# -------------------------------------------------------------------------

test_that("sw_pipeline_run(backend = 'targets') delegates correctly", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times2", .swt_times2),
    sw_step("plus1",  .swt_plus1)
  ))

  res <- sw_pipeline_run(
    pip, input = 5,
    backend = "targets",
    callr_function = NULL,
    trace = FALSE
  )

  expect_equal(res$result, 11)
  expect_equal(names(res$intermediates), c("times2", "plus1"))
  expect_equal(res$steps_completed, 2L)
})

# -------------------------------------------------------------------------
# Error propagation
# -------------------------------------------------------------------------

test_that("a failing step surfaces a wrapped error from the targets backend", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times2", .swt_times2),
    sw_step("boom",   .swt_boom)
  ))

  expect_error(
    sw_pipeline_run_targets(pip, input = 1,
                            callr_function = NULL, trace = FALSE),
    "failed under targets backend"
  )
})

# -------------------------------------------------------------------------
# Step-name sanitisation
# -------------------------------------------------------------------------

test_that("step names that collide with reserved internal targets are rejected", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("sw_input", .swt_times2)
  ))

  expect_error(
    sw_pipeline_to_targets(pip, input = 1),
    "collide with internal target names"
  )
})

test_that("step names with spaces are sanitised but kept as intermediates keys", {
  skip_targets_unavailable()
  local_targets_dir()

  pip <- sw_pipeline("toy", steps = list(
    sw_step("times 2", .swt_times2),
    sw_step("plus 1",  .swt_plus1)
  ))

  res <- sw_pipeline_run_targets(pip, input = 5,
                                 callr_function = NULL, trace = FALSE)

  expect_equal(res$result, 11)
  expect_equal(names(res$intermediates), c("times 2", "plus 1"))
})

# -------------------------------------------------------------------------
# Port typing round-trips through the targets backend
# -------------------------------------------------------------------------

test_that("typed ProcessingStep round-trips through targets backend + type check fires", {
  skip_targets_unavailable()
  local_targets_dir()

  # Use a typed pipeline — the new input_type property must survive the
  # saveRDS/readRDS round-trip that the targets backend performs for each
  # step definition.
  pip <- sw_step("double", .swt_times2, input_type = "numeric") %>>%
         sw_step("plus1",  .swt_plus1,  input_type = "numeric")

  # Happy path: numeric input passes type check in child session.
  res <- sw_pipeline_run_targets(pip, input = 8,
                                 callr_function = NULL, trace = FALSE)
  expect_equal(res$result, 17)

  # Wrong-type input: the type check in the child session should fire and
  # the wrapped error should surface through the targets backend.
  expect_error(
    sw_pipeline_run_targets(pip, input = "eight",
                            callr_function = NULL, trace = FALSE),
    "type check failed"
  )
})
