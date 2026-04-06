# tests/testthat/test-composable_pipeline.R
# Unit tests for R/composable_pipeline.R — Composable Pipeline Framework (R7)

# =========================================================================
# ProcessingStep tests
# =========================================================================

test_that("sw_step rejects empty name", {
  skip_if_not_installed("S7")
  expect_error(sw_step("", identity), "non-empty")
})

test_that("sw_step rejects non-character name", {
  skip_if_not_installed("S7")
  expect_error(sw_step(42, identity), "non-empty character string")
})

test_that("sw_step rejects non-function FUN", {
  skip_if_not_installed("S7")
  expect_error(sw_step("test", 42), "must be a function")
})

test_that("sw_step rejects unknown function name", {
  skip_if_not_installed("S7")
  expect_error(sw_step("test", "nonexistent_function_xyz123"),
               "Cannot find function")
})

test_that("sw_step resolves character function names", {
  skip_if_not_installed("S7")
  step <- sw_step("test", "identity")
  expect_equal(step@name, "test")
  expect_true(is.function(step@FUN))
})

test_that("sw_step creates valid ProcessingStep with function", {
  skip_if_not_installed("S7")
  step <- sw_step("double", function(x) x * 2)
  expect_equal(step@name, "double")
  expect_true(is.function(step@FUN))
  expect_equal(step@ARGS, list())
})

test_that("sw_step stores ARGS correctly", {
  skip_if_not_installed("S7")
  step <- sw_step("scale", function(x, factor) x * factor,
                   list(factor = 10))
  expect_equal(step@ARGS, list(factor = 10))
})

test_that("sw_step rejects non-list ARGS", {
  skip_if_not_installed("S7")
  expect_error(sw_step("test", identity, "not_a_list"),
               "must be a list")
})

# =========================================================================
# execute_step tests
# =========================================================================

test_that("execute_step runs function with input", {
  skip_if_not_installed("S7")
  step <- sw_step("double", function(x) x * 2)
  result <- execute_step(step, 5)
  expect_equal(result, 10)
})

test_that("execute_step passes ARGS correctly", {
  skip_if_not_installed("S7")
  step <- sw_step("add", function(x, y) x + y, list(y = 3))
  result <- execute_step(step, 7)
  expect_equal(result, 10)
})

test_that("execute_step rejects non-ProcessingStep", {
  skip_if_not_installed("S7")
  expect_error(execute_step("not_a_step", 5), "ProcessingStep")
})

test_that("execute_step propagates function errors", {
  skip_if_not_installed("S7")
  step <- sw_step("fail", function(x) stop("intentional error"))
  expect_error(execute_step(step, 5), "intentional error")
})

# =========================================================================
# Pipeline creation tests
# =========================================================================

test_that("sw_pipeline creates empty pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test_pipe")
  expect_equal(pip@name, "test_pipe")
  expect_equal(length(pip@steps), 0)
})

test_that("sw_pipeline rejects empty name", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline(""), "non-empty")
})

test_that("sw_pipeline rejects non-character name", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline(42), "non-empty character string")
})

test_that("sw_pipeline creates pipeline with steps", {
  skip_if_not_installed("S7")
  steps <- list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  )
  pip <- sw_pipeline("test_pipe", steps = steps)
  expect_equal(length(pip@steps), 2)
})

test_that("sw_pipeline rejects duplicate step names", {
  skip_if_not_installed("S7")
  steps <- list(
    sw_step("dup", function(x) x + 1),
    sw_step("dup", function(x) x * 2)
  )
  expect_error(sw_pipeline("test_pipe", steps = steps), "Duplicate")
})

test_that("sw_pipeline rejects non-list steps", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline("test_pipe", steps = "not_a_list"),
               "must be a list")
})

# =========================================================================
# Pipeline management tests
# =========================================================================

test_that("sw_pipeline_add appends step to end", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test")
  pip <- sw_pipeline_add(pip, sw_step("s1", function(x) x + 1))
  pip <- sw_pipeline_add(pip, sw_step("s2", function(x) x * 2))
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2"))
})

test_that("sw_pipeline_add inserts at position 0", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  pip <- sw_pipeline_add(pip, sw_step("s0", function(x) x - 1), after = 0)
  expect_equal(sw_pipeline_step_names(pip), c("s0", "s1"))
})

test_that("sw_pipeline_add inserts at specific position", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s3", function(x) x + 3)
  ))
  pip <- sw_pipeline_add(pip, sw_step("s2", function(x) x + 2), after = 1)
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2", "s3"))
})

test_that("sw_pipeline_add rejects duplicate step names", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  expect_error(
    sw_pipeline_add(pip, sw_step("s1", function(x) x * 2)),
    "already exists"
  )
})

test_that("sw_pipeline_add rejects invalid after position", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  expect_error(sw_pipeline_add(pip, sw_step("s2", identity), after = 5),
               "must be between")
})

test_that("sw_pipeline_remove by name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_remove(pip, "s1")
  expect_equal(sw_pipeline_step_names(pip), "s2")
})

test_that("sw_pipeline_remove by index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_remove(pip, 2)
  expect_equal(sw_pipeline_step_names(pip), "s1")
})

test_that("sw_pipeline_remove rejects unknown name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  expect_error(sw_pipeline_remove(pip, "nonexistent"), "No step named")
})

test_that("sw_pipeline_remove rejects out-of-range index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  expect_error(sw_pipeline_remove(pip, 5), "must be between")
})

test_that("sw_pipeline_remove rejects empty pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test")
  expect_error(sw_pipeline_remove(pip, "s1"), "no steps to remove")
})

test_that("sw_pipeline_replace by name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_replace(pip, "s1",
                              sw_step("s1_new", function(x) x + 10))
  expect_equal(sw_pipeline_step_names(pip), c("s1_new", "s2"))
})

test_that("sw_pipeline_replace by index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_replace(pip, 2, sw_step("s2_new", function(x) x * 3))
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2_new"))
})

test_that("sw_pipeline_replace rejects duplicate new name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("s1", function(x) x + 1),
    sw_step("s2", function(x) x * 2)
  ))
  expect_error(
    sw_pipeline_replace(pip, 2, sw_step("s1", function(x) x * 3)),
    "already exists"
  )
})

test_that("sw_pipeline_length returns correct count", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test")
  expect_equal(sw_pipeline_length(pip), 0)

  pip <- sw_pipeline_add(pip, sw_step("s1", identity))
  expect_equal(sw_pipeline_length(pip), 1)

  pip <- sw_pipeline_add(pip, sw_step("s2", identity))
  expect_equal(sw_pipeline_length(pip), 2)
})

test_that("sw_pipeline_step_names returns character vector", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_step("alpha", identity),
    sw_step("beta", identity),
    sw_step("gamma", identity)
  ))
  expect_equal(sw_pipeline_step_names(pip), c("alpha", "beta", "gamma"))
})

# =========================================================================
# Pipeline execution tests
# =========================================================================

test_that("sw_pipeline_run executes linear pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("math", steps = list(
    sw_step("double", function(x) x * 2),
    sw_step("add_one", function(x) x + 1),
    sw_step("square", function(x) x^2)
  ))
  result <- sw_pipeline_run(pip, input = 5, trace = FALSE)
  # 5 * 2 = 10, 10 + 1 = 11, 11^2 = 121
  expect_equal(result$result, 121)
  expect_equal(result$steps_completed, 3)
})

test_that("sw_pipeline_run stores intermediates", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("math", steps = list(
    sw_step("double", function(x) x * 2),
    sw_step("add_one", function(x) x + 1)
  ))
  result <- sw_pipeline_run(pip, input = 5, trace = FALSE)
  expect_equal(result$intermediates$double, 10)
  expect_equal(result$intermediates$add_one, 11)
})

test_that("sw_pipeline_run reports step failure with context", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("fail_pipe", steps = list(
    sw_step("ok", function(x) x + 1),
    sw_step("bad", function(x) stop("boom"))
  ))
  expect_error(
    sw_pipeline_run(pip, input = 5, trace = FALSE),
    "failed at step 2.*bad.*boom"
  )
})

test_that("sw_pipeline_run rejects empty pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("empty")
  expect_error(sw_pipeline_run(pip, input = 5), "no steps to execute")
})

test_that("sw_pipeline_run works with data.frame input", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("df_pipe", steps = list(
    sw_step("add_col", function(df) {
      df$y <- df$x * 2
      df
    }),
    sw_step("filter", function(df) df[df$x > 2, , drop = FALSE])
  ))
  input_df <- data.frame(x = 1:5)
  result <- sw_pipeline_run(pip, input = input_df, trace = FALSE)
  expect_equal(nrow(result$result), 3)
  expect_true("y" %in% names(result$result))
})

test_that("sw_pipeline_run with trace = TRUE prints messages", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("verbose", steps = list(
    sw_step("step1", identity)
  ))
  expect_message(
    sw_pipeline_run(pip, input = 1, trace = TRUE),
    "step1"
  )
})

# =========================================================================
# Pipeline composition tests
# =========================================================================

test_that("sw_pipeline_concat merges two pipelines", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("first", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline("second", steps = list(
    sw_step("s2", function(x) x * 2)
  ))
  combined <- sw_pipeline_concat(pip1, pip2)
  expect_equal(sw_pipeline_step_names(combined), c("s1", "s2"))
  expect_equal(sw_pipeline_length(combined), 2)
})

test_that("sw_pipeline_concat uses custom name", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("a", steps = list(sw_step("s1", identity)))
  pip2 <- sw_pipeline("b", steps = list(sw_step("s2", identity)))
  combined <- sw_pipeline_concat(pip1, pip2, name = "custom_name")
  expect_equal(combined@name, "custom_name")
})

test_that("sw_pipeline_concat default name combines originals", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("first")
  pip2 <- sw_pipeline("second")
  combined <- sw_pipeline_concat(pip1, pip2)
  expect_equal(combined@name, "first + second")
})

test_that("sw_pipeline_concat rejects duplicate step names", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("a", steps = list(sw_step("dup", identity)))
  pip2 <- sw_pipeline("b", steps = list(sw_step("dup", identity)))
  expect_error(sw_pipeline_concat(pip1, pip2), "Duplicate step name")
})

test_that("sw_pipeline_concat result is executable", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("p1", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline("p2", steps = list(
    sw_step("s2", function(x) x * 3)
  ))
  combined <- sw_pipeline_concat(pip1, pip2)
  result <- sw_pipeline_run(combined, input = 4, trace = FALSE)
  # (4 + 1) * 3 = 15
  expect_equal(result$result, 15)
})

# =========================================================================
# sw_pipeline_show tests
# =========================================================================

test_that("sw_pipeline_show prints pipeline info", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("demo", steps = list(
    sw_step("s1", identity),
    sw_step("s2", function(x) x * 2, list(extra = TRUE))
  ))
  out <- capture.output(sw_pipeline_show(pip))
  expect_true(any(grepl("Pipeline: demo", out)))
  expect_true(any(grepl("Steps: 2", out)))
  expect_true(any(grepl("s1", out)))
  expect_true(any(grepl("s2", out)))
})

test_that("sw_pipeline_show handles empty pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("empty")
  out <- capture.output(sw_pipeline_show(pip))
  expect_true(any(grepl("Steps: 0", out)))
})

# =========================================================================
# Convenience step constructors
# =========================================================================

test_that("sw_step_read_fcs creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_step_read_fcs()
  expect_equal(step@name, "read_fcs")
  expect_true(is.function(step@FUN))
})

test_that("sw_step_remove_margins creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_step_remove_margins()
  expect_equal(step@name, "remove_margins")
})

test_that("sw_step_signal_qc creates valid step with args", {
  skip_if_not_installed("S7")
  step <- sw_step_signal_qc(IT_limit = 0.6)
  expect_equal(step@name, "signal_qc")
  expect_equal(step@ARGS$IT_limit, 0.6)
})

test_that("sw_step_batch_correct creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_step_batch_correct(markers = c("CD3", "CD4"))
  expect_equal(step@name, "batch_correct")
  expect_equal(step@ARGS$markers, c("CD3", "CD4"))
})

test_that("sw_step_cluster creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_step_cluster(n_metaclusters = 25)
  expect_equal(step@name, "cluster")
  expect_equal(step@ARGS$n_metaclusters, 25)
})

# =========================================================================
# Integration / composability tests
# =========================================================================

test_that("pipeline supports functional composition pattern", {
  skip_if_not_installed("S7")
  # Build up a pipeline incrementally
  pip <- sw_pipeline("composed")
  pip <- sw_pipeline_add(pip, sw_step("normalize", function(x) x / max(x)))
  pip <- sw_pipeline_add(pip, sw_step("shift", function(x) x - mean(x)))
  pip <- sw_pipeline_add(pip, sw_step("abs", function(x) abs(x)))

  input <- c(2, 4, 6, 8, 10)
  result <- sw_pipeline_run(pip, input = input, trace = FALSE)

  # Verify: normalize → shift → abs
  normalized <- input / max(input)
  shifted <- normalized - mean(normalized)
  expected <- abs(shifted)
  expect_equal(result$result, expected)
})

test_that("pipeline modification preserves immutability", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("original", steps = list(
    sw_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline_add(pip1, sw_step("s2", function(x) x * 2))

  # Original pipeline should be unchanged
  expect_equal(sw_pipeline_length(pip1), 1)
  expect_equal(sw_pipeline_step_names(pip1), "s1")
  expect_equal(sw_pipeline_length(pip2), 2)
  expect_equal(sw_pipeline_step_names(pip2), c("s1", "s2"))
})

test_that("pipeline handles list inputs (multi-sample pattern)", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("batch", steps = list(
    sw_step("process_each", function(samples) {
      lapply(samples, function(s) s * 2)
    }),
    sw_step("combine", function(samples) {
      do.call(c, samples)
    })
  ))
  input <- list(a = 1:3, b = 4:6)
  result <- sw_pipeline_run(pip, input = input, trace = FALSE)
  expect_equal(result$result, c(2, 4, 6, 8, 10, 12))
})
