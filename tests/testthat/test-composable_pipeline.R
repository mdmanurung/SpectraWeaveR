# tests/testthat/test-composable_pipeline.R
# Unit tests for R/composable_pipeline.R — Composable Pipeline Framework (R7)

# =========================================================================
# ProcessingStep tests
# =========================================================================

test_that("sw_pipeline_step rejects empty name", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step("", identity), "non-empty")
})

test_that("sw_pipeline_step rejects non-character name", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step(42, identity), "non-empty character string")
})

test_that("sw_pipeline_step rejects non-function FUN", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step("test", 42), "must be a function")
})

test_that("sw_pipeline_step rejects unknown function name", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step("test", "nonexistent_function_xyz123"),
               "Cannot find function")
})

test_that("sw_pipeline_step resolves character function names", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("test", "identity")
  expect_equal(step@name, "test")
  expect_true(is.function(step@FUN))
})

test_that("sw_pipeline_step creates valid ProcessingStep with function", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("double", function(x) x * 2)
  expect_equal(step@name, "double")
  expect_true(is.function(step@FUN))
  expect_equal(step@ARGS, list())
})

test_that("sw_pipeline_step stores ARGS correctly", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("scale", function(x, factor) x * factor,
                   list(factor = 10))
  expect_equal(step@ARGS, list(factor = 10))
})

test_that("sw_pipeline_step rejects non-list ARGS", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step("test", identity, "not_a_list"),
               "must be a list")
})

# =========================================================================
# sw_pipeline_step_run tests
# =========================================================================

test_that("sw_pipeline_step_run runs function with input", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("double", function(x) x * 2)
  result <- sw_pipeline_step_run(step, 5)
  expect_equal(result, 10)
})

test_that("sw_pipeline_step_run passes ARGS correctly", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("add", function(x, y) x + y, list(y = 3))
  result <- sw_pipeline_step_run(step, 7)
  expect_equal(result, 10)
})

test_that("sw_pipeline_step_run rejects non-ProcessingStep", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step_run("not_a_step", 5), "ProcessingStep")
})

test_that("sw_pipeline_step_run propagates function errors", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step("fail", function(x) stop("intentional error"))
  expect_error(sw_pipeline_step_run(step, 5), "intentional error")
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
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  )
  pip <- sw_pipeline("test_pipe", steps = steps)
  expect_equal(length(pip@steps), 2)
})

test_that("sw_pipeline rejects duplicate step names", {
  skip_if_not_installed("S7")
  steps <- list(
    sw_pipeline_step("dup", function(x) x + 1),
    sw_pipeline_step("dup", function(x) x * 2)
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
  pip <- sw_pipeline_add(pip, sw_pipeline_step("s1", function(x) x + 1))
  pip <- sw_pipeline_add(pip, sw_pipeline_step("s2", function(x) x * 2))
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2"))
})

test_that("sw_pipeline_add inserts at position 0", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  pip <- sw_pipeline_add(pip, sw_pipeline_step("s0", function(x) x - 1), after = 0)
  expect_equal(sw_pipeline_step_names(pip), c("s0", "s1"))
})

test_that("sw_pipeline_add inserts at specific position", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s3", function(x) x + 3)
  ))
  pip <- sw_pipeline_add(pip, sw_pipeline_step("s2", function(x) x + 2), after = 1)
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2", "s3"))
})

test_that("sw_pipeline_add rejects duplicate step names", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  expect_error(
    sw_pipeline_add(pip, sw_pipeline_step("s1", function(x) x * 2)),
    "already exists"
  )
})

test_that("sw_pipeline_add rejects invalid after position", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  expect_error(sw_pipeline_add(pip, sw_pipeline_step("s2", identity), after = 5),
               "must be between")
})

test_that("sw_pipeline_remove by name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_remove(pip, "s1")
  expect_equal(sw_pipeline_step_names(pip), "s2")
})

test_that("sw_pipeline_remove by index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_remove(pip, 2)
  expect_equal(sw_pipeline_step_names(pip), "s1")
})

test_that("sw_pipeline_remove rejects unknown name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  expect_error(sw_pipeline_remove(pip, "nonexistent"), "No step named")
})

test_that("sw_pipeline_remove rejects out-of-range index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
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
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_replace(pip, "s1",
                              sw_pipeline_step("s1_new", function(x) x + 10))
  expect_equal(sw_pipeline_step_names(pip), c("s1_new", "s2"))
})

test_that("sw_pipeline_replace by index", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  pip <- sw_pipeline_replace(pip, 2, sw_pipeline_step("s2_new", function(x) x * 3))
  expect_equal(sw_pipeline_step_names(pip), c("s1", "s2_new"))
})

test_that("sw_pipeline_replace rejects duplicate new name", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("s1", function(x) x + 1),
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  expect_error(
    sw_pipeline_replace(pip, 2, sw_pipeline_step("s1", function(x) x * 3)),
    "already exists"
  )
})

test_that("sw_pipeline_length returns correct count", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test")
  expect_equal(sw_pipeline_length(pip), 0)

  pip <- sw_pipeline_add(pip, sw_pipeline_step("s1", identity))
  expect_equal(sw_pipeline_length(pip), 1)

  pip <- sw_pipeline_add(pip, sw_pipeline_step("s2", identity))
  expect_equal(sw_pipeline_length(pip), 2)
})

test_that("sw_pipeline_step_names returns character vector", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("test", steps = list(
    sw_pipeline_step("alpha", identity),
    sw_pipeline_step("beta", identity),
    sw_pipeline_step("gamma", identity)
  ))
  expect_equal(sw_pipeline_step_names(pip), c("alpha", "beta", "gamma"))
})

# =========================================================================
# Pipeline execution tests
# =========================================================================

test_that("sw_pipeline_run executes linear pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("math", steps = list(
    sw_pipeline_step("double", function(x) x * 2),
    sw_pipeline_step("add_one", function(x) x + 1),
    sw_pipeline_step("square", function(x) x^2)
  ))
  result <- sw_pipeline_run(pip, input = 5, trace = FALSE)
  # 5 * 2 = 10, 10 + 1 = 11, 11^2 = 121
  expect_equal(result$result, 121)
  expect_equal(result$steps_completed, 3)
})

test_that("sw_pipeline_run stores intermediates", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("math", steps = list(
    sw_pipeline_step("double", function(x) x * 2),
    sw_pipeline_step("add_one", function(x) x + 1)
  ))
  result <- sw_pipeline_run(pip, input = 5, trace = FALSE)
  expect_equal(result$intermediates$double, 10)
  expect_equal(result$intermediates$add_one, 11)
})

test_that("sw_pipeline_run reports step failure with context", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("fail_pipe", steps = list(
    sw_pipeline_step("ok", function(x) x + 1),
    sw_pipeline_step("bad", function(x) stop("boom"))
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
    sw_pipeline_step("add_col", function(df) {
      df$y <- df$x * 2
      df
    }),
    sw_pipeline_step("filter", function(df) df[df$x > 2, , drop = FALSE])
  ))
  input_df <- data.frame(x = 1:5)
  result <- sw_pipeline_run(pip, input = input_df, trace = FALSE)
  expect_equal(nrow(result$result), 3)
  expect_true("y" %in% names(result$result))
})

test_that("sw_pipeline_run with trace = TRUE prints messages", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("verbose", steps = list(
    sw_pipeline_step("step1", identity)
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
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline("second", steps = list(
    sw_pipeline_step("s2", function(x) x * 2)
  ))
  combined <- sw_pipeline_concat(pip1, pip2)
  expect_equal(sw_pipeline_step_names(combined), c("s1", "s2"))
  expect_equal(sw_pipeline_length(combined), 2)
})

test_that("sw_pipeline_concat uses custom name", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("a", steps = list(sw_pipeline_step("s1", identity)))
  pip2 <- sw_pipeline("b", steps = list(sw_pipeline_step("s2", identity)))
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
  pip1 <- sw_pipeline("a", steps = list(sw_pipeline_step("dup", identity)))
  pip2 <- sw_pipeline("b", steps = list(sw_pipeline_step("dup", identity)))
  expect_error(sw_pipeline_concat(pip1, pip2), "Duplicate step name")
})

test_that("sw_pipeline_concat result is executable", {
  skip_if_not_installed("S7")
  pip1 <- sw_pipeline("p1", steps = list(
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline("p2", steps = list(
    sw_pipeline_step("s2", function(x) x * 3)
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
    sw_pipeline_step("s1", identity),
    sw_pipeline_step("s2", function(x) x * 2, list(extra = TRUE))
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

test_that("sw_pipeline_step_read_fcs creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step_read_fcs()
  expect_equal(step@name, "read_fcs")
  expect_true(is.function(step@FUN))
})

test_that("sw_pipeline_step_filter_margins creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step_filter_margins()
  expect_equal(step@name, "remove_margins")
})

test_that("sw_pipeline_step_qc creates valid step with args", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step_qc(IT_limit = 0.6)
  expect_equal(step@name, "signal_qc")
  expect_equal(step@ARGS$IT_limit, 0.6)
})

test_that("sw_pipeline_step_correct creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step_correct(markers = c("CD3", "CD4"))
  expect_equal(step@name, "batch_correct")
  expect_equal(step@ARGS$markers, c("CD3", "CD4"))
})

test_that("sw_pipeline_step_cluster creates valid step", {
  skip_if_not_installed("S7")
  step <- sw_pipeline_step_cluster(n_metaclusters = 25)
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
  pip <- sw_pipeline_add(pip, sw_pipeline_step("normalize", function(x) x / max(x)))
  pip <- sw_pipeline_add(pip, sw_pipeline_step("shift", function(x) x - mean(x)))
  pip <- sw_pipeline_add(pip, sw_pipeline_step("abs", function(x) abs(x)))

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
    sw_pipeline_step("s1", function(x) x + 1)
  ))
  pip2 <- sw_pipeline_add(pip1, sw_pipeline_step("s2", function(x) x * 2))

  # Original pipeline should be unchanged
  expect_equal(sw_pipeline_length(pip1), 1)
  expect_equal(sw_pipeline_step_names(pip1), "s1")
  expect_equal(sw_pipeline_length(pip2), 2)
  expect_equal(sw_pipeline_step_names(pip2), c("s1", "s2"))
})

test_that("pipeline handles list inputs (multi-sample pattern)", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("batch", steps = list(
    sw_pipeline_step("process_each", function(samples) {
      lapply(samples, function(s) s * 2)
    }),
    sw_pipeline_step("combine", function(samples) {
      do.call(c, samples)
    })
  ))
  input <- list(a = 1:3, b = 4:6)
  result <- sw_pipeline_run(pip, input = input, trace = FALSE)
  expect_equal(result$result, c(2, 4, 6, 8, 10, 12))
})

# =========================================================================
# sw_plot_pipeline tests
# =========================================================================

test_that("sw_plot_pipeline text output works", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("demo", steps = list(
    sw_pipeline_step("load", identity),
    sw_pipeline_step("filter", function(x) x, list(threshold = 0.5)),
    sw_pipeline_step("cluster", function(x) x)
  ))
  out <- capture.output(sw_plot_pipeline(pip))
  expect_true(any(grepl("Pipeline: demo", out)))
  expect_true(any(grepl("load", out)))
  expect_true(any(grepl("filter", out)))
  expect_true(any(grepl("cluster", out)))
})

test_that("sw_plot_pipeline data style returns data.frame", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("demo", steps = list(
    sw_pipeline_step("s1", identity),
    sw_pipeline_step("s2", function(x) x, list(a = 1, b = 2))
  ))
  result <- sw_plot_pipeline(pip, style = "data")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_equal(result$name, c("s1", "s2"))
  expect_equal(result$n_args, c(0L, 2L))
})

test_that("sw_plot_pipeline handles empty pipeline", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("empty")
  out <- capture.output(sw_plot_pipeline(pip))
  expect_true(any(grepl("empty pipeline", out)))
})

test_that("sw_plot_pipeline rejects non-Pipeline", {
  skip_if_not_installed("S7")
  expect_error(sw_plot_pipeline("not_a_pipeline"), "Pipeline")
})


# =========================================================================
# %>>% infix operator (DSL)
# =========================================================================

test_that("%>>% composes two steps into a 2-step pipeline", {
  skip_if_not_installed("S7")
  s1 <- sw_pipeline_step("double", function(x) x * 2)
  s2 <- sw_pipeline_step("plus1",  function(x) x + 1)
  pip <- s1 %>>% s2
  expect_equal(sw_pipeline_step_names(pip), c("double", "plus1"))
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 11)
})

test_that("%>>% prepends a step to a pipeline (step %>>% pipe)", {
  skip_if_not_installed("S7")
  tail_pip <- sw_pipeline("tail", steps = list(
    sw_pipeline_step("plus1", function(x) x + 1)
  ))
  pip <- sw_pipeline_step("double", function(x) x * 2) %>>% tail_pip
  expect_equal(sw_pipeline_step_names(pip), c("double", "plus1"))
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 11)
})

test_that("%>>% appends a step to a pipeline (pipe %>>% step)", {
  skip_if_not_installed("S7")
  head_pip <- sw_pipeline("head", steps = list(
    sw_pipeline_step("double", function(x) x * 2)
  ))
  pip <- head_pip %>>% sw_pipeline_step("plus1", function(x) x + 1)
  expect_equal(sw_pipeline_step_names(pip), c("double", "plus1"))
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 11)
})

test_that("%>>% concatenates two pipelines (pipe %>>% pipe)", {
  skip_if_not_installed("S7")
  p1 <- sw_pipeline("p1", list(sw_pipeline_step("a", function(x) x * 2)))
  p2 <- sw_pipeline("p2", list(sw_pipeline_step("b", function(x) x + 1)))
  pip <- p1 %>>% p2
  expect_equal(sw_pipeline_step_names(pip), c("a", "b"))
  expect_equal(pip@name, "p1 >> p2")
})

test_that("%>>% is left-associative: a %>>% b %>>% c chains three steps", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline_step("a", function(x) x * 2) %>>%
         sw_pipeline_step("b", function(x) x + 1) %>>%
         sw_pipeline_step("c", function(x) x - 3)
  expect_equal(sw_pipeline_length(pip), 3L)
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 8)
})

test_that("%>>% rejects NULL operands with a clear message", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("a", identity)
  expect_error(NULL %>>% s, "not NULL")
  expect_error(s %>>% NULL, "not NULL")
})

test_that("%>>% rejects non-step/non-pipeline operands", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("a", identity)
  expect_error("foo" %>>% s, "ProcessingStep or Pipeline")
  expect_error(s %>>% 42,    "ProcessingStep or Pipeline")
})

test_that("%>>% rewraps duplicate-name error with DSL phrasing", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("dup", function(x) x)
  expect_error(s %>>% s, "^'%>>%':")
})

test_that("%>>% is purely functional: neither operand is mutated", {
  skip_if_not_installed("S7")
  s1 <- sw_pipeline_step("a", function(x) x * 2)
  s2 <- sw_pipeline_step("b", function(x) x + 1)
  p1 <- sw_pipeline("p1", list(s1))
  p2 <- sw_pipeline("p2", list(s2))
  out <- p1 %>>% p2
  expect_equal(sw_pipeline_step_names(p1), "a")
  expect_equal(sw_pipeline_step_names(p2), "b")
  expect_equal(p1@name, "p1")
  expect_equal(p2@name, "p2")
  expect_equal(sw_pipeline_step_names(out), c("a", "b"))
})

test_that("%>>% with an empty pipeline still works", {
  skip_if_not_installed("S7")
  empty <- sw_pipeline("empty")
  s <- sw_pipeline_step("a", function(x) x * 2)
  pip <- empty %>>% s
  expect_equal(sw_pipeline_step_names(pip), "a")
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 10)
})


# =========================================================================
# Port typing (input_type / output_type)
# =========================================================================

test_that("sw_pipeline_step accepts and stores input_type and output_type", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("double", function(x) x * 2,
               input_type = "numeric", output_type = "numeric")
  expect_equal(s@input_type,  "numeric")
  expect_equal(s@output_type, "numeric")
})

test_that("sw_pipeline_step defaults port types to character(0) (opt out)", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("x", identity)
  expect_identical(s@input_type,  character(0))
  expect_identical(s@output_type, character(0))
})

test_that("sw_pipeline_step rejects non-character port types", {
  skip_if_not_installed("S7")
  expect_error(sw_pipeline_step("x", identity, input_type = 1), "input_type")
  expect_error(sw_pipeline_step("x", identity, output_type = list()), "output_type")
})

test_that("sw_pipeline_run enforces declared input_type at each step", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline_step("double", function(x) x * 2, input_type = "numeric") %>>%
         sw_pipeline_step("plus1",  function(x) x + 1,  input_type = "numeric")
  expect_equal(sw_pipeline_run(pip, input = 5, trace = FALSE)$result, 11)
  expect_error(
    sw_pipeline_run(pip, input = "five", trace = FALSE),
    "type check failed at step 1"
  )
})

test_that("sw_pipeline_run skips type check when input_type is empty", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("any", list(
    sw_pipeline_step("first",  function(x) as.character(x)),  # numeric -> character
    sw_pipeline_step("second", function(x) paste0(x, "!"))    # no input_type declared
  ))
  res <- sw_pipeline_run(pip, input = 42, trace = FALSE)
  expect_equal(res$result, "42!")
})

test_that("sw_pipeline_run type-check accepts inherited classes via methods::is", {
  skip_if_not_installed("S7")
  # tibble inherits from data.frame
  skip_if_not_installed("tibble")
  pip <- sw_pipeline("df", list(
    sw_pipeline_step("noop", function(x) x, input_type = "data.frame")
  ))
  res <- sw_pipeline_run(pip, input = tibble::tibble(x = 1:3),
                         trace = FALSE)
  expect_s3_class(res$result, "tbl_df")
})

test_that("static pre-run walk warns on declared-type mismatch but still runs", {
  skip_if_not_installed("S7")
  # Scenario: step 1 declares it outputs "character", but the actual
  # function is identity — so at runtime the value is still numeric and
  # step 2's runtime `is(x, "numeric")` check passes. The declared-type
  # mismatch is caught only by the static pre-run walk.
  pip <- sw_pipeline("lied", list(
    sw_pipeline_step("lies",   function(x) x,
            input_type = "numeric", output_type = "character"),
    sw_pipeline_step("trusts", function(x) x + 1,
            input_type = "numeric", output_type = "numeric")
  ))
  res <- expect_warning(
    sw_pipeline_run(pip, input = 5, trace = FALSE),
    "declared types incompatible"
  )
  expect_equal(res$result, 6)
})

test_that("built-in convenience constructors declare sensible default port types", {
  skip_if_not_installed("S7")
  s_read   <- sw_pipeline_step_read_fcs()
  s_margin <- sw_pipeline_step_filter_margins()
  s_norm   <- sw_pipeline_step_normalize()
  expect_equal(s_read@input_type,    "character")
  expect_equal(s_read@output_type,   c("flowFrame", "flowSet"))
  expect_equal(s_margin@input_type,  "flowFrame")
  expect_equal(s_margin@output_type, "flowFrame")
  expect_equal(s_norm@input_type,    c("tbl_df", "data.frame"))
  expect_equal(s_norm@output_type,   c("tbl_df", "data.frame"))
})

test_that("convenience constructors allow overriding default port types", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step_read_fcs(input_type = "list", output_type = "list")
  expect_equal(s@input_type,  "list")
  expect_equal(s@output_type, "list")
})

test_that("safe accessors tolerate objects missing the new properties", {
  skip_if_not_installed("S7")
  s <- sw_pipeline_step("x", identity)
  # Normal happy path.
  expect_identical(
    SpectraWeaveR:::.step_input_type(s), character(0)
  )
  # Simulate a legacy object whose `@input_type` throws: a plain list that
  # is not a ProcessingStep triggers the tryCatch in the accessor.
  expect_identical(
    SpectraWeaveR:::.step_input_type(list(name = "legacy")),
    character(0)
  )
  expect_identical(
    SpectraWeaveR:::.step_output_type(list(name = "legacy")),
    character(0)
  )
})


# =========================================================================
# Auto-print via S7 method registration
# =========================================================================

test_that("print(pipeline) produces the same output as sw_pipeline_show()", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("demo", list(
    sw_pipeline_step("a", function(x) x * 2, input_type = "numeric",
            output_type = "numeric"),
    sw_pipeline_step("b", function(x) x + 1)
  ))
  show_out  <- capture.output(sw_pipeline_show(pip))
  print_out <- capture.output(print(pip))
  expect_equal(print_out, show_out)
})

test_that("sw_pipeline_show renders declared port types in its output", {
  skip_if_not_installed("S7")
  pip <- sw_pipeline("typed", list(
    sw_pipeline_step("double", function(x) x * 2,
            input_type = "numeric", output_type = "numeric")
  ))
  out <- capture.output(sw_pipeline_show(pip))
  expect_true(any(grepl("numeric -> numeric", out)))
})
