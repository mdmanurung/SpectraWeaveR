#' @title Composable Pipeline Framework (R7/S7)
#'
#' @description
#' A composable, step-based pipeline framework for flow cytometry analysis,
#' inspired by CytoPipeline (UCLouvain-CBIO). Uses the S7 OOP system instead
#' of S4 for modern, clean class definitions with properties and generics.
#'
#' The framework provides two main classes:
#' \itemize{
#'   \item \code{ProcessingStep}: A single processing operation wrapping a
#'     function and its arguments.
#'   \item \code{Pipeline}: An ordered collection of steps that can be
#'     composed, modified, and executed.
#' }
#'
#' @section Key design choices inspired by CytoPipeline:
#' \itemize{
#'   \item Lazy evaluation: steps are defined but not executed until
#'     \code{execute()} is called
#'   \item Composability: steps can wrap any R function, enabling user-defined
#'     custom functions alongside built-in SpectraWeaveR functions
#'   \item Data threading: the output of each step is piped as the first
#'     argument to the next step
#' }
#'
#' @section Functions that could be ported from CytoPipeline:
#' The following CytoPipeline utilities would complement SpectraWeaveR:
#' \itemize{
#'   \item \strong{estimateScaleTransforms()}: Compute logicle/linear scale
#'     transforms from pooled samples. SpectraWeaveR currently only uses
#'     arcsinh(cofactor=6000); CytoPipeline supports estimateLogicle and
#'     linearQuantile per-channel estimation.
#'   \item \strong{applyScaleTransforms()}: Apply a pre-computed transformList
#'     to flowFrames. Could unify SpectraWeaveR's scattered transformation
#'     logic.
#'   \item \strong{compensateFromMatrix()}: Apply spillover compensation from
#'     FCS keywords or external CSV. SpectraWeaveR currently delegates this
#'     to AutoSpectral; a standalone function would help non-AutoSpectral
#'     workflows.
#'   \item \strong{aggregateAndSample()}: Pool events from multiple FCS files
#'     and downsample. Useful for computing shared transforms or quick
#'     exploratory analysis.
#'   \item \strong{removeDoubletsCytoPipeline()}: Doublet removal using
#'     FSC-A/FSC-H ratio with MAD-based gating. Simpler than building a
#'     full openCyto gating template just for singlet gating.
#'   \item \strong{qualityControlFlowAI()}: FlowAI-based anomaly detection as
#'     an alternative to PeacoQC.
#'   \item \strong{collectNbOfRetainedEvents()}: Track event counts at each
#'     pipeline step. SpectraWeaveR's QC summary covers one step; this would
#'     provide a full audit trail.
#'   \item \strong{singletsGate()}: Standalone singlet gate via parallelogram
#'     method. Lightweight alternative to full openCyto gating.
#'   \item \strong{Channel classification (areSignalCols/areFluoCols)}:
#'     SpectraWeaveR's sw_get_fluor_channels uses regex; CytoPipeline's
#'     approach inspects flowFrame metadata for more robust detection.
#'   \item \strong{Pipeline visualization (plotCytoPipelineProcessingQueue)}:
#'     Render a flowchart of processing steps. Would be valuable for
#'     documenting and presenting SpectraWeaveR pipelines.
#' }
#'
#' @name composable_pipeline
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Internal helper: ensure S7 is available
# ---------------------------------------------------------------------------
.check_s7 <- function() {
  if (!requireNamespace("S7", quietly = TRUE)) {
    stop(
      "Package 'S7' is required for the composable pipeline framework.\n",
      "Install it with: install.packages('S7')",
      call. = FALSE
    )
  }
}

# ---------------------------------------------------------------------------
# ProcessingStep — R7 class
# ---------------------------------------------------------------------------

#' Create an R7 ProcessingStep Class
#'
#' Returns the S7 class object for \code{ProcessingStep}. The class is created
#' once and cached. This function ensures S7 is loaded before attempting class
#' creation.
#'
#' @return The S7 class object for ProcessingStep.
#' @keywords internal
.get_ProcessingStep_class <- function() {
  # Return cached class if already created
  ns <- environment(.get_ProcessingStep_class)
  if (!is.null(ns$.ProcessingStep_class)) {
    return(ns$.ProcessingStep_class)
  }

  .check_s7()

  cls <- S7::new_class("ProcessingStep",
    properties = list(
      name = S7::class_character,
      FUN  = S7::class_function,
      ARGS = S7::class_list
    ),
    validator = function(self) {
      errors <- character()
      if (length(self@name) != 1 || nchar(self@name) == 0) {
        errors <- c(errors, "'name' must be a single non-empty string")
      }
      if (length(errors) > 0) errors else NULL
    }
  )

  ns$.ProcessingStep_class <- cls
  cls
}

# Cache slot
.ProcessingStep_class <- NULL

#' Create a Processing Step
#'
#' Constructs a \code{ProcessingStep} object representing a single operation
#' in the pipeline. Each step wraps a function and its arguments.
#'
#' @param name Character; a unique name for this step (e.g., "load_fcs",
#'   "signal_qc").
#' @param FUN A function to execute, or a character string naming a function.
#'   If a character string, it is resolved via \code{match.fun()} at creation
#'   time.
#' @param ARGS A named list of additional arguments to pass to \code{FUN}.
#'   The data from the previous step is always passed as the first argument;
#'   \code{ARGS} provides the remaining arguments. Default: \code{list()}.
#'
#' @return A \code{ProcessingStep} S7 object.
#'
#' @details
#' When the step is executed, \code{FUN} is called as:
#' \code{do.call(FUN, c(list(previous_result), ARGS))}
#'
#' For the first step in a pipeline, \code{previous_result} comes from the
#' initial input provided to \code{execute()}.
#'
#' @examples
#' \dontrun{
#' # Wrap an existing SpectraWeaveR function
#' step_qc <- sw_step("signal_qc", sw_signal_qc, list(IT_limit = 0.55))
#'
#' # Wrap a custom function
#' step_filter <- sw_step("filter_low", function(ff) ff[1:100, ], list())
#'
#' # Use a function name (resolved via match.fun)
#' step_log <- sw_step("log_transform", "log1p", list())
#' }
#'
#' @export
sw_step <- function(name, FUN, ARGS = list()) {
  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("'name' must be a single non-empty character string.", call. = FALSE)
  }

  # Resolve character function names

  if (is.character(FUN)) {
    if (length(FUN) != 1 || nchar(FUN) == 0) {
      stop("'FUN' as character must be a single non-empty function name.",
           call. = FALSE)
    }
    resolved <- tryCatch(match.fun(FUN), error = function(e) NULL)
    if (is.null(resolved)) {
      stop("Cannot find function '", FUN, "'.", call. = FALSE)
    }
    FUN <- resolved
  }

  if (!is.function(FUN)) {
    stop("'FUN' must be a function or a character string naming a function.",
         call. = FALSE)
  }

  if (!is.list(ARGS)) {
    stop("'ARGS' must be a list.", call. = FALSE)
  }

  cls <- .get_ProcessingStep_class()
  cls(name = name, FUN = FUN, ARGS = ARGS)
}

#' Execute a Processing Step
#'
#' Runs the function wrapped by a \code{ProcessingStep}, passing \code{input}
#' as the first argument followed by the step's \code{ARGS}.
#'
#' @param step A \code{ProcessingStep} object.
#' @param input The data to pass as the first argument to the step's function.
#'
#' @return The result of calling \code{step@@FUN(input, ...)}.
#'
#' @export
execute_step <- function(step, input) {
  .check_s7()
  cls <- .get_ProcessingStep_class()
  if (!S7::S7_inherits(step, cls)) {
    stop("'step' must be a ProcessingStep object.", call. = FALSE)
  }
  do.call(step@FUN, c(list(input), step@ARGS))
}


# ---------------------------------------------------------------------------
# Pipeline — R7 class
# ---------------------------------------------------------------------------

#' Create an R7 Pipeline Class
#'
#' Returns the S7 class object for \code{Pipeline}. Created once and cached.
#'
#' @return The S7 class object for Pipeline.
#' @keywords internal
.get_Pipeline_class <- function() {
  ns <- environment(.get_Pipeline_class)
  if (!is.null(ns$.Pipeline_class)) {
    return(ns$.Pipeline_class)
  }

  .check_s7()
  step_cls <- .get_ProcessingStep_class()

  cls <- S7::new_class("Pipeline",
    properties = list(
      name  = S7::class_character,
      steps = S7::class_list
    ),
    validator = function(self) {
      errors <- character()
      if (length(self@name) != 1 || nchar(self@name) == 0) {
        errors <- c(errors, "'name' must be a single non-empty string")
      }

      # Validate all items in steps are ProcessingStep objects
      step_cls <- .get_ProcessingStep_class()
      for (i in seq_along(self@steps)) {
        if (!S7::S7_inherits(self@steps[[i]], step_cls)) {
          errors <- c(errors,
            paste0("Element ", i, " of 'steps' is not a ProcessingStep"))
        }
      }

      # Check for duplicate step names
      step_names <- vapply(self@steps, function(s) s@name, character(1))
      dups <- step_names[duplicated(step_names)]
      if (length(dups) > 0) {
        errors <- c(errors,
          paste0("Duplicate step name(s): ",
                 paste(unique(dups), collapse = ", ")))
      }

      if (length(errors) > 0) errors else NULL
    }
  )

  ns$.Pipeline_class <- cls
  cls
}

# Cache slot
.Pipeline_class <- NULL

#' Create a Pipeline
#'
#' Constructs a \code{Pipeline} object — an ordered collection of
#' \code{ProcessingStep} objects that can be composed, modified, and executed.
#'
#' @param name Character; a descriptive name for the pipeline
#'   (e.g., "spectral_flow_analysis").
#' @param steps A list of \code{ProcessingStep} objects, or \code{NULL} for
#'   an empty pipeline. Steps are executed in order. Default: \code{NULL}.
#'
#' @return A \code{Pipeline} S7 object.
#'
#' @examples
#' \dontrun{
#' # Create an empty pipeline and add steps later
#' pip <- sw_pipeline("my_analysis")
#' pip <- sw_pipeline_add(pip, sw_step("load", sw_read_fcs))
#'
#' # Create a pipeline with steps
#' pip <- sw_pipeline("qc_pipeline", steps = list(
#'   sw_step("qc", sw_signal_qc, list(IT_limit = 0.55)),
#'   sw_step("filter", my_custom_filter)
#' ))
#' }
#'
#' @export
sw_pipeline <- function(name, steps = NULL) {
  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("'name' must be a single non-empty character string.", call. = FALSE)
  }

  if (is.null(steps)) {
    steps <- list()
  }

  if (!is.list(steps)) {
    stop("'steps' must be a list of ProcessingStep objects or NULL.",
         call. = FALSE)
  }

  cls <- .get_Pipeline_class()
  cls(name = name, steps = steps)
}


# ---------------------------------------------------------------------------
# Pipeline management functions
# ---------------------------------------------------------------------------

#' Add a Step to a Pipeline
#'
#' Returns a new pipeline with the step appended at the given position.
#'
#' @param pipeline A \code{Pipeline} object.
#' @param step A \code{ProcessingStep} object to add.
#' @param after Integer; position after which to insert. If \code{NULL}
#'   (default), the step is appended at the end.
#'
#' @return A new \code{Pipeline} object with the step added.
#'
#' @export
sw_pipeline_add <- function(pipeline, step, after = NULL) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()
  step_cls <- .get_ProcessingStep_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }
  if (!S7::S7_inherits(step, step_cls)) {
    stop("'step' must be a ProcessingStep object.", call. = FALSE)
  }

  steps <- pipeline@steps

  # Check for duplicate names
  existing_names <- vapply(steps, function(s) s@name, character(1))
  if (step@name %in% existing_names) {
    stop("A step named '", step@name, "' already exists in the pipeline.",
         call. = FALSE)
  }

  if (is.null(after)) {
    steps <- c(steps, list(step))
  } else {
    if (!is.numeric(after) || after < 0 || after > length(steps)) {
      stop("'after' must be between 0 and ", length(steps), ".",
           call. = FALSE)
    }
    after <- as.integer(after)
    if (after == 0) {
      steps <- c(list(step), steps)
    } else {
      steps <- c(steps[seq_len(after)], list(step),
                 steps[seq_along(steps) > after])
    }
  }

  cls <- .get_Pipeline_class()
  cls(name = pipeline@name, steps = steps)
}

#' Remove a Step from a Pipeline
#'
#' Returns a new pipeline with the specified step removed.
#'
#' @param pipeline A \code{Pipeline} object.
#' @param name_or_index Character name or integer index of the step to remove.
#'
#' @return A new \code{Pipeline} object with the step removed.
#'
#' @export
sw_pipeline_remove <- function(pipeline, name_or_index) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }

  steps <- pipeline@steps
  if (length(steps) == 0) {
    stop("Pipeline has no steps to remove.", call. = FALSE)
  }

  if (is.character(name_or_index)) {
    step_names <- vapply(steps, function(s) s@name, character(1))
    idx <- match(name_or_index, step_names)
    if (is.na(idx)) {
      stop("No step named '", name_or_index, "' found in pipeline.",
           call. = FALSE)
    }
  } else if (is.numeric(name_or_index)) {
    idx <- as.integer(name_or_index)
    if (idx < 1 || idx > length(steps)) {
      stop("'name_or_index' must be between 1 and ", length(steps), ".",
           call. = FALSE)
    }
  } else {
    stop("'name_or_index' must be a character name or integer index.",
         call. = FALSE)
  }

  steps <- steps[-idx]
  cls <- .get_Pipeline_class()
  cls(name = pipeline@name, steps = steps)
}

#' Replace a Step in a Pipeline
#'
#' Returns a new pipeline with the specified step replaced by a new step.
#'
#' @param pipeline A \code{Pipeline} object.
#' @param name_or_index Character name or integer index of the step to replace.
#' @param new_step A \code{ProcessingStep} object to insert in its place.
#'
#' @return A new \code{Pipeline} object with the step replaced.
#'
#' @export
sw_pipeline_replace <- function(pipeline, name_or_index, new_step) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()
  step_cls <- .get_ProcessingStep_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }
  if (!S7::S7_inherits(new_step, step_cls)) {
    stop("'new_step' must be a ProcessingStep object.", call. = FALSE)
  }

  steps <- pipeline@steps
  if (length(steps) == 0) {
    stop("Pipeline has no steps to replace.", call. = FALSE)
  }

  if (is.character(name_or_index)) {
    step_names <- vapply(steps, function(s) s@name, character(1))
    idx <- match(name_or_index, step_names)
    if (is.na(idx)) {
      stop("No step named '", name_or_index, "' found in pipeline.",
           call. = FALSE)
    }
  } else if (is.numeric(name_or_index)) {
    idx <- as.integer(name_or_index)
    if (idx < 1 || idx > length(steps)) {
      stop("'name_or_index' must be between 1 and ", length(steps), ".",
           call. = FALSE)
    }
  } else {
    stop("'name_or_index' must be a character name or integer index.",
         call. = FALSE)
  }

  # Check duplicate names (excluding the step being replaced)
  other_names <- vapply(steps[-idx], function(s) s@name, character(1))
  if (new_step@name %in% other_names) {
    stop("A step named '", new_step@name, "' already exists in the pipeline.",
         call. = FALSE)
  }

  steps[[idx]] <- new_step
  cls <- .get_Pipeline_class()
  cls(name = pipeline@name, steps = steps)
}

#' Get Step Names from a Pipeline
#'
#' @param pipeline A \code{Pipeline} object.
#'
#' @return Character vector of step names.
#'
#' @export
sw_pipeline_step_names <- function(pipeline) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }
  vapply(pipeline@steps, function(s) s@name, character(1))
}

#' Get Number of Steps in a Pipeline
#'
#' @param pipeline A \code{Pipeline} object.
#'
#' @return Integer.
#'
#' @export
sw_pipeline_length <- function(pipeline) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }
  length(pipeline@steps)
}


# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------

#' Execute a Pipeline
#'
#' Runs all steps in order, threading data from one step to the next.
#' The output of each step becomes the first argument of the next step.
#'
#' @param pipeline A \code{Pipeline} object.
#' @param input The initial input data (e.g., file paths, a flowFrame, etc.).
#' @param trace Logical; if \code{TRUE}, prints step-by-step progress
#'   messages. Default: \code{TRUE}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{result}}{The final output from the last step}
#'     \item{\code{intermediates}}{A named list of intermediate results keyed
#'       by step name (only if \code{trace = TRUE})}
#'     \item{\code{steps_completed}}{Integer count of steps executed}
#'   }
#'
#' @details
#' Execution is sequential: each step receives the output of the previous step
#' as its first argument. This follows the data-threading pattern from
#' CytoPipeline's \code{execute()} method.
#'
#' If a step fails, execution stops and the error is reported with the step
#' name for easy debugging.
#'
#' @examples
#' \dontrun{
#' pip <- sw_pipeline("transform", steps = list(
#'   sw_step("double", function(x) x * 2),
#'   sw_step("add_one", function(x) x + 1)
#' ))
#' result <- sw_pipeline_run(pip, input = 5)
#' # result$result == 11 (5 * 2 + 1)
#' }
#'
#' @export
sw_pipeline_run <- function(pipeline, input, trace = TRUE) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }

  steps <- pipeline@steps
  if (length(steps) == 0) {
    stop("Pipeline '", pipeline@name, "' has no steps to execute.",
         call. = FALSE)
  }

  intermediates <- list()
  current <- input
  n_steps <- length(steps)

  for (i in seq_along(steps)) {
    step <- steps[[i]]
    step_name <- step@name

    if (trace) {
      message(sprintf("=== [%d/%d] %s: %s ===",
                      i, n_steps, pipeline@name, step_name))
    }

    current <- tryCatch(
      execute_step(step, current),
      error = function(e) {
        stop(sprintf("Pipeline '%s' failed at step %d ('%s'): %s",
                     pipeline@name, i, step_name, conditionMessage(e)),
             call. = FALSE)
      }
    )

    intermediates[[step_name]] <- current
  }

  if (trace) {
    message(sprintf("=== Pipeline '%s' complete (%d steps) ===",
                    pipeline@name, n_steps))
  }

  list(
    result = current,
    intermediates = intermediates,
    steps_completed = n_steps
  )
}


# ---------------------------------------------------------------------------
# Pipeline composition
# ---------------------------------------------------------------------------

#' Concatenate Two Pipelines
#'
#' Returns a new pipeline containing all steps from both pipelines in order.
#'
#' @param pipeline1 A \code{Pipeline} object.
#' @param pipeline2 A \code{Pipeline} object.
#' @param name Character; name for the new combined pipeline. If \code{NULL},
#'   uses \code{"pipeline1_name + pipeline2_name"}.
#'
#' @return A new \code{Pipeline} object.
#'
#' @export
sw_pipeline_concat <- function(pipeline1, pipeline2, name = NULL) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline1, pip_cls)) {
    stop("'pipeline1' must be a Pipeline object.", call. = FALSE)
  }
  if (!S7::S7_inherits(pipeline2, pip_cls)) {
    stop("'pipeline2' must be a Pipeline object.", call. = FALSE)
  }

  if (is.null(name)) {
    name <- paste0(pipeline1@name, " + ", pipeline2@name)
  }

  all_steps <- c(pipeline1@steps, pipeline2@steps)

  # Check for duplicate step names across the two pipelines
  all_names <- vapply(all_steps, function(s) s@name, character(1))
  dups <- all_names[duplicated(all_names)]
  if (length(dups) > 0) {
    stop("Duplicate step name(s) across pipelines: ",
         paste(unique(dups), collapse = ", "), call. = FALSE)
  }

  cls <- .get_Pipeline_class()
  cls(name = name, steps = all_steps)
}


# ---------------------------------------------------------------------------
# Show / print helpers
# ---------------------------------------------------------------------------

#' Print a Pipeline Summary
#'
#' Displays a human-readable summary of a pipeline and its steps.
#'
#' @param pipeline A \code{Pipeline} object.
#'
#' @return Invisible \code{NULL}. Called for its side effect of printing.
#'
#' @export
sw_pipeline_show <- function(pipeline) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }

  cat(sprintf("Pipeline: %s\n", pipeline@name))
  cat(sprintf("Steps: %d\n", length(pipeline@steps)))

  if (length(pipeline@steps) > 0) {
    cat("---\n")
    for (i in seq_along(pipeline@steps)) {
      step <- pipeline@steps[[i]]
      fun_name <- tryCatch(
        deparse(substitute(step@FUN)),
        error = function(e) "<function>"
      )
      # Try to get a more informative function name
      fun_label <- tryCatch({
        env <- environment(step@FUN)
        if (!is.null(env)) {
          # Check if the function matches a known name in the parent env
          body_str <- paste(deparse(body(step@FUN)), collapse = " ")
          if (nchar(body_str) > 60) {
            body_str <- paste0(substr(body_str, 1, 57), "...")
          }
          body_str
        } else {
          "<function>"
        }
      }, error = function(e) "<function>")

      n_args <- length(step@ARGS)
      args_str <- if (n_args > 0) {
        paste0(" (", n_args, " arg", if (n_args > 1) "s", ")")
      } else {
        ""
      }
      cat(sprintf("  %d. %s%s\n", i, step@name, args_str))
    }
  }

  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Convenience: pre-built step constructors for SpectraWeaveR functions
# ---------------------------------------------------------------------------

#' Create a Step for Reading FCS Files
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_read_fcs}}.
#'
#' @param ... Additional arguments passed to \code{sw_read_fcs}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_read_fcs <- function(...) {
  sw_step("read_fcs", sw_read_fcs, list(...))
}

#' Create a Step for Margin Removal
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_remove_margins}}.
#'
#' @param ... Additional arguments passed to \code{sw_remove_margins}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_remove_margins <- function(...) {
  sw_step("remove_margins", sw_remove_margins, list(...))
}

#' Create a Step for Signal QC
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_signal_qc}}.
#'
#' @param ... Additional arguments passed to \code{sw_signal_qc}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_signal_qc <- function(...) {
  sw_step("signal_qc", sw_signal_qc, list(...))
}

#' Create a Step for Batch Correction
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_batch_correct}}.
#'
#' @param ... Additional arguments passed to \code{sw_batch_correct}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_batch_correct <- function(...) {
  sw_step("batch_correct", sw_batch_correct, list(...))
}

#' Create a Step for Data Normalisation
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_normalize}}.
#'
#' @param ... Additional arguments passed to \code{sw_normalize}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_normalize <- function(...) {
  sw_step("normalize", sw_normalize, list(...))
}

#' Create a Step for SOM Clustering
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_create_som}}.
#'
#' @param ... Additional arguments passed to \code{sw_create_som}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_create_som <- function(...) {
  sw_step("create_som", sw_create_som, list(...))
}

#' Create a Step for ComBat Correction
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_correct_data}}.
#'
#' @param ... Additional arguments passed to \code{sw_correct_data}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_correct_data <- function(...) {
  sw_step("correct_data", sw_correct_data, list(...))
}

#' Create a Step for Clustering
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_cluster}}.
#'
#' @param ... Additional arguments passed to \code{sw_cluster}.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_step_cluster <- function(...) {
  sw_step("cluster", sw_cluster, list(...))
}


# ---------------------------------------------------------------------------
# Pipeline visualization
# ---------------------------------------------------------------------------

#' Plot Pipeline Processing Queue
#'
#' Renders a text-based flowchart of the pipeline steps, showing the
#' sequential flow of data through each processing step. Inspired by
#' CytoPipeline's \code{plotCytoPipelineProcessingQueue()}.
#'
#' @param pipeline A \code{Pipeline} object.
#' @param style Character; output style. \code{"text"} (default) prints a
#'   text-based flowchart to the console. \code{"data"} returns a data.frame
#'   of step information without printing.
#'
#' @return For \code{style = "text"}, invisible \code{NULL} (called for side
#'   effects). For \code{style = "data"}, a \code{data.frame} with columns
#'   \code{step_number}, \code{name}, and \code{n_args}.
#'
#' @export
sw_plot_pipeline <- function(pipeline, style = c("text", "data")) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }

  style <- match.arg(style)

  steps <- pipeline@steps
  n_steps <- length(steps)

  step_info <- data.frame(
    step_number = seq_len(n_steps),
    name = vapply(steps, function(s) s@name, character(1)),
    n_args = vapply(steps, function(s) length(s@ARGS), integer(1)),
    stringsAsFactors = FALSE
  )

  if (style == "data") {
    return(step_info)
  }

  # Text-based flowchart
  header <- sprintf(" Pipeline: %s (%d steps)", pipeline@name, n_steps)
  width <- max(nchar(header) + 4, 40)

  cat(strrep("=", width), "\n")
  cat(header, "\n")
  cat(strrep("=", width), "\n")

  if (n_steps == 0) {
    cat("  (empty pipeline)\n")
    cat(strrep("=", width), "\n")
    return(invisible(NULL))
  }

  for (i in seq_len(n_steps)) {
    step <- steps[[i]]
    args_str <- if (length(step@ARGS) > 0) {
      arg_names <- names(step@ARGS)
      if (!is.null(arg_names) && any(arg_names != "")) {
        paste0("(", paste(arg_names[arg_names != ""], collapse = ", "), ")")
      } else {
        paste0("(", length(step@ARGS), " args)")
      }
    } else {
      ""
    }

    box_content <- sprintf("[%d] %s %s", i, step@name, args_str)
    box_width <- nchar(box_content) + 4

    cat("  ", strrep("-", box_width), "\n", sep = "")
    cat("  | ", box_content, " |\n", sep = "")
    cat("  ", strrep("-", box_width), "\n", sep = "")

    if (i < n_steps) {
      cat("       |\n")
      cat("       v\n")
    }
  }

  cat(strrep("=", width), "\n")
  invisible(NULL)
}
