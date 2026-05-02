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
#'     SpectraWeaveR's sw_channel_get_fluor uses regex; CytoPipeline's
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
      name        = S7::class_character,
      FUN         = S7::class_function,
      ARGS        = S7::class_list,
      input_type  = S7::class_character,
      output_type = S7::class_character
    ),
    validator = function(self) {
      errors <- character()
      if (length(self@name) != 1 || nchar(self@name) == 0) {
        errors <- c(errors, "'name' must be a single non-empty string")
      }
      if (length(errors) > 0) errors else NULL
    }
  )

  # Register an S7 print method on the class so top-level auto-print at the
  # REPL shows the same flowchart that sw_pipeline_show() has always shown.
  # Requires R/zzz.R to call `S7::methods_register()` at package load.
  S7::method(print, cls) <- function(x, ...) {
    cat(sprintf("<ProcessingStep> %s", x@name))
    itypes <- .step_input_type(x)
    otypes <- .step_output_type(x)
    if (length(itypes) > 0 || length(otypes) > 0) {
      cat(sprintf("  [%s -> %s]",
                  if (length(itypes) > 0) paste(itypes, collapse = "|") else "*",
                  if (length(otypes) > 0) paste(otypes, collapse = "|") else "*"))
    }
    if (length(x@ARGS) > 0) {
      cat(sprintf("  (%d arg%s)", length(x@ARGS),
                  if (length(x@ARGS) > 1) "s" else ""))
    }
    cat("\n")
    invisible(x)
  }

  ns$.ProcessingStep_class <- cls
  cls
}

# Cache slot
.ProcessingStep_class <- NULL

# ---------------------------------------------------------------------------
# Safe accessors for optional type properties.
#
# Wrapped in tryCatch so that ProcessingStep objects serialised against an
# OLDER schema (e.g. RDS files already in a user's `_targets/.sw_defs/` store
# from a prior SpectraWeaveR version) degrade gracefully to "opt out" rather
# than throwing "no property called input_type". Legacy pipelines therefore
# keep working, and the next call to `sw_pipeline_run_targets()` rewrites
# the RDS files with the new schema.
# ---------------------------------------------------------------------------
.step_input_type <- function(step) {
  tryCatch(step@input_type,  error = function(e) character(0))
}

.step_output_type <- function(step) {
  tryCatch(step@output_type, error = function(e) character(0))
}

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
#' @param input_type Optional character vector of class names that this
#'   step's input is expected to match. Checked at run time via
#'   \code{methods::is()}, so inherited, virtual, and S4 classes all work.
#'   Pass \code{character(0)} (the default) to opt out of the check.
#' @param output_type Optional character vector of class names declaring
#'   what this step produces. Used for the static pre-run compatibility
#'   walk in \code{\link{sw_pipeline_run}}. \code{character(0)} (default)
#'   means "unknown / don't check".
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
#' step_qc <- sw_pipeline_step("signal_qc", sw_qc_run, list(IT_limit = 0.55))
#'
#' # Wrap a custom function
#' step_filter <- sw_pipeline_step("filter_low", function(ff) ff[1:100, ], list())
#'
#' # Use a function name (resolved via match.fun)
#' step_log <- sw_pipeline_step("log_transform", "log1p", list())
#'
#' # Declare port types for run-time validation
#' step_typed <- sw_pipeline_step("double", function(x) x * 2,
#'                       input_type = "numeric", output_type = "numeric")
#' }
#'
#' @export
sw_pipeline_step <- function(name, FUN, ARGS = list(),
                    input_type  = character(0),
                    output_type = character(0)) {
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

  if (!is.character(input_type)) {
    stop("'input_type' must be a character vector.", call. = FALSE)
  }
  if (!is.character(output_type)) {
    stop("'output_type' must be a character vector.", call. = FALSE)
  }

  cls <- .get_ProcessingStep_class()
  cls(name = name, FUN = FUN, ARGS = ARGS,
      input_type = input_type, output_type = output_type)
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
sw_pipeline_step_run <- function(step, input) {
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

  # Register an S7 print method so typing `pip` at the REPL shows the
  # flowchart automatically (same output as sw_pipeline_show()).
  S7::method(print, cls) <- function(x, ...) {
    sw_pipeline_show(x)
    invisible(x)
  }

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
#' pip <- sw_pipeline_add(pip, sw_pipeline_step("load", sw_io_read_fcs))
#'
#' # Create a pipeline with steps
#' pip <- sw_pipeline("qc_pipeline", steps = list(
#'   sw_pipeline_step("qc", sw_qc_run, list(IT_limit = 0.55)),
#'   sw_pipeline_step("filter", my_custom_filter)
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
#' @param backend Execution backend. Either \code{"sequential"} (the default)
#'   to run steps in-process as a simple for-loop, or \code{"targets"} to
#'   dispatch each step as a \code{targets::tar_target()} (enables
#'   incremental reruns and parallelism). The \code{targets} package must be
#'   installed for the latter.
#' @param store When \code{backend = "targets"}, the path to the targets
#'   data store. Use a stable, persistent path to benefit from incremental
#'   reruns. Default: \code{"_targets"}.
#' @param script When \code{backend = "targets"}, the path to the
#'   \code{_targets.R} script that will be (re)written. Default:
#'   \code{"_targets.R"}.
#' @param ... When \code{backend = "targets"}, additional arguments forwarded
#'   to \code{\link{sw_pipeline_run_targets}} (and from there to
#'   \code{targets::tar_make()}).
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
#' When \code{backend = "targets"}, execution is delegated to
#' \code{\link{sw_pipeline_run_targets}}, which writes a \code{_targets.R}
#' script and invokes \code{targets::tar_make()}. The returned list has the
#' same shape as the sequential backend.
#'
#' @examples
#' \dontrun{
#' pip <- sw_pipeline("transform", steps = list(
#'   sw_pipeline_step("double", function(x) x * 2),
#'   sw_pipeline_step("add_one", function(x) x + 1)
#' ))
#' result <- sw_pipeline_run(pip, input = 5)
#' # result$result == 11 (5 * 2 + 1)
#'
#' # Same pipeline through the targets backend (cached, incremental reruns):
#' result2 <- sw_pipeline_run(pip, input = 5,
#'                            backend = "targets",
#'                            store   = tempfile("targets_store_"))
#' }
#'
#' @seealso \code{\link{sw_pipeline_run_targets}},
#'   \code{\link{sw_pipeline_to_targets}}
#' @export
sw_pipeline_run <- function(pipeline, input, trace = TRUE,
                            backend = c("sequential", "targets"),
                            store = "_targets",
                            script = "_targets.R",
                            ...) {
  .check_s7()
  pip_cls <- .get_Pipeline_class()

  if (!S7::S7_inherits(pipeline, pip_cls)) {
    stop("'pipeline' must be a Pipeline object.", call. = FALSE)
  }

  backend <- match.arg(backend)
  if (backend == "targets") {
    return(sw_pipeline_run_targets(
      pipeline = pipeline,
      input    = input,
      script   = script,
      store    = store,
      trace    = trace,
      ...
    ))
  }

  steps <- pipeline@steps
  if (length(steps) == 0) {
    stop("Pipeline '", pipeline@name, "' has no steps to execute.",
         call. = FALSE)
  }

  # --- Static pre-run type-compatibility walk -----------------------------
  # Consumes `output_type` declarations so they are not dead metadata.
  # Uses string intersection (we only have class names here, no live
  # object), so we emit warnings rather than errors — the actual guarantee
  # comes from the runtime `methods::is()` check in the loop below.
  if (length(steps) >= 2) {
    for (i in seq_len(length(steps) - 1)) {
      out_t <- .step_output_type(steps[[i]])
      in_t  <- .step_input_type(steps[[i + 1]])
      if (length(out_t) > 0 && length(in_t) > 0 &&
          length(intersect(out_t, in_t)) == 0) {
        warning(sprintf(
          paste0("Pipeline '%s': declared types incompatible between ",
                 "step %d ('%s', outputs %s) and step %d ('%s', expects %s). ",
                 "Execution will proceed but may fail."),
          pipeline@name,
          i,     steps[[i]]@name,     paste(out_t, collapse = "|"),
          i + 1, steps[[i + 1]]@name, paste(in_t,  collapse = "|")),
          call. = FALSE)
      }
    }
  }

  intermediates <- list()
  current <- input
  n_steps <- length(steps)

  # Pre-compute whether the pipeline is "typed at all" so that fully
  # untyped legacy pipelines don't flood trace output with "skipped" notes,
  # while mixed pipelines still surface opt-outs to the user.
  any_typed <- any(vapply(steps,
                          function(s) length(.step_input_type(s)) > 0,
                          logical(1)))

  for (i in seq_along(steps)) {
    step <- steps[[i]]
    step_name <- step@name

    if (trace) {
      message(sprintf("=== [%d/%d] %s: %s ===",
                      i, n_steps, pipeline@name, step_name))
    }

    # --- Runtime input-type check ---------------------------------------
    declared_in <- .step_input_type(step)
    if (length(declared_in) > 0) {
      ok <- any(vapply(declared_in,
                       function(t) methods::is(current, t),
                       logical(1)))
      if (!isTRUE(ok)) {
        stop(sprintf(
          paste0("Pipeline '%s' type check failed at step %d ('%s'): ",
                 "expected input of type %s, got %s."),
          pipeline@name, i, step_name,
          paste(declared_in, collapse = "|"),
          paste(class(current), collapse = "/")),
          call. = FALSE)
      }
    } else if (trace && any_typed) {
      message(sprintf("    (step '%s': type check skipped — no input_type declared)",
                      step_name))
    }

    current <- tryCatch(
      sw_pipeline_step_run(step, current),
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
# Pipeline composition DSL: `%>>%`
# ---------------------------------------------------------------------------

# Internal helper: wrap a bare ProcessingStep as a singleton Pipeline.
.sw_step_as_pipeline <- function(step) {
  cls <- .get_Pipeline_class()
  cls(name = step@name, steps = list(step))
}

#' Compose Pipeline Steps with `\%>>\%`
#'
#' Infix operator that composes any two of \code{ProcessingStep} /
#' \code{Pipeline} into a new \code{Pipeline}. Inspired by
#' \href{https://mlr3pipelines.mlr-org.com/}{mlr3pipelines}' \code{\%>>\%},
#' but adapted for SpectraWeaveR's functional, linear pipeline model.
#'
#' The four supported cases all return a fresh \code{Pipeline} and are
#' purely functional — the left- and right-hand operands are never
#' modified:
#'
#' \itemize{
#'   \item \code{step \%>>\% step}  — build a 2-step pipeline
#'   \item \code{step \%>>\% pipe}  — prepend the step to the pipeline
#'   \item \code{pipe \%>>\% step}  — append the step to the pipeline
#'   \item \code{pipe \%>>\% pipe}  — concatenate the two pipelines
#' }
#'
#' Duplicate step names across the operands raise an error, inherited from
#' \code{\link{sw_pipeline_concat}}. \code{NULL} operands are rejected
#' up-front with a clear message. R's custom infix operators are
#' left-associative, so \code{a \%>>\% b \%>>\% c} reads as
#' \code{(a \%>>\% b) \%>>\% c}, matching user intent.
#'
#' @param lhs A \code{ProcessingStep} or \code{Pipeline}.
#' @param rhs A \code{ProcessingStep} or \code{Pipeline}.
#'
#' @return A new \code{Pipeline} object.
#'
#' @examples
#' \dontrun{
#' pip <- sw_pipeline_step("double", function(x) x * 2) \%>>\%
#'        sw_pipeline_step("plus1",  function(x) x + 1)
#' sw_pipeline_run(pip, input = 5)  # 11
#' }
#'
#' @seealso \code{\link{sw_pipeline_concat}},
#'   \code{\link{sw_pipeline_add}}
#' @name grapes-greater-greater-grapes
#' @export
`%>>%` <- function(lhs, rhs) {
  .check_s7()

  if (is.null(lhs) || is.null(rhs)) {
    stop("'%>>%' operands must be ProcessingStep or Pipeline objects, ",
         "not NULL.", call. = FALSE)
  }

  step_cls <- .get_ProcessingStep_class()
  pip_cls  <- .get_Pipeline_class()

  lhs_is_step <- S7::S7_inherits(lhs, step_cls)
  lhs_is_pipe <- S7::S7_inherits(lhs, pip_cls)
  rhs_is_step <- S7::S7_inherits(rhs, step_cls)
  rhs_is_pipe <- S7::S7_inherits(rhs, pip_cls)

  if (!(lhs_is_step || lhs_is_pipe)) {
    stop("'%>>%' left operand must be a ProcessingStep or Pipeline.",
         call. = FALSE)
  }
  if (!(rhs_is_step || rhs_is_pipe)) {
    stop("'%>>%' right operand must be a ProcessingStep or Pipeline.",
         call. = FALSE)
  }

  lhs_pip <- if (lhs_is_step) .sw_step_as_pipeline(lhs) else lhs
  rhs_pip <- if (rhs_is_step) .sw_step_as_pipeline(rhs) else rhs

  combined_name <- paste0(lhs_pip@name, " >> ", rhs_pip@name)

  tryCatch(
    sw_pipeline_concat(lhs_pip, rhs_pip, name = combined_name),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("^Duplicate step name", msg)) {
        stop("'%>>%': ", tolower(substring(msg, 1, 1)), substring(msg, 2),
             call. = FALSE)
      }
      stop("'%>>%': ", msg, call. = FALSE)
    }
  )
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

      n_args <- length(step@ARGS)
      args_str <- if (n_args > 0) {
        paste0(" (", n_args, " arg", if (n_args > 1) "s" else "", ")")
      } else {
        ""
      }

      # Include declared port types when present, formatted as [in -> out].
      itypes <- .step_input_type(step)
      otypes <- .step_output_type(step)
      type_str <- if (length(itypes) > 0 || length(otypes) > 0) {
        sprintf(" [%s -> %s]",
                if (length(itypes) > 0) paste(itypes, collapse = "|") else "*",
                if (length(otypes) > 0) paste(otypes, collapse = "|") else "*")
      } else {
        ""
      }

      cat(sprintf("  %d. %s%s%s\n", i, step@name, type_str, args_str))
    }
  }

  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Convenience: pre-built step constructors for SpectraWeaveR functions
# ---------------------------------------------------------------------------

# Sensible default port types for the built-in convenience constructors.
# These reflect the real first-argument / return types of the underlying
# SpectraWeaveR functions (verified by reading their roxygen blocks).
# A `character(0)` entry means "output is a compound list / heterogeneous
# structure that cannot be cleanly type-matched — opt out of the static
# and runtime checks for that edge".
.SW_STEP_DEFAULT_TYPES <- list(
  read_fcs       = list(input = "character",
                        output = c("flowFrame", "flowSet")),
  remove_margins = list(input = "flowFrame",
                        output = "flowFrame"),
  signal_qc      = list(input = "flowFrame",
                        output = character(0)),
  batch_correct  = list(input = c("tbl_df", "data.frame"),
                        output = c("tbl_df", "data.frame")),
  normalize      = list(input = c("tbl_df", "data.frame"),
                        output = c("tbl_df", "data.frame")),
  create_som     = list(input = c("tbl_df", "data.frame"),
                        output = character(0)),
  correct_data   = list(input = c("tbl_df", "data.frame"),
                        output = c("tbl_df", "data.frame")),
  cluster        = list(input = c("tbl_df", "data.frame", "matrix"),
                        output = character(0)),
  annotate       = list(input = character(0),
                        output = c("tbl_df", "data.frame")),
  differential   = list(input = c("tbl_df", "data.frame"),
                        output = c("tbl_df", "data.frame"))
)

#' Create a Step for Reading FCS Files
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_io_read_fcs}}. Declares sensible port types by default
#' (input: \code{"character"}, output: \code{c("flowFrame","flowSet")});
#' override via the \code{input_type}/\code{output_type} parameters.
#'
#' @param ... Additional arguments passed to \code{sw_io_read_fcs}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_read_fcs <- function(...,
                             input_type  = .SW_STEP_DEFAULT_TYPES$read_fcs$input,
                             output_type = .SW_STEP_DEFAULT_TYPES$read_fcs$output) {
  sw_pipeline_step("read_fcs", sw_io_read_fcs, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Margin Removal
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_filter_margins}}. Defaults: input/output
#' \code{"flowFrame"}.
#'
#' @param ... Additional arguments passed to \code{sw_filter_margins}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_filter_margins <- function(...,
                                   input_type  = .SW_STEP_DEFAULT_TYPES$remove_margins$input,
                                   output_type = .SW_STEP_DEFAULT_TYPES$remove_margins$output) {
  sw_pipeline_step("remove_margins", sw_filter_margins, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Signal QC
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_qc_run}}. Default input \code{"flowFrame"}; output is
#' left unspecified because \code{sw_qc_run} returns a compound list
#' (\code{$FinalFF}, \code{$PlotPath}, \code{$Summary}). Downstream steps
#' that need to thread the flowFrame should declare
#' \code{input_type = "flowFrame"} and extract \code{$FinalFF} themselves.
#'
#' @param ... Additional arguments passed to \code{sw_qc_run}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_qc <- function(...,
                              input_type  = .SW_STEP_DEFAULT_TYPES$signal_qc$input,
                              output_type = .SW_STEP_DEFAULT_TYPES$signal_qc$output) {
  sw_pipeline_step("signal_qc", sw_qc_run, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Batch Correction
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_correct_run}}. Defaults: tibble/data.frame in and out.
#'
#' @param ... Additional arguments passed to \code{sw_correct_run}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_correct <- function(...,
                                  input_type  = .SW_STEP_DEFAULT_TYPES$batch_correct$input,
                                  output_type = .SW_STEP_DEFAULT_TYPES$batch_correct$output) {
  sw_pipeline_step("batch_correct", sw_correct_run, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Data Normalisation
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_correct_normalize}}. Defaults: tibble/data.frame in and out.
#'
#' @param ... Additional arguments passed to \code{sw_correct_normalize}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_normalize <- function(...,
                              input_type  = .SW_STEP_DEFAULT_TYPES$normalize$input,
                              output_type = .SW_STEP_DEFAULT_TYPES$normalize$output) {
  sw_pipeline_step("normalize", sw_correct_normalize, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for SOM Clustering
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_correct_som}}. Default input: tibble/data.frame; output
#' left unspecified (returns an integer vector of cluster labels).
#'
#' @param ... Additional arguments passed to \code{sw_correct_som}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_som <- function(...,
                               input_type  = .SW_STEP_DEFAULT_TYPES$create_som$input,
                               output_type = .SW_STEP_DEFAULT_TYPES$create_som$output) {
  sw_pipeline_step("create_som", sw_correct_som, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for ComBat Correction
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_correct_apply}}. Defaults: tibble/data.frame in and out.
#'
#' @param ... Additional arguments passed to \code{sw_correct_apply}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_correct_apply <- function(...,
                                 input_type  = .SW_STEP_DEFAULT_TYPES$correct_data$input,
                                 output_type = .SW_STEP_DEFAULT_TYPES$correct_data$output) {
  sw_pipeline_step("correct_data", sw_correct_apply, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Clustering
#'
#' Convenience constructor for a \code{ProcessingStep} that wraps
#' \code{\link{sw_cluster_run}}. Default input: tibble/data.frame/matrix;
#' output left unspecified (returns an \code{sw_cluster_result} list).
#'
#' @param ... Additional arguments passed to \code{sw_cluster_run}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_cluster <- function(...,
                            input_type  = .SW_STEP_DEFAULT_TYPES$cluster$input,
                            output_type = .SW_STEP_DEFAULT_TYPES$cluster$output) {
  sw_pipeline_step("cluster", sw_cluster_run, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Cell Type Annotation
#'
#' Convenience constructor for a \code{\link{ProcessingStep}} that wraps
#' \code{\link{sw_annotate_run}}. Accepts either an
#' \code{sw_cluster_result} or a MFI tibble from \code{\link{sw_cluster_mfi}}
#' as input and returns an annotation tibble.
#'
#' @param ... Additional arguments passed to \code{sw_annotate_run}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_annotate <- function(...,
                             input_type  = .SW_STEP_DEFAULT_TYPES$annotate$input,
                             output_type = .SW_STEP_DEFAULT_TYPES$annotate$output) {
  sw_pipeline_step("annotate", sw_annotate_run, list(...),
          input_type = input_type, output_type = output_type)
}

#' Create a Step for Differential Expression
#'
#' Convenience constructor for a \code{\link{ProcessingStep}} that wraps
#' \code{\link{sw_diff_expression}}. Input should be a
#' \code{data.frame}/\code{tibble} with a \code{cluster} column; output is a
#' long-format DE result tibble.
#'
#' @param ... Additional arguments passed to \code{sw_diff_expression}.
#' @param input_type Character vector overriding the default input port type.
#' @param output_type Character vector overriding the default output port type.
#'
#' @return A \code{ProcessingStep} object.
#'
#' @export
sw_pipeline_step_diff <- function(...,
                                  input_type  = .SW_STEP_DEFAULT_TYPES$differential$input,
                                  output_type = .SW_STEP_DEFAULT_TYPES$differential$output) {
  sw_pipeline_step("differential", sw_diff_expression, list(...),
          input_type = input_type, output_type = output_type)
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
