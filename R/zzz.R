#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # S7 is in Suggests, so the composable pipeline framework is optional.
  # When it IS available, we force creation of the ProcessingStep and
  # Pipeline classes at package load. This is necessary for two reasons:
  #
  #   1. Each lazy-class getter registers an S7 `print` method on the class
  #      the first time it is created (see R/composable_pipeline.R). If we
  #      skip this at load and let the class be created later, the method
  #      would register too late to be picked up by `S7::methods_register()`.
  #
  #   2. `S7::methods_register()` publishes S7 methods to base R's S3
  #      dispatch table so that typing `pip` at the REPL auto-prints the
  #      flowchart (rather than the default S7 generic object printer).
  #
  # Without this, the auto-print feature would silently not work.
  if (requireNamespace("S7", quietly = TRUE)) {
    tryCatch({
      .get_ProcessingStep_class()
      .get_Pipeline_class()
      S7::methods_register()
    }, error = function(e) {
      # Deliberate no-op: if class creation or method registration fails
      # (e.g. in an unusual S7 upgrade scenario), fall back to the lazy
      # on-demand path. Users can still use sw_pipeline_show() explicitly.
      NULL
    })
  }
}
