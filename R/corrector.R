# Corrector contract (C0): one interface for every batch-correction engine.
#
# Design ideas (frozen model, arrival cases, abstention) follow the cyTether
# implementation bundle v1 (MIT), docs/02_interfaces_and_frozen_contract.md.
# The code is new; defects listed in the roadmap review are not carried over.

# Engine registry ---------------------------------------------------------------

# Registered engines, by name. Engine files add themselves when sourced.
.sw_corrector_engines <- new.env(parent = emptyenv())

# Version tag of the frozen model layout; bump when fields change.
.sw_frozen_schema_version <- "sw_frozen_model/1"

# Columns every effects table must have.
.sw_effect_cols <- c(
  "batch",
  "population",
  "marker",
  "estimate",
  "se_donor",
  "n_donors",
  "n_aliquots",
  "n_aliquots_anchor",
  "support",
  "source"
)

#' Register a batch-correction engine
#'
#' Adds an engine to the registry that \code{\link{sw_corrector_fit}} and
#' \code{\link{sw_corrector_integrate}} dispatch on. The paired reference
#' mean shift engine (\code{"paired_shift"}) registers itself this way.
#'
#' An engine is two functions. \code{fit(xm, cells, anchor, min_cells,
#' hierarchy, seed, ...)} gets the marker matrix \code{xm} (cells by markers,
#' already transformed) and the per-cell table \code{cells} (columns
#' \code{row}, \code{sample_id}, \code{batch}, \code{population},
#' \code{is_reference}, \code{reference_id}, \code{aliquot}). It returns a
#' list with \code{effects} (columns \code{batch}, \code{population},
#' \code{marker}, \code{estimate}, \code{se_donor}, \code{n_donors},
#' \code{n_aliquots}, \code{n_aliquots_anchor}, \code{support},
#' \code{source}), \code{targets} (anything the engine needs later to
#' estimate a new batch) and \code{excluded} (summaries left out, with a
#' reason).
#'
#' \code{apply(model, xm, cells, new_batches)} runs at integrate time for
#' batches the model has not seen. It returns a list with \code{effects} for
#' those batches only and \code{excluded}. It must not change effects already
#' in \code{model}.
#'
#' @param name One string: the engine name used as \code{engine = }.
#' @param fit,apply Functions, as described above.
#' @param can_freeze \code{TRUE} when the engine keeps earlier corrected
#'   values unchanged as new batches arrive. A refitting method sets
#'   \code{FALSE}, so the invariance gate fails it visibly. **No default.**
#' @param description One string shown in diagnostics.
#'
#' @return The engine name, invisibly.
#'
#' @examples
#' zero_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, ...) {
#'   list(effects = NULL, targets = NULL, excluded = NULL)
#' }
#' zero_apply <- function(model, xm, cells, new_batches) {
#'   list(effects = NULL, excluded = NULL)
#' }
#' sw_corrector_register_engine("example_zero", zero_fit, zero_apply,
#'   can_freeze = TRUE, description = "Example engine that never corrects."
#' )
#' "example_zero" %in% sw_corrector_engines()
#'
#' @export
sw_corrector_register_engine <- function(
  name,
  fit,
  apply,
  can_freeze,
  description = NA_character_
) {
  .sw_stop_missing(c(
    name = missing(name),
    fit = missing(fit),
    apply = missing(apply),
    can_freeze = missing(can_freeze)
  ))
  if (
    !is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)
  ) {
    stop("'name' must be one non-empty string.", call. = FALSE)
  }
  if (!is.function(fit) || !is.function(apply)) {
    stop("'fit' and 'apply' must be functions.", call. = FALSE)
  }
  if (!is.logical(can_freeze) || length(can_freeze) != 1 || is.na(can_freeze)) {
    stop("'can_freeze' must be TRUE or FALSE.", call. = FALSE)
  }
  assign(
    name,
    list(
      name = name,
      fit = fit,
      apply = apply,
      can_freeze = can_freeze,
      description = as.character(description)
    ),
    envir = .sw_corrector_engines
  )
  invisible(name)
}

#' List registered correction engines
#'
#' @return A sorted character vector of engine names.
#'
#' @examples
#' sw_corrector_engines()
#'
#' @export
sw_corrector_engines <- function() {
  sort(ls(.sw_corrector_engines))
}

# Look up an engine or stop with the list of known names.
.sw_corrector_engine <- function(engine) {
  if (!exists(engine, envir = .sw_corrector_engines, inherits = FALSE)) {
    stop(
      sprintf(
        "Unknown engine '%s'. Registered: %s.",
        engine,
        paste(sw_corrector_engines(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  get(engine, envir = .sw_corrector_engines, inherits = FALSE)
}

# Shared checks ------------------------------------------------------------------

# Stop naming every required argument the caller left out.
.sw_stop_missing <- function(flags) {
  absent <- names(flags)[flags]
  if (length(absent) > 0) {
    stop(
      "Required argument(s) with no default: ",
      paste0("'", absent, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Normalise and check list(name = "asinh", cofactor = <number>).
.sw_check_transform <- function(transform) {
  if (
    !is.list(transform) ||
      !is.character(transform[["name"]]) ||
      length(transform[["name"]]) != 1 ||
      is.na(transform[["name"]])
  ) {
    stop(
      "'transform' must be a list with a 'name' string, ",
      "for example list(name = \"asinh\", cofactor = 6000).",
      call. = FALSE
    )
  }
  out <- list(name = transform[["name"]])
  if (transform[["name"]] == "asinh") {
    cf <- transform[["cofactor"]]
    if (!is.numeric(cf) || length(cf) != 1 || is.na(cf) || cf <= 0) {
      stop(
        "An asinh 'transform' needs a positive 'cofactor'; it has no default.",
        call. = FALSE
      )
    }
    out$cofactor <- as.double(cf)
  } else if (!is.null(transform[["cofactor"]])) {
    out$cofactor <- as.double(transform[["cofactor"]])
  }
  out
}

# Build the marker matrix and the per-cell table from x, manifest and labels.
.sw_corrector_cells <- function(x, manifest, population, markers, sample_ids) {
  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("A matrix 'x' must be numeric.", call. = FALSE)
    }
    if (is.null(sample_ids)) {
      stop("A matrix 'x' needs 'sample_ids', one per row.", call. = FALSE)
    }
    ids <- as.character(sample_ids)
  } else if (is.data.frame(x)) {
    if (!is.null(sample_ids)) {
      ids <- as.character(sample_ids)
    } else if ("sample_id" %in% names(x)) {
      ids <- as.character(x[["sample_id"]])
    } else {
      stop(
        "'x' needs a 'sample_id' column, or pass 'sample_ids'.",
        call. = FALSE
      )
    }
  } else {
    stop("'x' must be a numeric matrix or a data frame.", call. = FALSE)
  }
  if (length(ids) != nrow(x) || anyNA(ids)) {
    stop(
      "Sample ids must be one non-missing value per row of 'x'.",
      call. = FALSE
    )
  }
  absent <- setdiff(markers, colnames(x))
  if (length(absent) > 0) {
    stop(
      "Marker(s) not in 'x': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  xm <- if (is.matrix(x)) {
    x[, markers, drop = FALSE]
  } else {
    as.matrix(as.data.frame(x)[markers])
  }
  if (!is.numeric(xm)) {
    stop("Marker columns of 'x' must be numeric.", call. = FALSE)
  }
  storage.mode(xm) <- "double"
  if (any(!is.finite(xm))) {
    stop("Marker values must be finite; found NA, NaN or Inf.", call. = FALSE)
  }
  if (
    !is.atomic(population) || length(population) != nrow(x) || anyNA(population)
  ) {
    stop(
      "'population' must be one non-missing label per row of 'x'.",
      call. = FALSE
    )
  }

  # Every sample in x must have one manifest row.
  if (
    !is.data.frame(manifest) ||
      !all(c("sample_id", "batch") %in% names(manifest))
  ) {
    stop(
      "'manifest' must come from sw_sample_manifest() (needs sample_id, batch).",
      call. = FALSE
    )
  }
  uncovered <- setdiff(unique(ids), as.character(manifest$sample_id))
  if (length(uncovered) > 0) {
    stop(
      "Sample(s) in 'x' missing from 'manifest': ",
      paste(utils::head(uncovered, 10), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  m <- match(ids, as.character(manifest$sample_id))
  role <- if ("reference_role" %in% names(manifest)) {
    as.character(manifest$reference_role[m])
  } else {
    rep("none", length(m))
  }
  ref <- if ("reference_id" %in% names(manifest)) {
    as.character(manifest$reference_id[m])
  } else {
    rep(NA_character_, length(m))
  }
  # Only measured references (technical replicate or anchor) are evidence,
  # the same rule sw_design_audit() uses for a shared reference.
  is_ref <- role %in% c("technical_replicate", "anchor") & !is.na(ref)

  has_aliquot <- "aliquot" %in% names(manifest)
  aliquot <- if (has_aliquot) as.character(manifest$aliquot[m]) else ref
  if (has_aliquot && any(is_ref & is.na(aliquot))) {
    stop(
      "Reference sample(s) have a missing 'aliquot': ",
      paste(unique(ids[is_ref & is.na(aliquot)]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  aliquot_source <- if (has_aliquot) {
    "manifest aliquot column"
  } else {
    paste0(
      "defaulted to reference_id: the manifest has no aliquot column, so all ",
      "cells of one reference in one batch form one summary"
    )
  }

  cells <- tibble::tibble(
    row = seq_len(nrow(x)),
    sample_id = ids,
    batch = as.character(manifest$batch[m]),
    population = as.character(population),
    is_reference = is_ref,
    reference_id = ifelse(is_ref, ref, NA_character_),
    aliquot = ifelse(is_ref, aliquot, NA_character_)
  )
  list(xm = xm, cells = cells, aliquot_source = aliquot_source)
}

# The audit must hold an eligible row for this engine and cover the same batches.
.sw_check_audit <- function(design_audit, engine, batches) {
  if (!inherits(design_audit, "sw_design_audit")) {
    stop("'design_audit' must come from sw_design_audit().", call. = FALSE)
  }
  elig <- design_audit$eligibility
  row <- elig[!is.na(elig$name) & elig$name == engine, , drop = FALSE]
  if (nrow(row) == 0) {
    stop(sprintf("audit has no %s eligibility row", engine), call. = FALSE)
  }
  if (!isTRUE(row$eligible[1])) {
    stop(
      sprintf(
        "The design audit marks engine '%s' ineligible: %s",
        engine,
        row$reason[1]
      ),
      call. = FALSE
    )
  }
  if (!setequal(as.character(design_audit$batches), batches)) {
    stop(
      "The design audit covers batches ",
      paste(sort(design_audit$batches), collapse = ", "),
      " but 'x' has ",
      paste(batches, collapse = ", "),
      ". Audit the same data you fit.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Package version, or NA when the code is sourced rather than installed.
.sw_pkg_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SpectraWeaveR")),
    error = function(e) NA_character_
  )
}

# Empty tables in the shared layout.
.sw_empty_effects <- function() {
  tibble::tibble(
    batch = character(),
    population = character(),
    marker = character(),
    estimate = double(),
    se_donor = double(),
    n_donors = integer(),
    n_aliquots = integer(),
    n_aliquots_anchor = integer(),
    support = character(),
    source = character()
  )
}

.sw_empty_excluded <- function() {
  tibble::tibble(
    reference_id = character(),
    aliquot = character(),
    batch = character(),
    population = character(),
    n_cells = integer(),
    reason = character()
  )
}

# Digest and schema --------------------------------------------------------------

# MD5 of a canonical text form of every field except the digest itself.
# hexNumeric deparse keeps doubles exact, so the hash moves with any change.
.sw_model_digest <- function(model) {
  body <- unclass(model)
  body$digest <- NULL
  txt <- deparse(
    body,
    control = c(
      "keepNA",
      "keepInteger",
      "niceNames",
      "showAttributes",
      "hexNumeric"
    )
  )
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  writeLines(txt, f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

# Fields every frozen model carries.
.sw_frozen_fields <- c(
  "schema_version",
  "engine",
  "can_freeze",
  "markers",
  "transform",
  "population_levels",
  "anchor",
  "min_cells",
  "hierarchy",
  "known_batches",
  "effects",
  "targets",
  "excluded",
  "aliquot_source",
  "seed",
  "software",
  "created",
  "parent_digest",
  "digest"
)

# Stop unless `model` has the current schema and an intact digest.
.sw_validate_model <- function(model) {
  if (!inherits(model, "sw_frozen_model")) {
    stop("'model' must be an sw_frozen_model.", call. = FALSE)
  }
  absent <- setdiff(.sw_frozen_fields, names(model))
  if (length(absent) > 0) {
    stop(
      "Frozen model is missing field(s): ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!identical(model$schema_version, .sw_frozen_schema_version)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(model$schema_version),
        .sw_frozen_schema_version
      ),
      call. = FALSE
    )
  }
  if (!all(.sw_effect_cols %in% names(model$effects))) {
    stop("Frozen model 'effects' lacks required columns.", call. = FALSE)
  }
  if (!identical(model$digest, .sw_model_digest(model))) {
    stop(
      "Digest mismatch: the model changed after it was frozen. ",
      "Refit, or load the original file.",
      call. = FALSE
    )
  }
  invisible(model)
}

# Fit ----------------------------------------------------------------------------

#' Fit a frozen batch-correction model
#'
#' Estimates technical batch effects with one registered engine and freezes
#' them, with everything needed to apply them later. The first engine is
#' \code{"paired_shift"}, a paired reference mean shift; see
#' \code{\link{sw_corrector_paired_shift}}.
#'
#' \code{x} must already be on the scale named in \code{transform}. This
#' function records \code{transform} in the model and never applies it.
#' \code{\link{sw_corrector_integrate}} rejects data declared on another
#' scale.
#'
#' Checks, all of which stop with a message: every required argument is set;
#' the manifest covers every sample in \code{x}; the markers exist and hold
#' no \code{NA}; the design audit has an eligible row named after
#' \code{engine} and covers the same batches as \code{x}; the anchor batch is
#' present.
#'
#' Reference samples are those with \code{reference_role}
#' \code{"technical_replicate"} or \code{"anchor"} and a
#' \code{reference_id}, the same rule \code{\link{sw_design_audit}} uses.
#' The replicate unit is the reference aliquot, never the cell.
#'
#' @param x A numeric matrix (cells by markers) or a data frame with one row
#'   per cell, already transformed. A data frame needs a \code{sample_id}
#'   column unless \code{sample_ids} is given.
#' @param manifest An \code{sw_manifest} from \code{\link{sw_sample_manifest}}.
#' @param design_audit An \code{sw_design_audit} from
#'   \code{\link{sw_design_audit}} for the same batches.
#' @param population Character vector of frozen per-cell population labels,
#'   one per row of \code{x}, made before correction (gating or a frozen
#'   classifier). **Required.**
#' @param markers Character vector of marker columns to correct. **Required.**
#' @param transform \code{list(name = "asinh", cofactor = <number>)}: the
#'   scale \code{x} is on. The cofactor has no default; 6000 is common for
#'   spectral flow, but it is your choice. **Required.**
#' @param anchor The batch that defines the target scale. Its effect is zero.
#'   **Required, no default**: a default (for example the first sorted batch)
#'   would turn that batch's noise into bias everywhere.
#' @param engine Name of a registered engine (\code{\link{sw_corrector_engines}}).
#'   **Required.**
#' @param min_cells Minimum cells per reference summary (reference, aliquot,
#'   batch, population). Smaller summaries are excluded and listed with a
#'   reason. **Required.**
#' @param hierarchy Technical hierarchy (study, run, instrument) for engines
#'   that model one. \code{"paired_shift"} accepts only \code{NULL}.
#' @param seed Seed recorded in the model and passed to the engine.
#'   \code{"paired_shift"} draws no random numbers. **Required.**
#' @param ... Passed to the engine's fit function.
#' @param sample_ids Optional sample id per row; required for a matrix
#'   \code{x}.
#'
#' @return An \code{sw_frozen_model}: a list with \code{schema_version},
#'   \code{engine}, \code{can_freeze}, \code{markers}, \code{transform},
#'   \code{population_levels}, \code{anchor}, \code{min_cells},
#'   \code{hierarchy}, \code{known_batches}, \code{effects} (tibble: batch,
#'   population, marker, estimate, se_donor, n_donors, n_aliquots,
#'   n_aliquots_anchor, support, source), \code{targets}, \code{excluded},
#'   \code{aliquot_source}, \code{seed}, \code{software}, \code{created},
#'   \code{parent_digest} and \code{digest}.
#'
#' @seealso \code{\link{sw_corrector_integrate}},
#'   \code{\link{sw_corrector_save}}, \code{\link{sw_corrector_compare}}
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   sample_id = paste0("s", 1:6),
#'   batch = rep(c("b1", "b2", "b3"), each = 2),
#'   reference_role = rep(c("technical_replicate", "none"), 3),
#'   reference_id = rep(c("R1", NA), 3)
#' )
#' manifest <- sw_sample_manifest(meta,
#'   sample_id = "sample_id", batch = "batch",
#'   reference_role = "reference_role", reference_id = "reference_id"
#' )
#' shift <- c(b1 = 0, b2 = 0.5, b3 = -0.3)
#' cell_batch <- rep(meta$batch, each = 50)
#' x <- tibble::tibble(
#'   sample_id = rep(meta$sample_id, each = 50),
#'   CD3 = rnorm(300, 2) + shift[cell_batch],
#'   CD4 = rnorm(300, 1) + shift[cell_batch]
#' )
#' # Hand-built audit with a paired_shift row; sw_design_audit() makes the real one.
#' audit <- structure(list(
#'   batches = c("b1", "b2", "b3"),
#'   eligibility = tibble::tibble(
#'     name = "paired_shift", sw_function = "sw_corrector_fit",
#'     kind = "backend", eligible = TRUE, reason_code = "eligible",
#'     reason = "example", override_allowed = TRUE
#'   )
#' ), class = "sw_design_audit")
#' model <- sw_corrector_fit(x, manifest, audit,
#'   population = rep("T", 300), markers = c("CD3", "CD4"),
#'   transform = list(name = "asinh", cofactor = 6000), anchor = "b1",
#'   engine = "paired_shift", min_cells = 10, seed = 1
#' )
#' model$effects
#'
#' @export
sw_corrector_fit <- function(
  x,
  manifest,
  design_audit,
  population,
  markers,
  transform,
  anchor,
  engine,
  min_cells,
  hierarchy = NULL,
  seed,
  ...,
  sample_ids = NULL
) {
  .sw_stop_missing(c(
    x = missing(x),
    manifest = missing(manifest),
    design_audit = missing(design_audit),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    engine = missing(engine),
    min_cells = missing(min_cells),
    seed = missing(seed)
  ))
  transform <- .sw_check_transform(transform)
  if (
    !is.character(markers) ||
      length(markers) == 0 ||
      anyNA(markers) ||
      anyDuplicated(markers) > 0
  ) {
    stop("'markers' must be unique, non-missing strings.", call. = FALSE)
  }
  if (!is.character(anchor) || length(anchor) != 1 || is.na(anchor)) {
    stop("'anchor' must be one batch name.", call. = FALSE)
  }
  if (!is.character(engine) || length(engine) != 1 || is.na(engine)) {
    stop("'engine' must be one engine name.", call. = FALSE)
  }
  if (
    !is.numeric(min_cells) ||
      length(min_cells) != 1 ||
      is.na(min_cells) ||
      min_cells < 1 ||
      min_cells != round(min_cells)
  ) {
    stop("'min_cells' must be one whole number of at least 1.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("'seed' must be one number.", call. = FALSE)
  }

  d <- .sw_corrector_cells(x, manifest, population, markers, sample_ids)
  batches <- sort(unique(d$cells$batch))
  if (!anchor %in% batches) {
    stop(
      sprintf(
        "Anchor batch '%s' is not in 'x'. Batches: %s.",
        anchor,
        paste(batches, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  eng <- .sw_corrector_engine(engine)
  .sw_check_audit(design_audit, engine, batches)

  fitted <- eng$fit(
    xm = d$xm,
    cells = d$cells,
    anchor = anchor,
    min_cells = as.integer(min_cells),
    hierarchy = hierarchy,
    seed = seed,
    ...
  )
  effects <- if (is.null(fitted$effects)) {
    .sw_empty_effects()
  } else {
    fitted$effects
  }
  if (!all(.sw_effect_cols %in% names(effects))) {
    stop(
      sprintf(
        "Engine '%s' returned effects without the required columns.",
        engine
      ),
      call. = FALSE
    )
  }
  excluded <- if (is.null(fitted$excluded)) {
    .sw_empty_excluded()
  } else {
    fitted$excluded
  }

  model <- structure(
    list(
      schema_version = .sw_frozen_schema_version,
      engine = engine,
      can_freeze = eng$can_freeze,
      markers = markers,
      transform = transform,
      population_levels = sort(unique(d$cells$population)),
      anchor = anchor,
      min_cells = as.integer(min_cells),
      hierarchy = hierarchy,
      known_batches = batches,
      effects = effects,
      targets = fitted$targets,
      excluded = excluded,
      aliquot_source = d$aliquot_source,
      seed = seed,
      software = list(
        R = R.version.string,
        SpectraWeaveR = .sw_pkg_version()
      ),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      parent_digest = NA_character_,
      digest = NA_character_
    ),
    class = "sw_frozen_model"
  )
  model$digest <- .sw_model_digest(model)
  model
}

#' @export
print.sw_frozen_model <- function(x, ...) {
  cat(
    "<sw_frozen_model> engine:",
    x$engine,
    "| anchor:",
    x$anchor,
    "| batches:",
    length(x$known_batches),
    "| markers:",
    length(x$markers),
    "| populations:",
    length(x$population_levels),
    "\n"
  )
  cat(
    "Transform:",
    x$transform[["name"]],
    if (!is.null(x$transform[["cofactor"]])) {
      paste0("(cofactor ", x$transform[["cofactor"]], ")")
    },
    "| min_cells:",
    x$min_cells,
    "| digest:",
    x$digest,
    "\n"
  )
  cat(
    nrow(x$effects),
    "effect(s);",
    nrow(x$excluded),
    "reference summary(ies) excluded\n"
  )
  invisible(x)
}

# Integrate ----------------------------------------------------------------------

#' Apply a frozen model to new or old data
#'
#' Corrects \code{x} with the frozen effects and reports one action per
#' (batch, population, marker). Earlier effects never change, so values
#' corrected before stay the same when new batches arrive.
#'
#' Arrival cases:
#' \describe{
#'   \item{Known batch}{Use the frozen effect. Reference samples in
#'     \code{x} are not used again.}
#'   \item{New batch with a bridging reference}{A reference
#'     (\code{reference_id}) in the new batch also has a frozen summary in
#'     the anchor batch. The engine estimates that batch's effect from its
#'     own incoming reference cells against the frozen anchor summaries. The
#'     effect is added to a copy of the model (returned as \code{model});
#'     old effects are untouched.}
#'   \item{New batch without a bridging reference}{Abstain.}
#'   \item{Known batch or population without an effect}{Abstain. A batch
#'     seen at fit time is never re-estimated; refit to get a new model.}
#' }
#' Abstain keeps the original value exactly. A population label not in the
#' model is abstained too.
#'
#' @param model An \code{sw_frozen_model}.
#' @param x Data in the same form as for \code{\link{sw_corrector_fit}}.
#' @param manifest An \code{sw_manifest} covering every sample in \code{x}.
#' @param population Per-cell population labels for \code{x}.
#' @param transform The scale \code{x} is on, as
#'   \code{list(name, cofactor)}. It must equal the model's transform.
#'   **Required.**
#' @param unsupported \code{"abstain"} keeps uncorrectable values and
#'   reports them; \code{"error"} stops when any action is an abstention.
#'   **Required, no default.**
#' @param sample_ids Optional sample id per row; required for a matrix
#'   \code{x}.
#'
#' @return An \code{sw_corrector_result}: a list with \code{x} (corrected
#'   data, same shape, row order and columns as the input; marker values
#'   as double), \code{actions} (tibble, one row per batch, population and
#'   marker: \code{action} \code{"apply"} or \code{"abstain"},
#'   \code{reason_code}, \code{reason}, estimate and support columns,
#'   \code{n_cells}), \code{model} (the input model, or a copy extended
#'   with new-batch effects, with a new digest and \code{parent_digest}),
#'   \code{excluded} (reference summaries left out in this call) and
#'   \code{aliquot_source}.
#'
#' @seealso \code{\link{sw_corrector_fit}}, \code{\link{sw_corrector_diagnose}}
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   sample_id = paste0("s", 1:6),
#'   batch = rep(c("b1", "b2", "b3"), each = 2),
#'   reference_role = rep(c("technical_replicate", "none"), 3),
#'   reference_id = rep(c("R1", NA), 3)
#' )
#' manifest <- sw_sample_manifest(meta,
#'   sample_id = "sample_id", batch = "batch",
#'   reference_role = "reference_role", reference_id = "reference_id"
#' )
#' shift <- c(b1 = 0, b2 = 0.5, b3 = -0.3)
#' cell_batch <- rep(meta$batch, each = 50)
#' x <- tibble::tibble(
#'   sample_id = rep(meta$sample_id, each = 50),
#'   CD3 = rnorm(300, 2) + shift[cell_batch]
#' )
#' audit <- structure(list(
#'   batches = c("b1", "b2"),
#'   eligibility = tibble::tibble(
#'     name = "paired_shift", sw_function = "sw_corrector_fit",
#'     kind = "backend", eligible = TRUE, reason_code = "eligible",
#'     reason = "example", override_allowed = TRUE
#'   )
#' ), class = "sw_design_audit")
#' tf <- list(name = "asinh", cofactor = 6000)
#' old <- cell_batch != "b3"
#' model <- sw_corrector_fit(x[old, ], manifest, audit,
#'   population = rep("T", sum(old)), markers = "CD3", transform = tf,
#'   anchor = "b1", engine = "paired_shift", min_cells = 10, seed = 1
#' )
#' # b3 is new; its reference R1 bridges it to the anchor.
#' res <- sw_corrector_integrate(model, x, manifest,
#'   population = rep("T", 300), transform = tf, unsupported = "abstain"
#' )
#' res$actions
#'
#' @export
sw_corrector_integrate <- function(
  model,
  x,
  manifest,
  population,
  transform,
  unsupported,
  sample_ids = NULL
) {
  .sw_stop_missing(c(
    model = missing(model),
    x = missing(x),
    manifest = missing(manifest),
    population = missing(population),
    transform = missing(transform),
    unsupported = missing(unsupported)
  ))
  if (
    !is.character(unsupported) ||
      length(unsupported) != 1 ||
      !unsupported %in% c("abstain", "error")
  ) {
    stop("'unsupported' must be \"abstain\" or \"error\".", call. = FALSE)
  }
  .sw_validate_model(model)
  transform <- .sw_check_transform(transform)
  if (!identical(transform, model$transform)) {
    stop(
      "Transform mismatch: the model was fitted on ",
      .sw_transform_label(model$transform),
      " but 'x' is declared as ",
      .sw_transform_label(transform),
      ".",
      call. = FALSE
    )
  }
  absent <- setdiff(model$markers, colnames(x))
  if (length(absent) > 0) {
    stop(
      "Marker set mismatch: model marker(s) not in 'x': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  d <- .sw_corrector_cells(x, manifest, population, model$markers, sample_ids)
  eng <- .sw_corrector_engine(model$engine)

  # Only batches the model has never seen may get new effects.
  new_batches <- setdiff(sort(unique(d$cells$batch)), model$known_batches)
  ext <- if (length(new_batches) > 0) {
    eng$apply(model, d$xm, d$cells, new_batches)
  } else {
    list(effects = NULL, excluded = NULL)
  }
  new_eff <- if (is.null(ext$effects)) .sw_empty_effects() else ext$effects
  new_exc <- if (is.null(ext$excluded)) .sw_empty_excluded() else ext$excluded
  new_eff <- new_eff[new_eff$batch %in% new_batches, , drop = FALSE]

  # Extend a copy; the input model and its effects stay as they were.
  updated <- model
  if (nrow(new_eff) > 0) {
    added <- sort(unique(new_eff$batch))
    updated$effects <- dplyr::bind_rows(model$effects, new_eff)
    updated$known_batches <- sort(c(model$known_batches, added))
    updated$excluded <- dplyr::bind_rows(
      model$excluded,
      new_exc[new_exc$batch %in% added, , drop = FALSE]
    )
    updated$parent_digest <- model$digest
    updated$digest <- .sw_model_digest(updated)
  }

  applied <- .sw_apply_effects(
    updated,
    d$xm,
    d$cells,
    new_batches = new_batches,
    estimated_now = unique(new_eff$batch)
  )
  if (unsupported == "error" && any(applied$actions$action == "abstain")) {
    ab <- applied$actions[applied$actions$action == "abstain", ]
    stop(
      sprintf(
        "%d (batch, population, marker) combination(s) cannot be corrected, e.g. %s / %s / %s: %s",
        nrow(ab),
        ab$batch[1],
        ab$population[1],
        ab$marker[1],
        ab$reason[1]
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      x = .sw_corrector_restore(x, applied$xm, model$markers),
      actions = applied$actions,
      model = updated,
      excluded = new_exc,
      aliquot_source = d$aliquot_source
    ),
    class = "sw_corrector_result"
  )
}

# Short text for a transform, used in messages.
.sw_transform_label <- function(transform) {
  if (is.null(transform[["cofactor"]])) {
    transform[["name"]]
  } else {
    sprintf(
      "%s (cofactor %s)",
      transform[["name"]],
      format(transform[["cofactor"]])
    )
  }
}

# Write corrected marker columns back into the input's shape.
.sw_corrector_restore <- function(x, xm, markers) {
  out <- x
  if (is.matrix(x)) {
    out[, markers] <- xm
  } else {
    for (m in markers) {
      out[[m]] <- unname(xm[, m])
    }
  }
  out
}

# Subtract each (batch, population, marker) effect; abstain where none exists.
.sw_apply_effects <- function(model, xm, cells, new_batches, estimated_now) {
  markers <- model$markers
  groups <- unique(cells[, c("batch", "population")])
  groups <- groups[order(groups$batch, groups$population), , drop = FALSE]
  group_key <- paste(groups$batch, groups$population, sep = "\r")
  n_cells <- as.integer(table(factor(
    paste(cells$batch, cells$population, sep = "\r"),
    levels = group_key
  )))

  # One action row per (batch, population, marker) present in x.
  idx <- rep(seq_len(nrow(groups)), each = length(markers))
  actions <- tibble::tibble(
    batch = groups$batch[idx],
    population = groups$population[idx],
    marker = rep(markers, nrow(groups)),
    n_cells = n_cells[idx]
  )
  eff <- model$effects[, .sw_effect_cols]
  actions <- dplyr::left_join(
    actions,
    eff,
    by = c("batch", "population", "marker"),
    relationship = "many-to-one"
  )

  has_eff <- !is.na(actions$estimate)
  pop_unknown <- !actions$population %in% model$population_levels
  is_new <- actions$batch %in% new_batches
  now <- actions$batch %in% estimated_now
  actions$action <- ifelse(has_eff & !pop_unknown, "apply", "abstain")
  actions$reason_code <- dplyr::case_when(
    pop_unknown ~ "population_unknown",
    has_eff & now ~ "estimated_new_batch",
    has_eff ~ "frozen_effect",
    is_new ~ "no_bridging_reference",
    .default = "known_batch_no_effect"
  )
  actions$reason <- dplyr::case_when(
    pop_unknown ~ "Population label is not in the frozen model.",
    has_eff & now ~ sprintf(
      "New batch: effect estimated from %d reference donor(s) in this batch against frozen anchor summaries.",
      actions$n_donors
    ),
    has_eff ~ sprintf(
      "Frozen effect from %d reference donor(s) (%s).",
      actions$n_donors,
      actions$support
    ),
    is_new ~ "New batch has no reference shared with the anchor (after min_cells); original values kept.",
    .default = "Batch is in the model but has no effect for this population; original values kept. Refit to change this."
  )
  actions <- actions[, c(
    "batch",
    "population",
    "marker",
    "action",
    "reason_code",
    "reason",
    "estimate",
    "se_donor",
    "n_donors",
    "n_aliquots",
    "n_aliquots_anchor",
    "support",
    "n_cells"
  )]

  # Subtract per cell; cells under an abstain action keep their exact value.
  out <- xm
  act_key <- paste(
    actions$batch,
    actions$population,
    actions$marker,
    sep = "\r"
  )
  shift <- ifelse(actions$action == "apply", actions$estimate, NA_real_)
  for (m in markers) {
    s <- shift[match(
      paste(cells$batch, cells$population, m, sep = "\r"),
      act_key
    )]
    ok <- !is.na(s)
    out[ok, m] <- xm[ok, m] - s[ok]
  }
  list(xm = out, actions = actions)
}

#' @export
print.sw_corrector_result <- function(x, ...) {
  tab <- table(x$actions$action)
  cat(
    "<sw_corrector_result>",
    nrow(x$actions),
    "(batch, population, marker) action(s):",
    paste(names(tab), tab, sep = " = ", collapse = ", "),
    "\n"
  )
  invisible(x)
}

# Diagnose -----------------------------------------------------------------------

#' Diagnostics for an integrated result
#'
#' One row per (batch, population, marker) present in the integrated data:
#' the action and its reason, the estimate, its donor-level standard error
#' (\code{NA} with fewer than two donors), the number of reference donors and
#' aliquots, support, and where the aliquot key came from.
#'
#' @param result An \code{sw_corrector_result} from
#'   \code{\link{sw_corrector_integrate}}.
#'
#' @return A tibble: the \code{actions} columns plus \code{engine},
#'   \code{anchor}, \code{batch_status} (\code{"anchor"}, \code{"known"},
#'   \code{"new_estimated"}, \code{"new_unbridged"}),
#'   \code{n_excluded_summaries} (reference summaries below
#'   \code{min_cells} for that batch and population) and
#'   \code{aliquot_source}.
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   sample_id = c("s1", "s2"), batch = c("b1", "b2"),
#'   reference_role = "technical_replicate", reference_id = "R1"
#' )
#' manifest <- sw_sample_manifest(meta,
#'   sample_id = "sample_id", batch = "batch",
#'   reference_role = "reference_role", reference_id = "reference_id"
#' )
#' x <- tibble::tibble(
#'   sample_id = rep(c("s1", "s2"), each = 40),
#'   CD3 = rnorm(80) + rep(c(0, 0.4), each = 40)
#' )
#' audit <- structure(list(
#'   batches = c("b1", "b2"),
#'   eligibility = tibble::tibble(
#'     name = "paired_shift", sw_function = "sw_corrector_fit",
#'     kind = "backend", eligible = TRUE, reason_code = "eligible",
#'     reason = "example", override_allowed = TRUE
#'   )
#' ), class = "sw_design_audit")
#' tf <- list(name = "asinh", cofactor = 6000)
#' model <- sw_corrector_fit(x, manifest, audit,
#'   population = rep("T", 80), markers = "CD3", transform = tf,
#'   anchor = "b1", engine = "paired_shift", min_cells = 10, seed = 1
#' )
#' res <- sw_corrector_integrate(model, x, manifest, rep("T", 80),
#'   transform = tf, unsupported = "abstain"
#' )
#' sw_corrector_diagnose(res)
#'
#' @export
sw_corrector_diagnose <- function(result) {
  if (!inherits(result, "sw_corrector_result")) {
    stop("'result' must come from sw_corrector_integrate().", call. = FALSE)
  }
  model <- result$model
  out <- result$actions
  out$engine <- model$engine
  out$anchor <- model$anchor
  out$batch_status <- dplyr::case_when(
    out$batch == model$anchor ~ "anchor",
    out$reason_code == "estimated_new_batch" ~ "new_estimated",
    out$reason_code == "no_bridging_reference" ~ "new_unbridged",
    .default = "known"
  )
  exc <- dplyr::bind_rows(model$excluded, result$excluded)
  exc <- unique(exc)
  exc_key <- paste(exc$batch, exc$population, sep = "\r")
  out$n_excluded_summaries <- vapply(
    paste(out$batch, out$population, sep = "\r"),
    function(k) sum(exc_key == k),
    integer(1),
    USE.NAMES = FALSE
  )
  out$aliquot_source <- result$aliquot_source
  out
}

# Save and load ------------------------------------------------------------------

#' Save and load a frozen model
#'
#' Stores the model as an RDS file and checks its digest on load. RDS keeps
#' every double exactly and needs no new dependency; JSON would need
#' \code{jsonlite}, which is not a package import. The digest is an MD5 of a
#' canonical text form of every field (doubles written in hex), computed
#' with base R (\code{tools::md5sum}). Load fails when the file was changed
#' after saving, or when its schema version is not the one this package
#' reads. Only load files you trust: \code{readRDS} can run code stored in
#' a crafted file.
#'
#' @param model An \code{sw_frozen_model}.
#' @param path File path.
#'
#' @return \code{sw_corrector_save()} returns \code{path} invisibly;
#'   \code{sw_corrector_load()} returns the \code{sw_frozen_model}.
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   sample_id = c("s1", "s2"), batch = c("b1", "b2"),
#'   reference_role = "technical_replicate", reference_id = "R1"
#' )
#' manifest <- sw_sample_manifest(meta,
#'   sample_id = "sample_id", batch = "batch",
#'   reference_role = "reference_role", reference_id = "reference_id"
#' )
#' x <- tibble::tibble(
#'   sample_id = rep(c("s1", "s2"), each = 40),
#'   CD3 = rnorm(80) + rep(c(0, 0.4), each = 40)
#' )
#' audit <- structure(list(
#'   batches = c("b1", "b2"),
#'   eligibility = tibble::tibble(
#'     name = "paired_shift", sw_function = "sw_corrector_fit",
#'     kind = "backend", eligible = TRUE, reason_code = "eligible",
#'     reason = "example", override_allowed = TRUE
#'   )
#' ), class = "sw_design_audit")
#' model <- sw_corrector_fit(x, manifest, audit,
#'   population = rep("T", 80), markers = "CD3",
#'   transform = list(name = "asinh", cofactor = 6000),
#'   anchor = "b1", engine = "paired_shift", min_cells = 10, seed = 1
#' )
#' path <- tempfile(fileext = ".rds")
#' sw_corrector_save(model, path)
#' identical(sw_corrector_load(path)$digest, model$digest)
#'
#' @export
sw_corrector_save <- function(model, path) {
  .sw_stop_missing(c(model = missing(model), path = missing(path)))
  .sw_validate_model(model)
  saveRDS(model, path)
  invisible(path)
}

#' @rdname sw_corrector_save
#' @export
sw_corrector_load <- function(path) {
  .sw_stop_missing(c(path = missing(path)))
  model <- readRDS(path)
  .sw_validate_model(model)
  model
}

# Compare ------------------------------------------------------------------------

#' Metric records for correction candidates, including the uncorrected data
#'
#' Computes per-comparison metric records for each candidate and for the
#' uncorrected data (\code{candidate_id = "uncorrected"}, always included).
#' The result feeds \code{\link{sw_batch_candidates_from_metrics}}, then
#' \code{\link{sw_batch_decide}}.
#'
#' Comparisons are the population-free batch pairs of
#' \code{design_audit$comparisons}. Metric definitions (new code, so every
#' record carries \code{calibration = "uncalibrated"}):
#' \describe{
#'   \item{\code{emd_within_population}}{The 1-Wasserstein distance between
#'     the two batches, per (population, marker), on matched reference
#'     samples. For each reference id measured in both batches, the cells of
#'     that reference in batch A and in batch B are compared by
#'     \eqn{\frac{1}{K}\sum_k |Q_A(p_k) - Q_B(p_k)|}, with
#'     \eqn{p_k = (k - 0.5)/K}, \eqn{K} = \code{n_quantiles} and
#'     \code{stats::quantile(type = 7)}. The value is the mean over
#'     reference ids; the weight is the number of reference ids used.}
#'   \item{\code{mad_compression_fraction}}{Within each sample and population,
#'     per marker: \eqn{\max(0, 1 - MAD_{after}/MAD_{before})}, with
#'     \code{stats::mad}. It is the fraction by which a candidate shrank
#'     within-sample spread; an additive shift gives 0. The value is the
#'     median over samples in the two batches; the weight is that count.}
#'   \item{\code{joint_structure_distance}}{Within each sample and
#'     population: the mean absolute difference between the off-diagonal
#'     Pearson correlations of the markers before and after correction. The
#'     value is the mean over samples; \code{marker} is \code{NA}.}
#' }
#' Within-sample metrics are used for spread and structure because pooled
#' spread drops legitimately when batch offsets are removed.
#'
#' Status per record: \code{"ok"}, \code{"no_shared_reference"} (EMD only),
#' \code{"too_few_cells"} (fewer than \code{min_cells} on a side, or no
#' sample with enough cells), \code{"constant_marker"} (a marker without
#' spread, so its MAD or correlation is undefined), or \code{"single_marker"}
#' (structure needs at least two markers).
#'
#' @param candidates Named list of corrected data, each with the same rows
#'   and row order as \code{uncorrected}: a matrix, data frame, or an
#'   \code{sw_corrector_result}. The name \code{"uncorrected"} is reserved.
#' @param design_audit An \code{sw_design_audit} for the same batches.
#' @param uncorrected The data before correction, already transformed.
#' @param manifest An \code{sw_manifest} covering every sample.
#' @param population Per-cell population labels.
#' @param markers Marker columns to compare.
#' @param min_cells Minimum cells per reference per batch (EMD) or per
#'   sample (MAD, structure). **Required.**
#' @param sample_ids Optional sample id per row; required for a matrix.
#' @param n_quantiles Number of quantile points for the 1-Wasserstein
#'   approximation. More points resolve tails better.
#'
#' @return A long tibble with columns \code{candidate_id}, \code{metric},
#'   \code{population}, \code{marker}, \code{comparison_id}, \code{value},
#'   \code{status}, \code{weight} and \code{calibration}.
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   sample_id = paste0("s", 1:4), batch = rep(c("b1", "b2"), each = 2),
#'   reference_role = rep(c("technical_replicate", "none"), 2),
#'   reference_id = rep(c("R1", NA), 2)
#' )
#' manifest <- sw_sample_manifest(meta,
#'   sample_id = "sample_id", batch = "batch",
#'   reference_role = "reference_role", reference_id = "reference_id"
#' )
#' cell_batch <- rep(meta$batch, each = 50)
#' x <- tibble::tibble(
#'   sample_id = rep(meta$sample_id, each = 50),
#'   CD3 = rnorm(200) + (cell_batch == "b2") * 0.5,
#'   CD4 = rnorm(200)
#' )
#' shifted <- x
#' shifted$CD3 <- x$CD3 - (cell_batch == "b2") * 0.5
#' audit <- structure(list(
#'   batches = c("b1", "b2"),
#'   comparisons = tibble::tibble(
#'     comparison_id = "b1__vs__b2", batch_a = "b1", batch_b = "b2",
#'     population = NA_character_, shared_reference_ids = list("R1"),
#'     support = "reference_supported"
#'   )
#' ), class = "sw_design_audit")
#' sw_corrector_compare(list(shift = shifted), audit,
#'   uncorrected = x, manifest = manifest, population = rep("T", 200),
#'   markers = c("CD3", "CD4"), min_cells = 10
#' )
#'
#' @export
sw_corrector_compare <- function(
  candidates,
  design_audit,
  uncorrected,
  manifest,
  population,
  markers,
  min_cells,
  sample_ids = NULL,
  n_quantiles = 100
) {
  .sw_stop_missing(c(
    candidates = missing(candidates),
    design_audit = missing(design_audit),
    uncorrected = missing(uncorrected),
    manifest = missing(manifest),
    population = missing(population),
    markers = missing(markers),
    min_cells = missing(min_cells)
  ))
  if (!inherits(design_audit, "sw_design_audit")) {
    stop("'design_audit' must come from sw_design_audit().", call. = FALSE)
  }
  if (
    !is.list(candidates) ||
      is.data.frame(candidates) ||
      is.null(names(candidates)) ||
      any(!nzchar(names(candidates))) ||
      anyDuplicated(names(candidates)) > 0
  ) {
    stop(
      "'candidates' must be a list with unique, non-empty names.",
      call. = FALSE
    )
  }
  if ("uncorrected" %in% names(candidates)) {
    stop(
      "The name 'uncorrected' is reserved for the identity candidate.",
      call. = FALSE
    )
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 2) {
    stop("'min_cells' must be one number of at least 2.", call. = FALSE)
  }

  base <- .sw_corrector_cells(
    uncorrected,
    manifest,
    population,
    markers,
    sample_ids
  )
  cells <- base$cells
  # Each candidate as a marker matrix aligned to the uncorrected rows.
  mats <- c(
    list(uncorrected = base$xm),
    lapply(candidates, function(cand) {
      if (inherits(cand, "sw_corrector_result")) {
        cand <- cand$x
      }
      if (nrow(cand) != nrow(base$xm)) {
        stop(
          "Every candidate must have the same rows as 'uncorrected'.",
          call. = FALSE
        )
      }
      if (is.data.frame(cand) && "sample_id" %in% names(cand)) {
        if (!identical(as.character(cand$sample_id), cells$sample_id)) {
          stop(
            "Candidate rows are not in the order of 'uncorrected'.",
            call. = FALSE
          )
        }
      }
      .sw_corrector_cells(
        cand,
        manifest,
        population,
        markers,
        cells$sample_id
      )$xm
    })
  )

  comps <- design_audit$comparisons
  if (
    !is.data.frame(comps) ||
      !all(
        c("comparison_id", "batch_a", "batch_b", "population") %in% names(comps)
      )
  ) {
    stop(
      "'design_audit$comparisons' must be a table with comparison_id, ",
      "batch_a, batch_b and population, as sw_design_audit() returns.",
      call. = FALSE
    )
  }
  comps <- comps[is.na(comps$population), , drop = FALSE]
  pops <- sort(unique(cells$population))
  probs <- (seq_len(n_quantiles) - 0.5) / n_quantiles

  records <- list()
  for (i in seq_len(nrow(comps))) {
    pair <- c(comps$batch_a[i], comps$batch_b[i])
    cid <- comps$comparison_id[i]
    in_pair <- cells$batch %in% pair
    refs_a <- unique(cells$reference_id[
      cells$is_reference & cells$batch == pair[1]
    ])
    refs_b <- unique(cells$reference_id[
      cells$is_reference & cells$batch == pair[2]
    ])
    shared <- sort(intersect(refs_a, refs_b))
    for (pop in pops) {
      in_pop <- in_pair & cells$population == pop
      samples <- sort(unique(cells$sample_id[in_pop]))
      sample_rows <- lapply(samples, function(s) {
        which(in_pop & cells$sample_id == s)
      })
      sample_rows <- sample_rows[lengths(sample_rows) >= min_cells]
      for (cand in names(mats)) {
        xm <- mats[[cand]]
        records[[length(records) + 1]] <- .sw_compare_emd(
          xm,
          cells,
          in_pop,
          pair,
          shared,
          markers,
          min_cells,
          probs,
          cand,
          pop,
          cid
        )
        records[[length(records) + 1]] <- .sw_compare_mad(
          xm,
          base$xm,
          sample_rows,
          markers,
          cand,
          pop,
          cid
        )
        records[[length(records) + 1]] <- .sw_compare_structure(
          xm,
          base$xm,
          sample_rows,
          markers,
          cand,
          pop,
          cid
        )
      }
    }
  }
  out <- dplyr::bind_rows(records)
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      candidate_id = character(),
      metric = character(),
      population = character(),
      marker = character(),
      comparison_id = character(),
      value = double(),
      status = character(),
      weight = double()
    )
  }
  out$calibration <- rep("uncalibrated", nrow(out))
  out
}

# One record per marker: mean 1-Wasserstein distance over shared references.
.sw_compare_emd <- function(
  xm,
  cells,
  in_pop,
  pair,
  shared,
  markers,
  min_cells,
  probs,
  cand,
  pop,
  cid
) {
  rows <- lapply(markers, function(m) {
    w <- vapply(
      shared,
      function(r) {
        is_r <- in_pop & cells$is_reference & cells$reference_id %in% r
        va <- xm[is_r & cells$batch == pair[1], m]
        vb <- xm[is_r & cells$batch == pair[2], m]
        if (length(va) < min_cells || length(vb) < min_cells) {
          return(NA_real_)
        }
        qa <- stats::quantile(va, probs, type = 7, names = FALSE)
        qb <- stats::quantile(vb, probs, type = 7, names = FALSE)
        mean(abs(qa - qb))
      },
      numeric(1)
    )
    w <- w[!is.na(w)]
    status <- if (length(shared) == 0) {
      "no_shared_reference"
    } else if (length(w) == 0) {
      "too_few_cells"
    } else {
      "ok"
    }
    tibble::tibble(
      candidate_id = cand,
      metric = "emd_within_population",
      population = pop,
      marker = m,
      comparison_id = cid,
      value = if (length(w) > 0) mean(w) else NA_real_,
      status = status,
      weight = as.double(length(w))
    )
  })
  dplyr::bind_rows(rows)
}

# One record per marker: median over samples of max(0, 1 - MAD after / before).
.sw_compare_mad <- function(xm, base, sample_rows, markers, cand, pop, cid) {
  rows <- lapply(markers, function(m) {
    comp <- vapply(
      sample_rows,
      function(ii) {
        before <- stats::mad(base[ii, m])
        if (before <= 0) {
          return(NA_real_)
        }
        max(0, 1 - stats::mad(xm[ii, m]) / before)
      },
      numeric(1)
    )
    comp <- comp[!is.na(comp)]
    status <- if (length(sample_rows) == 0) {
      "too_few_cells"
    } else if (length(comp) == 0) {
      "constant_marker"
    } else {
      "ok"
    }
    tibble::tibble(
      candidate_id = cand,
      metric = "mad_compression_fraction",
      population = pop,
      marker = m,
      comparison_id = cid,
      value = if (length(comp) > 0) stats::median(comp) else NA_real_,
      status = status,
      weight = as.double(length(comp))
    )
  })
  dplyr::bind_rows(rows)
}

# One record: mean over samples of the mean |cor before - cor after|.
.sw_compare_structure <- function(
  xm,
  base,
  sample_rows,
  markers,
  cand,
  pop,
  cid
) {
  d <- vapply(
    sample_rows,
    function(ii) {
      b <- base[ii, markers, drop = FALSE]
      a <- xm[ii, markers, drop = FALSE]
      sds <- c(apply(b, 2, stats::sd), apply(a, 2, stats::sd))
      if (length(markers) < 2 || any(sds <= 0)) {
        return(NA_real_)
      }
      cb <- stats::cor(b)
      ca <- stats::cor(a)
      off <- upper.tri(cb)
      mean(abs(cb[off] - ca[off]))
    },
    numeric(1)
  )
  d <- d[!is.na(d)]
  status <- if (length(markers) < 2) {
    "single_marker"
  } else if (length(sample_rows) == 0) {
    "too_few_cells"
  } else if (length(d) == 0) {
    "constant_marker"
  } else {
    "ok"
  }
  tibble::tibble(
    candidate_id = cand,
    metric = "joint_structure_distance",
    population = pop,
    marker = NA_character_,
    comparison_id = cid,
    value = if (length(d) > 0) mean(d) else NA_real_,
    status = status,
    weight = as.double(length(d))
  )
}
