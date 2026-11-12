# Tests for the corrector contract, using a toy engine so they do not
# depend on any real correction method.

# Toy engine: a fixed shift of 1 for every batch except the anchor.
toy_effects <- function(cells, markers, anchor, source) {
  g <- unique(cells[cells$is_reference, c("batch", "population")])
  idx <- rep(seq_len(nrow(g)), each = length(markers))
  tibble::tibble(
    batch = g$batch[idx],
    population = g$population[idx],
    marker = rep(markers, nrow(g)),
    estimate = ifelse(g$batch[idx] == anchor, 0, 1),
    se_donor = NA_real_,
    n_donors = 1L,
    n_aliquots = 1L,
    n_aliquots_anchor = 1L,
    support = "toy",
    source = source
  )
}
sw_corrector_register_engine(
  "toy_shift",
  fit = function(xm, cells, anchor, min_cells, hierarchy, seed, ...) {
    list(
      effects = toy_effects(cells, colnames(xm), anchor, "fit"),
      targets = list(anchor = anchor),
      excluded = NULL
    )
  },
  apply = function(model, xm, cells, new_batches) {
    sub <- cells[cells$batch %in% new_batches, ]
    list(
      effects = toy_effects(sub, colnames(xm), model$anchor, "integrate"),
      excluded = NULL
    )
  },
  can_freeze = TRUE,
  description = "Test engine."
)

toy_data <- function(batches = c("b1", "b2"), n = 20, seed = 2) {
  set.seed(seed)
  meta <- data.frame(
    sample_id = paste0("s_", rep(batches, each = 2), "_", 1:2),
    batch = rep(batches, each = 2),
    reference_role = rep(c("technical_replicate", "none"), length(batches)),
    reference_id = rep(c("R1", NA), length(batches))
  )
  manifest <- sw_sample_manifest(
    meta,
    sample_id = "sample_id",
    batch = "batch",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  ids <- rep(meta$sample_id, each = n)
  x <- tibble::tibble(
    sample_id = ids,
    CD3 = stats::rnorm(length(ids)),
    CD4 = stats::rnorm(length(ids)),
    note = "kept"
  )
  list(
    x = x,
    manifest = manifest,
    population = rep("T", length(ids)),
    batch = meta$batch[match(ids, meta$sample_id)]
  )
}

toy_audit <- function(batches, name = "toy_shift") {
  structure(
    list(
      batches = sort(batches),
      eligibility = tibble::tibble(
        name = name,
        sw_function = "sw_corrector_fit",
        kind = "backend",
        eligible = TRUE,
        reason_code = "eligible",
        reason = "test",
        override_allowed = TRUE
      ),
      comparisons = tibble::tibble(
        comparison_id = "b1__vs__b2",
        batch_a = "b1",
        batch_b = "b2",
        population = NA_character_,
        support = "reference_supported"
      )
    ),
    class = "sw_design_audit"
  )
}

tf <- list(name = "asinh", cofactor = 6000)

toy_args <- function(d) {
  list(
    x = d$x,
    manifest = d$manifest,
    design_audit = toy_audit(unique(d$batch)),
    population = d$population,
    markers = c("CD3", "CD4"),
    transform = tf,
    anchor = "b1",
    engine = "toy_shift",
    min_cells = 5,
    seed = 1
  )
}

test_that("the registry lists and rejects engines", {
  expect_true("toy_shift" %in% sw_corrector_engines())
  d <- toy_data()
  args <- toy_args(d)
  args$engine <- "nope"
  expect_error(do.call(sw_corrector_fit, args), "Unknown engine 'nope'")
  expect_error(
    sw_corrector_register_engine("bad", fit = identity, apply = identity),
    "can_freeze"
  )
})

test_that("fit returns a frozen model with the full schema", {
  d <- toy_data()
  model <- do.call(sw_corrector_fit, toy_args(d))
  expect_s3_class(model, "sw_frozen_model")
  expect_true(all(.sw_frozen_fields %in% names(model)))
  expect_identical(model$schema_version, "sw_frozen_model/1")
  expect_identical(model$transform, list(name = "asinh", cofactor = 6000))
  expect_identical(model$population_levels, "T")
  expect_identical(model$known_batches, c("b1", "b2"))
  expect_identical(model$digest, .sw_model_digest(model))
  expect_output(print(model), "sw_frozen_model")
})

test_that("missing required arguments are named", {
  d <- toy_data()
  for (a in c(
    "x",
    "manifest",
    "design_audit",
    "population",
    "markers",
    "transform",
    "anchor",
    "engine",
    "min_cells",
    "seed"
  )) {
    args <- toy_args(d)
    args[[a]] <- NULL
    expect_error(do.call(sw_corrector_fit, args), sprintf("'%s'", a))
  }
  model <- do.call(sw_corrector_fit, toy_args(d))
  full <- list(
    model = model,
    x = d$x,
    manifest = d$manifest,
    population = d$population,
    transform = tf,
    unsupported = "abstain"
  )
  for (a in names(full)) {
    args <- full
    args[[a]] <- NULL
    expect_error(do.call(sw_corrector_integrate, args), sprintf("'%s'", a))
  }
})

test_that("fit checks transform, markers, manifest and audit", {
  d <- toy_data()
  args <- toy_args(d)
  args$transform <- list(name = "asinh")
  expect_error(do.call(sw_corrector_fit, args), "cofactor")
  args <- toy_args(d)
  args$markers <- c("CD3", "CD8")
  expect_error(do.call(sw_corrector_fit, args), "CD8")
  args <- toy_args(d)
  args$x$CD3[1] <- NA
  expect_error(do.call(sw_corrector_fit, args), "finite")
  args <- toy_args(d)
  args$manifest <- d$manifest[-1, ]
  expect_error(do.call(sw_corrector_fit, args), "missing from 'manifest'")
  args <- toy_args(d)
  args$design_audit <- toy_audit(c("b1", "b2"), name = "other")
  expect_error(
    do.call(sw_corrector_fit, args),
    "audit has no toy_shift eligibility row"
  )
  args <- toy_args(d)
  args$population <- args$population[-1]
  expect_error(do.call(sw_corrector_fit, args), "population")
})

test_that("integrate keeps shape, order and non-marker columns", {
  d <- toy_data()
  model <- do.call(sw_corrector_fit, toy_args(d))
  res <- sw_corrector_integrate(
    model,
    d$x,
    d$manifest,
    d$population,
    transform = tf,
    unsupported = "abstain"
  )
  expect_s3_class(res, "sw_corrector_result")
  expect_identical(names(res$x), names(d$x))
  expect_identical(res$x$sample_id, d$x$sample_id)
  expect_identical(res$x$note, d$x$note)
  b2 <- d$batch == "b2"
  expect_identical(res$x$CD3[b2], d$x$CD3[b2] - 1)
  expect_identical(res$x$CD3[!b2], d$x$CD3[!b2])
  expect_true(all(res$actions$action == "apply"))
  expect_identical(res$model, model)

  # A matrix in gives a matrix out.
  xm <- as.matrix(d$x[, c("CD3", "CD4")])
  res_m <- sw_corrector_integrate(
    model,
    xm,
    d$manifest,
    d$population,
    transform = tf,
    unsupported = "abstain",
    sample_ids = d$x$sample_id
  )
  expect_true(is.matrix(res_m$x))
  expect_identical(unname(res_m$x[, "CD3"]), res$x$CD3)
})

test_that("integrate rejects a mismatched transform or marker set", {
  d <- toy_data()
  model <- do.call(sw_corrector_fit, toy_args(d))
  expect_error(
    sw_corrector_integrate(
      model,
      d$x,
      d$manifest,
      d$population,
      transform = list(name = "asinh", cofactor = 150),
      unsupported = "abstain"
    ),
    "Transform mismatch"
  )
  expect_error(
    sw_corrector_integrate(
      model,
      d$x[, c("sample_id", "CD3")],
      d$manifest,
      d$population,
      transform = tf,
      unsupported = "abstain"
    ),
    "Marker set mismatch"
  )
  expect_error(
    sw_corrector_integrate(
      model,
      d$x,
      d$manifest,
      d$population,
      transform = tf,
      unsupported = "maybe"
    ),
    "unsupported"
  )
})

test_that("a new batch extends a copy of the model only", {
  d3 <- toy_data(batches = c("b1", "b2", "b3"))
  old <- d3$batch != "b3"
  d <- list(
    x = d3$x[old, ],
    manifest = d3$manifest,
    population = d3$population[old],
    batch = d3$batch[old]
  )
  model <- do.call(sw_corrector_fit, toy_args(d))
  res <- sw_corrector_integrate(
    model,
    d3$x,
    d3$manifest,
    d3$population,
    transform = tf,
    unsupported = "abstain"
  )
  expect_identical(model$known_batches, c("b1", "b2"))
  expect_identical(res$model$known_batches, c("b1", "b2", "b3"))
  expect_identical(res$model$parent_digest, model$digest)
  act <- res$actions[res$actions$batch == "b3", ]
  expect_true(all(act$reason_code == "estimated_new_batch"))
})

test_that("save and load round trip; a tampered file fails", {
  d <- toy_data()
  model <- do.call(sw_corrector_fit, toy_args(d))
  path <- tempfile(fileext = ".rds")
  expect_identical(sw_corrector_save(model, path), path)
  loaded <- sw_corrector_load(path)
  expect_identical(loaded$digest, model$digest)
  expect_identical(.sw_model_digest(loaded), model$digest)

  tampered <- readRDS(path)
  tampered$effects$estimate[2] <- tampered$effects$estimate[2] + 1e-9
  saveRDS(tampered, path)
  expect_error(sw_corrector_load(path), "Digest mismatch")

  old_schema <- model
  old_schema$schema_version <- "sw_frozen_model/0"
  saveRDS(old_schema, path)
  expect_error(sw_corrector_load(path), "schema_version")

  partial <- unclass(model)
  partial$effects <- NULL
  saveRDS(structure(partial, class = "sw_frozen_model"), path)
  expect_error(sw_corrector_load(path), "missing field")

  # Integrate refuses a model changed in memory, too.
  model$anchor <- "b2"
  expect_error(
    sw_corrector_integrate(
      model,
      d$x,
      d$manifest,
      d$population,
      transform = tf,
      unsupported = "abstain"
    ),
    "Digest mismatch"
  )
})

test_that("compare always includes the uncorrected candidate", {
  d <- toy_data(n = 60)
  shifted <- d$x
  shifted$CD3 <- d$x$CD3 + (d$batch == "b2") * 2
  rec <- sw_corrector_compare(
    list(shifted = shifted),
    toy_audit(c("b1", "b2")),
    uncorrected = d$x,
    manifest = d$manifest,
    population = d$population,
    markers = c("CD3", "CD4"),
    min_cells = 10
  )
  expect_identical(
    names(rec),
    c(
      "candidate_id",
      "metric",
      "population",
      "marker",
      "comparison_id",
      "value",
      "status",
      "weight",
      "calibration"
    )
  )
  expect_setequal(unique(rec$candidate_id), c("uncorrected", "shifted"))
  expect_setequal(
    unique(rec$metric),
    c(
      "emd_within_population",
      "mad_compression_fraction",
      "joint_structure_distance"
    )
  )
  emd <- rec[rec$metric == "emd_within_population" & rec$marker == "CD3", ]
  expect_gt(emd$value[emd$candidate_id == "shifted"], 1.5)
  expect_true(all(
    rec$value[
      rec$candidate_id == "uncorrected" &
        rec$metric != "emd_within_population"
    ] ==
      0
  ))
  expect_error(
    sw_corrector_compare(
      list(uncorrected = d$x),
      toy_audit(c("b1", "b2")),
      uncorrected = d$x,
      manifest = d$manifest,
      population = d$population,
      markers = "CD3",
      min_cells = 10
    ),
    "reserved"
  )
  expect_error(
    sw_corrector_compare(
      list(a = d$x),
      toy_audit(c("b1", "b2")),
      uncorrected = d$x,
      manifest = d$manifest,
      population = d$population,
      markers = "CD3"
    ),
    "'min_cells'"
  )
})

test_that("compare marks pairs without a shared reference", {
  d <- toy_data(n = 30)
  man <- d$manifest
  man$reference_id[man$batch == "b2" & !is.na(man$reference_id)] <- "R2"
  rec <- sw_corrector_compare(
    list(same = d$x),
    toy_audit(c("b1", "b2")),
    uncorrected = d$x,
    manifest = man,
    population = d$population,
    markers = c("CD3", "CD4"),
    min_cells = 10
  )
  emd <- rec[rec$metric == "emd_within_population", ]
  expect_true(all(emd$status == "no_shared_reference"))
  expect_true(all(is.na(emd$value)))
})

test_that("synthetic_reference rows are never reference evidence", {
  meta <- data.frame(
    sample_id = c("a1", "a2", "b1", "b2", "c1", "c2"),
    batch = c("b1", "b1", "b2", "b2", "b3", "b3"),
    reference_role = c(
      "technical_replicate", "none", "technical_replicate", "none",
      "synthetic_reference", "synthetic_reference"
    ),
    reference_id = c("R1", NA, "R1", NA, "S1", "S1"),
    source_batch = c(NA, NA, NA, NA, "b1", "b2")
  )
  manifest <- sw_sample_manifest(
    meta, "sample_id", "batch",
    reference_role = "reference_role", reference_id = "reference_id",
    source_batch = "source_batch"
  )
  ids <- rep(meta$sample_id, each = 20)
  set.seed(3)
  x <- tibble::tibble(
    sample_id = ids,
    CD3 = stats::rnorm(length(ids)),
    CD4 = stats::rnorm(length(ids))
  )
  batch <- meta$batch[match(ids, meta$sample_id)]
  args <- list(
    x = x, manifest = manifest, design_audit = toy_audit(unique(batch)),
    population = rep("T", length(ids)), markers = c("CD3", "CD4"),
    transform = tf, anchor = "b1", engine = "toy_shift", min_cells = 5, seed = 1
  )
  built <- SpectraWeaveR:::.sw_corrector_cells(
    x, manifest, args$population, args$markers, ids
  )
  is_syn <- built$cells$batch == "b3"
  expect_false(any(built$cells$is_reference[is_syn]))
  expect_true(all(is.na(built$cells$reference_id[is_syn])))
  expect_true(any(built$cells$is_reference[!is_syn]))
  model <- do.call(sw_corrector_fit, args)
  # The toy engine takes effects from is_reference rows; b3 has only
  # synthetic rows, so it gets none.
  expect_setequal(unique(model$effects$batch), c("b1", "b2"))
  expect_false("b3" %in% model$effects$batch)
})
