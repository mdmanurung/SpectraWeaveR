# Tests for the paired reference mean shift engine and its invariants.

# Simulated data: additive batch shifts per (batch, population, marker),
# repeated reference donors, and participants with a case effect on m1 in T.
ps_sim <- function(
  batches = paste0("b", 1:3),
  ref_batches = batches,
  refs = c("R1", "R2"),
  aliquots = 1,
  n_part = 6,
  n_cells = 100,
  seed = 11
) {
  set.seed(seed)
  markers <- c("m1", "m2", "m3")
  pops <- c("T", "B")
  base <- rbind(T = c(2, 1, 0.5), B = c(0.5, 2, 1))
  colnames(base) <- markers
  shift <- array(
    round(stats::runif(length(batches) * 2 * 3, -0.8, 0.8), 2),
    dim = c(length(batches), 2, 3),
    dimnames = list(batches, pops, markers)
  )
  shift[1, , ] <- 0
  ref_offset <- matrix(
    stats::rnorm(length(refs) * 3, 0, 0.3),
    nrow = length(refs),
    dimnames = list(refs, markers)
  )

  # One manifest row per sample: reference aliquots, then participants.
  meta <- list()
  for (b in batches) {
    if (b %in% ref_batches) {
      for (r in refs) {
        for (a in seq_len(aliquots)) {
          meta[[length(meta) + 1]] <- data.frame(
            sample_id = paste(r, b, a, sep = "_"),
            batch = b,
            donor = r,
            aliquot = paste(r, b, a, sep = "_"),
            condition = "reference",
            reference_role = "technical_replicate",
            reference_id = r
          )
        }
      }
    }
    for (p in seq_len(n_part)) {
      meta[[length(meta) + 1]] <- data.frame(
        sample_id = paste0("P", b, "_", p),
        batch = b,
        donor = paste0("P", b, "_", p),
        aliquot = NA_character_,
        condition = if (p %% 2 == 0) "case" else "ctrl",
        reference_role = "none",
        reference_id = NA_character_
      )
    }
  }
  meta <- do.call(rbind, meta)
  part_offset <- matrix(
    stats::rnorm(nrow(meta) * 3, 0, 0.1),
    nrow = nrow(meta),
    dimnames = list(meta$sample_id, markers)
  )

  # Cells: biology + donor offset + batch shift + measurement noise.
  rows <- list()
  for (i in seq_len(nrow(meta))) {
    for (k in pops) {
      off <- if (meta$reference_role[i] == "none") {
        part_offset[meta$sample_id[i], ]
      } else {
        ref_offset[meta$reference_id[i], ]
      }
      mu <- base[k, ] +
        off +
        shift[meta$batch[i], k, ] +
        (meta$condition[i] == "case" & k == "T") * c(1, 0, 0)
      vals <- matrix(
        stats::rnorm(n_cells * 3, rep(mu, each = n_cells), 0.3),
        ncol = 3,
        dimnames = list(NULL, markers)
      )
      rows[[length(rows) + 1]] <- tibble::tibble(
        sample_id = meta$sample_id[i],
        population = k,
        tibble::as_tibble(vals)
      )
    }
  }
  cells <- dplyr::bind_rows(rows)
  manifest <- sw_sample_manifest(
    meta,
    sample_id = "sample_id",
    batch = "batch",
    donor = "donor",
    aliquot = "aliquot",
    condition = "condition",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  list(
    x = cells[, c("sample_id", markers)],
    population = cells$population,
    batch = meta$batch[match(cells$sample_id, meta$sample_id)],
    manifest = manifest,
    meta = meta,
    shift = shift,
    markers = markers
  )
}

# Hand-built audit with one eligibility row, so tests do not depend on
# sw_design_audit() having a paired_shift row yet.
ps_audit <- function(batches, eligible = TRUE, name = "paired_shift") {
  batches <- sort(unique(batches))
  pairs <- if (length(batches) >= 2) {
    utils::combn(batches, 2, simplify = FALSE)
  } else {
    list()
  }
  structure(
    list(
      batches = batches,
      n_batches = length(batches),
      eligibility = tibble::tibble(
        name = name,
        sw_function = "sw_corrector_fit",
        kind = "backend",
        eligible = eligible,
        reason_code = if (eligible) {
          "eligible"
        } else {
          "missing_reference_in_some_batch"
        },
        reason = if (eligible) "test" else "a batch has no shared reference",
        override_allowed = TRUE
      ),
      comparisons = tibble::tibble(
        comparison_id = vapply(pairs, paste, character(1), collapse = "__vs__"),
        batch_a = vapply(pairs, `[`, character(1), 1),
        batch_b = vapply(pairs, `[`, character(1), 2),
        population = NA_character_,
        support = "reference_supported"
      )
    ),
    class = "sw_design_audit"
  )
}

ps_tf <- list(name = "asinh", cofactor = 6000)

ps_fit <- function(s, keep = rep(TRUE, nrow(s$x)), min_cells = 10, ...) {
  sw_corrector_fit(
    s$x[keep, ],
    s$manifest,
    ps_audit(s$batch[keep]),
    population = s$population[keep],
    markers = s$markers,
    transform = ps_tf,
    anchor = "b1",
    engine = "paired_shift",
    min_cells = min_cells,
    seed = 1,
    ...
  )
}

ps_integrate <- function(model, s, keep = rep(TRUE, nrow(s$x))) {
  sw_corrector_integrate(
    model,
    s$x[keep, ],
    s$manifest,
    population = s$population[keep],
    transform = ps_tf,
    unsupported = "abstain"
  )
}

test_that("paired_shift is registered and can freeze", {
  expect_true("paired_shift" %in% sw_corrector_engines())
  expect_true(.sw_corrector_engine("paired_shift")$can_freeze)
})

test_that("known additive shifts are recovered on reference samples", {
  s <- ps_sim()
  model <- ps_fit(s)
  eff <- model$effects
  true <- s$shift[cbind(eff$batch, eff$population, eff$marker)]
  expect_lt(max(abs(eff$estimate - true)), 0.1)
  expect_true(all(eff$n_donors == 2L))
  expect_true(all(eff$n_aliquots == 2L))
})

test_that("participants get the reference effect and keep their biology", {
  s <- ps_sim()
  model <- ps_fit(s)
  res <- ps_integrate(model, s)
  is_part <- !grepl("^R", s$x$sample_id)

  # Corrected participant values carry no residual batch shift.
  true <- s$shift[cbind(s$batch, s$population, "m1")]
  resid <- res$x$m1 - (s$x$m1 - true)
  expect_lt(max(abs(tapply(resid[is_part], s$batch[is_part], mean))), 0.1)

  # Case-control difference in m1 (T) per donor survives correction.
  in_t <- is_part & s$population == "T"
  donor_mean <- tapply(res$x$m1[in_t], s$x$sample_id[in_t], mean)
  cond <- s$meta$condition[match(names(donor_mean), s$meta$sample_id)]
  diff <- mean(donor_mean[cond == "case"]) - mean(donor_mean[cond == "ctrl"])
  expect_lt(abs(diff - 1), 0.2)

  # Within a batch the correction is a constant, so contrasts are unchanged.
  b2 <- in_t & s$batch == "b2"
  before <- tapply(s$x$m1[b2], s$x$sample_id[b2], mean)
  after <- tapply(res$x$m1[b2], s$x$sample_id[b2], mean)
  expect_equal(
    unname(after - after[1]),
    unname(before - before[1]),
    tolerance = 1e-12
  )
})

test_that("the anchor batch has zero effect and unchanged values", {
  s <- ps_sim()
  model <- ps_fit(s)
  anc <- model$effects[model$effects$batch == "b1", ]
  expect_true(nrow(anc) > 0)
  expect_identical(anc$estimate, rep(0, nrow(anc)))
  expect_true(all(anc$support == "anchor"))
  res <- ps_integrate(model, s)
  in_b1 <- s$batch == "b1"
  expect_identical(res$x[in_b1, ], s$x[in_b1, ])
})

test_that("new batches do not change rows corrected earlier", {
  s <- ps_sim(batches = paste0("b", 1:6), ref_batches = paste0("b", 1:5))
  model <- ps_fit(s, keep = s$batch %in% paste0("b", 1:3))
  ia <- which(s$batch %in% paste0("b", 1:4))
  ib <- which(s$batch %in% paste0("b", 1:6))
  res_a <- ps_integrate(model, s, keep = ia)
  res_b <- ps_integrate(model, s, keep = ib)

  # Rows of b1..b4 are bitwise identical whether or not b5 and b6 arrive.
  expect_identical(res_b$x[match(ia, ib), ], res_a$x)

  # Old effects in the extended copy are the frozen ones; the input is untouched.
  old_b <- res_b$model$effects[res_b$model$effects$source == "fit", ]
  expect_identical(old_b, model$effects)
  expect_identical(res_b$model$parent_digest, model$digest)
  expect_false(identical(res_b$model$digest, model$digest))
  expect_setequal(res_b$model$known_batches, paste0("b", 1:5))

  # The new batch b4 is estimated from its own references.
  b4 <- res_a$model$effects[res_a$model$effects$batch == "b4", ]
  true <- s$shift[cbind(b4$batch, b4$population, b4$marker)]
  expect_lt(max(abs(b4$estimate - true)), 0.1)
  expect_true(all(b4$source == "integrate"))

  # b6 has no bridging reference: abstain and keep the original values.
  act6 <- res_b$actions[res_b$actions$batch == "b6", ]
  expect_true(all(act6$action == "abstain"))
  expect_true(all(act6$reason_code == "no_bridging_reference"))
  in6 <- s$batch[ib] == "b6"
  expect_identical(res_b$x[in6, ], s$x[ib[in6], ])
})

test_that("row order does not change per-row results", {
  s <- ps_sim(batches = paste0("b", 1:4))
  keep_fit <- s$batch %in% paste0("b", 1:3)
  model <- ps_fit(s, keep = keep_fit)
  res <- ps_integrate(model, s)

  set.seed(3)
  perm <- sample(nrow(s$x))
  sp <- s
  sp$x <- s$x[perm, ]
  sp$population <- s$population[perm]
  sp$batch <- s$batch[perm]
  res_p <- ps_integrate(model, sp)
  expect_identical(res_p$x[order(perm), ], res$x)
  expect_identical(res_p$model$effects, res$model$effects)

  # Fitting on permuted rows gives the same effects.
  model_p <- ps_fit(sp, keep = sp$batch %in% paste0("b", 1:3))
  expect_identical(model_p$effects, model$effects)
  expect_identical(model_p$targets, model$targets)
})

test_that("integrating twice with the same model gives identical results", {
  s <- ps_sim(batches = paste0("b", 1:4))
  model <- ps_fit(s, keep = s$batch %in% paste0("b", 1:3))
  expect_identical(ps_integrate(model, s), ps_integrate(model, s))
})

test_that("abstain keeps original values exactly", {
  s <- ps_sim()
  model <- ps_fit(s)
  # A population the model never saw, inside a known batch.
  s$population[s$batch == "b2" & s$population == "B"] <- "NK"
  res <- ps_integrate(model, s)
  nk <- s$population == "NK"
  act <- res$actions[res$actions$population == "NK", ]
  expect_true(all(act$action == "abstain"))
  expect_true(all(act$reason_code == "population_unknown"))
  expect_identical(res$x[nk, ], s$x[nk, ])
  expect_error(
    sw_corrector_integrate(
      model,
      s$x,
      s$manifest,
      s$population,
      transform = ps_tf,
      unsupported = "error"
    ),
    "cannot be corrected"
  )
})

test_that("a fit-time batch without a bridging reference is never re-estimated", {
  s <- ps_sim(batches = paste0("b", 1:3), ref_batches = c("b1", "b2"))
  model <- ps_fit(s)
  expect_false("b3" %in% model$effects$batch)
  res <- ps_integrate(model, s)
  act <- res$actions[res$actions$batch == "b3", ]
  expect_true(all(act$reason_code == "known_batch_no_effect"))
})

test_that("min_cells excludes small summaries with a stated reason", {
  s <- ps_sim(aliquots = 2)
  # Keep only 5 T cells of one aliquot of R1 in b2.
  small <- s$x$sample_id == "R1_b2_2" & s$population == "T"
  drop <- which(small)[-(1:5)]
  s$x <- s$x[-drop, ]
  s$population <- s$population[-drop]
  s$batch <- s$batch[-drop]
  model <- ps_fit(s, min_cells = 10)

  exc <- model$excluded
  expect_equal(nrow(exc), 1)
  expect_identical(exc$aliquot, "R1_b2_2")
  expect_identical(exc$n_cells, 5L)
  expect_match(exc$reason, "below min_cells = 10")
  eff <- model$effects[model$effects$batch == "b2", ]
  expect_true(all(eff$n_aliquots[eff$population == "T"] == 3L))
  expect_true(all(eff$n_aliquots[eff$population == "B"] == 4L))

  diag <- sw_corrector_diagnose(ps_integrate(model, s))
  row <- diag[diag$batch == "b2" & diag$population == "T", ]
  expect_true(all(row$n_excluded_summaries == 1L))
})

test_that("aliquots are summarised separately and weighted equally", {
  meta <- data.frame(
    sample_id = c("a1", "a2", "a3"),
    batch = c("b1", "b2", "b2"),
    aliquot = c("x1", "y1", "y2"),
    reference_role = "technical_replicate",
    reference_id = "R1"
  )
  manifest <- sw_sample_manifest(
    meta,
    sample_id = "sample_id",
    batch = "batch",
    aliquot = "aliquot",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  # Anchor at 0; in b2 one aliquot of 30 cells at 1 and one of 10 at 3.
  x <- tibble::tibble(
    sample_id = rep(c("a1", "a2", "a3"), c(20, 30, 10)),
    m1 = rep(c(0, 1, 3), c(20, 30, 10))
  )
  fit <- function(man) {
    sw_corrector_fit(
      x,
      man,
      ps_audit(c("b1", "b2")),
      population = rep("T", 60),
      markers = "m1",
      transform = ps_tf,
      anchor = "b1",
      engine = "paired_shift",
      min_cells = 5,
      seed = 1
    )
  }
  model <- fit(manifest)
  eff <- model$effects[model$effects$batch == "b2", ]
  # Equal aliquot weights give 2; a pooled cell mean would give 1.5.
  expect_equal(eff$estimate, 2)
  expect_identical(eff$n_aliquots, 2L)
  expect_identical(model$aliquot_source, "manifest aliquot column")

  # Without an aliquot column the key falls back to the reference id.
  manifest2 <- sw_sample_manifest(
    meta[, c("sample_id", "batch", "reference_role", "reference_id")],
    sample_id = "sample_id",
    batch = "batch",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  model2 <- fit(manifest2)
  expect_match(model2$aliquot_source, "defaulted to reference_id")
  expect_equal(model2$effects$estimate[model2$effects$batch == "b2"], 1.5)
  res <- sw_corrector_integrate(
    model2,
    x,
    manifest2,
    rep("T", 60),
    transform = ps_tf,
    unsupported = "abstain"
  )
  expect_match(sw_corrector_diagnose(res)$aliquot_source[1], "defaulted")
})

test_that("paired_shift rejects a hierarchy and needs its audit row", {
  s <- ps_sim()
  expect_error(ps_fit(s, hierarchy = ~study), "paired_shift has no hierarchy")
  args <- list(
    x = s$x,
    manifest = s$manifest,
    population = s$population,
    markers = s$markers,
    transform = ps_tf,
    anchor = "b1",
    engine = "paired_shift",
    min_cells = 10,
    seed = 1
  )
  expect_error(
    do.call(
      sw_corrector_fit,
      c(args, list(design_audit = ps_audit(s$batch, name = "cycombine")))
    ),
    "audit has no paired_shift eligibility row"
  )
  expect_error(
    do.call(
      sw_corrector_fit,
      c(args, list(design_audit = ps_audit(s$batch, eligible = FALSE)))
    ),
    "ineligible"
  )
  expect_error(
    do.call(
      sw_corrector_fit,
      c(args, list(design_audit = ps_audit(c("b1", "b2"))))
    ),
    "same data"
  )
  args$anchor <- "b9"
  expect_error(
    do.call(sw_corrector_fit, c(args, list(design_audit = ps_audit(s$batch)))),
    "Anchor batch 'b9'"
  )
})

test_that("diagnose reports one row per batch, population and marker", {
  s <- ps_sim(batches = paste0("b", 1:4), ref_batches = paste0("b", 1:3))
  model <- ps_fit(s, keep = s$batch %in% paste0("b", 1:3))
  diag <- sw_corrector_diagnose(ps_integrate(model, s))
  expect_equal(nrow(diag), 4 * 2 * 3)
  expect_setequal(
    unique(diag$batch_status),
    c("anchor", "known", "new_unbridged")
  )
  expect_true(all(
    c("n_donors", "n_aliquots", "support", "reason") %in% names(diag)
  ))
})

test_that("save and load keep the paired_shift model exactly", {
  s <- ps_sim()
  model <- ps_fit(s)
  path <- tempfile(fileext = ".rds")
  sw_corrector_save(model, path)
  loaded <- sw_corrector_load(path)
  expect_identical(loaded, model)
  expect_identical(ps_integrate(loaded, s)$x, ps_integrate(model, s)$x)
})

test_that("compare rewards the paired shift without distorting spread", {
  s <- ps_sim()
  model <- ps_fit(s)
  res <- ps_integrate(model, s)
  rec <- sw_corrector_compare(
    list(paired_shift = res),
    ps_audit(s$batch),
    uncorrected = s$x,
    manifest = s$manifest,
    population = s$population,
    markers = s$markers,
    min_cells = 10
  )
  expect_setequal(unique(rec$candidate_id), c("uncorrected", "paired_shift"))
  expect_true(all(rec$calibration == "uncalibrated"))
  emd <- rec[rec$metric == "emd_within_population" & rec$status == "ok", ]
  expect_lt(
    sum(emd$value[emd$candidate_id == "paired_shift"]),
    sum(emd$value[emd$candidate_id == "uncorrected"]) / 2
  )
  mad <- rec[rec$metric == "mad_compression_fraction", ]
  expect_true(all(mad$value == 0))
  jsd <- rec[rec$metric == "joint_structure_distance", ]
  expect_lt(max(jsd$value), 1e-10)

  # The records feed sw_batch_candidates_from_metrics().
  if (!exists("sw_batch_candidates_from_metrics", mode = "function")) {
    f <- testthat::test_path("..", "..", "R", "batch_decide.R")
    skip_if_not(file.exists(f), "batch_decide.R not available")
    env <- new.env()
    sys.source(f, envir = env)
    from_metrics <- env$sw_batch_candidates_from_metrics
  } else {
    from_metrics <- sw_batch_candidates_from_metrics
  }
  cand <- from_metrics(rec, identity_id = "uncorrected")
  expect_equal(nrow(cand), 2)
  expect_gt(
    cand$emd_relative_reduction[cand$candidate_id == "paired_shift"],
    0.5
  )
})
