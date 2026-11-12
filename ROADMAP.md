# SpectraWeaveR feature roadmap

Goal: one R package that takes spectral flow data from raw FCS to batch-checked, clustered output, with every correction step gated by what the design can support.

Status key: **done** (code and tests in the repo), **next** (designed, not written), **planned**, **gated** (starts only after a stated evidence gate), **parked** (set aside; the row names its return condition), **frozen** (code exists; no new work until a stated condition holds).

Decisions of 2026-11-09. Track B items B3 to B10 are parked. They return when independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. The flowstate I/O items (D1, D2) and the batch-side items (A1 to A6, C0 to C5: design audit, decision layer, IQR flags, manifest, correctors) are frozen: no new work until a real multi-batch dataset exists. Update 2026-11-11 (user, plan grill): a real multi-batch dataset is on disk (`cytobatchflagr_zenodo`: 18 batches, 4 repeated controls per batch), so the new batch-diagnostic, reference-check and synthetic-reference functions in `plans/2026-11-11_reference-and-batch-diagnostics.md` may be built once that plan is approved, and A2c (`sw_correct_metrics`) is unfrozen. A2b and A3 stay frozen.

## Source material

| Source | Role | License | Integration |
|---|---|---|---|
| AutoSpectral | Unmixing engine | AGPL | Wrap from Suggests. No code copied. |
| spectracle, flowstate | Spectra from controls; data.table FCS I/O | AGPL | Wrap from Suggests. No code copied. |
| CytoBatchFlagR | Control-based batch flags (IQR, EMD, cluster) | GPL-3 | Clean rewrite from the method description only. Function bodies not read. |
| batchdx-cyto | Design audit, decision layer | MIT | Ported to R. |
| `plans/AutoSpectral_Implementation_Bundle_v1.1.zip` | Study-level calibration: pick AutoSpectral settings per study, then audit | `[gap]` no license file in bundle | Track B. Confirm the license before copying `prototype/*.R`; otherwise reimplement from the spec. |
| `plans/cyTether_implementation_benchmark_bundle_v1.zip` | Reference-anchored correction onto a frozen calibration scale | MIT | Track C. |
| Prism-style figure theme | `.sw_theme_prism()` in `R/theme_prism.R` | MIT (package) | Theme code supplied by the package author, who states they hold its copyright (user decision 2026-11-11, grill R3.3). The code matches `theme_prism()` in ggprism (github.com/csdaw/ggprism, GPL >= 3); the copyright statement is the author's and was not verified here. Made standalone: no ggprism dependency. |

Both bundles are design plus prototypes. Neither is validated on real FCS data. Their own validation notes say so, and this roadmap keeps that wording.

## Rules for every track

1. **Wrap AGPL and GPL code; do not copy it.** Check with `requireNamespace(..., quietly = TRUE)` and call `pkg::fn()`.
2. **Impactful choices have no default.** A threshold, limit, tolerance, grid or seed that can change a ranking or a selected setting is a required argument. Its roxygen text keeps the "uncalibrated" label until calibration evidence exists. See the register at the end.
3. **Design before correction.** `sw_design_audit()` runs first. A confounded or non-identifiable contrast blocks automatic correction; abstaining or keeping the uncorrected data is a valid result.
4. **Abstain, never invent.** Missing input gives a stated status and reason, never a zero or a fabricated flag.
5. **Plots need a grill plan** (plot-first rule). No plotting function is added until the user approves a plan row for it.
6. **Held-out evidence stays closed** during selection. Truth files and audit blocks are never read while fitting or tuning.

## Gold standard: latest AutoSpectral

AutoSpectral is the reference for unmixing. As of 2026-11-08 the latest is **1.8.8, commit `3c504c4560cf1f0fd5dce298eaabaafaebfae9be`** (HEAD of the default branch, same commit the study-calibration bundle pins). Re-check the latest release before each milestone; when it moves, update this line and the pin together.

How extensions are weighed against it:

1. **Baseline first.** Every extension (spectracle spectra, study-calibration candidates, secondary correction, GLS noise, AF subspaces) is scored against a run of the unmodified AutoSpectral default on the same events, controls and panel. No extension is reported alone.
2. **Non-inferior or dropped.** An extension is kept only if it is no worse than the baseline on the critical guardrails (dim-positive recovery, negative width, rare-AF coverage) and better on at least one ranking criterion. The non-inferiority margin is the user's choice (required argument, no default).
3. **Call, don't rebuild.** Where AutoSpectral already has a routine (extraction, AF, variants, joint unmixing, noise model, mismatch and variability checks), the package wraps it. A reimplementation needs a written reason and a parity test against the AutoSpectral output.
4. **Same units.** Compare in AutoSpectral's units: spectra in rows, detectors in columns, linear intensity. Any other convention is converted at the function boundary with a name check.
5. **Record the baseline.** Each result stores the AutoSpectral version, commit, AutoSpectralRcpp version and whether the R or compiled backend ran.
6. **Version drift is a finding.** If the installed AutoSpectral differs from the pin, results carry a warning; archived OMIP outputs come from earlier code and serve as regression cases only after a rerun on the pin.
7. **Suggests floor.** `DESCRIPTION` has `AutoSpectral (>= 1.8.8)`, the version the wrappers were read from.

Where an extension is already outside AutoSpectral's scope (batch detection, batch correction, cyTether), the baseline is "no correction" plus the corrector with the same inputs, and unmixed AutoSpectral output is the fixed upstream input.

## Track A: batch diagnosis and control QC

| # | Item | Function(s) | Status |
|---|---|---|---|
| A1 | Design audit: pair support classes, confounding, reference graph, method eligibility | `sw_design_audit()` in `R/design_audit.R` | **done**; frozen: no new work until a real multi-batch dataset exists. Batchdx port plus `study_col`, `run_col`, `donor_col`, `aliquot_col`, `anchor`, an `anchor_bridge` table and a `paired_shift` eligibility row (`tests/testthat/test-design_audit.R`). |
| A2 | Decision layer with identity candidate always kept and the states of batchdx | `sw_batch_policy()`, `sw_batch_decide()`, `sw_batch_candidates_from_metrics()` | **done**, 68 expectations pass. Open: it has no donor-level protected-contrast gate (A2b), and no bridge from `sw_correct_evaluate_*()` output (A2c). Frozen: no new work until a real multi-batch dataset exists. |
| A2b | Add a gate on donor-level protected contrasts. Without it a method that mixes batches well but erases condition differences can pass. Needs donor and condition from the manifest. | `sw_batch_decide()` | **next** Frozen: no new work until a real multi-batch dataset exists. |
| A2c | Bridge: turn `sw_correct_*` results into the long metric tibble `sw_batch_candidates_from_metrics()` reads. cyCombine's `evaluate_*` lists lack per-comparison EMD, joint-structure distance and compression fraction, so this needs new metric code. | `sw_correct_metrics()` | **next** Frozen: no new work until a real multi-batch dataset exists. |
| A3 | Gate `sw_correct_run()` and `sw_pipeline_run_all()` on the audit | argument `design_audit =` | **next**, after C0 (needs the shared manifest as join key) Frozen: no new work until a real multi-batch dataset exists. |
| A4 | Control IQR flags (clean rewrite). Thresholds, multiplier, minimum cells and quantile type are required arguments. | `sw_flag_iqr()` in `R/batch_flag.R` | **done**, 54 expectations pass. `min_effect` (required named vector per metric, `effect` column, status `below_min_effect`) added in P0c. Open: no value is chosen; DR9 needs values from a null-arm run and user approval. Frozen: no new work until a real multi-batch dataset exists. |
| A5 | EMD and cluster-based flags | reuse `sw_correct_emd()` and `sw_cluster_run()` | **partly done** (2026-11-12, P3b): per-population EMD rows through `sw_batch_scores(population =)` and `sw_batch_interaction()`, ranked rows, no flags. Gate or cluster labels; cluster labels cap at `eligible_for_review`. Tests: `tests/testthat/test-batch_scores.R`. Unfrozen by the user 2026-11-12. |
| A6 | Drift and composition detectors from batchdx | `sw_flag_drift()`, `sw_flag_composition()` | **partly done** (2026-11-12, P3b) as score rows, not flags: `batch_order_rho` (Spearman of batch median against batch order) and per-batch population share rows `pop_share`, `share_rest_z` in `sw_batch_scores(population =)`. No `sw_flag_*` functions. Unfrozen by the user 2026-11-12. |

## Track B: unmixing and study calibration (AutoSpectral bundle)

Pin: the latest AutoSpectral (see "Gold standard"). Currently 1.8.8, commit `3c504c4560cf1f0fd5dce298eaabaafaebfae9be`. The `Suggests` floor is 1.8.8. Record the R, AutoSpectral, AutoSpectralRcpp versions and backend in every result. B2's spectracle path is an extension: score it against the AutoSpectral-extracted spectra on the same controls before recommending it.

| # | Item | Function(s) | Status |
|---|---|---|---|
| B1 | Expose AutoSpectral diagnostics (mismatch, variability, condition number, AF and unmix comparison) | `sw_unmix_assess_*()`, `sw_unmix_compare*()` in `R/unmix_assess.R` | **done**. Two stay exported (`sw_unmix_check_controls()`, `sw_unmix_condition_number()`). Seven became internal on 2026-11-09 (see the register note). Tests: `tests/testthat/test-unmix_assess.R`; the real-AutoSpectral call test needs AutoSpectral. Also `sw_unmix_diagnose()` has `max_sif` (check id `spectra_sif`): the diagonal of the hotspot matrix (the spreading inflation factor), with the largest off-diagonal partner named in the reason. It uses the formula of `AutoSpectral::calculate.hotspot.matrix()` in 1.8.8, in base R; a parity test is in `tests/testthat/test-unmix_report.R`. |
| B2 | spectracle spectra to AutoSpectral spectra, with an explicit rule for the appended AF row | `sw_io_spectracle_to_spectra()`, `sw_unmix_spectracle()` in `R/spectracle.R` | **done** (`tests/testthat/test-spectracle.R`); the real-spectracle smoke test is skipped where spectracle is not installed |
| B3 | Study manifest intake and preflight: assay role vs analysis use, detector map, control compatibility, checksums, lineage. JSON Schema from the bundle. | `sw_study_manifest()`, `sw_study_preflight()` | **parked** Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B4 | Candidate runner: fixed OLS and the pinned AutoSpectral default first, one change per block after that. Each candidate is a typed settings list with a stable ID. | `sw_study_candidate()`, `sw_study_run()` | **parked** Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B5 | Score table keyed by `candidate_id, settings_id, phase, block_id, criterion_id`; losses where larger is worse; the full expected grid recorded so absent scores cannot vanish | `sw_study_scores()` | **parked** Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B6 | Selection: critical guardrail rejects; rank by worst criterion mean of `(loss - limit)/scale`; near-best tolerance then lower complexity | `sw_study_select()` | **parked**. Reimplement from the spec: the bundle has no license file and its prototype is standalone base R. Note `study_policy.R` rejects on any breach, while the spec says critical breaches only; follow the spec. Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B7 | Lock and audit: write the lock atomically before final scores open; observed mode, or exact binomial block-risk mode (`stats::binom.test()`). With a 0.20 risk limit and zero failures, about 17 blocks are needed to qualify, so the bundle's 4-block minimum can only return inconclusive. The block count is a required argument. | `sw_study_lock()`, `sw_study_audit()` | **parked** Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B8 | Run report (JSON, CSV, readable) with provenance and event counts | `sw_study_report()` | **parked** Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B9 | Secondary correction, GLS noise model, AF subspaces | candidates in B4 | **parked** on B4 to B7 passing on one real regime Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |
| B10 | Ideas from public Ozette Resolve descriptions, built from open methods only. Resolve's model, positivity formula and output scale are not published, so none of it is copied or inferred. Candidates: non-negative least squares (Lawson-Hanson) as an unmixing candidate; a control event filter by cosine similarity to an empirical median spectrum; paired-unstained provenance. A positivity score would be our own design (for example a tail probability against matched unstained events) and is not "Resolve's score". Review: `scratchpad/research/ozette_resolve_public_review.md`. | candidates in B4 | **parked** like B9; each is scored against unmodified AutoSpectral Return condition: independent evaluation tubes (matched raw controls plus held-out stained tubes) exist. |

B3 to B10 are parked. Open data need: complete matched raw controls and independent evaluation tubes. The bundle records that the OMIP-102 and OMIP-095 FlowRepository retrievals returned service errors, so B3 to B7 are tested on synthetic fixtures first.

## Track C: reference-anchored integration (cyTether bundle)

The bundle's R scaffold implements paired-reference mean shift only. Its Python hierarchical empirical Bayes (HEB) prototype uses moment estimates and has not beaten simple paired shift: in the bundle's two-seed pilot, paired shift gave lower protected-contrast error. So HEB is gated, and paired shift is the default until it is not.

| # | Item | Function(s) | Status |
|---|---|---|---|
| C0 | **Corrector contract, shared manifest and population model.** This comes before any engine. See "Corrector contract" below. | `sw_sample_manifest()`, `sw_corrector_*()` in `R/manifest.R`, `R/corrector.R` | **done**; frozen: no new work until a real multi-batch dataset exists. Tests: `test-manifest.R`, `test-corrector.R`. Open: `sw_pipeline_run_all()` still has no population step. P0b (2026-11-12): role `synthetic_reference` with required `source_batch`, never a bridge; audit lists it in `synthetic_references`. |
| C1 | Calibration graph: a (donor, batch, population) reference graph, `study_col` and `run_col`, a check that each batch maps to one study, donor and aliquot counts, and an eligibility rule "every batch shares a reference donor with the anchor batch". The current audit counts any shared reference, but both cyTether engines can only use donors measured in the anchor batch. | extend `sw_design_audit()` | **done**; frozen: no new work until a real multi-batch dataset exists. `study_col`, `run_col`, `aliquot_col`, `anchor`, `anchor_bridge` and the `paired_shift` eligibility row exist and are tested. |
| C2 | Paired-reference mean shift as the first engine, `engine = "paired_shift"`. Fit on reference summaries keyed by (donor, aliquot, batch, population); cells are not the replicate unit. Verify against the bundle's `S*/reference_shift.npz` outputs. Anchor choice is a required argument: defaulting to the first sorted batch turns its noise into bias everywhere. | `sw_corrector_fit(engine = "paired_shift")` in `R/corrector_paired_shift.R` | **done**; frozen: no new work until a real multi-batch dataset exists. Anchor has no default. Tests: `test-corrector_paired_shift.R` on synthetic data. Missing: the check against the bundle's `S*/reference_shift.npz` outputs. |
| C3 | Frozen model: marker names, transform, schema, digest; reject incompatible marker sets; abstain for unsupported phenotype or marker and keep the original value | `sw_corrector_save()`, `sw_corrector_load()` | **done**; frozen: no new work until a real multi-batch dataset exists. Round trip, tampered-file failure, marker and transform mismatch at integrate, and abstain are tested. |
| C4 | Invariants as tests: earlier corrected values bitwise unchanged when new data arrive; order invariance; idempotent re-integration | tests for C2, C3 | **done**: new batches leave earlier rows unchanged, order invariance and repeat integration are tested in `test-corrector_paired_shift.R`; frozen: no new work until a real multi-batch dataset exists. |
| C5 | Bridge to benchmark scenarios S0 to S12 (export from the bundle's `.npz`) to compare `sw_corrector_*` against `sw_correct_run()` | `inst/benchmarks/` | **planned**; needs Python, kept out of the package build |
| C6 | Joint HEB by restricted marginal likelihood, validated against `lme4` and a small Stan fit. Fix the prototype gaps first: no pooled variance across markers and populations, no hard-coded two runs per study, no cells-as-independent variance with a 0.02 floor, aliquot level estimable, a new run conditional on the frozen study posterior. Fit raw summaries with donor × population fixed effects and the anchor fixed at zero, instead of modelling covariance of paired differences. | `sw_corrector_fit(engine = "heb")` | **gated**: must beat C2 across at least 100 seeds on held-out scenarios, with an ablation. Preregister the metric set and the two-stage rule now. |
| C7 | Location-scale, robust tails, low-rank signatures, instrument and stain crossed effects | extensions | **gated**: each needs a predeclared safety gate and a failure mode |
| C8 | **Synthetic references from real study events** (plan `plans/2026-11-11_reference-and-batch-diagnostics.md`, DR21, DR23, DR25). Diagnostic comparator only, never a bridge. `n_events`, `n_draws`, `seed`, `strata`, `n_parts`, `events_per_part` and `max_overlap_share` are required, no defaults. | `sw_synref_downsample()`, `sw_synref_mix()`, `sw_synref_split()`, `sw_synref_spread()`, `sw_synref_fidelity()`, `sw_synref_goal()`, `sw_synref_check_overlap()`, `sw_synref_as_manifest()`, `sw_synref_save()`, `sw_synref_load()` in `R/synref_draw.R`, `R/synref_checks.R`, `R/synref_io.R` | **done (P1)**: class, draws, split, overlap guard, manifest rows, save and load with digest. **done (P2b)**: `sw_synref_mix()`, `sw_synref_spread()`, `sw_synref_fidelity()`, `sw_synref_goal()` (tests in `test-synref_draw.R`, `test-synref_checks.R`). Tests: `test-synref_draw.R`, `test-synref_checks.R`, `test-synref_io.R`. `sw_synref_as_manifest()` rows pass `sw_sample_manifest()` and `sw_design_audit()` (tested). |
| C9 | **Reference-sample consistency checks** (plan `plans/2026-11-11_reference-and-batch-diagnostics.md`, DR2, DR15, DR16, DR37). Does a repeated reference look the same in every batch, read against a chance floor. Numbers only, ranked within check and marker, no per-reference verdict. `min_events`, `roles`, `n_events`, `seed`, `order`, `min_batches_per_reference`, `transform`, `cofactor`, `n_splits`, `n_quantiles` are required, no defaults. | `sw_ref_inventory()`, `sw_ref_consistency()`, `sw_ref_noise_floor()` in `R/ref_inventory.R`, `R/ref_consistency.R`, `R/ref_noise_floor.R`; internal `.sw_batch_scores_groups()` in `R/batch_scores.R` | **done (P2c)**. Tests: `test-ref_inventory.R`, `test-ref_consistency.R`, `test-ref_noise_floor.R`. Open: `n_quantiles = NULL` means the exact EMD (user to confirm the quantile-grid meaning); no real-data run yet. |

**Evidence note.** The bundle's pilot has two seeds, four stages and two reference donors. The direction flips by stage and seed (S6 seed 7001: HEB 0.165 vs shift 0.185; S4: shift wins), and the S4 to S7 truths are additive shifts, which is the case where paired shift is correctly specified. So the pilot cannot separate the two methods. Paired shift is the default because it is simpler, not because it was shown better. Add scenarios with many runs and few reference cells, where shrinkage can matter.

**Scale limit to document.** An additive shift on `asinh(x/cofactor)` is a gain change only for bright populations. Dim and negative unmixed values see a scale change instead.

Do not merge C6 and C7 into one untestable model. The bundle's own rule applies: each change needs an ablation, a simulated truth and a named failure mode.

## Corrector contract (C0)

cyTether is the archetype: other correctors (cyCombine through `sw_correct_run()`, ComBat, CytoNorm) sit behind the same functions and are compared on the same audit and metrics.

```r
sw_sample_manifest(df, sample_id, batch, study = NULL, run = NULL, instrument = NULL,
                   donor = NULL, aliquot = NULL, condition = NULL,
                   reference_role = NULL, reference_id = NULL)
sw_corrector_fit(x, manifest, design_audit, population, markers, transform,
                 hierarchy, anchor, engine, min_cells, ...)      # -> sw_frozen_model
sw_corrector_integrate(model, x, manifest, population,
                       unsupported = c("abstain", "error"))       # -> list(x, actions)
sw_corrector_diagnose(result)
sw_corrector_save(model, path); sw_corrector_load(path)
sw_corrector_compare(candidates, design_audit, metrics, gates)   # identity always included
```

- **One manifest.** `sw_correct_prepare()` hard-codes `sample, batch, condition`; the audit uses `sample_id` and a `batch_col`. The manifest replaces both and is the join key for the audit, the corrector and the study-calibration intake (B3).
- **Population model.** cyTether fits per (batch, population, marker) and needs frozen per-cell labels before correcting. `sw_pipeline_run_all()` clusters after correction, so it has no population step. C0 adds one: gate-derived labels or a frozen classifier.
- **Frozen model schema.** schema version, markers (name, clone, fluorochrome, detector), transform and cofactor, population model, hierarchy, anchor, hyperparameters, effects and covariance per (batch, population, marker), support, abstention rules, software versions, seed, digest. Record the transform and cofactor in the model and in the B7 lock, and check both at integrate.
- **Arrival cases.** Known batch: use the frozen effect. New run in a known study: condition on the frozen study posterior. New study: estimate its effects only, leave old ones untouched. No bridging reference: abstain and keep the original value.
- **`diagnose` returns** one row per (batch, population, marker): estimate, SE, donors, aliquots, support, connectivity, action (apply, shrink, abstain) and reason.
- **`sw_correct_run()` becomes an adapter** with `can_freeze = FALSE`. A method that refits on new data then fails the invariance gate (C4) visibly, instead of silently.
- **Selection** must include a donor-level protected-contrast metric. In the bundle's pilot, naive batch centering has the best held-out reference MAD (0.084) and the second-worst protected-contrast error (0.167).

## Track D: I/O and package engineering

| # | Item | Status |
|---|---|---|
| D1 | flowstate as an alternative FCS reader, with bridges to flowFrame and tibble | **done** (`R/flowstate_io.R`, `test-flowstate_io.R`); frozen: no new work until a real multi-batch dataset exists. The real-flowstate test is skipped where flowstate is not installed. |
| D2 | `DESCRIPTION`: add Suggests `flowstate`, `spectracle`; Remotes `nlaniewski/flowstate`, `nlaniewski/spectracle`; extend `inst/scripts/install_dependencies.R` | **done**; frozen: no new work until a real multi-batch dataset exists. |
| D3 | `NAMESPACE`, man pages, `CLAUDE.md` pipeline table, `PLAN.md`, `PROGRESS.md` | after each item; run `devtools::document()` |
| D4 | `R CMD check` on a machine with the Bioconductor stack | **planned**; not possible in the current shell (no flowCore) |
| D5 | Panel object (marker, clone, fluorochrome, detector). `sw_correct_prepare()` matches `$PnN` only. | **planned**, before C3 |
| D6 | Provenance block in every result: package versions, seeds, `sessionInfo()`. `sw_io_event_audit()` exists but nothing calls it. | **planned** |
| D7 | Exporters: FCS write, `SingleCellExperiment`, AnnData; chunked application for millions of cells | **planned** |
| D8 | Live/dead, debris and doublet gates inside `sw_pipeline_run_all()` (the helpers exist in `gating_utils.R` but the pipeline does not call them) | **planned** |
| D9 | `SingleCellExperiment` API, pkgdown, vignette for the audit-then-correct workflow | **planned**, last |
| D10 | Folder structure: package vs `experiments/` (datasets registry, benchmarks, integration, slurm, runs). Done 2026-11-09. `inst/scripts/setup_pixi_env.sh` moved to `tools/setup_pixi_env.sh` on 2026-11-09; the `setup` task in `pixi.toml` and `.Rbuildignore` are updated. | `experiments/README.md` | **done** |
| D11 | Cleaning stage in one call: margins (raw), arcsinh, PeacoQC, then nonDebris, singlets and live gates through the openCyto path. Output is the live singlet events plus a per-sample audit table, the gate values and the PeacoQC summary. Nothing after live. Covers the debris, doublet and live/dead part of D8 for stand-alone use; D8 (wiring into `sw_pipeline_run_all()`) stays open. `sw_gate_run()` now also takes a `flowSet`. `sw_transform_estimate()` is logicle only and has no cofactor, so the arcsinh list is built in the function and applied with `sw_transform_apply()`. Margins take a required `margin_range = c(min, max)` for the fluorescence channels, passed to PeacoQC `RemoveMargins()` through `sw_filter_margins_peacoqc(channel_specifications = list(AllFluoChannels = margin_range))`; scatter keeps the FCS ranges. A lower limit above 0 is read by PeacoQC as 0 (`min(minRange, 0)`); `sw_clean_run()` warns, naming `margin_range` and the value, and continues (a warning rather than an error). Measurements: `experiments/runs/2026-11-11_clean-run-smoke/warnings.txt` (run folder not tracked by git) lists PeacoQC margin warnings for A3 99.34, A6 99.81, B11 99.48, B2 99.16, D1 99.86 percent; benchmark 05 passed c(-21761, 1e7) (`12_cleanup.R:28`, user 2026-11-09); margin_range required (user, 2026-11-11 plan grill). A sample emptied by the margin step, and a PeacoQC error on a sample, stop the run with the sample named; no minimum event count is set. `pool_samples = TRUE` writes `groupBy` = number of samples (openCyto pools only with a non-empty `groupBy`), so it gives one shared cut. `sw_gate_run()` now does the same for any template row with `collapseDataForGating` TRUE and an empty `groupBy` (message, temporary copy of the template; rows with an explicit `groupBy` are not touched), so direct callers with `pool_samples = TRUE` in `sw_gate_template()` now pool. | `sw_clean_run()` in `R/clean_gate.R` | **done** (2026-11-11; margin_range, error naming and pooled groupBy fix added the same day); tests: `tests/testthat/test-clean_gate.R`. Every impactful value is a required argument (10); a call with missing ones stops once with a message that lists all of them (meaning, how to choose, benchmark 05 value). See the register. |

## Integration audit (2026-11-08)

Checked: dead code, overlap, and wiring between old and new modules.

| Finding | Decision |
|---|---|
| No dead helpers; every export is defined. Four `sw_pipeline_step_*` constructors (som, correct_apply, annotate, diff) had no caller and no test. | `sw_pipeline_step_annotate()` and `sw_pipeline_step_diff()` were removed (unused, 2026-11-09). som and correct_apply stay: add tests, or cut them if no pipeline needs them. |
| `sw_io_flowstate_transform()` vs `sw_transform_estimate()`/`sw_transform_apply()`; `sw_io_flowstate_to_tibble()` vs `sw_io_ff_to_tibble()`. | Keep: they act on a different class (flowstate vs flowFrame) and flowstate transforms in place. Docs point to the flowFrame equivalents. |
| `sw_unmix_quality()` (CV per channel) and `sw_unmix_spillover_matrix()` (cosine similarity of reference spectra) vs the new `sw_unmix_assess_*()`. | Keep. Different questions. `sw_unmix_spillover_matrix()` is not a spillover spreading matrix from unmixed controls; documented 2026-11-11 (name kept, roxygen says cosine similarity). |
| `sw_correct_evaluate_*()` output cannot feed `sw_batch_decide()`. | A2c. |
| `sw_batch_decide()`, `sw_flag_iqr()`, `sw_design_audit()` are not composable-pipeline steps and are unknown to the assistant tool list (`R/assistant.R`). | Add step constructors and tool entries after A3, not before. |
| `sw_pipeline_run_all()` hands PeacoQC and gating untransformed frames, never calls `sw_transform_*`, and runs no audit. Its defaults (cofactor 6000, seed 42, SOM 8×8, 20 metaclusters, IT limit 0.55, MAD 6) are impactful choices without register rows. The user confirmed cofactor 6000, n_metaclusters 20, seed 42, gate_node 'lymphocytes' and unmix_method 'AutoSpectral' as they are (2026-11-09). `sw_correct_prepare()` drops markers missing in a sample without a message. In `sw_correct_evaluate_quick()` output, `emd_reduction_pct` is the MAD reduction. | Add register rows or grandfather them in writing; fix the label; stop the silent drop. These pre-date this roadmap. |
| `sw_design_audit()` eligibility table maps to `sw_correct_detect_batch`, `sw_flag_iqr`, `sw_correct_run`. It now has a `paired_shift` row (C1). | Done. |

Fixed in this pass: `sw_design_audit()` now needs `minimum_batches_origin`, type-checks its floor arguments, rejects `reference_role` without `reference_id`, warns on an all-NA donor column, checks nesting without a rank warning, and stays quiet about disconnected references when none exist.

## Suggested order

Batch-side steps 1 to 3 are frozen and B3 to B8 are parked (see the decisions at the top). The list stays as the plan for when they return.

1. C0 and C1: manifest, corrector contract, population model, extended audit.
2. A2b, A2c, A3: protected-contrast gate, metric bridge, gating on the audit.
3. C2 to C4: paired-shift engine with a frozen model, checked against the bundle's outputs. Pull C5's export forward; it is the cheapest oracle.
4. D1, D2, B1, B2, D5, D6: independent of the above; run alongside.
5. B3 to B8 on synthetic fixtures.
6. C6, C7, B9 once their gates are met.

### B11: unmix diagnostics suite (from benchmark 05)

| # | Item | Function(s) | Status |
|---|---|---|---|
| B11 | Spike-in recovery with an independence check, negative shape and marker screen, scatter and spike figures (F6 to F8), Nature theme for all unmix figures, as diagnostics on AutoSpectral outputs (from benchmark 05). Unmixing stays with the user. Limits unsourced: `not_judged`. Plan: `plans/unmix_diagnostics_plan_v2.md`; decisions: `plans/unmix_diagnostics_decisions.md`. | `sw_unmix_spike_*()`, `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()`, `sw_unmix_figure_scatter_grid()`, `sw_unmix_figure_spike_*()` | **done** (2026-11-11). All four PRs are merged into main: theme and F1 to F4, spike build and score, screen and F6 to F8, diagnose, report, dashboard and quarto. `baseline_tol` and `min_negatives`: caller chooses; no recommended value (benchmark 05 used 0.05 and 1 as test values only). |

### Reference and batch diagnostics

| # | Item | Function(s) | Status |
|---|---|---|---|
| P2a | Planted shifts for synthetic references (DR24, DR37): four families, truth table apart from the data, invertible. Plan: `plans/2026-11-11_reference-and-batch-diagnostics.md`. | `sw_synref_shift()` | **done** on branch `feat/p2a-synref-shift`, not merged. |

Module 1 is not a rebuild of AutoSpectral (rule 3 of "Gold standard"): AutoSpectral 1.8.8 has `sim.flow.data()` for simulated events but no spike into real cells, so spike build and score are new code that calls no AutoSpectral routine.

## Register of choices the user must set

Each item below changes a result, so none has a package default.

| Choice | Where | Source of the starting value |
|---|---|---|
| Shift family, shift size per target, cofactor, seed | `sw_synref_shift()` (P2a) | user, DR24 and DR37, 2026-11-12; no values chosen yet |
| Minimum number of batches; require a shared reference | `sw_design_audit()` | batchdx uses 3 and `TRUE`; "configuration default" |
| Decision gates: minimum coverage, EMD reduction, MAD compression, joint-structure change, no-effect EMD | `sw_batch_decide()` | batchdx: 0.6, 0.2, 0.25, 0.1, 0.05; all "uncalibrated engineering default" |
| IQR multiplier; positivity thresholds; minimum cells per control side | `sw_flag_iqr()` | CytoBatchFlagR describes 1.5; batchdx adds 100 cells; thresholds are analyst-set |
| Candidate grid, criterion limits and scales, near-best tolerance, event minima | `sw_study_*()` | Bundle examples are synthetic and illustrative only |
| Audit mode and risk limit; number of independent blocks | `sw_study_audit()` | Bundle demo: 24 blocks |
| Reference-shift options: anchor batch, minimum cells per reference summary | `sw_corrector_fit()` | cyTether scaffold: `min_cells = 3` |
| Seeds and replicate counts for benchmarks | C5, C6 | Bundle: at least 100 seeds for publication-scale claims |
| Origin of the batch-count floor | `sw_design_audit()` | Required string; stored in the result |
| Population labels (`population`), `min_events` per file x population, `group_summary` (`median` or `mean`) | `sw_batch_scores()` (P3b), `sw_batch_interaction()` (`which_check`, `min_populations`) | None yet; the benchmark label source (gate or k-means k) is open until the user answers |
| Minimum effect size for a control flag | `sw_flag_iqr()` (A4 follow-up) | Argument added (P0c); values per metric open (DR9, null-arm run, user approval) |
| Limit on donor-level protected-contrast change | `sw_batch_decide()` (A2b) | None yet |
| Reference checks: `min_events`, `min_batches_per_reference` (plan proposes 3), `n_events`, `n_splits` (plan proposes 20), `n_quantiles`, `roles`, `seed`, `order` | `sw_ref_inventory()`, `sw_ref_consistency()`, `sw_ref_noise_floor()` (C9) | Plan proposals; benchmark values per the 2026-11-12 plan: `min_events = 500`, `n_splits = 20`, `min_batches_per_reference = 3`. Reasoning, not sourced numbers |
| `n_events`, `seed`, reference roles, minimum donors per batch, optional noise floor | `sw_ref_vs_study()` (P3a, plan DR17, DR18) | Plan proposes 20,000 events (reported, not verified) and 4 donors (dataset rule D5, source not opened); both stay uncalibrated, the user sets them |
| Minimum events per population, label source | `sw_ref_composition()` (P3a, DR19) | None; floor is the analytic multinomial standard error |
| Quantile probabilities, minimum events per population | `sw_ref_coverage()` (P3a, DR20) | Plan proposes 0.01 and 0.99 (weak source: CytoBatchNorm percentile range, as reported by a survey) |
| Anchor batch; hierarchy (study, run, instrument); population model | `sw_corrector_fit()` | The bundle defaults to the first sorted batch, which this package will not copy |
| Number of audit blocks | `sw_study_audit()` | Binomial bound, about 17 at a 0.20 limit with zero failures |
| Cofactor, seed, SOM grid, metaclusters, IT limit, MAD limit in `sw_pipeline_run_all()` | existing pipeline | Pre-existing defaults (6000, 42, 8×8, 20, 0.55, 6). Confirmed by the user 2026-11-09: cofactor 6000, seed 42, n_metaclusters 20, gate_node 'lymphocytes', unmix_method 'AutoSpectral'. SOM grid, IT limit and MAD limit still need a confirm. |
| Verdict limits for unmix diagnostics | `sw_unmix_limits_autospectral()` | AutoSpectral article 01: `max_sif` above 6 most likely causes problems, 4 to 6 is worth reviewing (Mage, Konecny, Mair 2026, Cytometry A 109(6):436-455, doi:10.1002/cyto.a.70044; preprint doi 10.1101/2025.04.17.649396). Mage's own diagonal cutoffs are 3 and 8, rough guidelines. User chose warn 4, fail 6 on 2026-11-09. Condition number, max cosine and negative spread have no sourced cutoff and stay not judged. |
| `margin_range`: Lower and upper limit (raw scale) of the margin step, applied to every fluorescence channel. PeacoQC removes an event at or beyond a limit in any fluorescence channel, and the event at the observed minimum or maximum of a channel. Scatter channels keep the range stored in the FCS file. PeacoQC reads a lower limit above 0 as 0; `sw_clean_run()` warns and continues. | `sw_clean_run()` (D11) | caller chooses; no recommended value. How to choose: look at the raw fluorescence densities and at the share of events the margin step flags per channel. Measurements: `experiments/runs/2026-11-11_clean-run-smoke/warnings.txt` (run folder not tracked by git) lists PeacoQC margin warnings for A3 99.34, A6 99.81, B11 99.48, B2 99.16, D1 99.86 percent; benchmark 05 passed c(-21761, 1e7) (`12_cleanup.R:28`, user 2026-11-09); margin_range required (user, 2026-11-11 plan grill). |
| `debris_window`: Lower and upper FSC-A limit (raw scale) of the events the nonDebris density is built from. openCyto drops events outside [min, max] before it looks for peaks and valleys, so the cut is a valley of the density of the events inside the window and lies inside it. Events below the cut are debris. Events above max are kept. If the window holds one peak or no valley, the cut falls on the lowest event in the window. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used c(2e5, 1.5e6) with mindensity on 40-colour Aurora PBMC, for reference only. How to choose: Plot raw FSC-A as a density and FSC-A against SSC-A (ggplot2::geom_hex() on sw_io_ff_to_tibble(ff)). Debris is the low-FSC cloud beside the main cells. Look at one sample per batch, because the valley can move. The window must hold the valley you want. Use Inf as the upper bound to keep every event above min in the density. |
| `singlet_wider_gate`: Width of the openCyto singletGate on FSC-A against FSC-H. TRUE widens it and keeps more events near the edge. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used TRUE on 40-colour Aurora PBMC, for reference only. How to choose: Plot FSC-A against FSC-H of non-debris events (raw scale) and look at the width of the diagonal. Run both settings on a trial sample and compare the singlets rows of the audit table. |
| `pool_samples`: TRUE fits each gate once on the pooled events of all samples and applies it to every sample. FALSE fits one gate per sample. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used TRUE (gates pooled over donors) on 40-colour Aurora PBMC, for reference only. How to choose: Plot per-sample densities of FSC-A and of the viability channel on one panel. If the valleys line up across donors and batches, pooling is defensible. If they shift by batch, do not pool. Compare the $gates table of a trial run with each setting. |
| `viability_channel`: Channel of the viability dye. Dead cells are the high tail. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used "LIVE DEAD NIR-A" on 40-colour Aurora PBMC, for reference only. How to choose: List channels with sw_channel_get_fluor(ff) and check the panel sheet for the amine-reactive dye. Plot its density after arcsinh and look for a dead peak or tail. |
| `viability_method`: How the live cut is found: "mindensity" takes the valley of the density, "boundary" takes the fixed cut in viability_cut. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used "mindensity" on 40-colour Aurora PBMC, for reference only. How to choose: Plot the density of the arcsinh-transformed viability channel per sample (sw_gate_density_cut() shows one valley cut per vector). Use "mindensity" when every sample shows two peaks. Use "boundary" when a sample has no valley. |
| `viability_cut`: Live cut on the arcsinh scale for viability_method = "boundary" (events below it are live). Must be NULL for "mindensity". | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used NULL (mindensity finds the valley) on 40-colour Aurora PBMC, for reference only. How to choose: For "boundary", read the cut from the same density plot as viability_method, or from an unstained or fully live control. For "mindensity", pass NULL on purpose. |
| `IT_limit`: Isolation Tree limit of PeacoQC. It trades stringency against the number of events kept. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used 0.55 on 40-colour Aurora PBMC, for reference only. How to choose: Run sw_qc_run() on one trial sample per batch with a few values, read pct_removed, and inspect the QC plot with sw_qc_plot(). Summarise several samples with sw_qc_summary() on sw_qc_batch() output. |
| `MAD`: Median Absolute Deviation limit of PeacoQC. It trades stringency against the number of events kept. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used 6 on 40-colour Aurora PBMC, for reference only. How to choose: Same trial runs as IT_limit: vary it with sw_qc_run(), read pct_removed, inspect with sw_qc_plot(), compare samples with sw_qc_summary(). |
| `cofactor`: Divisor of the arcsinh transformation of all fluorescence channels, asinh(x / cofactor). Scatter stays raw. It shifts the QC and viability cuts. | `sw_clean_run()` (D11) | caller chooses; no recommended value; benchmark 05 used 6000 on 40-colour Aurora PBMC, for reference only. How to choose: Transform a few channels, including the viability channel, with several cofactors (flowCore::arcsinhTransform(a = 0, b = 1 / cofactor)) and plot their densities. Pick the value where the negative peak is compact and the positive peak is not squeezed. |

| Spike levels; cofactor; spiked fraction; seed; populations | `sw_unmix_spike_build()` (B11) | user, 2026-11-11, grill R3, R4, R5. Benchmark 05 used levels 0 to 3, cofactor 6000, fraction 0.75 (implied by levels 0 to 3), seed 42, CD4T, CD8T and CD3- |
| Loss cut; baseline tolerance | `sw_unmix_spike_score()` (B11) | user, 2026-11-11, grill R7, R2. Benchmark 05 used loss cut 0.5; `baseline_tol`: caller chooses; no recommended value (benchmark 05 used 0.05 and 1 as test values only) |
| Co-spike level, second dye and its levels | `sw_unmix_spike_cospike()` (B11) | user, 2026-11-11, grill R2.9, R8. Benchmark 05 used level 1, GLUT1 at 0, 1, 2 |
| Positive cut, squeeze fraction, minimum negatives, reference setup | `sw_unmix_negative_shape()` (B11) | user, 2026-11-11, grill R9, R10. Benchmark 05 used arcsinh 1 and `as188_joint`; `min_negatives`: caller chooses; no recommended value (benchmark 05 used 0.05 and 1 as test values only) |
| Minimum median rise; minimum donors | `sw_unmix_marker_screen()` (B11) | user, 2026-11-11, grill R11. Benchmark 05 used 0.10 and 5 |
| Verdict limits for spike and screen checks | `sw_unmix_diagnose()` (B11) | user, 2026-11-11, grill Q16: none sourced, so `not_judged` unless the caller passes a limit |
| `n_events` and `seed`: events drawn per file and the seed of the draw, for figures F3 and F4 | `sw_unmix_sample_events()` (B11) | caller chooses; no recommended value (benchmark 01, Mendeley 7-colour, used 50,000 events and seed 1 as test values only; `experiments/benchmarks/01_mendeley-7color/README.md`, user, 2026-11-09). Both are required arguments with no default. |

Note, 2026-11-09: `sw_unmix_assess_mismatch()`, `sw_unmix_assess_mismatch_clusters()`, `sw_unmix_assess_variability()`, `sw_unmix_assess_alignment()`, `sw_unmix_compare()`, `sw_unmix_compare_data()` and `sw_unmix_compare_af()` became internal (`@keywords internal`, no export) because no run used them. Code, help text and tests are kept.

An earlier stage's choice is not settled by its code. When a later step inherits one, it gets its own row here and a one-line confirm from the user.

## Pre-existing defaults on impactful choices, unconfirmed

These defaults pre-date this roadmap. Each can change a result. The user confirmed only cofactor 6000, n_metaclusters 20, seed 42, gate_node 'lymphocytes' and unmix_method 'AutoSpectral' for `sw_pipeline_run_all()` (2026-11-09). Nothing in this table is changed here. Values were read from the function signatures in `R/`.

| Function | Default |
|---|---|
| `sw_correct_run()`, `sw_correct_som()` | `xdim = 8`, `ydim = 8`, `rlen = 10`, `seed = 42` |
| `sw_correct_detect_batch()` | `seed = 42` |
| `sw_plot_batch_dimred()` | `seed = 42` |
| `sw_correct_prepare()` | `cofactor = 6000` |
| `sw_gate_run()` | `cofactor = 6000` |
| `sw_correct_emd()`, `sw_correct_evaluate_emd()` | `binSize = 0.1` |
| `sw_cluster_run()` | `xdim = 10`, `ydim = 10`, `rlen = 10`, `n_metaclusters = 20`, `k = 30`, `seed = 42` |
| `sw_corrector_compare()` | `n_quantiles = 100` |
| `sw_annotate_run()` | `min_score = 0.3` |
| `sw_diff_abundance()` | `min_prop = 0.01`, `fdr_threshold = 0.05` |
| `sw_diff_expression()` | `min_cells = 3L`, `fdr_threshold = 0.05` |
| `sw_plot_volcano()` | `fdr_threshold = 0.05`, `fc_threshold = 0.5` |
| `sw_dimred_run()` | `n_neighbors = 15L`, `min_dist = 0.2`, `max_cells = 50000L`, `seed = 42L` |
| `sw_gate_singlets()` | `nmad = 4` |
| `sw_qc_run()` | `IT_limit = 0.55`, `MAD = 6` |
| `sw_qc_summary()` | `threshold = 30` |
| `sw_unmix_quality()` | `cv_threshold = 0.5` |
