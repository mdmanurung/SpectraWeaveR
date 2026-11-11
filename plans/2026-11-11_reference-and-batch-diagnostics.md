# Plan: batch-effect diagnostics, reference-sample consistency and synthetic references (2026-11-11)

Status: **draft, 24 of 35 decision rows settled (DR1 to DR23 and DR25); open: DR24, DR26 to DR33, DR37, DR38.** No code written. Nothing is built before the user approves the plan.

Plot-first note: no figure was opened for this plan. Every number quoted from plans, handoffs or papers is a reported statement and carries `unsupported: no inspected figure`. Each benchmark below must save and inspect its figure before any conclusion cites it.

## Goal

User request (2026-11-11): plan all functions for (1) batch-effect diagnostics, (2) checking whether the batch reference samples are consistent across batches, and (3) utilities that build synthetic reference samples for batch correction, for example by downsampling real study cells.

Three limits shape every row:

- Impactful choices are required arguments with no default. Each needs a plan row the user approves (decision rows `DR1` to `DR33`).
- Diagnostics rank per check and marker. No function returns a per-file, per-reference or per-batch pass or fail (Q14).
- AGPL and GPL code is wrapped, never copied. Every new function here is own code or a clean rewrite from a public definition.

## What exists

Opened by me this step: `CLAUDE.md`, `ROADMAP.md` (lines 1 to 30), `R/manifest.R` (lines 118 to 130), `R/utils.R` (lines 390 to 404), `R/corrector.R` (lines 1320 to 1420, plus a grep), `R/design_audit.R` (lines 18 to 28, 262), `plans/2026-11-10_diagnostics-next.md` (lines 1 to 150), `experiments/README.md`, `experiments/datasets/registry.csv`, `experiments/benchmarks/04_mendeley-correction/README.md`. Everything else below is from the surveys and is marked "survey".

Batch diagnostics (survey, signatures re-checked where noted):

- `sw_batch_scores()` gives pair_emd, rest_emd, rest_emd_z, median_range, eta2, spread_ratio and order_rho per marker and batch. `data, batch, markers, n_events, seed, order, reference_group` are required. It has no population argument and no role argument.
- `sw_flag_iqr()` is a Tukey fence per batch, marker and metric. It has no minimum effect (ROADMAP A4).
- `sw_otd_scores()` scores files against a pooled reference. It does not compare batches.
- `sw_cluster_scores()` already returns per-file cluster shares as robust z, with `reference_group`.
- `sw_diagnose()` stacks nine modules. It has no reference, flag_iqr or design module. `sw_diagnosis_top()` ranks per module, check and item.
- No function turns batch scores into a batch-by-marker table. No plot function exists.

Reference consistency (survey):

- `sw_design_audit()` only checks that a reference id exists in two batches. It never reads event data.
- `sw_corrector_compare()` compares each shared reference id between batch pairs by quantile distance. I read its signature (`n_quantiles = 100` is a default) and grepped `R/corrector.R` for hold-out text. I found none, but I did not read the body at lines 1420 to 1550, so "no hold-out" is unverified.
- `sw_corrector_diagnose()` gives `se_donor`, which is NA with one reference donor.
- The cyTether bundle logs an "aliquot noise floor" as a diagnostic. Whether the repo implements it was not checked.

Synthetic references: none exist as package code.

- `sw_io_subsample()` has `seed = NULL` and `setup = "forceNEvent"` as defaults (opened by me). It returns one pooled `flowFrame` without per-sample metadata.
- Allowed `reference_role` values are `none`, `technical_replicate`, `anchor`, `biological_pseudo_control` (opened by me, `R/manifest.R` 124 to 126, also `R/design_audit.R` 262). The audit text says a `biological_pseudo_control` "never counts as a measured technical reference" (opened, line 23).
- The spike builder `sw_unmix_spike_build()` is the closest pattern: assignment table, realised fraction, stop on overlap, design written to FCS keywords (survey).
- Fault injectors exist only as script-local functions in `experiments/` (survey).

Data (survey, plus registry opened by me):

- CytoBatchFlagR extract (Zenodo 15388817): 4 control samples across 18 batches, one Aurora. It is the only real set with the same material in every batch.
- Aurora ex vivo extract: 149 files in 4 batches (one document says 160), one `BMV_REF` file per batch. Whether the four `BMV_REF` files are the same material is unverified (survey). Donor is nested in batch.
- Mendeley 7-colour file (`autospectral_mendeley` in the registry): one file, no batches. Benchmark 04 uses it for unmixing correction (opened).
- FR-FCM-Z247 stays sealed as the only out-of-sample check (survey).

Freeze status needs one user line (`DR1`). `ROADMAP.md` line 7 (opened) freezes A1 to A6 and C0 to C5 until a real multi-batch dataset exists. The survey reports that the 2026-11-09 goal update and row G32 lift the freeze for diagnostics, but `ROADMAP.md` was never edited to say so. This plan assumes the diagnostics goal covers the new functions and keeps A2b, A2c and A3 frozen unless the user says otherwise.

## Design rules

1. Reuse first. Batch and reference scores call the internal `.sw_batch_scores_core` and `sw_emd_1d`. Synthetic builders reuse the RNG-restore pattern of `sw_io_subsample()` and the assignment-table pattern of the spike builder.
2. New arguments never get a default that changes a number. NULL is allowed only as an explicit "skip", as in `sw_diagnose()`.
3. Output is a long tibble per check, marker and, where relevant, batch, reference id or population. Rankings use `sw_diagnosis_top()` or `sw_batch_matrix()`. A "met" column exists only in the benchmark scorecard, per check and fault class.
4. A synthetic reference is a diagnostic comparator, null arm, pseudo-batch base or benchmark truth. It is never a measured technical reference and never feeds the `paired_shift` or HEB bridge path (`DR3`, `DR4`).
5. Cells are not replicates. Donors and aliquots are (Q11). Every split, draw and floor function names its unit, and cell-level floors are labelled "sampling only".
6. Every random draw takes a seed. Every synthetic object stores provenance as source rows, so a draw rebuilds without storing cells.
7. No plot function is built before the user approves its plan row (`DR30`).

## Function inventory

Naming: `sw_batch_*` for batch diagnostics, `sw_ref_*` for reference checks on measured references, `sw_synref_*` for synthetic references. Required choices cite decision rows. No values appear in these tables.

### A. Changes to existing exports

| Name | Purpose | Inputs | Outputs | Required choices | Builds on | Source |
|---|---|---|---|---|---|---|
| `sw_io_subsample` (change) | Remove an unplanned impactful default | existing, with `seed` and `setup` now required | unchanged pooled `flowFrame` | setup, seed (`DR10`) | itself; `tests/testthat/test-utils.R` | `R/utils.R:397-403` (opened) |
| `sw_sample_manifest`, `sw_design_audit` (change) | Record a synthetic reference without counting it as measured | `reference_role` accepts `synthetic_reference`; a `source_batch` note | same classes; audit lists synthetic references apart from measured ones | role name, source_batch required or not (`DR3`) | `.sw_manifest_check_nesting`, `.sw_reference_components`, `.sw_anchor_bridge` | `R/manifest.R:124-126`, `R/design_audit.R:23, 262` (opened) |
| `sw_flag_iqr` (change) | Close ROADMAP A4: stop flagging noise | existing plus `min_effect` per metric | same tibble plus `effect`; status `below_min_effect` | min_effect per metric and its unit (`DR9`) | `.sw_flag_metrics`, `.sw_flag_panel` | ROADMAP A4 (survey); CytoBatchFlagR over-flags when spread is tiny (survey, full text); clean rewrite, tool is GPL-3 |
| `sw_batch_scores` (change, optional) | Role split from plan row D7 | existing plus `role` | same tibble plus `role` | approve or decline (`DR6`) | `.sw_batch_scores_core`, `reference_group` | `plans/2026-11-10_diagnostics-next.md` lines 65, 84 (opened): D7 still "needs user approval" |
| `sw_corrector_compare` (change) | Make held-out scoring explicit; drop an impactful default | existing, `n_quantiles` required, plus a table of reference samples excluded from scoring | same long metric table plus `n_excluded` | n_quantiles, exclusion unit (`DR11`) | `emd_within_population` | signature opened; body at 1420-1550 not read, read it in PR 0 |
| `sw_diagnose` (change) | Stack the new outputs | existing nine modules plus `reference`, `batch_pop`, `flag_iqr`, each result or explicit NULL | same `sw_diagnosis` table, new modules | which modules run (`DR33`) | `.sw_diag_key`, `sw_diagnosis_top` | `R/diagnose.R` (survey) |

### B. Batch-effect diagnostics (new)

| Name | Purpose | Inputs | Outputs | Required choices | Builds on | Source |
|---|---|---|---|---|---|---|
| `sw_batch_scores_by_population` | Is the effect marker-specific or population-specific: the `sw_batch_scores` checks inside each population | data, labels (per-event, per file), batch, markers, n_events, seed, order, reference_group, min_events_per_population, label_source (gate or clustering), checks | long tibble: check, population, marker, batch, value, n_events, n_files; NA when too few events | n_events, seed, min events, label source, checks (`DR7`, `DR8`, `DR12`) | `.sw_batch_scores_core`, `sw_emd_1d` | own code; cyCombine and CytoNorm compute EMD per population (survey, full text) |
| `sw_batch_matrix` | Ranked batch-by-marker view of any batch score or flag table; the missing table | scores (from `sw_batch_scores`, `_by_population`, `sw_ref_consistency`, `sw_flag_iqr`), check, direction, rank scope | long tibble: check, batch, marker, value, rank; wide view as attribute | check, direction per check, rank scope (`DR2`, `DR38`) | `sw_diagnosis_top` ranking | own design; gap (a) in survey |
| `sw_batch_interaction` | Does a marker's batch shift differ between populations; flags a misspecified additive per-population model | per-population scores, batch, markers, min_populations, scale unit, contrast base | tibble: marker, batch, shift_pooled, shift_sd_across_populations, ratio, n_populations | min_populations, scale unit, contrast base (`DR7`) | `sw_batch_scores_by_population` | own design, no literature source; ROADMAP C0 and C2 note on additive shifts (survey) |
| `sw_batch_otd` | Joint multivariate batch distance, pairwise or batch against the rest | data, batch, markers, n_events per batch, blur, seed, max_iter, tol, mode | tibble: batch (pair), otd, converged | all (`DR12`; phase 3) | `sw_sinkhorn_divergence`, `.sw_ot_entropic` | Glehr et al. 2026 definitions (survey); clean rewrite exists, code licence unknown |
| `sw_correct_metrics` | ROADMAP A2c: bridge from `sw_correct_*` output to `sw_batch_candidates_from_metrics()` | before, after, batch, markers, candidate_id, labels, n_quantiles, bin size, min cells | long tibble: candidate_id, metric, marker, comparison_id, value, status, weight | bin size, labels, weights, identity_id; gated by `DR1` | `sw_correct_emd`, `sw_correct_evaluate_mad`, `sw_batch_candidates_from_metrics` | ROADMAP A2c (survey); cyCombine EMD and MAD (survey, MIT + file LICENSE) |
| `sw_plot_batch_matrix` | Dot-matrix figure of one batch-by-marker check, source data returned | `sw_batch_matrix` result, check, order, scale, n_markers | list(plot, source_data); needs ggplot2 (Suggests) | all; needs an approved plot row (`DR30`) | figure pattern of `sw_unmix_figure_*` | plot-first rule, ROADMAP rule 5 |

### C. Reference consistency (new, measured references)

These answer "is the reference the same material in every batch?" with rows per reference id, marker and batch. None returns a verdict per reference.

| Name | Purpose | Inputs | Outputs | Required choices | Builds on | Source |
|---|---|---|---|---|---|---|
| `sw_ref_inventory` | List each reference by batch with event counts, markers present, populations present | data, manifest, markers, labels (optional), min_events | tibble: reference_id, batch, sample_id, n_events, marker, present, population, n_population_events, under_min | markers, min_events, label source (`DR15`) | `sw_sample_manifest`, `sw_design_audit` reference graph, `sw_io_event_audit` | own code; gap: the audit never reads events (survey) |
| `sw_ref_consistency` | Per reference id, marker and batch: pair_emd, rest_emd, rest_emd_z, median shift z, order_rho, event-count ratio | data, manifest, design_audit, markers, n_events, seed, order, min_batches_per_reference, min_events, roles | long tibble class `sw_ref_consistency` with status `not_assessed` for thin rows | markers, n_events, seed, order, min batches, min events, roles, transform and cofactor (`DR12`, `DR15`, `DR37`) | `.sw_batch_scores_core` with `reference_group = reference_id` and per-group rows kept; run `sw_flag_iqr` on the same files for the Tukey view | own composition; EMD per marker per batch pair from CytoNorm and cyCombine, median pairwise EMD from CytoBatchFlagR (survey, full text, GPL-3 so rewrite only) |
| `sw_ref_noise_floor` | Per marker, how different two draws of the same reference are by chance, so between-batch EMD can be read against it | data, manifest, markers, n_events, n_splits, seed, n_quantiles, split unit | tibble: reference_id, marker, floor_type (sampling or aliquot), value | n_events, n_splits, seed, split unit, quantile (`DR16`) | `sw_emd_1d`, `sw_synref_split` | cyTether "aliquot noise floor" (survey); AutoSpectral bundle: more events do not add independent runs (survey) |
| `sw_ref_vs_study` | Does the reference move between batches by less than donors differ inside one batch | data, manifest, markers, n_events, seed, spread scope, min_donors_per_batch, reference roles, optional noise floor | tibble: marker, reference_id, batch_a, batch_b, ref_emd, donor_emd_median, floor_emd, ratio | spread scope, min donors, comparator statistic (`DR17`, `DR18`) | `sw_emd_1d`, manifest donor and reference columns | `experiments/integration/aurora_exvivo/03_followups.R` header (survey); own design, no literature source |
| `sw_ref_composition` | Population shares of one reference across batches with a noise floor from event counts | data, manifest, labels, min_population_events | tibble: reference_id, population, batch, share, floor, cv_across_batches | labels, min events, noise model (`DR19`) | `sw_cluster_scores` when k-means labels suffice | Ivison 2018, CytofRUV Hellinger idea (survey, text read; GPL-3 so definitions only) |
| `sw_ref_coverage` | Does the reference span the study range per marker and population | reference data, study data, batch, markers, probs, labels, min_events | tibble: batch, marker, population, share_below, share_above, n_reference_cells | probs, markers, labels (`DR20`) | `sw_ref_inventory` | CytoNorm: references must span the full expression range (survey, full text) |
| `sw_ref_split` | Fit and evaluation split of reference samples; unit is donor or aliquot, never cell | manifest, unit, scheme, n_eval, seed | tibble: fold, sample_id, reference_id, role (fit or eval); digest | unit, scheme, n_eval, seed (`DR11`) | spike builder assignment-table pattern | gap (d) in survey; cyTether Q1 held-out donor endpoint (survey) |
| `sw_ref_influence` | Leave one reference out and refit; show which reference drives a correction | all `sw_corrector_fit` required arguments plus `drop_ids` and `drop_unit` | tibble: dropped_id, batch, population, marker, effect_full, effect_dropped, delta, status | engine, anchor, min_cells, seed, drop unit (`DR28`; phase 3) | `sw_corrector_fit`, `sw_corrector_diagnose` | gap (f) in survey; own design |

### D. Synthetic references (new)

`sw_synref` is an S3 class: `events` (named list of matrices), `provenance` (synref_id, draw, source_sample, source_row, batch, donor, population, stratum), `recipe` (every argument), `digest`, `scope` (`within_batch` or `cross_batch`), `certainty = "diagnostic_only"`. Draw functions stop when a stratum is short. They never filter silently.

| Name | Purpose | Inputs | Outputs | Required choices | Builds on | Source |
|---|---|---|---|---|---|---|
| `sw_synref_downsample` | Core builder: draw events from real study cells into one reference per group. Covers a random draw, a per-batch aggregate (equal events per file), a donor-balanced pool and a population-balanced draw through `groups`, `allocation` and `strata` | data, manifest, groups (which files feed which reference), source_samples, allocation, strata, labels, n_events, replace, short_stratum, n_draws, seed | `sw_synref` with realised counts per stratum | all of these (`DR21`, `DR22`, `DR23`) | RNG-restore pattern of `sw_io_subsample`, spike builder realised-fraction report | precedents: CytoBatchFlagR pooled control-free samples, CytoNorm without-controls vignette, iMUBAC, Glehr pooled reference (survey, as reported); recipe choices are own design |
| `sw_synref_mix` | Mix real populations at set proportions; report requested against realised | data, labels, proportions, n_events, source_samples, replace, short_stratum, seed | `sw_synref` plus requested and realised proportions | all (`DR27`) | `sw_synref_downsample`, spike builder disjoint-population rule | gap (c) in survey; cyTether rare-population failure row (survey) |
| `sw_synref_split` | Cut one real sample into disjoint pseudo-aliquots or pseudo-batches with identical biology | one source, number of parts, events per part, replace, unit, seed | `sw_synref` labelled `floor = "sampling only"` | parts, events per part, replace, seed (`DR21`) | `sw_synref_downsample` | own design; cells are not replicates (Q11, survey) |
| `sw_synref_shift` | Plant a known batch effect into a synthetic or real sample; returns the truth table apart from the data | x, target batches, family, params table (batch, marker, population, size), cofactor, order, seed | list(data, truth) | family, sizes, cofactor, targets, seed (`DR24`) | `sw_io_ff_to_tibble`, `sw_transform_apply` | ROADMAP note: additive shifts on arcsinh act as a scale change for dim and negative values (survey); fault injectors in `experiments/` (survey) |
| `sw_synref_goal` | Quantile goal table across batches for quantile alignment; role "goal", never a reference | data, batch, markers, quantile set, goal type, n_events per batch, seed | `sw_synref_goal` table with per-batch quantiles and spread | all (`DR26`; phase 3) | `sw_anomaly_features` aggregate pattern | CytoNorm goal options (survey; README, vignette and paper disagree on mean or median) |
| `sw_synref_spread` | Repeat a recipe and report how far each marker quantile moves between draws | recipe, data, n_draws, seeds, markers, quantiles | tibble: marker, quantile, sd_across_draws, range_across_draws, n_draws | n_draws, seeds, quantiles (`DR22`) | `sw_synref_downsample` | no source gives a replicate count (survey); own design |
| `sw_synref_fidelity` | What a draw changed against its source: marker EMD to source, population share change | synref, data, markers, labels, n_events, seed | tibble: marker, emd_to_source, population, share_source, share_synref, share_ratio | n_events, seed, labels | `sw_emd_1d` | own code |
| `sw_synref_check_overlap` | Stop when a synthetic reference shares cells, files or donors with data it calibrates or is scored on | synref, calibrate_ids, eval ids, overlap unit, manifest, max overlap share | tibble: synref_id, eval_file, unit, n_shared, share, status; error on overlap | unit, max overlap share (`DR25`) | provenance table, spike builder independence check | AutoSpectral bundle: injections from one parent stay in one dependence group (survey) |
| `sw_synref_as_manifest` | Rows for `sw_sample_manifest` with role `synthetic_reference` and a parent note | synref, manifest, batch assignment | manifest rows | role name, batch assignment (`DR3`) | `sw_sample_manifest` | `R/manifest.R` (opened) |
| `sw_synref_save`, `sw_synref_load` | Frozen object with digest; optional FCS with provenance keywords | synref, run_dir (must not exist), format | files in a new dir; load fails on digest mismatch | run_dir, format | `sw_signal_reference_save/load` pattern, spike builder keywords | survey |

New exports: 4 batch functions, 1 bridge, 1 plot proposal, 8 reference functions, 11 synthetic functions (25), against 136 today.

## Validation strategy

Two layers, as the user asked: simulated ground truth, then real data.

### 1. Unit tests (mocks, no real data, `tests/testthat/`)

Each function gets a truth check that needs no flowCore:

- Draw functions: realised count per stratum equals the request; every drawn row is a source row; the same seed gives the same draw; a short stratum stops or caps exactly as asked; planted overlaps are found by `sw_synref_check_overlap`.
- `sw_synref_shift`: applying then inverting the shift returns the input within tolerance; the truth table matches the mean difference of shifted and unshifted cells.
- Reference and batch checks: on tiny simulated matrices, a planted shift in one batch of one reference ranks first in its marker column; identical pseudo-aliquots stay near the noise floor; thin rows get `not_assessed`.
- Role change: `synthetic_reference` never creates an anchor bridge in `sw_design_audit()`; existing eligibility rows are unchanged (regression test).

### 2. Simulated ground truth (`experiments/benchmarks/12_reference-checks`, `13_synref-truth`)

Generators (population mixtures on arcsinh, cofactor 6000, with dim and negative tails) stay script-local in `experiments/benchmarks/_common/`, as `make_file()` does in benchmark 07 (`DR5`). Scenarios:

- Reference faults: marker drift in one batch of one reference, population removed or rescaled, material swap, event loss, dim-marker degradation.
- Batch faults: pure marker shift, pure population shift, both, in each of the four shift families, plus drift over run order.
- Null arms with the same draws and no injection, rerun at more seeds where the earlier null arms were (G40, G41, survey).
- Sizes include the smallest size already tested on real files (D2 answer, opened).

Scoring: AUROC per check, fault class and size, from rows ranked within check and marker, never from a per-file summary. The bar is the D1 answer (opened): mean AUROC at least 0.90 over at least 5 seeds, lower 95% CI at least 0.80, shuffled-label null at most 0.55 (`DR14`). "Met" is a column of the scorecard that benchmark 11 reads, per check and fault class. Dev seeds only. Sealed test seeds stay closed (`DR13`).

Synthetic-reference truth: `sw_synref_fidelity` and `sw_synref_spread` against known mixtures. Spread must shrink as events grow. A reference drawn within one batch must not recover a planted shift of that batch, and the benchmark reports this on purpose.

### 3. Real data

| Dataset | Use | What it can and cannot show |
|---|---|---|
| CytoBatchFlagR, 4 controls x 18 batches | Reference checks: observed values per marker and batch; seeded fault injection into a copy of one control in chosen batches (shift, population removal, event loss); optional material swap (an unstimulated control placed where a stimulated one belongs, `DR29`) | The only same-material multi-batch set. Four controls, so AUROC has wide intervals. Frozen and bridge trials on this set met no arm of the D1 bar (survey, text only). Expect some checks to land "observed only". |
| Aurora ex vivo extract, 4 `BMV_REF` files | `sw_ref_vs_study`, observed values per marker next to same-batch donor EMD | Identity of the four files is unverified, so no AUROC and no claim that they are the same material. |
| Mendeley 7-colour file (`autospectral_mendeley`) | Real cells as the source for `sw_synref_downsample`, `_mix`, `_split`: fidelity, overlap guard, draw spread, pseudo-aliquot null, planted shifts and recovery by the batch scores | One file, no batches. It tests the builders and the null, not batch detection. Content of the file not opened by me. |
| BMV unmixing controls (`05_bmv-unmix-compare`) | Not used | Unmixing controls, not a batch design. |
| FR-FCM-Z247 | Not touched | Sealed as the out-of-sample check. |

Integration test `experiments/integration/synref_mendeley.R` (full stack): split the Mendeley file into pseudo-batches, plant a known shift, run `sw_correct_run()` and check that `sw_batch_scores()` before and after reports the planted marker. It needs cyCombine and flowCore.

Every benchmark writes to a new dated run folder with `report/provenance.csv` and `sessionInfo.txt`, saves PNG figures with source CSV, and cites only figures it has opened. Benchmarks call exported functions only (`experiments/README.md` rule 2, opened).

## File layout

R files (one topic each, roxygen with the "uncalibrated" label on every impactful argument):

- `R/batch_scores_population.R`: `sw_batch_scores_by_population`, `sw_batch_interaction`
- `R/batch_matrix.R`: `sw_batch_matrix`; `R/batch_plots.R`: `sw_plot_batch_matrix` (only after `DR30`)
- `R/batch_otd.R`, `R/batch_correct_metrics.R`
- `R/ref_inventory.R`, `R/ref_consistency.R`, `R/ref_noise_floor.R`, `R/ref_vs_study.R`, `R/ref_composition.R`, `R/ref_coverage.R`, `R/ref_split.R`, `R/ref_influence.R`
- `R/synref_draw.R` (class, downsample, mix, split, goal), `R/synref_shift.R`, `R/synref_checks.R` (spread, fidelity, overlap), `R/synref_io.R` (manifest rows, save, load)
- Edits: `R/utils.R`, `R/manifest.R`, `R/design_audit.R`, `R/batch_flag.R`, `R/batch_scores.R` (only if `DR6`), `R/corrector.R`, `R/diagnose.R`

Tests: one `tests/testthat/test-<file>.R` per R file above. Add `NAMESPACE` exports through roxygen. Update `tests/testthat/test-utils.R` for `sw_io_subsample`.

Experiments: `experiments/benchmarks/12_reference-checks/`, `13_synref-truth/`, `experiments/integration/synref_mendeley.R`, `experiments/integration/cytobatchflagr/09_reference_checks.R`, generators under `experiments/benchmarks/_common/`. Benchmark 11 gains the new rows.

Docs: `vignettes/reference-consistency.qmd`, `vignettes/synthetic-references.qmd`, `reference/synthetic-references.qmd`, new sections in `reference/diagnostics.qmd`, an entry in `_quarto.yml`. Update the `CLAUDE.md` table and the ROADMAP register in the same PR as each function, and add a ledger line, benchmark row and coverage-table row for every new check.

## Build order

Small PRs. Each one is tests plus docs plus roxygen.

0. No code. Answer `DR1` to `DR5`. Read `R/corrector.R` lines 1420 to 1550 to settle the hold-out question. Save this plan file.
1. `sw_io_subsample` with required seed and setup, callers and tests migrated (`DR10`).
2. Role `synthetic_reference` in manifest and audit, with regression tests (`DR3`).
3. `sw_synref` class, `sw_synref_downsample`, `sw_synref_split`, `sw_synref_check_overlap`, `sw_synref_as_manifest`.
4. `sw_synref_mix`, `sw_synref_shift`, `sw_synref_spread`, `sw_synref_fidelity`.
5. `sw_ref_inventory`, `sw_ref_consistency`, `sw_ref_noise_floor`, on a shared internal that keeps per-group rows.
6. `sw_ref_vs_study`, `sw_ref_composition`, `sw_ref_coverage`.
7. `sw_flag_iqr` with `min_effect` (`DR9`).
8. `sw_batch_scores_by_population`, `sw_batch_matrix`, `sw_batch_interaction`, then `sw_diagnose` modules (`DR33`).
9. Benchmarks 12 and 13, the Mendeley integration test, CytoBatchFlagR checks, figures inspected, reference page and vignettes.
10. `sw_ref_split`, hold-out in `sw_corrector_compare`, then `sw_ref_influence`.
11. Gated: `sw_batch_otd`, `sw_synref_goal`, `sw_correct_metrics`, `sw_plot_batch_matrix`, and the optional `role` split.
12. Docs sync: `CLAUDE.md`, ROADMAP drift (C6 "gated", A5, A6), scorecard rows.

Synthetic builders come before reference checks because the null arms and the noise floor need pseudo-aliquots.

## Out of scope

- Any correction that uses a synthetic reference as a bridge. A reference-free, lower-certainty workflow would need its own plan and gate (`DR4`).
- Density-dependent downsampling (SPADE), geometric sketching, k-means coverage draws. Both flatten composition. Licences and the sketching paper were not opened (`DR32`).
- kBET and LISI mixing, variance components and ICC, classifier two-sample score, CytoMDS, Levey-Jennings charts (`DR32`).
- Wrappers of GPL code: cytoFlagR, harmony, CytoNorm (`DR31`).
- Exported simulators (`DR5`), `sw_batch_power`, a per-batch drift slope function (`order_rho` already exists), a batch composition function (`sw_cluster_scores` covers it).
- A2b, A3 and gating `sw_correct_run()` on the audit (ROADMAP, frozen).
- New real data. The dataset question G42 is open and unchanged.

## Risks to keep in view

- One real same-material set. Real-data claims for the new reference checks rest on 4 controls x 18 batches. The identity of the `BMV_REF` files is unverified.
- Synthetic references drawn from a batch's own cells carry that batch's effect with its biology. In the cyTether pilot, naive batch centering gave low technical error and 1.70 times the protected-contrast error of `paired_shift` (survey, README text, `unsupported: no inspected figure`). Closed dead ends: internal anchors, rolling reference, single-file bridge (survey).
- Several new required arguments change signatures (`sw_diagnose`, `sw_flag_iqr`, `sw_io_subsample`, `sw_corrector_compare`). Callers in `experiments/` need a migration pass. Benchmarks 07 and 08 `source()` files from `R/` (survey) and break rule 2.
- Plan drift: counts and status differ between `CLAUDE.md`, `ROADMAP.md` and the handoffs. Treat the handoffs as newer, and restate status in the new plan file.

## Settled decisions (user, 2026-11-11)

| ID | Decision | Note |
|---|---|---|
| DR1 | The diagnostics goal covers the new functions, and A2c (`sw_correct_metrics`) is unfrozen too. A2b and A3 stay frozen. | `cytobatchflagr_zenodo` (18 batches, 4 repeated controls each) is on disk, which meets the ROADMAP freeze condition. Edit `ROADMAP.md` line 7 when the plan is committed. Departs from the pick (A2c stayed frozen). |
| DR2 | Ban per-batch and per-reference roll-ups over markers, and any per-reference "consistent" flag. Rank only within check and marker. | Follows Q14. |
| DR3 | New `reference_role` value `synthetic_reference`, treated like `biological_pseudo_control`, never a bridge. `source_batch` note required. | Edits `R/manifest.R` and `R/design_audit.R`. |
| DR4 | A reference-free workflow that uses a synthetic reference is allowed only as a separate, lower-certainty workflow with its own gate, planned later. Nothing in this plan feeds a synthetic reference to `paired_shift` or the HEB bridge. | Departs from the pick (never). |
| DR5 | Simulators and fault injectors are script-local in `experiments/benchmarks/_common/`. Only `sw_synref_shift` is exported. | |
| DR6 | Decline the role split in `sw_batch_scores()` for now. `sw_ref_vs_study` carries the contrast. | |
| DR7 | Extend `sw_batch_scores()` with `population` (NULL means skip). `group_summary` is required whenever `population` is set. | Departs from the pick (new function). About 20 call sites; no caller breaks. |
| DR8 | Population-level checks accept gate or cluster labels. The output records `label_source`. Cluster-label results are flagged batch-dependent and capped at `eligible_for_review`. | |
| DR9 | `sw_flag_iqr()` gets a required named `min_effect` vector, one value per metric, with the unit stated. Values come from a null-arm run on the CytoBatchFlagR extract and need your approval. No number is set. | |
| DR10 | `seed` and `setup` become required in `sw_io_subsample()`. Migrate callers and `test-utils.R`. | |
| DR11 | `n_quantiles` becomes required in `sw_corrector_compare()`. Add the hold-out exclusion only if `R/corrector.R` lines 1420-1550 show fit samples are reused for scoring. | |
| DR12 | EMD-based checks take `n_events` as a required argument. The proposal for the plan is 20,000 events per file and arcsinh with cofactor 6000. OTD takes its own lower cap. | The 20,000 comes from a reported CytoBatchFlagR test, not verified. |
| DR13 | New benchmarks use the dev seeds 441 to 460 (at least the D1 floor of 5). Test seeds 431 to 440 are never read. | |
| DR14 | A check counts as "met" only under D1: mean AUROC at least 0.90 over at least 5 seeds, lower 95% CI at least 0.80, null at most 0.55, per fault class at the smallest real-tested size. Real-data reference checks with 4 controls may be "observed only"; the scorecard says so and the bar is not loosened. | Source: `plans/2026-11-10_diagnostics-next.md` lines 111-118. |
| DR15 | `sw_ref_inventory` and `sw_ref_consistency` take required `min_batches_per_reference` and `min_events`. The plan proposes 3 batches for z-based rows; `pair_emd` alone works at 2. `min_events` is named by the user later. | The 3 is reasoning, not a sourced number. |
| DR16 | `sw_ref_noise_floor` always reports a cell split labelled "sampling only" and adds an aliquot floor when a reference has at least 2 aliquots in a batch. `n_splits` is a required argument the user sets; the plan proposes 20. | |
| DR17 | `sw_ref_vs_study` summarises donor spread as the median of same-batch donor pairwise EMD, with at least 4 donors per batch (dataset rule D5). | Source not opened by me (`03_followups.R`, G42). |
| DR18 | Donor spread comes from same-batch donors only. | Donor is nested in batch (`diagnostics-next.md`). |
| DR19 | `sw_ref_composition` uses an analytic multinomial floor from event counts. No bootstrap, no extra seed. | Own design, no source. |
| DR20 | `sw_ref_coverage` takes a required vector of probabilities. The plan proposes 0.01 and 0.99. | Weak source (CytoBatchNorm percentile range, as reported by the survey). The value is the user's. |
| DR21 | `sw_synref_downsample` and `sw_synref_split` draw without replacement, equal events per file, with strata chosen from none, file, donor or population, and stop on shortfall. Every part is a required argument with no default. The output states that cells are not replicates. | Equal-n and population strata are reasoning, not sourced. |
| DR22 | Synthetic references use the same `n_events` as the EMD checks (proposal 20,000, DR12). `n_draws` is a required argument the user sets; the plan proposes 20 for `sw_synref_spread`. | No source gives a replicate count. |
| DR23 | Synthetic references are within-batch only (`sw_synref_downsample`). Cross-batch blends exist only through `sw_synref_goal` with role `goal`. | Reasoning, not sourced. |
| DR25 | Any event overlap between a synthetic reference and calibration or scoring data stops. Donor is the dependence unit when scoring. `max_overlap_share` is required; the plan proposes 0. | Dependence-group rule from the AutoSpectral bundle (not opened by me). |
