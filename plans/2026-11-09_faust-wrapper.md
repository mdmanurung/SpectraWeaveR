# Plan: FAUST wrapper for annotation before and after batch correction

Branch: `worktree-cytoscan-diagnostics`. Date: 2026-11-09.
Status: **awaiting approval**. No package code written. Route: wrap now, clean-room later.

## Purpose

Annotate cells with FAUST (Greene et al. 2021, doi:10.1016/j.patter.2021.100372) on the same events
before and after batch correction. The output answers one question: did correction keep the biology?
It reports label agreement (ARI), per-phenotype proportions and per-phenotype marker medians, and
whether a phenotype splits by batch. It is not a fault detector and not a differential-abundance tool.

## Source, licence and install

- faust 0.5.5 is GPL (>= 3) (`DESCRIPTION`). It stays in `Suggests` with `Remotes: RGLab/FAUST`,
  `RGLab/scamp`. Wrapper code calls `faust::` through a lookup helper; nothing is copied or paraphrased.
- Install check, 2026-11-09, `/exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript`:
  `faust: FALSE`, `scamp: FALSE`; `flowWorkspace: TRUE`, `flowCore: TRUE`. The worktree has no `.pixi/`.
  Nothing was installed. The real-data check below is blocked until faust installs.
- API facts and where they were read (`gh api repos/RGLab/FAUST/contents/...`):
  - `man/faust.Rd`: `faust(gatingSet, startingCellPop, activeChannels, channelBounds, experimentalUnit,
    imputationHierarchy, projectPath, depthScoreThreshold, selectionQuantile, nameOccuranceNum,
    supervisedList, debugFlag, threadNum, seedValue, drawAnnotationHistograms, annotationsApproved,
    densitySubSampleThreshold, densitySubSampleSize, densitySubSampleIterations, archDescriptionList,
    plottingDevice)`. Input is a `GatingSet`; `activeChannels` must match the `desc` field of the
    flowFrame parameters exactly. Returns `NULL`; output is `projectPath/faustData/faustCountMatrix.rds`
    (rows = `sampleNames(gatingSet)`, columns = phenotype strings, `0_0_0_0_0` = unclassified). It
    check-points intermediates; `experimentalUnit` and `imputationHierarchy` name `pData` columns.
  - `vignettes/faustTuning.Rmd`: `faust()` wraps `generateAnnotationThresholds()` then
    `discoverPhenotypes()`; `plotData/scoreLines.*` shows depth scores; unlink `faustData` before a
    re-run with new settings; `startingCellPop = "root"` works on an ungated GatingSet;
    `faustData/sampleData/<sample>/exprsMat.rds` holds per-sample expression.
  - `man/makeAnnotationEmbedding.Rd`: `makeAnnotationEmbedding(projectPath, sampleNameVec, vizType,
    embeddingSelectionQuantile, embeddingDepthScoreThreshold)` returns one row per observation with
    `faustLabels`, `sampleOfOrigin`, `umapX`, `umapY` and winsorised, unit-scaled marker columns.

## Steps (phase 1)

1. `R/faust_wrap.R`: `.sw_faust_fn(name)` mirrors `.sw_cytoscan_fn()` (requireNamespace, install hint).
2. `sw_faust_gatingset(data, markers, sample_table, n_events, subsample_seed)`: per-file matrices
   (the `events.rds` layout), subsampled to `n_events` rows with `subsample_seed`, to `flowFrame` with
   `desc` = marker names, then `flowSet`, then `flowWorkspace::GatingSet()`; `pData` gets `name`,
   `batch`, `donor`, `role` from `sample_table`. One identity check: `names(data)` equals
   `sample_table$sample`. Markers are aligned by name, never by position. The kept row index is
   returned as an attribute so `sw_faust_labels()` joins labels to the same rows.
3. `sw_faust_run(gs, project_path, starting_pop, markers, depth_score_threshold, selection_quantile,
   experimental_unit, imputation_hierarchy, channel_bounds, name_occurrence_num, supervised_list,
   seed, threads, density_subsample_threshold, density_subsample_size, density_subsample_iterations,
   annotations_approved, draw_histograms)`: stops if `project_path/faustData` exists (run dir never
   overwritten), calls `faust()`, returns `project_path` with class `sw_faust_run`, a `status`
   (`complete`, or `thresholds_only` when `annotations_approved = FALSE` stops before discovery) and an
   attribute echoing every argument plus `packageVersion("faust")`. `draw_histograms` passes through to
   `drawAnnotationHistograms`; the `plotData/hist_<marker>_ab_1.*` plots show thresholds per
   experimental unit and change no numbers. Fixed inside: `plottingDevice = "png"`,
   `debugFlag = FALSE`, `archDescriptionList = list(targetArch = "singleCPU")`.
4. `sw_faust_labels(run, data)`: long tibble `sample, cell, phenotype` from per-cell labels, plus
   `sample, phenotype, n, proportion` from `faustCountMatrix.rds`. Per-cell route: the documented one is
   `makeAnnotationEmbedding()` (`faustLabels`, `sampleOfOrigin`); medians use the caller's `data`,
   joined to labels by within-sample row position. Row-order preservation is unverified, so the first
   real run adds one check: per-sample label counts equal the count-matrix row. A run with status
   `thresholds_only` stops here with a message naming the `scoreLines` review and the re-run.
5. `sw_faust_compare(before, after, data_before, data_after, batch, role, markers, match_by,
   ari_scope, median_scope)`: long tibble in the `sw_diagnose()` shape (`module = "faust"`, `subject`,
   `check`, `item`, `value`, `text`). Checks: `ari` (per file or pooled), `prop_before`, `prop_after`,
   `prop_diff`, `median_before`, `median_after`, `median_diff` (item = marker), `prop_by_batch` (per
   phenotype, batch range of proportions, split by `role` ref vs donor). Phenotypes present on one side
   only get `NA` on the other, never a drop. The unclassified label is matched as all-zero tokens, not
   as the literal `0_0_0_0_0` (its length follows the marker count, unverified), and stays as a row.
   One identity check: `before` and `after` carry equal kept-row attributes per sample, else stop;
   `data_after` must keep the row order of `data_before` (a property of `sw_correct_run()` output
   that is unverified). No pass or fail column.
6. Docs: roxygen, `NAMESPACE`, a CLAUDE.md table row, `reference/faust.qmd`; `DESCRIPTION` Suggests
   and Remotes; `inst/scripts/install_dependencies.R` entry. Then atomic commits.

Phase-1 "optimised" means: `n_events` per file (required), `project_path` per run reused as the
cache (FAUST check-points), `threads` passed to `threadNum`, both runs sharing one marker vector.

## Impactful choices (one row each; required argument, no default)

| # | Choice | Argument | Source |
|---|---|---|---|
| F1 | Marker set given to FAUST, identical for before and after | `markers` | needs user approval |
| F2 | Starting population (`"root"` or a gate node) | `starting_pop` | needs user approval |
| F3 | Depth-score threshold (FAUST convention 0.01) | `depth_score_threshold` | needs user approval |
| F4 | Selection quantile (FAUST convention 0.5) | `selection_quantile` | needs user approval |
| F5 | Experimental unit (`""` = per sample; donor; batch). Batch standardises thresholds within batch and hides part of the effect under test | `experimental_unit` | needs user approval |
| F6 | Imputation hierarchy (`""` or a `pData` column); same caveat as F5 | `imputation_hierarchy` | needs user approval |
| F7 | Channel bounds: `""` empirical, or a 2 x markers matrix (vignette: Low at the image of 0 under the transform) | `channel_bounds` | needs user approval |
| F8 | Name occurrence count for gating out (FAUST convention 0) | `name_occurrence_num` | needs user approval |
| F9 | Supervised boundary preferences per marker, or `NA` | `supervised_list` | needs user approval |
| F10 | Seed | `seed` | needs user approval |
| F11 | Density subsample threshold, size, iterations (they move gate locations) | three `density_subsample_*` | needs user approval |
| F12 | Skip FAUST's threshold review (`annotationsApproved = TRUE`) or stop after thresholds | `annotations_approved` | needs user approval |
| F13 | Events per file and the subsample seed | `n_events`, `subsample_seed` | needs user approval |
| F14 | Threads; whether `threadNum` changes results is unverified | `threads` | needs user approval |
| F15 | Input scale: caller passes transformed data (arcsinh cofactor 6000); wrapper never transforms | inherited: `CLAUDE.md` | needs user approval |
| F16 | Correction that makes the "after" data, with its settings | inherited: `sw_correct_run()` | needs user approval |
| F17 | Phenotype match rule between runs (exact string; marker order fixed by F1) | `match_by` | needs user approval |
| F18 | ARI per file or pooled | `ari_scope` | needs user approval |
| F19 | Median: per phenotype per file, or pooled per phenotype | `median_scope` | needs user approval |

## Tests

- `tests/testthat/test-faust_wrap.R`, mocked with `local_mocked_bindings(.sw_faust_fn = ...)`. The
  mock `faust` writes a fake `faustData/faustCountMatrix.rds` and a fake embedding to `project_path`,
  because the wrapper reads from disk. Cases: missing argument stops; existing `faustData` stops;
  marker and `desc` mismatch stops; labels tibble shape; compare table shape; one-sided phenotype gives
  `NA`; all-zero label kept; kept-row mismatch stops; `thresholds_only` stops in `sw_faust_labels()`;
  ARI 1 when labels are identical.
- Real-data check (after faust installs): `experiments/integration/aurora_exvivo/05_faust_before_after.R`
  on `experiments/runs/2026-11-09_aurora-exvivo-real/events.rds` (149 files, 4 batches, 5000 events
  each, one `BMV_REF` file per batch). Writes a new dated run folder with the compare table and
  `fig_faust_prop_by_batch.png`: per-phenotype proportion by batch, one point per file, REF files
  marked, donors as points over the batch median. Conclusions cite that figure after inspection. An ARI
  figure (`fig_faust_ari.png`) needs F16 approved first.

## Deferred to phase 2

A clean-room rewrite of one step from the paper only, never from the GPL source. Trigger: on the Aurora
run with approved settings, one profiled step takes more than the user's stated budget (to be set in
hours per run) or fails on memory. Profile with `system.time()` per phase (`generateAnnotationThresholds`,
`discoverPhenotypes`) before any rewrite. Also deferred: FAUST's slurm `archDescriptionList`, PFDA.

## Risks

- Donor is nested in batch (37, 37, 37, 38 files). A phenotype that splits by batch among donor
  files is confounded with donor; only the four `BMV_REF` files give an unconfounded batch contrast.
  The `role` split in `prop_by_batch` keeps the two apart; the plan draws no causal conclusion.
- Install: the README pins BiocManager 3.13 and the `NAMESPACE` imports `flowWorkspace::getData`;
  the env runs R 4.4.3 with flowWorkspace 4.x. Check `getData` exists before relying on the install.
- Runtime and memory: the docs give no figures. Dip tests run on all cells (faust.Rd), so `n_events`
  is the lever. Unverified until the first run.
- Licence: GPL-3 through Suggests only; no code in this repo derives from faust.
- Two runs can select different markers; F1 fixes the vector, but FAUST may still drop a marker by
  depth score on one side, which changes phenotype strings. `match_by` and the `NA` rule surface this.
