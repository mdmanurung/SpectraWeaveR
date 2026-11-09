# Plan: file-anomaly, header and batch-score diagnostics (CytoScan family)

Branch: `worktree-cytoscan-diagnostics` (based on local `main` 6845d1d, not `origin/main`).
Status: **awaiting approval**. No package code written yet.

## Sources and licence

- CytoScan, CytoTools: GPL (>= 2). CytoScan-manuscript: no licence.
- Wrappers call the installed packages (Suggests). They copy nothing.
- Light versions are clean-room: written from the method text (CytoScan README,
  vignette, manuscript), not from their code. Disclosure: the author of this
  plan read CytoScan's `Flag()`, `generateFeatures.R` and CytoTools'
  `getFCSmetadata()` source during scoping. To keep the light code independent:
  header audit returns a long tibble (file x parameter x keyword), not
  concatenated strings; 1D EMD is written as the mean absolute difference of
  matched quantiles, not a `transport` call.
- File headers say "clean-room from written method description; no upstream code".

## Modules (this branch)

| # | File | Exports | Needs |
|---|------|---------|-------|
| 1 | `R/fcs_header_audit.R` | `sw_fcs_header_audit()` | flowCore |
| 2 | `R/anomaly_features.R` | `sw_anomaly_features()` | base R |
| 3 | `R/anomaly_flag.R` | `sw_anomaly_flag_outlier()`, `sw_anomaly_flag_novelty()` | isotree, mclust (Suggests) |
| 4 | `R/batch_scores.R` | `sw_batch_emd_scores()`, adapter to `sw_batch_candidates_from_metrics()` | base R |
| 5 | `R/cytoscan_wrap.R` | `sw_cytoscan_run()`, `sw_cytotools_header()` | CytoScan, CytoTools (Suggests, `Remotes:`) |

Tests: one file per module. Wrappers use `local_mocked_bindings`; isotree tests use
`skip_if_not_installed()`. Real data checks go in `experiments/`.
Docs: one row per module in the CLAUDE.md table, one `reference/*.qmd`.
Env: add `r-isotree` to `pixi.toml`. CytoScan, CytoTools, flowFP stay GitHub/Bioc Suggests.

## Impactful choices (each its own row; every one is a required argument, no default)

| # | Choice | Source | Your confirm |
|---|--------|--------|--------------|
| C1 | Feature set: summary (mean, sd, median per channel), quantiles, EMD vs aggregate | CytoScan convention | ? |
| C2 | Quantile grid step | CytoScan convention (0.1) | ? |
| C3 | Features deferred: `binning`, `fingerprint` (need flowFP, not installed) | scope cut | ? |
| C4 | Aggregate slot (reference or test), `agg_size`, per-file `n_events`, seed | CytoScan convention (10000, 1000, 42) | ? |
| C5 | Input scale: caller supplies already-transformed data (arcsinh cofactor 6000 per CLAUDE.md); the module never transforms | repo rule | ? |
| C6 | Channels used for features | analyst | ? |
| C7 | Outlier flag: `ntrees`, `ndim`, score cutoff | CytoScan convention (1000, 1, 0.5) | ? |
| C8 | Novelty flag: GMM fit on reference features, log-density quantile cutoff | CytoScan convention (0.05) | ? |
| C9 | Header audit: keywords compared (`$PnN`, `$PnS`, `$PnV`, cytometer id key, export-time key) and what counts as mismatch | analyst | ? |
| C10 | Batch scores: per-marker EMD between batch pairs; mapping to `records` (`comparison_id`, `weight`, `status`); how pair EMDs aggregate | batchdx-cyto adapter | ? |
| C11 | EMD: own implementation (no `transport`) | clean-room | ? |

Values in parentheses are the sources' conventions, shown for reference only.
They are not defaults in the code, and examples mark them "illustrative".

## Deferred (separate plan)

Correction add-ons. None of the three repos has a corrector. Proposed method:
CytoNorm-style per-cluster quantile spline. It needs a `can_freeze` answer for
the corrector contract and a decision on whether it fits the shift-based
`effects` table. Not started.

## Order of work

1. Row confirms (you).
2. Header audit + its tests (smallest, no new deps).
3. Features + EMD + batch-score adapter.
4. Flags (after `pixi install` of isotree).
5. Wrappers, docs, CLAUDE.md rows, NAMESPACE via roxygen, `devtools::test()`.
6. Atomic commits on the worktree branch.

## Settled in grilling (2026-11-09)

| Q | Decision | Source |
|---|---|---|
| Q1 | ROADMAP Track A stays frozen. Module 4 (batch EMD scores) is dropped from this branch. Modules 1, 2, 3, 5 go ahead. | user |
| Q2 | Features: quantiles + EMD (departs from the recommendation, quantiles only). Needs aggregate rows (slot, size, per-file n, seed). Preprint: EMD detected 72% / 83% vs quantiles 97% / 100% (summary read, not verbatim). | user |
| Q3 | Caller supplies transformed data. Module never transforms; required `scale` label. | user |
| Q4 | Caller names header keywords (required). `$PnN` and `$PnS` always checked. | user |
| Q5 | Mismatch = differs from a named reference file (required). | user |
| Q6 | Copy iidx dashboard (Apache-2.0) and adapt to sw_* results. Needs LICENSE text, NOTICE entry, changed-file marks, and a data-contract row. | user |
| Q7 | PeacoQC (GPL >= 3): wrap, never copy. Add `sw_qc_plot()` (PlotPeacoQC) and `sw_qc_heatmap()` (PeacoQCHeatmap). | user |

## Added modules

| # | File | Exports | Needs |
|---|------|---------|-------|
| 6 | `R/qc_plots.R` | `sw_qc_plot()`, `sw_qc_heatmap()` | PeacoQC (already Imports) |
| 7 | `R/qc_app.R` (+ `NOTICE`, Apache-2.0 headers on adapted parts) | `sw_qc_app()`: tabs for header audit, anomaly flags with feature view, PeacoQC plot and heatmap | shiny, bslib, ggiraph, shinycssloaders (Suggests) |

## Open rows (not asked yet)

- None left to ask. C2 quantile step, agg_size, per-file events, seed, forest settings, cutoff and
  GMM quantile are required arguments with no default (ROADMAP rule); examples mark values illustrative.
- GMM fact (corrected): a first test on isotropic random data (always XII, G=1) was
  uninformative, because spherical is the true model there. On correlated synthetic data
  (3 latent factors, 60 files, mclust 6.1.2): 9 and 54 features gave XXX, G=1; 252 features
  gave VEI, G=9 (9 components for 60 points). Synthetic only; not your cohort.

## Settled, round 3

| Q | Decision | Source |
|---|---|---|
| Q8 | iidx: take only the Shiny app framework (layout, theme, widgets). No differential analysis; the package stays preprocessing. Replaces the earlier Q6 panel idea. | user |
| Q9 | EMD: exact quantile coupling, no bins, own code. | user |
| Q10 | Aggregate: required `aggregate_from` ("reference" or "all"), plus required `agg_size`, per-file events, seed. | user |
| Q11 | Novelty flag: stop with an error naming file and feature counts when reference files <= features. No PCA. | user |

| Q12 | New `sw_qc_app()` for file-level QC. `sw_unmix_app()` untouched. | user |
| Q13 | Widgets: cards (bslib) + ggiraph + shinycssloaders only. No shinyjs, shinyBS, bsplus, shinythemes, shinyWidgets. | user |

## Final scope (awaiting go)

Modules: 1 header audit, 2 features (quantiles + exact-coupling EMD), 3 flags (isotree outlier,
mclust novelty with n_ref > features check), 5 CytoScan/CytoTools wrappers, 6 `sw_qc_plot()` /
`sw_qc_heatmap()`, 7 `sw_qc_app()`. Dropped: batch EMD scores (Track A frozen), correction
add-ons (separate plan), iidx differential analysis. PeacoQC is never copied.

## Goal update (2026-11-09, user): consistent diagnosis of technical problems, batch effects extensively

User message: "work autonomously until we can consistently diagnose for any technical problems
including but not limited to batch effects extensively". This overrides the ROADMAP Track A
freeze for this branch (the user's later, explicit instruction). Impactful choices are still
required arguments with no default, so autonomy does not need new answers.

Done so far: header audit, quantile + exact-EMD features, outlier/novelty flags, CytoScan and
CytoTools wrappers, PeacoQC plot and heatmap wrappers, benchmark v1 (50 iterations).

| # | Module | Exports | Fault classes it covers |
|---|--------|---------|-------------------------|
| 8 | `R/fcs_signal_audit.R` | `sw_fcs_signal_audit()` | saturation, dead channel, NA/Inf, negative share, event-rate instability and gaps, channel-pair correlation change (spillover/unmixing error) |
| 9 | `R/batch_scores.R` | `sw_batch_scores()` | per-marker batch shift (pairwise exact EMD, median range), variance explained by batch, within-batch acquisition-order drift |
| 10 | `R/diagnose.R` | `sw_diagnose()` | runs 1, 3, 8, 9 and returns one findings table |
| 7 | `R/qc_app.R` | `sw_qc_app()` | shows a diagnosis (Shiny, bslib, ggiraph) |
| 11 | benchmark v2 | `experiments/benchmarks/07_fault-detection/` | fault injection, detection matrix per tool; gaps listed, not hidden |

"Consistently" is measured, not declared: benchmark v2 reports detection (AUROC) per fault class
and per tool over many seeds with intervals. The pass bar is the user's choice and is not set here.
Known gap from v1: a rare-population change (20% of events) is invisible at file level.
