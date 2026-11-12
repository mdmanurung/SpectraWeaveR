# SpectraWeaveR — Repository Summary

## What This Repository Is

**SpectraWeaveR** is an R package for end-to-end processing of 40-color spectral flow cytometry data. It orchestrates five established Bioconductor/GitHub packages into a single reproducible pipeline, handling all format conversions between their incompatible data structures. Additionally, it provides CytoPipeline-ported utility functions and a composable S7-based pipeline framework.

## Repository Structure

```
SpectraWeaveR/
├── DESCRIPTION          # Package metadata, dependencies, and version
├── NAMESPACE            # 136 exported functions and imported namespaces
├── LICENSE              # MIT license
├── README.md            # Project overview, installation, usage examples
├── CLAUDE.md            # This file — repo summary for AI assistants
├── PLAN.md              # Full implementation plan with completion status
├── PROGRESS.md          # Session-by-session progress log
├── R/
│   ├── annotate.R                # Cluster annotation from marker scores
│   ├── assistant.R               # LLM assistant (ellmer) with tool list
│   ├── batch_correct.R           # Step 4: batch correction (cyCombine), diagnostics, plots
│   ├── batch_decide.R            # Batch decision policy and states (batchdx port)
│   ├── batch_flag.R              # Control IQR flags (Tukey fence)
│   ├── batch_scores.R            # Per-marker batch scores: pairwise EMD, median range, eta2, order drift
│   ├── anomaly_features.R        # Per-file quantile + exact-EMD features; sw_emd_1d()
│   ├── anomaly_flag.R            # Outlier (isolation forest) and novelty (GMM) flags
│   ├── cluster_scores.R          # Cluster shares and per-cluster medians as robust z (rare-population checks)
│   ├── cytoscan_wrap.R           # CytoScan and CytoTools wrappers (GPL, Suggests)
│   ├── diagnose.R                # sw_diagnose(): stacks the diagnostic modules
│   ├── fcs_header_audit.R        # FCS header audit against a reference file
│   ├── fcs_signal_audit.R        # Per-file signal checks (saturation, dead channel, rate, correlation)
│   ├── qc_app.R                  # sw_qc_app(): Shiny explorer (layout adapted from iidx, Apache-2.0)
│   ├── relativise.R              # CD3 (reference-population) relativisation, Glehr et al. 2026
│   ├── otd.R                     # Debiased Sinkhorn divergence and per-file OTD
│   ├── negative_audit.R          # Negative-population checks (GMM negative per channel)
│   ├── cluster_stability.R       # Clustering stability across seeds and sizes
│   ├── qc_plots.R                # PeacoQC plot and heatmap wrappers
│   ├── cluster.R                 # Step 5: clustering (kohonen SOM / FastPG)
│   ├── composable_pipeline.R     # S7 composable step and pipeline framework
│   ├── corrector.R               # Corrector contract: fit, integrate, diagnose, save, load, compare
│   ├── corrector_paired_shift.R  # Paired-reference mean shift engine
│   ├── design_audit.R            # Design audit: pair support, anchor bridge, method eligibility
│   ├── differential.R            # Differential abundance and expression
│   ├── dimred.R                  # Dimensionality reduction (UMAP/PCA)
│   ├── flowstate_io.R            # flowstate reader, writer and bridges (wrapped)
│   ├── gate.R                    # Step 2: automated gating (openCyto)
│   ├── gating_utils.R            # Standalone gating utilities (from CytoPipeline)
│   ├── manifest.R                # Validated sample manifest
│   ├── pipeline.R                # End-to-end pipeline orchestrator
│   ├── qc.R                      # Step 3: signal stability QC (PeacoQC)
│   ├── spectracle.R              # Spectra from raw controls (spectracle, wrapped)
│   ├── targets_backend.R         # Optional targets backend for composable pipelines
│   ├── transforms.R              # Scale transforms (from CytoPipeline)
│   ├── unmix.R                   # Step 1: spectral unmixing (AutoSpectral)
│   ├── unmix_assess.R            # AutoSpectral diagnostic wrappers (7 internal, 2 exported)
│   ├── unmix_correct.R           # Spectra correction wrappers (AutoSpectral)
│   ├── unmix_dashboard.R         # Unmix dashboard and app (another worker)
│   ├── unmix_dashboard_quarto.R  # Quarto HTML dashboard (inst/quarto_dashboard template)
│   ├── unmix_dashboard_data.R    # Data model for the unmix dashboard
│   ├── unmix_diagnostics.R       # Unmixing quality: spillover matrix, CV
│   ├── unmix_figures.R           # Approved unmix figure panels F1 to F4
│   ├── unmix_report.R            # Unmix result, diagnose, report
│   ├── unmix_dashboard.R         # sw_unmix_dashboard(): one self-contained interactive HTML
│   ├── unmix_dashboard_data.R    # Run-folder reader and event binning shared by the dashboard and app
│   ├── unmix_app.R               # sw_unmix_app(): Shiny explorer over run folders
│   ├── utils.R                   # Format bridges, channel classes, subsample
│   └── zzz.R                     # Package load hooks (S7 classes)
├── tests/
│   └── testthat/        # 30 test files, ~11,100 lines (testthat edition 3)
├── inst/
│   └── scripts/
│       └── install_dependencies.R  # pak-based setup script
├── vignettes/
│   ├── spectral-unmixing.qmd       # Unmixing vignette
│   └── batch-correction.qmd       # Batch correction vignette
├── reference/           # 12 API reference .qmd pages
├── index.qmd            # Quarto website landing page
├── installation.qmd     # Installation guide
├── _quarto.yml          # Quarto website configuration
└── .github/
    └── workflows/
        └── quarto-publish.yml  # CI/CD for documentation site
```

## The Five-Step Pipeline

| Step | Function(s) | Underlying Package | Input → Output |
|------|-------------|-------------------|----------------|
| 1. Unmixing | `sw_unmix_pipeline()`, `sw_unmix_run()`, `sw_io_load_unmixed()` | AutoSpectral / flowCore | FCS files → unmixed FCS / flowSet |
| 1a. Setup | `sw_unmix_setup()` | AutoSpectral | control dir → sw_setup |
| 1b. Controls | `sw_unmix_prepare()` | AutoSpectral | sw_setup → spectra + flow.control |
| 1c. AF Spectra | `sw_unmix_extract_af()` | AutoSpectral | unstained FCS → AF spectra matrix |
| 1d. Variants | `sw_unmix_extract_variants()` | AutoSpectral | controls → spectral variants |
| 2. Gating | `sw_gate_run()`, `sw_gate_template()` | openCyto / flowWorkspace | flowSet → GatingSet |
| 3. Signal QC | `sw_qc_run()`, `sw_qc_batch()` | PeacoQC | flowFrame (gated) → flowFrame (cleaned) |
| 4. Batch correction | `sw_correct_run()` | cyCombine | flowSet → tibble (corrected) |
| 4a. Modular | `sw_correct_normalize()` → `sw_correct_som()` → `sw_correct_apply()` | cyCombine | step-by-step control |
| 4b. Diagnostics | `sw_correct_detect_batch()`, `sw_correct_emd()`, `sw_correct_evaluate_emd()`, `sw_correct_evaluate_mad()` | cyCombine | batch effect assessment |
| 4c. Visualization | `sw_plot_batch_densities()`, `sw_plot_batch_dimred()` | cyCombine / ggplot2 | density/UMAP/PCA plots |
| 5. Clustering | `sw_cluster_run()`, `sw_plot_clusters()`, `sw_cluster_predict()` | kohonen / FastPG | tibble/matrix → sw_cluster_result |
| End-to-end | `sw_pipeline_run_all()` | all of the above | FCS files → annotated clusters |
| Design audit | `sw_design_audit()` | batchdx-cyto (ported) | sample table → pair support classes, method eligibility |
| Batch decision | `sw_batch_policy()`, `sw_batch_decide()`, `sw_batch_candidates_from_metrics()` | batchdx-cyto (ported) | candidate metrics → state per candidate, identity always kept |
| Cleaning stage | `sw_clean_run()` | PeacoQC / openCyto / flowCore (wrapped) | FCS files, flowFrames or flowSet → margins (raw, `margin_range` on fluorescence channels), arcsinh, PeacoQC, nonDebris, singlets, live gates → live singlet flowFrames + per-sample audit, gate and QC tables. All 10 impactful values are required (`margin_range`, `debris_window`, `singlet_wider_gate`, `pool_samples`, `viability_channel`, `viability_method`, `viability_cut`, `IT_limit`, `MAD`, `cofactor`); a call with missing ones stops once and lists them with how to choose each. `pool_samples = TRUE` writes `groupBy` = number of samples, because openCyto pools only when `groupBy` is not empty. `sw_gate_run()` sets `groupBy` to the number of samples for any template row with `collapseDataForGating` TRUE and an empty `groupBy`; an explicit `groupBy` stays. A `margin_range` lower limit above 0 is read by PeacoQC as 0; the run warns and continues. A sample emptied by the margin step, or a PeacoQC error, stops the run naming the sample. |
| Control flags | `sw_flag_iqr()`, `sw_flag_iqr_summary()` | CytoBatchFlagR method (clean rewrite) | control events → Tukey-fence flags per batch; `min_effect` (required, one value per metric, uncalibrated) sets the smallest distance beyond the fence that counts: `effect` column, status `below_min_effect` (flag FALSE) |
| Unmix checks | `sw_unmix_check_controls()`, `sw_unmix_condition_number()` | AutoSpectral (wrapped) | spectra / FCS → diagnostics tibbles. The other seven wrappers in `R/unmix_assess.R` are internal since 2026-11-09. |
| Alt. I/O and spectra | `sw_io_*_flowstate()`, `sw_unmix_spectracle()`, `sw_io_spectracle_to_spectra()` | flowstate, spectracle (wrapped) | flowstate / controls → tibble, flowFrame, spectra matrix |
| Unmix report | `sw_unmix_result()`, `sw_unmix_diagnose()`, `sw_unmix_report()`, `sw_unmix_limits_autospectral()` | AutoSpectral 1.8.8 (reference) | unmix results → tables (incl. `max_sif`, check `spectra_sif`), provenance, HTML/CSV/JSON report; the limits helper gives `max_sif` warn 4, fail 6 |
| Unmix dashboard | `sw_unmix_dashboard()`, `sw_unmix_app()` | run folders (report/ and figures*/panels/source_data) | one offline HTML file with what-if SIF sliders (display only) and an AutoSpectral tab (its hotspot, similarity and spillover-spread matrices, AF spectra and per-stage spectra as interactive views; its JPG and PDF figures as a linked gallery), or a Shiny app that redraws events for any channel pair; `inst/dashboard/dashboard.html` is the page template |
| Unmix figures | `sw_unmix_sample_events()`, `sw_unmix_figure_*()`, `sw_unmix_save_figures()` | ggplot2 | approved panels F1 to F4 with source data |
| Unmix diagnostics | `sw_unmix_spike_build()`, `sw_unmix_spike_cospike()`, `sw_unmix_spike_score()`, `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()`, `sw_unmix_diagnose(spike=, screen=)`, `sw_unmix_dashboard_quarto()` | benchmark 05 scripts 75, 86b, 88-88h (ported) | spiked raw files and unmixed files → recovery, independence check, squeeze flags; check rows `spike_recovery` and `marker_screen` (item `setup|marker|pop`, not_judged unless `max_share_lost`, `min_median_recovery`, `max_share_moved` or `max_delta_squeezed` is passed); report CSVs `spike_summary.csv`, `baseline_check.csv`, `marker_screen.csv`; a Spike & screen dashboard tab; a Quarto dashboard from `inst/quarto_dashboard` (local file only). Real-data regression: `experiments/integration/spike_bmv.R`, `screen_bmv.R` |
| Unmix correction | `sw_unmix_correct()`, `sw_unmix_correct_signatures()`, `sw_unmix_correct_spillover()`, `sw_unmix_compare_spectra()` | AutoSpectral (wrapped) | spectra → candidate corrected spectra, change table |
| Manifest and correctors | `sw_sample_manifest()`, `sw_corrector_*()` | cyTether design (paired shift engine) | manifest + frozen model → corrected data with per-cell-group actions. Role `synthetic_reference` (needs `source_batch`) is listed in `sw_design_audit()$synthetic_references` and never bridges, like `biological_pseudo_control` |
| FCS header audit | `sw_fcs_header_audit()`, `sw_fcs_header_summary()` | CytoTools idea (clean rewrite) | FCS paths + required keywords + reference file → long tibble of channel and keyword differences |
| Signal audit | `sw_fcs_signal_audit()` | own code | per-file event matrices → numbers for saturation, dead channel, non-finite, event-rate and correlation shifts, and `bimodality`/`bimodality_z`; no cutoffs; `reference` (required, NULL = cohort) scores `median_z`/`sd_z` against a frozen reference |
| Frozen signal reference | `sw_signal_reference_fit()`, `sw_signal_reference_save()`, `sw_signal_reference_load()` | own code | clean reference files → per-channel centre and scale of file medians and SDs, S3 object with digest, for faults most files share |
| Synthetic references | `sw_synref_downsample()`, `sw_synref_mix()`, `sw_synref_split()`, `sw_synref_spread()`, `sw_synref_fidelity()`, `sw_synref_goal()`, `sw_synref_check_overlap()`, `sw_synref_as_manifest()`, `sw_synref_save()`, `sw_synref_load()` | own code | real events → within-batch draws without replacement (strata none, file, donor, population; stop on shortfall) or disjoint pseudo-aliquots (`floor = "sampling only"`) → S3 `sw_synref` with provenance, recipe and digest, `certainty = "diagnostic_only"`; overlap guard stops on shared event, file or donor above `max_overlap_share`; manifest rows with role `synthetic_reference`; `sw_synref_mix()` hits set population shares (requested vs realised, `short_stratum` stop or cap); `sw_synref_spread()` gives sd and range of marker quantiles across draws; `sw_synref_fidelity()` gives marker EMD and population share ratio against the source; `sw_synref_goal()` gives a batch quantile goal table with role `goal` (never a reference); `n_events`, `n_draws`, `seed`, `strata`, `n_parts`, `events_per_part`, `short_stratum`, `goal_statistic`, `probs` required |
| Anomaly features and flags | `sw_anomaly_features()`, `sw_emd_1d()`, `sw_anomaly_flag_outlier()`, `sw_anomaly_flag_novelty()` | CytoScan method (clean rewrite); isotree, mclust (Suggests) | transformed per-file data → quantile + exact-EMD features → isolation-forest score or GMM log-density. Every cutoff is a required argument. |
| Batch-effect scores | `sw_batch_scores()` | own code | per-file data + batch labels → pairwise EMD, batch-vs-rest EMD (`rest_emd`, `rest_emd_z`), median range, eta2, order drift per marker |
| Cluster scores | `sw_cluster_scores()` | own code (stats::kmeans) | per-file cluster shares and per-cluster marker medians as robust z against other files; `k`, `nstart`, events, seed, min events and `reference` (files to fit k-means on, or NULL for all) are required, no scaling; `sw_diagnose` ranks each marker and cluster as its own item |
| Negative-population audit | `sw_negative_audit()` | own code; symptoms from the Colibri blog; mclust (Suggests) | per file and channel, a 2-component 1D GMM (model "V") gives the negative: its median, robust SD, share, separation, hyper-negative share, and joint-negative correlation per channel pair, each with a robust z across files. `n_events`, `seed`, `hyper_k`, `pairs`, `min_separation` required; channels whose median separation is below `min_separation` get `NA` negative checks and `split_ok = 0` |
| Cluster stability | `sw_cluster_stability()` | own code (stats::kmeans) | ARI and moved-event share across seeds and subsample sizes; `k`, seeds, sizes, events, draw seed, nstart required |
| Optimal transport distance | `sw_sinkhorn_divergence()`, `sw_otd_scores()` | own code from the definitions (as used by Glehr et al. 2026) | per-file debiased Sinkhorn divergence to a pooled reference; blur, events, seed, max_iter, tol required |
| Diagnosis | `sw_diagnose()`, `sw_diagnosis_top()` | own code | stacks the header, signal, anomaly, batch, cluster, relativise, otd, negative and gate modules into one table; no pass or fail; files are ranked per check and marker, never by a per-file summary (user decision Q14) |
| Reference relativisation | `sw_relativise_split_kmeans()`, `sw_relativise_fit()`, `sw_relativise_diagnose()`, `sw_relativise_apply()`; maps `sw_relativise_batch_map()` (per-batch controls), `sw_relativise_internal_map()` (reference channel from the full stain) | Glehr et al. 2026 Box 1 (clean rewrite from the paper) | linear controls + control map → per-channel offset and gain per unit → gain diagnostics → relativised linear data. Split scale, seed, nstart, min events, targets and `input_scale = "linear"` are required. Spec: `plans/2026-11-09_reference-relativisation.md` |
| FAUST before and after correction | `sw_faust_gatingset()`, `sw_faust_run()`, `sw_faust_labels()`, `sw_faust_compare()` | FAUST, scamp (GPL, wrapped from Suggests) | per-file matrices → GatingSet → FAUST labels → label agreement (ARI), phenotype shares and medians before vs after correction; every FAUST knob required. Plan: `plans/2026-11-09_faust-wrapper.md` |
| CytoScan family wrappers | `sw_cytoscan_run()`, `sw_cytoscan_preprocess_none()`, `sw_cytotools_header()` | CytoScan, CytoTools (GPL, wrapped from Suggests) | FCS paths → CytoScan flags with the caller's preprocess function |
| QC plots and app | `sw_qc_plot()`, `sw_qc_heatmap()`, `sw_qc_app()` | PeacoQC (wrapped); layout adapted from iidx (Apache-2.0, see NOTICE) | QC result or diagnosis → PNG, heatmap, Shiny explorer |

## Extended Modules

### Composable Pipeline (`R/composable_pipeline.R`)
S7/R7-based step-and-pipeline framework inspired by CytoPipeline:
- `sw_pipeline_step()` / `sw_pipeline_step_run()` — define and run processing steps
- `sw_pipeline()` / `sw_pipeline_run()` — create and execute pipelines
- Pipeline manipulation: `add`, `remove`, `replace`, `concat`, `show`, `plot`
- 8 convenience step constructors for common operations
- **`targets` backend** (`R/targets_backend.R`): opt-in via
  `sw_pipeline_run(..., backend = "targets", store = "_targets")`, or use
  `sw_pipeline_run_targets()` / `sw_pipeline_to_targets()` directly. Each
  `sw_pipeline_step()` becomes one `targets::tar_target()` so reruns are incremental
  and steps can be parallelised. `targets` is in Suggests.

### Scale Transforms (`R/transforms.R`)
Ported from CytoPipeline (UCLouvain-CBIO):
- `sw_transform_estimate()` — logicle / linear quantile per-channel estimation
- `sw_transform_apply()` — apply pre-computed transformList

### Gating Utilities (`R/gating_utils.R`)
Ported from CytoPipeline / CytoPipelineUtils (UCLouvain-CBIO):
- `sw_gate_singlets()` — parallelogram singlet gate
- `sw_filter_doublets()` / `sw_filter_doublets_peacoqc()` — doublet removal
- `sw_filter_debris()` — manual polygon debris gate
- `sw_filter_margins_peacoqc()` — enhanced margin removal

### Utilities (`R/utils.R`)
- Format bridges: `sw_io_ff_to_tibble()`, `sw_io_tibble_to_ff()`, `sw_io_exprs_to_tibble()`
- Channel classification: `sw_channel_is_signal()`, `sw_channel_is_fluor()`
- Aggregation: `sw_io_subsample()` — pool and subsample events; `n_total_events`, `setup` and `seed` are required (`seed = NULL` means do not seed)
- Audit trail: `sw_io_event_audit()` — event count tracking

## Key Design Constraints

- **`truncate_max_range = FALSE`** must be set when reading Aurora FCS files via `flowCore::read.FCS()` — Aurora fluorescence intensities (~4×10⁶) exceed flowCore's default range cutoffs.
- **arcsinh cofactor = 6000** for spectral flow (not 5 as in CyTOF).
- **Marker naming consistency**: FCS channel names (e.g. `BV421-A`) must be mapped to marker names (e.g. `CD3`) before passing data to cyCombine.
- **RemoveMargins must run before transformation** — it relies on raw signal intensities.
- **Input validation runs before dependency checks** in unmixing functions for clear error messages.
- **AutoSpectral is in Suggests** (not Imports) — installed from GitHub: `remotes::install_github('DrCytometer/AutoSpectral')`.

## Dependencies

### Required (Imports)
- `flowCore` (≥ 2.6.0) — FCS I/O, flowFrame/flowSet classes
- `flowWorkspace` (≥ 4.6.0) — GatingSet, cytoset
- `openCyto` (≥ 2.6.0) — CSV-driven hierarchical gating
- `PeacoQC` (≥ 1.4.0) — Isolation Tree + MAD signal QC
- `kohonen` (≥ 3.0.0) — Self-organizing map clustering
- `cyCombine` (≥ 0.2.7) — ComBat batch correction on SOM clusters
- `dplyr` (≥ 1.0.0), `tibble` (≥ 3.1.0) — tibble manipulation
- `methods`, `stats` — base R

### Optional (Suggests)
- `AutoSpectral` (≥ 1.8.8) — spectral unmixing (GitHub)
- `AutoSpectralRcpp` (≥ 1.0.0) — faster unmixing (GitHub)
- `FastPG` — fast graph-based clustering
- `S7` (≥ 0.2.0) — R7 OOP for composable pipeline
- `ggcyto` (≥ 1.22.0) — gate visualization
- `ggplot2` — general plotting
- `uwot` (≥ 0.1.14) — UMAP dimensionality reduction
- `pheatmap` (≥ 1.0.12) — cluster annotation heatmaps
- `testthat` (≥ 3.0.0) — testing
- `knitr`, `rmarkdown`, `quarto` — documentation

## Validated Pipeline Order (den Braanker 2021, Quintelier 2021)

```
unmixing → import → RemoveMargins → transformation →
signal QC → gating → batch correction → clustering → visualization
```

## Build & Test Commands

```r
# Install all dependencies (first time)
source("inst/scripts/install_dependencies.R")

# Document (requires roxygen2)
devtools::document()

# Build
R CMD build .

# Check
R CMD check SpectraWeaveR_*.tar.gz

# Test only
devtools::test()
# or
testthat::test_dir("tests/testthat")
```

## Repository layout: package vs experiments

```
Package (what R CMD build ships)
  R/  man/  tests/  vignettes/  inst/  DESCRIPTION  NAMESPACE  LICENSE  README.md
Project docs and planning (not shipped)
  CLAUDE.md  PLAN.md  PROGRESS.md  ROADMAP.md  EVALUATION.md  plans/
Website (not shipped)
  _quarto.yml  index.qmd  installation.qmd  reference/
Dev environment (not shipped)
  pixi.toml  pixi.lock  .pixi/  tools/ (setup_pixi_env.sh)
Experiments, benchmarks, reproducibility (not shipped; see experiments/README.md)
  experiments/datasets/     registry + verified download script (data live in $SW_DATA_DIR, outside the repo)
  experiments/benchmarks/   one numbered folder per benchmark: README plan + run.R
  experiments/integration/  real-data smoke tests that need the full stack
  experiments/slurm/        sbatch wrappers
  experiments/runs/         dated output folders (git-ignored)
```

Rules: the package never reads `experiments/`; experiments call only exported package functions; unit tests in `tests/` use mocks; data are never committed; each run writes to a new dated folder.

## Roadmap and rules for new code

See `ROADMAP.md`. Latest AutoSpectral is the gold standard; extensions are scored against it. AGPL/GPL packages are wrapped, never copied. Impactful choices (thresholds, limits, tolerances, grids, seeds) are required arguments with no default.
