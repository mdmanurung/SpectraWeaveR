# SpectraWeaveR — Implementation Plan

## Goal

Build a complete, installable R package (`SpectraWeaveR`) that provides an end-to-end pipeline for 40-color spectral flow cytometry analysis. The package wraps five established tools—AutoSpectral, openCyto, PeacoQC, cyCombine, and kohonen/FastPG—into a coherent, format-safe workflow, with additional utilities ported from CytoPipeline.

---

## Phase 1: Package Skeleton ✅

- [x] `DESCRIPTION` — package name, version (0.1.0), license (MIT), all Imports/Suggests
- [x] `LICENSE` — MIT license file
- [x] `NAMESPACE` — export all 69 `sw_*` and `run_pipeline` functions; import from dependencies
- [x] `.Rbuildignore` — exclude non-package files (CLAUDE.md, PLAN.md, PROGRESS.md, etc.)
- [x] `README.md` — installation instructions and usage examples
- [x] `CLAUDE.md` — repository summary for AI assistants
- [x] `tests/testthat.R` — standard testthat runner

---

## Phase 2: Core R Source Files ✅

### `R/utils.R` — Format conversion utilities (581 lines)
- [x] `sw_read_fcs(files, ...)` — wrapper around `flowCore::read.FCS()` with `truncate_max_range = FALSE` enforced
- [x] `sw_flowframe_to_tibble(ff, sample, batch, condition)` — `exprs(ff)` → tibble with metadata columns
- [x] `sw_tibble_to_flowframe(df, markers)` — tibble → `flowCore::flowFrame`
- [x] `sw_exprs_to_tibble(mat, colnames, sample, batch, condition)` — matrix → tibble
- [x] `sw_get_fluor_channels(ff_or_fs)` — return channel names excluding FSC/SSC/Time
- [x] `sw_set_marker_names(ff, marker_map)` — rename channel descriptions to marker names
- [x] `sw_are_signal_cols(x)` — identify signal columns (ported from CytoPipeline)
- [x] `sw_are_fluor_cols(x)` — identify fluorochrome columns (ported from CytoPipeline)
- [x] `sw_aggregate_and_sample(fs, n_total_events, ...)` — pool and subsample events (ported from CytoPipeline)
- [x] `sw_collect_events_retained(intermediates)` — audit trail of event filtering

### `R/unmix.R` — Step 1: Spectral unmixing (1,152 lines)
- [x] `sw_autospectral_setup(control_dir, cytometer, ...)` — initialize AutoSpectral parameters
- [x] `sw_prepare_controls(setup, ...)` — gate, clean, and extract fluorophore spectra
- [x] `sw_extract_af_spectra(unstained_fcs, setup, spectra, ...)` — per-cell autofluorescence extraction
- [x] `sw_extract_spectral_variants(setup, spectra, af_spectra, ...)` — fluorophore emission variability mapping
- [x] `sw_unmix(input, spectra, setup, flow_control, ...)` — unmix fully-stained sample FCS files
- [x] `sw_unmix_pipeline(control_dir, sample_input, ...)` — one-call end-to-end unmixing orchestrator
- [x] `sw_load_unmixed(fcs_dir, pattern, ...)` — load pre-unmixed FCS files into a `flowSet`
- [x] `sw_remove_margins(ff, channel_specs)` — wrapper for `PeacoQC::RemoveMargins()`

### `R/gate.R` — Step 2: Automated gating (191 lines)
- [x] `sw_build_gating_template(output_file, template_type)` — write a default CSV gating template
- [x] `sw_gate(fcs_files, gating_template, transform_channels, ...)` — full openCyto gating workflow
- [x] `sw_extract_gated(gs, node)` — extract `flowFrame` list from a GatingSet node

### `R/qc.R` — Step 3: Signal quality control (200 lines)
- [x] `sw_signal_qc(ff, channels, IT_limit, MAD, output_dir, ...)` — wrapper for `PeacoQC::PeacoQC()`
- [x] `sw_signal_qc_batch(ff_list, ...)` — apply QC across a list of flowFrames
- [x] `sw_qc_summary(qc_results, threshold)` — tabulate and flag samples with high removal

### `R/batch_correct.R` — Step 4: Batch correction (802 lines)
#### All-in-one workflow
- [x] `sw_prepare_for_correction(ff_list, sample_meta, markers, cofactor)` — convert flowFrame list → cyCombine tibble with arcsinh(cofactor=6000)
- [x] `sw_batch_correct(uncorrected, markers, covar, ...)` — wrapper for `cyCombine::batch_correct()`
- [x] `sw_evaluate_correction(uncorrected, corrected, markers)` — quick MAD-based evaluation

#### Modular workflow
- [x] `sw_normalize(df, markers, norm_method)` — per-batch normalisation
- [x] `sw_create_som(df, markers, xdim, ydim, ...)` — SOM clustering for batch correction
- [x] `sw_correct_data(df, label, markers, covar, ...)` — ComBat correction with pre-computed labels

#### Diagnostics & evaluation
- [x] `sw_detect_batch_effect(df, markers, ...)` — density plots + EMD + MDS
- [x] `sw_compute_emd(df, markers, ...)` — Earth Mover's Distance computation
- [x] `sw_evaluate_emd(uncorrected, corrected, markers, ...)` — EMD-based correction evaluation
- [x] `sw_evaluate_mad(uncorrected, corrected, markers, ...)` — MAD-based correction evaluation

#### Visualization
- [x] `sw_plot_batch_densities(uncorrected, corrected, markers, ...)` — density ridgeline plots
- [x] `sw_plot_batch_dimred(df, markers, type, ...)` — UMAP/PCA projection by batch

### `R/cluster.R` — Step 5: Clustering (382 lines)
- [x] `sw_cluster(corrected, lineage_markers, method, ...)` — kohonen SOM (default) or FastPG clustering
- [x] `sw_get_cluster_assignments(cluster_result)` — per-cell cluster vector
- [x] `sw_cluster_mfis(cluster_result)` — MFI table per cluster
- [x] `sw_plot_clusters(cluster_result, plot_file)` — heatmap (pheatmap or base)
- [x] `sw_predict_clusters(cluster_result, newdata)` — project new data onto trained SOM

### `R/pipeline.R` — End-to-end orchestrator (258 lines)
- [x] `run_pipeline(fcs_dir, sample_meta, markers, lineage_markers, ...)` — 6-step pipeline with progress messages

---

## Phase 2b: Extended Modules ✅

### `R/composable_pipeline.R` — Composable pipeline framework (913 lines)
S7/R7-based step-and-pipeline framework ported/inspired by CytoPipeline.

- [x] `sw_step(name, FUN, ARGS)` — create a processing step
- [x] `execute_step(step, input)` — execute a single step
- [x] `sw_pipeline(name, steps)` — create a pipeline
- [x] `sw_pipeline_add/remove/replace/concat` — pipeline manipulation
- [x] `sw_pipeline_step_names/length/show` — pipeline introspection
- [x] `sw_pipeline_run(pipeline, input, trace)` — execute full pipeline
- [x] `sw_plot_pipeline(pipeline, style)` — text-based flowchart
- [x] 8 convenience step constructors (`sw_step_read_fcs`, `sw_step_remove_margins`, etc.)

### `R/transforms.R` — Scale transformation utilities (258 lines)
Ported from CytoPipeline (UCLouvain-CBIO).

- [x] `sw_estimate_scale_transforms(ff, ...)` — per-channel logicle / linear transforms
- [x] `sw_apply_scale_transforms(x, trans_list)` — apply pre-computed transformList

### `R/gating_utils.R` — Standalone gating utilities (455 lines)
Ported from CytoPipeline and CytoPipelineUtils (UCLouvain-CBIO).

- [x] `sw_singlets_gate(ff, ...)` — parallelogram singlet gate
- [x] `sw_remove_doublets(ff, ...)` — doublet removal via CytoPipeline algorithm
- [x] `sw_remove_doublets_peacoqc(ff, ...)` — doublet removal via PeacoQC
- [x] `sw_remove_debris_gate(ff, gate_data, ...)` — manual polygon debris gate
- [x] `sw_remove_margins_peacoqc(x, ...)` — enhanced margin removal via PeacoQC

---

## Phase 3: Documentation ✅

### Roxygen2 documentation
- [x] All 69 exported functions have roxygen2 documentation in source files
- [ ] `man/` pages need to be generated via `devtools::document()` (requires R with dependencies)

### Quarto documentation website
- [x] `_quarto.yml` — website configuration with navbar
- [x] `index.qmd` — landing page
- [x] `installation.qmd` — installation guide
- [x] `vignettes/spectral-unmixing.qmd` — unmixing vignette
- [x] `vignettes/batch-correction.qmd` — batch correction vignette
- [x] 11 API reference pages in `reference/` (one per module + overview)
- [x] `.github/workflows/quarto-publish.yml` — GitHub Actions deployment to Pages

---

## Phase 4: Tests ✅

Tests use synthetic in-memory data (no real FCS files required). External dependency tests use `skip_if_not_installed()`.

- [x] `test-utils.R` (378 lines) — format conversions, channel classification, aggregate/sample, event tracking
- [x] `test-unmix.R` (703 lines) — all 8 unmixing functions, input validation
- [x] `test-gate.R` (72 lines) — gating template, gating workflow, population extraction
- [x] `test-qc.R` (100 lines) — signal QC, batch QC, summary flagging
- [x] `test-batch_correct.R` (443 lines) — all-in-one workflow, modular workflow, diagnostics, evaluation, visualization
- [x] `test-cluster.R` (159 lines) — kohonen SOM and FastPG clustering, prediction, MFI, plotting
- [x] `test-pipeline.R` (98 lines) — pipeline input validation
- [x] `test-composable_pipeline.R` (560 lines) — step/pipeline creation, manipulation, execution, visualization
- [x] `test-transforms.R` (73 lines) — scale transform estimation and application
- [x] `test-gating_utils.R` (285 lines) — singlet gate, doublet removal, debris gate, margin removal

---

## Phase 5: Infrastructure ✅

- [x] `inst/scripts/install_dependencies.R` — pak-based setup script for all dependencies
- [x] `.github/workflows/quarto-publish.yml` — CI/CD for documentation site

---

## Phase 6: Remaining Work

- [ ] `R CMD build && R CMD check --no-manual` passes with 0 errors, 0 warnings
- [ ] All tests pass with real package installs
- [ ] Generate `man/` pages via `devtools::document()`
- [ ] Publish Quarto documentation site
- [ ] CRAN / Bioconductor submission preparation (if desired)
