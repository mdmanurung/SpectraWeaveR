# SpectraWeaveR — Implementation Plan

## Goal

Build a complete, installable R package (`SpectraWeaveR`) that provides an end-to-end pipeline for 40-color spectral flow cytometry analysis. The package wraps five established tools—AutoSpectral, openCyto, PeacoQC, cyCombine, and kohonen/FastPG—into a coherent, format-safe workflow, with additional utilities ported from CytoPipeline.

---

## Phase 1: Package Skeleton ✅

- [x] `DESCRIPTION` — package name, version (0.1.0), license (MIT), all Imports/Suggests
- [x] `LICENSE` — MIT license file
- [x] `NAMESPACE` — export all 69 `sw_*` and `sw_pipeline_run_all` functions; import from dependencies
- [x] `.Rbuildignore` — exclude non-package files (CLAUDE.md, PLAN.md, PROGRESS.md, etc.)
- [x] `README.md` — installation instructions and usage examples
- [x] `CLAUDE.md` — repository summary for AI assistants
- [x] `tests/testthat.R` — standard testthat runner

---

## Phase 2: Core R Source Files ✅

### `R/utils.R` — Format conversion utilities (581 lines)
- [x] `sw_io_read_fcs(files, ...)` — wrapper around `flowCore::read.FCS()` with `truncate_max_range = FALSE` enforced
- [x] `sw_io_ff_to_tibble(ff, sample, batch, condition)` — `exprs(ff)` → tibble with metadata columns
- [x] `sw_io_tibble_to_ff(df, markers)` — tibble → `flowCore::flowFrame`
- [x] `sw_io_exprs_to_tibble(mat, colnames, sample, batch, condition)` — matrix → tibble
- [x] `sw_channel_get_fluor(ff_or_fs)` — return channel names excluding FSC/SSC/Time
- [x] `sw_channel_set_markers(ff, marker_map)` — rename channel descriptions to marker names
- [x] `sw_channel_is_signal(x)` — identify signal columns (ported from CytoPipeline)
- [x] `sw_channel_is_fluor(x)` — identify fluorochrome columns (ported from CytoPipeline)
- [x] `sw_io_subsample(fs, n_total_events, ...)` — pool and subsample events (ported from CytoPipeline)
- [x] `sw_io_event_audit(intermediates)` — audit trail of event filtering

### `R/unmix.R` — Step 1: Spectral unmixing (1,152 lines)
- [x] `sw_unmix_setup(control_dir, cytometer, ...)` — initialize AutoSpectral parameters
- [x] `sw_unmix_prepare(setup, ...)` — gate, clean, and extract fluorophore spectra
- [x] `sw_unmix_extract_af(unstained_fcs, setup, spectra, ...)` — per-cell autofluorescence extraction
- [x] `sw_unmix_extract_variants(setup, spectra, af_spectra, ...)` — fluorophore emission variability mapping
- [x] `sw_unmix_run(input, spectra, setup, flow_control, ...)` — unmix fully-stained sample FCS files
- [x] `sw_unmix_pipeline(control_dir, sample_input, ...)` — one-call end-to-end unmixing orchestrator
- [x] `sw_io_load_unmixed(fcs_dir, pattern, ...)` — load pre-unmixed FCS files into a `flowSet`
- [x] `sw_filter_margins(ff, channel_specs)` — wrapper for `PeacoQC::RemoveMargins()`

### `R/gate.R` — Step 2: Automated gating (191 lines)
- [x] `sw_gate_template(output_file, template_type)` — write a default CSV gating template
- [x] `sw_gate_run(fcs_files, gating_template, transform_channels, ...)` — full openCyto gating workflow
- [x] `sw_gate_extract(gs, node)` — extract `flowFrame` list from a GatingSet node

### `R/qc.R` — Step 3: Signal quality control (200 lines)
- [x] `sw_qc_run(ff, channels, IT_limit, MAD, output_dir, ...)` — wrapper for `PeacoQC::PeacoQC()`
- [x] `sw_qc_batch(ff_list, ...)` — apply QC across a list of flowFrames
- [x] `sw_qc_summary(qc_results, threshold)` — tabulate and flag samples with high removal

### `R/batch_correct.R` — Step 4: Batch correction (802 lines)
#### All-in-one workflow
- [x] `sw_correct_prepare(ff_list, sample_meta, markers, cofactor)` — convert flowFrame list → cyCombine tibble with arcsinh(cofactor=6000)
- [x] `sw_correct_run(uncorrected, markers, covar, ...)` — wrapper for `cyCombine::batch_correct()`
- [x] `sw_correct_evaluate_quick(uncorrected, corrected, markers)` — quick MAD-based evaluation

#### Modular workflow
- [x] `sw_correct_normalize(df, markers, norm_method)` — per-batch normalisation
- [x] `sw_correct_som(df, markers, xdim, ydim, ...)` — SOM clustering for batch correction
- [x] `sw_correct_apply(df, label, markers, covar, ...)` — ComBat correction with pre-computed labels

#### Diagnostics & evaluation
- [x] `sw_correct_detect_batch(df, markers, ...)` — density plots + EMD + MDS
- [x] `sw_correct_emd(df, markers, ...)` — Earth Mover's Distance computation
- [x] `sw_correct_evaluate_emd(uncorrected, corrected, markers, ...)` — EMD-based correction evaluation
- [x] `sw_correct_evaluate_mad(uncorrected, corrected, markers, ...)` — MAD-based correction evaluation

#### Visualization
- [x] `sw_plot_batch_densities(uncorrected, corrected, markers, ...)` — density ridgeline plots
- [x] `sw_plot_batch_dimred(df, markers, type, ...)` — UMAP/PCA projection by batch

### `R/cluster.R` — Step 5: Clustering (382 lines)
- [x] `sw_cluster_run(corrected, lineage_markers, method, ...)` — kohonen SOM (default) or FastPG clustering
- [x] `sw_cluster_assignments(cluster_result)` — per-cell cluster vector
- [x] `sw_cluster_mfi(cluster_result)` — MFI table per cluster
- [x] `sw_plot_clusters(cluster_result, plot_file)` — heatmap (pheatmap or base)
- [x] `sw_cluster_predict(cluster_result, newdata)` — project new data onto trained SOM

### `R/pipeline.R` — End-to-end orchestrator (258 lines)
- [x] `sw_pipeline_run_all(fcs_dir, sample_meta, markers, lineage_markers, ...)` — 6-step pipeline with progress messages

---

## Phase 2b: Extended Modules ✅

### `R/composable_pipeline.R` — Composable pipeline framework (913 lines)
S7/R7-based step-and-pipeline framework ported/inspired by CytoPipeline.

- [x] `sw_pipeline_step(name, FUN, ARGS)` — create a processing step
- [x] `sw_pipeline_step_run(step, input)` — execute a single step
- [x] `sw_pipeline(name, steps)` — create a pipeline
- [x] `sw_pipeline_add/remove/replace/concat` — pipeline manipulation
- [x] `sw_pipeline_step_names/length/show` — pipeline introspection
- [x] `sw_pipeline_run(pipeline, input, trace)` — execute full pipeline
- [x] `sw_plot_pipeline(pipeline, style)` — text-based flowchart
- [x] 10 convenience step constructors (`sw_pipeline_step_read_fcs`, `sw_pipeline_step_filter_margins`, … `sw_pipeline_step_annotate`, `sw_pipeline_step_diff`)

### `R/transforms.R` — Scale transformation utilities (258 lines)
Ported from CytoPipeline (UCLouvain-CBIO).

- [x] `sw_transform_estimate(ff, ...)` — per-channel logicle / linear transforms
- [x] `sw_transform_apply(x, trans_list)` — apply pre-computed transformList

### `R/gating_utils.R` — Standalone gating utilities (455 lines)
Ported from CytoPipeline and CytoPipelineUtils (UCLouvain-CBIO).

- [x] `sw_gate_singlets(ff, ...)` — parallelogram singlet gate
- [x] `sw_filter_doublets(ff, ...)` — doublet removal via CytoPipeline algorithm
- [x] `sw_filter_doublets_peacoqc(ff, ...)` — doublet removal via PeacoQC
- [x] `sw_filter_debris(ff, gate_data, ...)` — manual polygon debris gate
- [x] `sw_filter_margins_peacoqc(x, ...)` — enhanced margin removal via PeacoQC

---

## Phase 3: Documentation ✅

### Roxygen2 documentation
- [x] All 90 exported functions have roxygen2 documentation in source files
- [x] `man/` Rd file names updated to match new function names (Session 12)
- [ ] Regenerate `man/` pages via `devtools::document()` (requires R with dependencies)

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
- [x] `test-gate.R` (~170 lines) — gating template (all 5 types), gating workflow integration, population extraction
- [x] `test-qc.R` (~145 lines) — signal QC validation + PeacoQC integration, batch QC, summary flagging
- [x] `test-batch_correct.R` (443 lines) — all-in-one workflow, modular workflow, diagnostics, evaluation, visualization
- [x] `test-cluster.R` (159 lines) — kohonen SOM and FastPG clustering, prediction, MFI, plotting
- [x] `test-pipeline.R` (~150 lines) — pipeline input validation + end-to-end smoke test with synthetic FCS
- [x] `test-composable_pipeline.R` (560 lines) — step/pipeline creation, manipulation, execution, visualization
- [x] `test-transforms.R` (~120 lines) — scale transform estimation + logicle integration, application verification
- [x] `test-gating_utils.R` (285 lines) — singlet gate, doublet removal, debris gate, margin removal
- [x] `test-dimred.R` (~120 lines) — dimensionality reduction (PCA + UMAP), subsampling, plotting
- [x] `test-unmix_diagnostics.R` (~160 lines) — SSM computation, orthogonality, empirical adjustment, quality metrics

---

## Phase 5: Infrastructure ✅

- [x] `inst/scripts/install_dependencies.R` — pak-based setup script for all dependencies
- [x] `.github/workflows/quarto-publish.yml` — CI/CD for documentation site

---

## Phase 6: Evaluation & Phase 1 Features ✅

### Bug Fixes (Session 9)
- [x] Fixed `sw_filter_margins()` channel selection — signal-only channels
- [x] Fixed `sw_pipeline_run_all()` gating bypass — uses margin-removed data
- [x] Completed EMD implementation in `sw_correct_evaluate_quick()`
- [x] Parameterized cofactor in `sw_gate_run()`

### Phase 1 Features (Session 10)
- [x] `R/dimred.R` — standalone UMAP/PCA module with `sw_dimred_run()` and `sw_plot_dimred()`
- [x] Extended `sw_gate_template()` — 4 new templates (myeloid, nk, treg, full_pbmc)
- [x] `R/unmix_diagnostics.R` — SSM diagnostics with `sw_unmix_spillover_matrix()`, `sw_plot_spillover_matrix()`, `sw_unmix_quality()`

---

## Phase 7: Remaining Work

### Phase 2 Features — HIGH Priority
- [x] **Automated cell type annotation** (`R/annotate.R`)
  - [x] `sw_annotate_run()` — cosine-similarity scoring of cluster MFI profiles against reference
  - [x] `sw_annotate_manual()` — accept user-supplied named vector of cluster → cell type mappings
  - [x] `sw_plot_annotation()` — heatmap or UMAP visualization of annotations
  - [x] Built-in reference profiles for 35 PBMC populations (`inst/extdata/pbmc_reference_matrix.csv` + `pbmc_reference_mask.csv`)
  - [x] `sw_annotate_load_ref()` — load built-in or custom reference CSV pair
  - [x] `sw_pipeline_step_annotate()` — composable pipeline step constructor
  - Leverages: `sw_cluster_mfi()` (already implemented), base R distance calculations
- [x] **Differential abundance & expression analysis** (`R/differential.R`)
  - [x] `sw_diff_abundance()` — edgeR QLF or voom test of cluster proportions between conditions
  - [x] `sw_diff_expression()` — limma per-cluster marker expression differences (per-sample medians)
  - [x] `sw_plot_volcano()` — volcano plot for DE results
  - [x] `sw_plot_boxplots()` — violin+boxplot per marker per group
  - [x] `sw_plot_abundance()` — bar or box chart of cluster proportions
  - [x] `sw_pipeline_step_diff()` — composable pipeline step constructor
  - New Suggests: `edgeR (>= 3.36.0)`, `limma (>= 3.50.0)`, `diffcyt (>= 1.3.0)`

### Phase 2 Features — MEDIUM Priority
- [ ] **Data export/interoperability** (`R/export.R`)
  - `sw_export_csv()`, `sw_export_fcs()`, `sw_export_h5ad()`
  - `sw_to_seurat()`, `sw_to_sce()` — convert to Seurat / SingleCellExperiment objects
  - Leverages: flowCore::write.FCS, anndata (R), SeuratObject, SingleCellExperiment
- [ ] **QC report generation** (`R/report.R`)
  - `sw_generate_report()` — parameterized Quarto/rmarkdown HTML report summarizing pipeline run
  - Event counts per step, QC removal rates, batch effect metrics, cluster compositions
  - Leverages: rmarkdown/quarto templates, plotly, DT
- [ ] **Panel design helper** (`R/panel.R`)
  - `sw_check_panel_compatibility()` — pairwise spectral overlap assessment
  - `sw_suggest_cofactors()` — per-fluorochrome optimal cofactor suggestion
  - `sw_fluorochrome_similarity()` — cosine similarity / spectral angle between reference spectra
  - Requires: reference spectra database in `inst/extdata/`

### Phase 2 Features — LOW Priority (Novel)
- [ ] Fluorochrome autofluorescence fingerprinting (`sw_af_fingerprint()`, `sw_classify_af_populations()`)
- [ ] Longitudinal batch monitoring (`sw_monitor_batch_drift()`, `sw_plot_batch_timeline()`, `sw_levy_jennings()`)

### Phase 8: Function Rename Overhaul ✅
- [x] Renamed all ~90 exported functions to `sw_{module}_{verb}()` hierarchy (Session 12)
- [x] All R source files updated (17 files)
- [x] NAMESPACE updated (90 exports with new names)
- [x] All 16 test files updated
- [x] All vignettes updated (8 files: 2 core + 5 tutorial + 1 llm-assistant)
- [x] All reference docs updated (12 .qmd files)
- [x] All LLM prompt templates updated (`inst/prompts/`)
- [x] `man/` Rd files renamed and content updated
- [x] CLAUDE.md, README.md, PLAN.md, EVALUATION.md, PROGRESS.md updated

### Documentation
- [ ] Add reference docs for `dimred.R` and `unmix_diagnostics.R` (new .qmd pages in `reference/`)
- [ ] Add vignette for composable pipeline framework
- [ ] Add vignette for clustering interpretation (MFI heatmaps, annotation)
- [ ] Convert vignette code blocks from `eval = FALSE` to `eval = TRUE` where possible
- [ ] Add troubleshooting / FAQ section to README

### Infrastructure
- [ ] `R CMD build && R CMD check --no-manual` passes with 0 errors, 0 warnings
- [ ] All tests pass with real package installs
- [ ] Regenerate `man/` pages via `devtools::document()` (old Rd content updated; filenames renamed)
- [ ] Publish Quarto documentation site
- [ ] CRAN / Bioconductor submission preparation (if desired)
