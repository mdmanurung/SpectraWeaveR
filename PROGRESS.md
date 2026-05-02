# SpectraWeaveR — Progress Log

---

## Session 1 — 2026-04-06 (Initial Setup)

### Completed

- [x] Explored repository: empty repo with only `README.md` and `.gitignore`
- [x] Confirmed no existing CI workflows, tests, or R source files
- [x] Created `DESCRIPTION` — package metadata, version 0.1.0, MIT license, all Imports/Suggests
- [x] Created `LICENSE` — MIT license text
- [x] Created `CLAUDE.md` — repository summary for AI assistants
- [x] Created `PLAN.md` — full phased implementation plan
- [x] Created `PROGRESS.md` — this file

---

## Session 2 — 2026-04-06 (Core Implementation)

### Completed

- [x] `NAMESPACE` — export declarations for all `sw_*` and `sw_pipeline_run_all` functions; imports from flowCore, dplyr, tibble, methods, stats
- [x] `.Rbuildignore` — exclude CLAUDE.md, PLAN.md, PROGRESS.md, .Rproj, .github
- [x] `tests/testthat.R` — standard testthat runner
- [x] `R/utils.R` — 6 format conversion utilities (sw_io_read_fcs, sw_io_ff_to_tibble, sw_io_tibble_to_ff, sw_io_exprs_to_tibble, sw_channel_get_fluor, sw_channel_set_markers)
- [x] `tests/testthat/test-utils.R` — unit tests for format conversions
- [x] `R/unmix.R` — 3 initial unmixing functions (sw_unmix_autospectral, sw_io_load_unmixed, sw_filter_margins)
- [x] `tests/testthat/test-unmix.R` — unit tests for unmixing wrappers
- [x] `R/gate.R` — 3 gating functions (sw_gate_template, sw_gate_run, sw_gate_extract)
- [x] `tests/testthat/test-gate.R` — unit tests for gating
- [x] `R/qc.R` — 3 QC functions (sw_qc_run, sw_qc_batch, sw_qc_summary)
- [x] `tests/testthat/test-qc.R` — unit tests for signal QC
- [x] `R/batch_correct.R` — 3 batch correction functions (sw_correct_prepare, sw_correct_run, sw_correct_evaluate_quick)
- [x] `tests/testthat/test-batch_correct.R` — unit tests for batch correction
- [x] `R/cluster.R` — 5 clustering functions using FlowSOM
- [x] `tests/testthat/test-cluster.R` — unit tests for clustering
- [x] `R/pipeline.R` — sw_pipeline_run_all() end-to-end orchestrator
- [x] `tests/testthat/test-pipeline.R` — unit tests for pipeline validation
- [x] `README.md` — updated with installation instructions, quick start, step-by-step usage examples

---

## Session 3 — Post-PR #10 (Replace FlowSOM with kohonen/FastPG)

### Completed

- [x] Replaced FlowSOM dependency with kohonen SOM (Imports) + FastPG (Suggests) in `R/cluster.R`
- [x] `sw_cluster_run()` now uses kohonen::som() with hierarchical metaclustering (ward.D2 + cutree), or FastPG::fastCluster()
- [x] Returns `sw_cluster_result` S3 class with assignments, data, model, etc.
- [x] Added `sw_cluster_predict()` — project new data onto trained SOM
- [x] Updated `DESCRIPTION` — replaced FlowSOM with kohonen (≥ 3.0.0) in Imports, FastPG in Suggests
- [x] Updated `NAMESPACE` — replaced sw_map_new_data with sw_cluster_predict
- [x] Updated `tests/testthat/test-cluster.R` — tests for both clustering methods

---

## Session 4 — Post-PR #10 (Unmixing Expansion)

### Completed

- [x] Expanded `R/unmix.R` from 3 to 8 functions covering the full AutoSpectral workflow:
  - `sw_unmix_setup()` — initialize parameters and control file
  - `sw_unmix_prepare()` — gate, clean, extract fluorophore spectra
  - `sw_unmix_extract_af()` — per-cell autofluorescence extraction (with multi-tissue support)
  - `sw_unmix_extract_variants()` — fluorophore emission variability mapping
  - `sw_unmix_run()` — unmix with OLS/WLS/Poisson/AutoSpectral methods
  - `sw_unmix_pipeline()` — one-call end-to-end unmixing orchestrator
- [x] AutoSpectral moved to Suggests (not Imports); all unmixing functions guard with `.check_autospectral()`
- [x] Input validation runs BEFORE dependency checks for clear error messages
- [x] Added AutoSpectralRcpp as optional Suggests for performance
- [x] Updated `NAMESPACE` with all new unmixing exports
- [x] Expanded `tests/testthat/test-unmix.R` to 703 lines covering all 8 functions

---

## Session 5 — Post-PR #10 (Batch Correction Expansion)

### Completed

- [x] Expanded `R/batch_correct.R` from 3 to 12 functions:
  - Modular workflow: `sw_correct_normalize()`, `sw_correct_som()`, `sw_correct_apply()`
  - Diagnostics: `sw_correct_detect_batch()`
  - Evaluation: `sw_correct_emd()`, `sw_correct_evaluate_emd()`, `sw_correct_evaluate_mad()`
  - Visualization: `sw_plot_batch_densities()`, `sw_plot_batch_dimred()`
- [x] Added DRY validation helpers: `.check_cycombine()`, `.validate_batch_df()`
- [x] Updated `NAMESPACE` with all 12 batch correction exports
- [x] Expanded `tests/testthat/test-batch_correct.R` to 443 lines

---

## Session 6 — Post-PR #10 (CytoPipeline Ports & Composable Pipeline)

### Completed

- [x] Created `R/composable_pipeline.R` (913 lines) — S7/R7 composable pipeline framework:
  - ProcessingStep and Pipeline S7 classes with lazy creation and caching
  - Full pipeline manipulation (add, remove, replace, concat, run, show, plot)
  - 8 convenience step constructors for common SpectraWeaveR operations
  - Text-based flowchart visualization
- [x] Created `R/transforms.R` (258 lines) — scale transformation utilities ported from CytoPipeline:
  - `sw_transform_estimate()` — per-channel logicle and linear quantile estimation
  - `sw_transform_apply()` — apply pre-computed transformList
- [x] Created `R/gating_utils.R` (455 lines) — standalone gating utilities ported from CytoPipeline/CytoPipelineUtils:
  - `sw_gate_singlets()` — parallelogram singlet gate
  - `sw_filter_doublets()` — CytoPipeline doublet removal algorithm
  - `sw_filter_doublets_peacoqc()` — PeacoQC-based doublet removal
  - `sw_filter_debris()` — manual polygon debris gate
  - `sw_filter_margins_peacoqc()` — enhanced margin removal with per-channel specs
- [x] Extended `R/utils.R` with:
  - `sw_channel_is_signal()` / `sw_channel_is_fluor()` — channel classification
  - `sw_io_subsample()` — aggregate and subsample flow frames
  - `sw_io_event_audit()` — audit trail of event filtering
- [x] Added S7 (≥ 0.2.0) and ggplot2 to Suggests in DESCRIPTION
- [x] Updated NAMESPACE with 31 additional exports (total: 69)
- [x] Created tests:
  - `tests/testthat/test-composable_pipeline.R` (560 lines)
  - `tests/testthat/test-transforms.R` (73 lines)
  - `tests/testthat/test-gating_utils.R` (285 lines)
  - Extended `tests/testthat/test-utils.R` to 378 lines

---

## Session 7 — Post-PR #11 (Setup Script & Documentation Site)

### Completed

- [x] Created `inst/scripts/install_dependencies.R` — pak-based setup script to install all dependencies
- [x] Created Quarto documentation website:
  - `_quarto.yml` — website configuration with navbar
  - `index.qmd` — landing page
  - `installation.qmd` — installation guide
  - `vignettes/spectral-unmixing.qmd` — unmixing vignette
  - `vignettes/batch-correction.qmd` — batch correction vignette
  - 11 API reference pages in `reference/` directory
- [x] Created `.github/workflows/quarto-publish.yml` — CI/CD for Quarto documentation site
- [x] Updated `.Rbuildignore` for Quarto files

---

## Session 8 — 2026-04-06 (Code Audit & Documentation Refresh)

### Completed

- [x] Audited all 10 R source files (5,192 lines total)
- [x] Audited all 10 test files (2,871 lines total)
- [x] Verified NAMESPACE has 69 exports, all matching implemented functions
- [x] Verified DESCRIPTION, .Rbuildignore, _quarto.yml are consistent
- [x] Updated `PLAN.md` — reflects actual implementation status with all phases checked
- [x] Updated `PROGRESS.md` — full session-by-session history (this file)
- [x] Updated `CLAUDE.md` — matches current repository structure, functions, and dependencies
- [x] Updated `README.md` — comprehensive feature listing and usage examples

---

## Session 9 — 2026-04-06 (Codebase Evaluation & Bug Fixes)

### Completed

- [x] Created `EVALUATION.md` — comprehensive codebase evaluation covering technical quality (8.5/10), bioinformatics correctness, test gaps, and 11 proposed new features
- [x] Fixed `sw_filter_margins()` (`R/unmix.R:1136-1139`) — now uses `sw_channel_is_signal(ff)` instead of `seq_len(ncol(ff))` to avoid margin removal on Time/Original_ID/scatter channels
- [x] Fixed `sw_pipeline_run_all()` gating step (`R/pipeline.R:173-176`) — now uses margin-removed `ff_list` instead of re-reading raw FCS files from disk
- [x] Completed EMD implementation in `sw_correct_evaluate_quick()` (`R/batch_correct.R:306`) — computes actual EMD values when cluster labels and cyCombine are available
- [x] Parameterized cofactor in `sw_gate_run()` (`R/gate.R:132`) — cofactor is now a function parameter with default 6000

---

## Session 10 — 2026-04-06 (Critical Test Gaps & Phase 1 Features)

### Completed

- [x] Added integration tests for `sw_qc_run()` and `sw_qc_batch()` in `test-qc.R`
- [x] Added integration tests for `sw_gate_run()` and `sw_gate_extract()` in `test-gate.R`
- [x] Added integration tests for `sw_transform_estimate()` and `sw_transform_apply()` in `test-transforms.R`
- [x] Added end-to-end smoke test for `sw_pipeline_run_all()` in `test-pipeline.R` (with synthetic FCS files)
- [x] Created `R/dimred.R` — standalone dimensionality reduction module:
  - `sw_dimred_run()` — UMAP (uwot) and PCA with auto-subsampling, metadata preservation
  - `sw_plot_dimred()` — ggplot2 scatter plots with continuous/discrete colour support
- [x] Created `tests/testthat/test-dimred.R` — 12 tests for dimred module
- [x] Extended `sw_gate_template()` in `R/gate.R`:
  - Added 4 new template types: `"myeloid"`, `"nk"`, `"treg"`, `"full_pbmc"`
  - Refactored template building with internal `.make_gate_row()` and `.build_template()` helpers
- [x] Added template validation tests in `test-gate.R` (parent hierarchy checks)
- [x] Created `R/unmix_diagnostics.R` — unmixing quality diagnostics module:
  - `sw_unmix_spillover_matrix()` — SSM via absolute cosine similarity of reference spectra, with optional empirical adjustment from unmixed data
  - `sw_plot_spillover_matrix()` — heatmap visualization (pheatmap with fallback to base)
  - `sw_unmix_quality()` — per-channel CV analysis with flagging
- [x] Created `tests/testthat/test-unmix_diagnostics.R` — 13 tests for SSM and quality metrics
- [x] Updated `NAMESPACE` — added 5 new exports
- [x] Updated `PROGRESS.md` and `PLAN.md`

---

## Session 11 — 2026-04-07 (Cell Type Annotation & Differential Analysis)

### Completed

- [x] Created `inst/extdata/pbmc_reference_matrix.csv` — built-in PBMC reference matrix (35 populations × 30 markers, values 0–3)
- [x] Created `inst/extdata/pbmc_reference_mask.csv` — binary marker relevance mask for informative markers per population
- [x] Created `R/annotate.R` (~430 lines) — cell type annotation module:
  - `sw_annotate_load_ref(name)` — load built-in or custom CSV reference pair
  - `sw_annotate_run(x, reference, markers, min_score, ...)` — cosine-similarity scoring of cluster MFI profiles vs reference; returns annotation tibble with top-1 and top-2 matches
  - `sw_annotate_manual(x, annotation_map)` — user-supplied cluster → cell type mapping
  - `sw_plot_annotation(annotation, cluster_result, dimred, type)` — heatmap or UMAP plot
  - Internal helpers: `.check_annotation_input()`, `.cosine_sim()`, `.score_cluster()`
- [x] Created `R/differential.R` (~520 lines) — differential abundance & expression module:
  - `sw_diff_abundance(x, meta, group_col, method, ...)` — edgeR QLF or voom test of cluster proportions; returns class `sw_da_result`
  - `sw_diff_expression(corrected, meta, group_col, ...)` — per-cluster limma on per-sample medians; returns long-format class `sw_de_result`
  - `sw_plot_volcano(de_result, cluster, fdr_threshold, fc_threshold, label_top)` — volcano plot
  - `sw_plot_boxplots(corrected, markers, group_col, ...)` — violin + boxplot per marker per group
  - `sw_plot_abundance(da_result, type, ...)` — bar or column chart of cluster proportions
  - Internal helpers: `.check_diffcyt_deps()`, `.check_limma()`, `.validate_meta_da()`, `.build_count_matrix()`, `.compute_cluster_medians()`
- [x] Updated `R/composable_pipeline.R`:
  - Added `annotate` and `differential` entries to `.SW_STEP_DEFAULT_TYPES`
  - Added `sw_pipeline_step_annotate()` convenience constructor (wraps `sw_annotate_run`)
  - Added `sw_pipeline_step_diff()` convenience constructor (wraps `sw_diff_expression`)
- [x] Updated `NAMESPACE` — added 11 new exports (total: 85)
- [x] Updated `DESCRIPTION` Suggests — added `edgeR (>= 3.36.0)`, `limma (>= 3.50.0)`, `diffcyt (>= 1.3.0)`
- [x] Created `tests/testthat/test-annotate.R` (~260 lines) — 34 tests for annotation module
- [x] Created `tests/testthat/test-differential.R` (~300 lines) — 36 tests for differential module
- [x] Updated `PLAN.md` Phase 7 — marked annotation + differential HIGH priority items as `[x]` completed

---

## Session 12 — 2026-05-02 (Full Function Rename Overhaul)

### Completed

Renamed all ~90 exported functions to a consistent `sw_{module}_{verb}()` hierarchy enabling autocompletion-driven discovery, mirroring Python submodule conventions. Clean break — no deprecated aliases.

**Rename summary by module:**

| Module prefix | Functions renamed | Examples |
|---|---|---|
| `sw_io_*` | 7 | `sw_read_fcs` → `sw_io_read_fcs`, `sw_flowframe_to_tibble` → `sw_io_ff_to_tibble`, `sw_aggregate_and_sample` → `sw_io_subsample` |
| `sw_channel_*` | 4 | `sw_are_signal_cols` → `sw_channel_is_signal`, `sw_set_marker_names` → `sw_channel_set_markers` |
| `sw_unmix_*` | 7 | `sw_autospectral_setup` → `sw_unmix_setup`, `sw_extract_af_spectra` → `sw_unmix_extract_af`, `sw_spillover_spreading_matrix` → `sw_unmix_spillover_matrix` |
| `sw_gate_*` | 4 | `sw_gate` → `sw_gate_run`, `sw_build_gating_template` → `sw_gate_template` |
| `sw_filter_*` | 5 | `sw_remove_margins` → `sw_filter_margins`, `sw_remove_doublets` → `sw_filter_doublets`, `sw_remove_debris_gate` → `sw_filter_debris` (consolidated from unmix.R and gating_utils.R) |
| `sw_qc_*` | 2 | `sw_signal_qc` → `sw_qc_run`, `sw_signal_qc_batch` → `sw_qc_batch` |
| `sw_transform_*` | 2 | `sw_estimate_scale_transforms` → `sw_transform_estimate`, `sw_apply_scale_transforms` → `sw_transform_apply` |
| `sw_correct_*` | 10 | `sw_batch_correct` → `sw_correct_run`, `sw_normalize` → `sw_correct_normalize`, `sw_evaluate_correction` → `sw_correct_evaluate_quick`, `sw_detect_batch_effect` → `sw_correct_detect_batch` |
| `sw_cluster_*` | 4 | `sw_cluster` → `sw_cluster_run`, `sw_cluster_mfis` → `sw_cluster_mfi`, `sw_predict_clusters` → `sw_cluster_predict` |
| `sw_dimred_*` | 1 | `sw_run_dimred` → `sw_dimred_run` |
| `sw_annotate_*` | 2 | `sw_annotate_clusters` → `sw_annotate_run`, `sw_load_reference` → `sw_annotate_load_ref` |
| `sw_diff_*` | 2 | `sw_differential_abundance` → `sw_diff_abundance`, `sw_differential_expression` → `sw_diff_expression` |
| `sw_pipeline_*` | 16 | `run_pipeline` → `sw_pipeline_run_all`, `sw_step` → `sw_pipeline_step`, `execute_step` → `sw_pipeline_step_run`, all `sw_step_*` → `sw_pipeline_step_*` |
| `sw_assistant_*` | 1 | `sw_mcp_server` → `sw_assistant_mcp` |

**Files modified:**
- [x] 17 R source files — function definitions and all internal cross-calls
- [x] `NAMESPACE` — all 90 export() entries updated
- [x] 16 test files — all function call sites
- [x] 8 vignette files (`vignettes/`) — code blocks and prose
- [x] 12 reference docs (`reference/`) — API pages
- [x] 4 LLM prompt templates (`inst/prompts/`) — code examples in prompts
- [x] 55 `man/` Rd files — content updated + files renamed to match new function names
- [x] `CLAUDE.md`, `README.md`, `PLAN.md`, `EVALUATION.md`

**Design decisions applied:**
- `sw_filter_*` is now a unified module — previously `sw_remove_margins()` was in `unmix.R` and `sw_remove_margins_peacoqc()` was in `gating_utils.R`; both now share `sw_filter_*` prefix
- `sw_evaluate_correction()` → `sw_correct_evaluate_quick()` to distinguish it from the cyCombine-powered `sw_correct_evaluate_emd()` and `sw_correct_evaluate_mad()`
- `sw_unmix_spillover_matrix()` replaces cryptic `sw_spillover_spreading_matrix()` acronym
- Internal step-label strings in `composable_pipeline.R` (e.g. `"batch_correct"`) were intentionally NOT renamed — they are opaque identifiers, not user-facing names

---

## Summary Statistics

| Category | Count |
|----------|-------|
| R source files | 17 |
| R source lines | ~7,500 |
| Exported functions | 90 |
| Test files | 16 |
| Test lines | ~4,700 |
| Documentation pages (.qmd) | 21 |
| Vignettes | 8 |
| CI/CD workflows | 1 |

### Notes

- All wrapper functions guard against missing packages via `requireNamespace()`
- Tests use `skip_if_not_installed()` for external dependencies
- roxygen2 documentation is embedded in source files; `man/` Rd file content and filenames are updated but `devtools::document()` must be re-run to fully regenerate from source

---

## What's Been Done (Complete)

| Area | Status | Details |
|------|--------|---------|
| Core pipeline (5 steps) | ✅ | unmix → gate → QC → batch correct → cluster |
| Composable pipeline (S7) | ✅ | Step/Pipeline classes, 10 convenience constructors |
| Scale transforms | ✅ | Logicle + linear quantile (ported from CytoPipeline) |
| Gating utilities | ✅ | Singlet gate, doublet removal, debris, margins |
| LLM assistant | ✅ | ellmer-based interactive pipeline builder |
| Bug fixes (4) | ✅ | Margin channels, pipeline gating, EMD, cofactor |
| Integration tests | ✅ | All 16 modules have integration tests |
| Dimensionality reduction | ✅ | UMAP + PCA with auto-subsampling |
| Gating templates (5) | ✅ | lymphocyte, myeloid, NK, Treg, full PBMC |
| Cell type annotation | ✅ | Cosine-similarity scoring + manual mapping + PBMC reference |
| Differential analysis | ✅ | edgeR DA + limma DE + volcano/boxplot/abundance plots |
| Function rename overhaul | ✅ | All 90 functions follow `sw_{module}_{verb}()` hierarchy |
| SSM diagnostics | ✅ | Spillover spreading matrix + per-channel QC |
| Documentation site | ✅ | Quarto website, 2 vignettes, 11 reference pages |
| Cell type annotation | ✅ | Cosine similarity vs 35-population PBMC reference |
| Differential analysis | ✅ | edgeR/voom DA + limma DE with volcano/boxplot/abundance plots |

## What Still Needs To Be Done

| Area | Priority | Effort | Details |
|------|----------|--------|---------|
| **Data export** | MEDIUM | Low-Medium | CSV, H5AD, FCS, Seurat, SCE interop |
| **QC report generation** | MEDIUM | Medium | Automated HTML pipeline summary reports |
| **Panel design helper** | MEDIUM | High | Spectral overlap checking, cofactor suggestion |
| **AF fingerprinting** | LOW | High | Tissue autofluorescence characterization |
| **Batch monitoring** | LOW | Medium | Levy-Jennings charts, drift detection |
| **Reference docs** | MEDIUM | Low | .qmd pages for dimred + unmix_diagnostics |
| **Composable pipeline vignette** | MEDIUM | Low | Tutorial for S7 pipeline framework |
| **Annotation vignette** | MEDIUM | Low | Tutorial for cluster annotation + DA/DE workflow |
| **R CMD check** | HIGH | Low | Requires R environment to run |
| **man/ pages** | HIGH | Low | `devtools::document()` needed |
