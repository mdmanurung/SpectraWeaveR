# SpectraWeaveR Comprehensive Codebase Evaluation

**Date**: 2026-04-06
**Scope**: Technical code quality, bioinformatics correctness, and feature proposals

---

## Part 1: Technical Evaluation

### Overall Code Quality: 8.5/10

#### Strengths

- **Excellent input validation** (~256 stop/warning calls across 11 files). Every exported function validates argument types, ranges, and column existence upfront. Consistent `call. = FALSE` pattern avoids confusing internal call stack pollution in error messages.
- **Well-organized code**: Clear `sw_*` prefix naming convention for all 80+ exported functions, internal helpers prefixed with `.` (e.g., `.check_cycombine()`, `.extract_marker_matrix()`), modular file structure grouping functions logically.
- **Comprehensive roxygen documentation**: All exported functions have `@param`, `@return`, `@details` sections with cross-references via `\code{\link{...}}`.
- **Dual-path workflows**: Both monolithic functions (`sw_pipeline_run_all()`, `sw_correct_run()`) and modular step-by-step alternatives, plus S7 composable pipeline framework for advanced users.
- **Modern R practices**: S7/R7 classes for composable pipeline, tibble-based data flow, proper `on.exit()` cleanup, explicit `.Random.seed` save/restore in `sw_io_subsample()`.
- **Robust dependency management**: `requireNamespace(..., quietly = TRUE)` checks with helpful installation instructions for all optional packages (AutoSpectral, FastPG, S7, ellmer).
- **Reproducibility**: Explicit `set.seed()` calls with configurable seed parameters throughout batch correction and clustering.

#### Bugs and Issues Found

| Issue | File:Line | Severity | Details |
|-------|-----------|----------|---------|
| **Margin removal uses ALL channels** | `unmix.R:1136-1139` | **Medium** | When `channel_specs = NULL`, passes `seq_len(ncol(ff))` to `PeacoQC::RemoveMargins()`, which includes Time, Original_ID, and scatter channels. Should filter to signal-only channels via `sw_channel_is_signal()`. |
| **Pipeline gating bypasses margin removal** | `pipeline.R:173-176` | **Medium** | When gating is enabled, `sw_gate_run()` re-reads FCS files from disk instead of using the already margin-removed `ff_list`. This means margin removal in Step 2 is effectively discarded. |
| **EMD evaluation returns empty tibble** | `batch_correct.R:306` | **Low** | `sw_correct_evaluate_quick()` always returns `emd = tibble::tibble()` — incomplete implementation. MAD-based evaluation works correctly, but the EMD component is a stub. |
| **Hardcoded cofactor in gating** | `gate.R:132` | **Low** | Cofactor 6000 is hardcoded inside `sw_gate_run()` function body rather than exposed as a function parameter. Users cannot override without providing a custom `transform_func`. |
| **PDF device error risk** | `cluster.R:306-320` | **Low** | If `grDevices::pdf()` fails (e.g., write-protected directory), the `dev.off()` call in the `finally` block could error on the wrong graphics device. |

### Test Suite Assessment

**310 tests across 11 files, 3371 lines of test code.**

| Module | Input Validation | Integration Tests | Overall Grade |
|--------|-----------------|-------------------|---------------|
| utils | A | A (real conversions) | **A** |
| unmix | A | B+ (needs AutoSpectral) | **B+** |
| batch_correct | A | A (real math, arcsinh verified) | **A** |
| cluster | A | A (real SOM + FastPG) | **A** |
| composable_pipeline | A | A (end-to-end execution) | **A** |
| gating_utils | A | A (real polygon gates) | **A** |
| gate | A | F (no workflow execution) | **C** |
| qc | B+ | F (parameter validation only) | **D** |
| transforms | B+ | F (no transform execution) | **D** |
| pipeline | A | F (structural/input only) | **D** |

**Critical test gaps**:
- `sw_pipeline_run_all()` has zero integration testing — only validates input arguments.
- `sw_qc_run()` never executes PeacoQC; only checks parameter ranges.
- `sw_gate_run()` never applies gates to a flowSet; only validates template CSV structure.
- `sw_transform_estimate()` / `sw_transform_apply()` never execute transformations.

---

## Part 2: Bioinformatics Evaluation

### Pipeline Order: CORRECT

The pipeline follows the validated order from den Braanker et al. (2021) and Quintelier et al. (2021, Nature Protocols):

```
Unmixing -> Load -> RemoveMargins -> [Transform -> Gating] -> Signal QC -> Batch Correction -> Clustering
```

Key validations:
- **Margin removal before transformation**: Correct — `PeacoQC::RemoveMargins` relies on raw fluorescence intensities to detect saturated events near detector limits.
- **QC after margin removal**: Correct — signal stability assessment runs on cleaned data.
- **Batch correction before clustering**: Correct — removes technical batch variation before biological analysis.

### Default Parameters: APPROPRIATE for Spectral Flow Cytometry

| Parameter | Default Value | Assessment |
|-----------|--------------|------------|
| arcsinh cofactor | 6000 | **Correct** for Aurora (~4M range). CyTOF uses 5; spectral flow requires ~6000. |
| `truncate_max_range` | FALSE | **Essential** for Aurora FCS files where fluorescence intensities exceed flowCore's default range cutoffs. |
| SOM grid (batch correction) | 8 x 8 | **Standard** for cyCombine workflows (64 nodes). |
| SOM grid (clustering) | 10 x 10 | **Appropriate** — larger grid allows finer population discrimination. |
| Metaclusters | 20 | **Reasonable** for 40-color panels (major + minor immune populations). |
| Hierarchical linkage | ward.D2 | **Standard** and robust for metaclustering. |
| QC IT_limit | 0.55 | **Conservative** Isolation Tree threshold (typical range 0.1–0.6). |
| QC MAD threshold | 6 | **Conservative** (typical range 4–8). |
| Normalization method | z-score | **Correct** for single-study data; rank-based available for multi-study meta-analysis. |
| Singlet gate geometry | Parallelogram | **More appropriate** for spectral flow than PeacoQC's semi-conic gate. |

### Bioinformatics Concerns

1. **QC on untransformed data when gating skipped**: If no gating template is provided in `sw_pipeline_run_all()`, `sw_qc_run()` receives untransformed data. PeacoQC still works on raw intensities, but this should be documented explicitly so users understand the behavior.

2. **Only "lymphocyte" gating template**: `sw_gate_template()` generates a single template type (nonDebris -> singlets -> lymphocytes via flowClust). Spectral panels frequently target myeloid, dendritic, or innate immune populations that need different gating hierarchies.

3. **Pipeline gating re-reads from disk** (see bug above): When gating is enabled, the margin-removed flowFrames from Step 2 are discarded because `sw_gate_run()` reads fresh FCS files. This means margin removal is bypassed for gated data — a bioinformatics concern since saturated events could persist through gating into downstream analysis.

---

## Part 3: Proposed New Features

### HIGH Priority

#### 1. Unmixing Quality Diagnostics
- **New file**: `R/unmix_diagnostics.R`
- **Functions**: `sw_unmix_quality()`, `sw_unmix_spillover_matrix()`, `sw_plot_spillover_matrix()`
- **Rationale**: The spillover spreading matrix (SSM) is the standard metric for evaluating spectral unmixing quality (Nguyen et al. 2013). With 40+ color panels, fluorochrome choice directly impacts data quality, yet the package provides no way to assess unmixing performance. This would calculate per-channel spreading coefficients from single-stain controls and visualize them as a heatmap, enabling users to identify problematic fluorochrome combinations.
- **Leverages**: flowCore expression matrices, custom SSM calculation
- **Complexity**: Medium

#### 2. Standalone Dimensionality Reduction Module
- **New file**: `R/dimred.R`
- **Functions**: `sw_run_umap()`, `sw_run_tsne()`, `sw_run_pca()`, `sw_plot_dimred()`
- **Rationale**: UMAP/tSNE are the standard visualization for high-dimensional cytometry data. Currently only available as a side effect in `sw_plot_batch_dimred()`. Users need standalone dimensionality reduction with consistent interface, marker coloring, sample/batch overlays, and population annotation support.
- **Leverages**: uwot (UMAP), Rtsne, stats::prcomp
- **Complexity**: Medium

#### 3. Differential Abundance and Expression Analysis
- **New file**: `R/differential.R`
- **Functions**: `sw_diff_abundance()`, `sw_diff_expression()`, `sw_plot_volcano()`, `sw_plot_boxplots()`
- **Rationale**: After clustering, the natural next question is "which populations differ between conditions?" This is the most common downstream analysis in spectral flow cytometry and currently requires leaving the package entirely. GLM-based abundance testing and per-cluster median expression comparison would complete the analysis workflow.
- **Leverages**: edgeR/diffcyt for GLM-based abundance testing, limma for expression, or simpler Wilcoxon/KS tests
- **Complexity**: Medium-High

#### 4. Automated Cell Type Annotation
- **New file**: `R/annotate.R`
- **Functions**: `sw_annotate_run()`, `sw_plot_annotation()`, `sw_marker_expression_summary()`
- **Rationale**: Manually annotating 20+ metaclusters by inspecting MFI heatmaps is tedious and error-prone. Reference-based annotation using known marker expression profiles (e.g., CD3+CD4+ = T helper cells) would dramatically improve usability. Particularly valuable for spectral panels where 40+ markers enable fine-grained phenotyping that's hard to do manually.
- **Leverages**: Marker-based scoring with positive/negative thresholds per lineage, or integration with scGate-style reference profiles
- **Complexity**: Medium-High

### MEDIUM Priority

#### 5. Export and Interoperability Module
- **New file**: `R/export.R`
- **Functions**: `sw_export_csv()`, `sw_export_h5ad()`, `sw_export_fcs()`, `sw_to_seurat()`, `sw_to_sce()`
- **Rationale**: Enables interop with Python (scanpy/AnnData), Seurat, and SingleCellExperiment ecosystems. Many downstream analyses (trajectory inference, RNA+protein integration) live in these ecosystems.
- **Leverages**: anndata (R), SeuratObject, SingleCellExperiment, flowCore::write.FCS
- **Complexity**: Medium

#### 6. Automated QC Report Generation
- **New file**: `R/report.R`
- **Functions**: `sw_generate_report()`, `sw_interactive_plot()`
- **Rationale**: Automated HTML/PDF reports summarizing each pipeline step (event counts per sample, QC removal rates, batch effect metrics before/after correction, cluster composition tables). Researchers need sharable, reproducible reports for lab meetings and publications.
- **Leverages**: rmarkdown/quarto parameterized templates, plotly for interactive plots, DT for tables
- **Complexity**: Medium

#### 7. Panel Design Helper
- **New file**: `R/panel.R`
- **Functions**: `sw_check_panel_compatibility()`, `sw_suggest_cofactors()`, `sw_fluorochrome_similarity()`
- **Rationale**: Spectral flow's unique challenge is fluorochrome selection. A tool that checks spectral overlap between fluorochromes, suggests per-fluorochrome optimal cofactors, and warns about problematic combinations would be novel. No existing R package provides comprehensive panel design assistance for spectral cytometry.
- **Leverages**: Reference spectra databases, spectral similarity metrics (cosine similarity, spectral angle)
- **Complexity**: High (requires reference spectra database)

#### 8. Complete EMD Implementation
- **Functions**: Fix `sw_correct_evaluate_quick()` to compute actual Earth Mover's Distance values; add `sw_plot_emd_heatmap()`
- **Rationale**: EMD is the gold-standard metric for batch effect quantification in cytometry. The current stub returns an empty tibble. Completing this would provide proper before/after batch correction comparison alongside the existing MAD metrics.
- **Leverages**: Already partially implemented via `sw_correct_emd()` / cyCombine's `evaluate_emd()`
- **Complexity**: Low

### LOW Priority (Novel Features)

#### 9. Fluorochrome Autofluorescence Fingerprinting
- **Functions**: `sw_af_fingerprint()`, `sw_classify_af_populations()`
- **Rationale**: Tissue autofluorescence is a major challenge in spectral flow, especially for tissue-derived samples. A module that characterizes tissue-specific AF spectral signatures and identifies AF-high subpopulations would be novel. Currently only basic AF extraction exists via AutoSpectral.
- **Complexity**: High

#### 10. Longitudinal Batch Monitoring
- **Functions**: `sw_monitor_batch_drift()`, `sw_plot_batch_timeline()`, `sw_levy_jennings()`
- **Rationale**: Core facilities running spectral instruments need to track instrument performance over time using QC beads. Levy-Jennings control charts and automated drift detection would serve this niche and complement the existing batch correction tools.
- **Complexity**: Medium

#### 11. Additional Gating Templates
- **Functions**: Extend `sw_gate_template()` with types: `"myeloid"`, `"dendritic"`, `"nk_cell"`, `"treg"`, `"full_immune"`
- **Rationale**: The single "lymphocyte" template limits usability. Pre-built templates for common immunophenotyping panels would help users get started faster and reduce gating errors.
- **Complexity**: Low per template (but requires domain expertise for each)

---

## Summary

SpectraWeaveR is a well-engineered package with production-quality code, correct bioinformatics methodology, and appropriate defaults for spectral flow cytometry.

### Completed (as of 2026-04-06)

1. **All 4 bugs fixed**: margin removal channel selection, pipeline gating bypass, EMD evaluation stub, cofactor parameterization
2. **Integration tests added** for all 4 previously-untested modules (QC, gating, transforms, pipeline)
3. **Phase 1 features implemented**:
   - Dimensionality reduction module (`R/dimred.R`) — UMAP + PCA
   - Additional gating templates (myeloid, NK, Treg, full PBMC)
   - SSM diagnostics module (`R/unmix_diagnostics.R`) — spillover spreading matrix + per-channel quality
4. See PR #17 for full details

### Still To Do

**High Priority (Phase 2 features)**:
- Automated cell type annotation (`R/annotate.R`)
- Differential abundance & expression analysis (`R/differential.R`)

**Medium Priority**:
- Data export/interoperability (`R/export.R`) — CSV, H5AD, FCS, Seurat, SCE
- Automated QC report generation (`R/report.R`)
- Panel design helper (`R/panel.R`)

**Low Priority (novel/niche)**:
- Fluorochrome autofluorescence fingerprinting
- Longitudinal batch monitoring / Levy-Jennings charts

**Infrastructure**:
- `R CMD check` with 0 errors/warnings
- Generate `man/` pages via `devtools::document()`
- Add reference docs for new modules (dimred, unmix_diagnostics)
- Composable pipeline vignette
- Publish Quarto documentation site
