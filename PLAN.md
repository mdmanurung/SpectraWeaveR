# SpectraWeaveR — Implementation Plan

## Goal

Build a complete, installable R package (`SpectraWeaveR`) that provides an end-to-end pipeline for 40-color spectral flow cytometry analysis. The package wraps five established tools—AutoSpectral, openCyto, PeacoQC, cyCombine, and FlowSOM—into a coherent, format-safe workflow.

---

## Phase 1: Package Skeleton

- [ ] `DESCRIPTION` — package name, version (0.1.0), license (MIT), all Imports/Suggests
- [ ] `LICENSE` — MIT license file
- [ ] `NAMESPACE` — export all `sw_*` and `run_pipeline` functions; import from dependencies
- [ ] `.Rbuildignore` — exclude non-package files (CLAUDE.md, PLAN.md, PROGRESS.md, etc.)
- [ ] `README.md` — update with installation instructions and usage example

---

## Phase 2: Core R Source Files

### `R/utils.R` — Format conversion utilities
These functions are the connective tissue of the pipeline. Every other module depends on them.

- [ ] `sw_read_fcs(files, ...)` — wrapper around `flowCore::read.FCS()` with `truncate_max_range = FALSE` enforced
- [ ] `sw_flowframe_to_tibble(ff, sample, batch, condition)` — `exprs(ff)` → tibble with metadata columns
- [ ] `sw_tibble_to_flowframe(df, markers)` — tibble → `flowCore::flowFrame`
- [ ] `sw_exprs_to_tibble(mat, colnames, sample, batch, condition)` — matrix → tibble
- [ ] `sw_get_fluor_channels(ff_or_cs)` — return channel names excluding FSC/SSC/Time
- [ ] `sw_set_marker_names(ff, marker_map)` — rename channel descriptions to marker names

### `R/unmix.R` — Step 1: Spectral unmixing
- [ ] `sw_unmix_autospectral(control_dir, asp_params, ...)` — wrapper for AutoSpectral workflow
- [ ] `sw_load_unmixed(fcs_dir, pattern, ...)` — load pre-unmixed FCS files from SpectroFlo/AutoSpectral into a `flowSet`
- [ ] `sw_remove_margins(ff, channel_specs)` — wrapper for `PeacoQC::RemoveMargins()` (must run pre-transformation)

### `R/gate.R` — Step 2: Automated gating (openCyto)
- [ ] `sw_build_gating_template(output_file, template_type)` — write a default CSV gating template for a standard lymphocyte panel
- [ ] `sw_gate(fcs_files, gating_template, transform_channels, ...)` — full openCyto gating workflow: load → GatingSet → transform → gt_gating
- [ ] `sw_extract_gated(gs, node)` — extract `flowFrame` list from a GatingSet node

### `R/qc.R` — Step 3: Signal quality control (PeacoQC)
- [ ] `sw_signal_qc(ff, channels, IT_limit, MAD, output_dir, ...)` — wrapper for `PeacoQC::PeacoQC()` with sensible spectral defaults
- [ ] `sw_signal_qc_batch(ff_list, ...)` — apply `sw_signal_qc` over a list of flowFrames; return cleaned list + summary stats
- [ ] `sw_qc_summary(qc_results)` — tabulate cells removed per sample; flag samples >30% removal

### `R/batch_correct.R` — Step 4: Batch correction (cyCombine)
- [ ] `sw_prepare_for_correction(ff_list, sample_meta, markers, cofactor)` — convert flowFrame list → cyCombine tibble with arcsinh transformation (cofactor=6000)
- [ ] `sw_batch_correct(uncorrected, markers, covar, xdim, ydim, norm_method, seed, ...)` — wrapper for `cyCombine::batch_correct()`
- [ ] `sw_evaluate_correction(uncorrected, corrected, markers)` — compute EMD and MAD metrics; return named list with interpretation flags

### `R/cluster.R` — Step 5: Clustering (FlowSOM)
- [ ] `sw_cluster(corrected, lineage_markers, xdim, ydim, n_metaclusters, seed, ...)` — convert tibble → flowFrame → `FlowSOM::FlowSOM()` wrapper
- [ ] `sw_get_cluster_assignments(fsom_result)` — return per-cell metacluster vector
- [ ] `sw_cluster_mfis(fsom_result)` — return median fluorescence intensity table per metacluster
- [ ] `sw_plot_clusters(fsom_result, plot_file)` — save `FlowSOMmary` PDF and star plot
- [ ] `sw_map_new_data(fsom_result, new_ff)` — project new samples onto trained FlowSOM via `NewData()`

### `R/pipeline.R` — End-to-end orchestrator
- [ ] `run_pipeline(fcs_dir, sample_meta, markers, lineage_markers, gating_template, output_dir, ...)` — calls all steps in validated order with progress messages; returns a named list with all intermediate and final results

---

## Phase 3: Documentation (man/)

One `.Rd` file per exported function, generated via roxygen2. Key parameters to document carefully:

- [ ] `sw_read_fcs.Rd` — emphasise `truncate_max_range`
- [ ] `sw_gate.Rd` — document gating template CSV format
- [ ] `sw_signal_qc.Rd` — document IT_limit / MAD defaults and when to adjust
- [ ] `sw_batch_correct.Rd` — document cofactor, covar requirement, EMD/MAD interpretation
- [ ] `sw_cluster.Rd` — document lineage vs functional marker split
- [ ] `run_pipeline.Rd` — full parameter table and example

---

## Phase 4: Tests (tests/testthat/)

Tests use synthetic in-memory data (no real FCS files required).

- [ ] `test-utils.R` — round-trip flowFrame ↔ tibble, channel extraction, marker renaming
- [ ] `test-unmix.R` — `sw_load_unmixed` file pattern matching, `sw_remove_margins` column check
- [ ] `test-gate.R` — `sw_build_gating_template` creates valid CSV; `sw_extract_gated` returns flowFrame
- [ ] `test-qc.R` — `sw_signal_qc` returns `FinalFF` and `GoodCells`; `sw_qc_summary` flags high-removal samples
- [ ] `test-batch_correct.R` — `sw_prepare_for_correction` column names; arcsinh cofactor applied; EMD/MAD metrics returned
- [ ] `test-cluster.R` — `sw_cluster` returns FlowSOM object; `sw_get_cluster_assignments` length matches input rows
- [ ] `test-pipeline.R` — `run_pipeline` returns named list with expected slots; smoke test end-to-end

---

## Phase 5: Final QC

- [ ] `R CMD build && R CMD check --no-manual` passes with 0 errors, 0 warnings
- [ ] All tests pass
- [ ] Parallel validation (code review + CodeQL) clean
- [ ] Update `README.md` with full usage example
- [ ] Update `PROGRESS.md` to mark all items complete
