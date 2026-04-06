# SpectraWeaveR — Repository Summary

## What This Repository Is

**SpectraWeaveR** is an R package for end-to-end processing of 40-color spectral flow cytometry data. It orchestrates five established Bioconductor/GitHub packages into a single reproducible pipeline, handling all format conversions between their incompatible data structures. Additionally, it provides CytoPipeline-ported utility functions and a composable S7-based pipeline framework.

## Repository Structure

```
SpectraWeaveR/
├── DESCRIPTION          # Package metadata, dependencies, and version
├── NAMESPACE            # 69 exported functions and imported namespaces
├── LICENSE              # MIT license
├── README.md            # Project overview, installation, usage examples
├── CLAUDE.md            # This file — repo summary for AI assistants
├── PLAN.md              # Full implementation plan with completion status
├── PROGRESS.md          # Session-by-session progress log
├── R/
│   ├── utils.R              # Format conversions + channel classification + aggregation
│   ├── unmix.R              # Step 1: Spectral unmixing (AutoSpectral)
│   ├── gate.R               # Step 2: Automated gating (openCyto)
│   ├── qc.R                 # Step 3: Signal stability QC (PeacoQC)
│   ├── batch_correct.R      # Step 4: Batch correction (cyCombine)
│   ├── cluster.R            # Step 5: Clustering (kohonen SOM / FastPG)
│   ├── pipeline.R           # Full end-to-end pipeline orchestrator
│   ├── composable_pipeline.R # S7-based composable step/pipeline framework
│   ├── transforms.R         # Scale transforms (ported from CytoPipeline)
│   └── gating_utils.R       # Standalone gating utilities (ported from CytoPipeline)
├── tests/
│   └── testthat/        # 10 test files, ~2,871 lines (testthat edition 3)
├── inst/
│   └── scripts/
│       └── install_dependencies.R  # pak-based setup script
├── vignettes/
│   ├── spectral-unmixing.qmd       # Unmixing vignette
│   └── batch-correction.qmd       # Batch correction vignette
├── reference/           # 11 API reference .qmd pages
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
| 1. Unmixing | `sw_unmix_pipeline()`, `sw_unmix()`, `sw_load_unmixed()` | AutoSpectral / flowCore | FCS files → unmixed FCS / flowSet |
| 1a. Setup | `sw_autospectral_setup()` | AutoSpectral | control dir → sw_setup |
| 1b. Controls | `sw_prepare_controls()` | AutoSpectral | sw_setup → spectra + flow.control |
| 1c. AF Spectra | `sw_extract_af_spectra()` | AutoSpectral | unstained FCS → AF spectra matrix |
| 1d. Variants | `sw_extract_spectral_variants()` | AutoSpectral | controls → spectral variants |
| 2. Gating | `sw_gate()`, `sw_build_gating_template()` | openCyto / flowWorkspace | flowSet → GatingSet |
| 3. Signal QC | `sw_signal_qc()`, `sw_signal_qc_batch()` | PeacoQC | flowFrame (gated) → flowFrame (cleaned) |
| 4. Batch correction | `sw_batch_correct()` | cyCombine | flowSet → tibble (corrected) |
| 4a. Modular | `sw_normalize()` → `sw_create_som()` → `sw_correct_data()` | cyCombine | step-by-step control |
| 4b. Diagnostics | `sw_detect_batch_effect()`, `sw_compute_emd()`, `sw_evaluate_emd()`, `sw_evaluate_mad()` | cyCombine | batch effect assessment |
| 4c. Visualization | `sw_plot_batch_densities()`, `sw_plot_batch_dimred()` | cyCombine / ggplot2 | density/UMAP/PCA plots |
| 5. Clustering | `sw_cluster()`, `sw_plot_clusters()`, `sw_predict_clusters()` | kohonen / FastPG | tibble/matrix → sw_cluster_result |
| End-to-end | `run_pipeline()` | all of the above | FCS files → annotated clusters |

## Extended Modules

### Composable Pipeline (`R/composable_pipeline.R`)
S7/R7-based step-and-pipeline framework inspired by CytoPipeline:
- `sw_step()` / `execute_step()` — define and run processing steps
- `sw_pipeline()` / `sw_pipeline_run()` — create and execute pipelines
- Pipeline manipulation: `add`, `remove`, `replace`, `concat`, `show`, `plot`
- 8 convenience step constructors for common operations

### Scale Transforms (`R/transforms.R`)
Ported from CytoPipeline (UCLouvain-CBIO):
- `sw_estimate_scale_transforms()` — logicle / linear quantile per-channel estimation
- `sw_apply_scale_transforms()` — apply pre-computed transformList

### Gating Utilities (`R/gating_utils.R`)
Ported from CytoPipeline / CytoPipelineUtils (UCLouvain-CBIO):
- `sw_singlets_gate()` — parallelogram singlet gate
- `sw_remove_doublets()` / `sw_remove_doublets_peacoqc()` — doublet removal
- `sw_remove_debris_gate()` — manual polygon debris gate
- `sw_remove_margins_peacoqc()` — enhanced margin removal

### Utilities (`R/utils.R`)
- Format bridges: `sw_flowframe_to_tibble()`, `sw_tibble_to_flowframe()`, `sw_exprs_to_tibble()`
- Channel classification: `sw_are_signal_cols()`, `sw_are_fluor_cols()`
- Aggregation: `sw_aggregate_and_sample()` — pool and subsample events
- Audit trail: `sw_collect_events_retained()` — event count tracking

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
- `AutoSpectral` (≥ 1.5.0) — spectral unmixing (GitHub)
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
