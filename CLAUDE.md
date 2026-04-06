# SpectraWeaveR — Repository Summary

## What This Repository Is

**SpectraWeaveR** is an R package for end-to-end processing of 40-color spectral flow cytometry data. It orchestrates five established Bioconductor/GitHub packages into a single reproducible pipeline, handling all format conversions between their incompatible data structures.

## Repository Structure

```
SpectraWeaveR/
├── DESCRIPTION          # Package metadata, dependencies, and version
├── NAMESPACE            # Exported functions and imported namespaces
├── LICENSE              # MIT license
├── README.md            # Project overview
├── CLAUDE.md            # This file — repo summary for AI assistants
├── PLAN.md              # Full implementation plan
├── PROGRESS.md          # Step-by-step progress log
├── R/
│   ├── unmix.R          # Step 1: Spectral unmixing (AutoSpectral / SpectroFlo)
│   ├── gate.R           # Step 2: Automated gating (openCyto)
│   ├── qc.R             # Step 3: Signal stability QC (PeacoQC)
│   ├── batch_correct.R  # Step 4: Batch correction (cyCombine)
│   ├── cluster.R        # Step 5: Clustering and visualization (FlowSOM)
│   ├── utils.R          # Format conversion utilities
│   └── pipeline.R       # Full end-to-end pipeline orchestrator
├── man/                 # Roxygen2-generated .Rd documentation
└── tests/
    └── testthat/        # Unit tests (testthat edition 3)
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
| 3. Signal QC | `sw_signal_qc()` | PeacoQC | flowFrame (gated) → flowFrame (cleaned) |
| 4. Batch correction | `sw_batch_correct()` | cyCombine | flowSet → tibble (corrected) |
| 5. Clustering | `sw_cluster()`, `sw_plot_clusters()` | FlowSOM | tibble/matrix → FlowSOM object |
| End-to-end | `run_pipeline()` | all of the above | FCS files → annotated clusters |

## Key Design Constraints

- **`truncate_max_range = FALSE`** must be set when reading Aurora FCS files via `flowCore::read.FCS()` — Aurora fluorescence intensities (~4×10⁶) exceed flowCore's default range cutoffs.
- **arcsinh cofactor = 6000** for spectral flow (not 5 as in CyTOF).
- **Marker naming consistency**: FCS channel names (e.g. `BV421-A`) must be mapped to marker names (e.g. `CD3`) before passing data to cyCombine.
- **Format bridge functions** (`sw_flowframe_to_tibble()`, `sw_tibble_to_flowframe()`, `sw_exprs_to_tibble()`) handle all inter-tool conversions.

## Dependencies

### Required (Imports)
- `flowCore` (≥ 2.6.0) — FCS I/O, flowFrame/flowSet classes
- `flowWorkspace` (≥ 4.6.0) — GatingSet, cytoset
- `openCyto` (≥ 2.6.0) — CSV-driven hierarchical gating
- `PeacoQC` (≥ 1.4.0) — Isolation Tree + MAD signal QC
- `FlowSOM` (≥ 2.6.0) — SOM clustering + MST visualization
- `cyCombine` (≥ 0.2.7) — ComBat batch correction on SOM clusters
- `dplyr`, `tibble` — tibble manipulation

### Optional (Suggests)
- `ggcyto` — gate visualization
- `uwot` — UMAP dimensionality reduction
- `pheatmap` — cluster annotation heatmaps

## Validated Pipeline Order (den Braanker 2021, Quintelier 2021)

```
unmixing → import → RemoveMargins → transformation →
signal QC → gating → batch correction → clustering → visualization
```

## Build & Test Commands

```r
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
