# SpectraWeaveR

<!-- badges: start -->
[![R-CMD-check](https://github.com/mdmanurung/SpectraWeaveR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdmanurung/SpectraWeaveR/actions/workflows/R-CMD-check.yaml)
[![R-universe status](https://mdmanurung.r-universe.dev/badges/SpectraWeaveR)](https://mdmanurung.r-universe.dev/SpectraWeaveR)
<!-- badges: end -->

End-to-end R pipeline for 40-color spectral flow cytometry analysis.

## Overview

SpectraWeaveR orchestrates five established Bioconductor/GitHub packages into a single reproducible pipeline for spectral flow cytometry data — from spectral unmixing through clustering and visualization. It also provides CytoPipeline-ported preprocessing utilities and a composable S7-based pipeline framework.

### Core Pipeline

| Step | Function | Package |
|------|----------|---------|
| 1. Unmixing | `sw_unmix_pipeline()`, `sw_unmix()`, `sw_load_unmixed()` | AutoSpectral / flowCore |
| 2. Gating | `sw_gate()`, `sw_build_gating_template()` | openCyto / flowWorkspace |
| 3. Signal QC | `sw_signal_qc()`, `sw_signal_qc_batch()` | PeacoQC |
| 4. Batch correction | `sw_batch_correct()`, `sw_prepare_for_correction()` | cyCombine |
| 5. Clustering | `sw_cluster()`, `sw_plot_clusters()` | kohonen / FastPG |
| End-to-end | `run_pipeline()` | all of the above |

### Extended Features

| Module | Functions | Description |
|--------|-----------|-------------|
| Batch diagnostics | `sw_detect_batch_effect()`, `sw_evaluate_emd()`, `sw_evaluate_mad()` | Batch effect detection and correction evaluation |
| Batch visualization | `sw_plot_batch_densities()`, `sw_plot_batch_dimred()` | Density and dimensionality reduction plots |
| Modular batch correction | `sw_normalize()` → `sw_create_som()` → `sw_correct_data()` | Fine-grained batch correction control |
| Scale transforms | `sw_estimate_scale_transforms()`, `sw_apply_scale_transforms()` | Logicle / linear transforms (CytoPipeline) |
| Gating utilities | `sw_singlets_gate()`, `sw_remove_doublets()`, `sw_remove_debris_gate()` | Standalone preprocessing gates (CytoPipeline) |
| Channel classification | `sw_are_signal_cols()`, `sw_are_fluor_cols()` | Channel type detection |
| Aggregation | `sw_aggregate_and_sample()` | Pool and subsample events |
| Audit trail | `sw_collect_events_retained()` | Track events across pipeline steps |
| Composable pipeline | `sw_step()`, `sw_pipeline()`, `sw_pipeline_run()` | S7-based step/pipeline framework |
| LLM assistant | `sw_assistant()`, `sw_quick_pipeline()` | LLM-powered pipeline builder (ellmer) |

## Installation

### From R-universe (recommended)

```r
install.packages("SpectraWeaveR",
  repos = c("https://mdmanurung.r-universe.dev",
            "https://bioc.r-universe.dev",
            "https://cloud.r-project.org"))
```

### From GitHub

```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowWorkspace", "openCyto", "PeacoQC"))

# Install kohonen from CRAN
install.packages("kohonen")

# Install cyCombine from GitHub
remotes::install_github("biosurf/cyCombine")

# Install SpectraWeaveR
remotes::install_github("mdmanurung/SpectraWeaveR")
```

Or use the included setup script to install all dependencies at once:

```r
source(system.file("scripts/install_dependencies.R", package = "SpectraWeaveR"))
```

### Optional Dependencies

```r
# ellmer for LLM-powered pipeline builder assistant
install.packages("ellmer")

# AutoSpectral for spectral unmixing (required only for Step 1)
remotes::install_github("DrCytometer/AutoSpectral")
remotes::install_github("DrCytometer/AutoSpectralRcpp")  # faster

# FastPG for graph-based clustering (alternative to kohonen SOM)
remotes::install_github("SamGG/FastPG")

# S7 for composable pipeline framework
install.packages("S7")

# Visualization extras
install.packages(c("pheatmap", "uwot"))
BiocManager::install("ggcyto")
```

## Quick Start

### Full Pipeline

```r
library(SpectraWeaveR)

# Define sample metadata
sample_meta <- data.frame(
  file     = c("sample1.fcs", "sample2.fcs", "sample3.fcs", "sample4.fcs"),
  sample   = c("S1", "S2", "S3", "S4"),
  batch    = c("B1", "B1", "B2", "B2"),
  condition = c("Ctrl", "Stim", "Ctrl", "Stim")
)

# Define markers
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14", "CD16",
             "CD45RA", "CD45RO", "HLA-DR")
lineage_markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")

# Run the full pipeline
results <- run_pipeline(
  fcs_dir         = "path/to/unmixed_fcs/",
  sample_meta     = sample_meta,
  markers         = markers,
  lineage_markers = lineage_markers,
  cofactor        = 6000,
  n_metaclusters  = 20,
  seed            = 42
)

# Access results
results$cluster_assignments  # per-cell metacluster labels
results$cluster_mfis         # MFI per metacluster
results$corrected            # batch-corrected expression data
```

### Step-by-Step Usage

```r
library(SpectraWeaveR)

# 1. Load unmixed FCS files
fs <- sw_load_unmixed("path/to/unmixed_fcs/")

# 2. Remove margin events (before transformation!)
ff_list <- lapply(seq_along(fs), function(i) sw_remove_margins(fs[[i]]))
names(ff_list) <- flowCore::sampleNames(fs)

# 3. Signal QC
qc_results <- sw_signal_qc_batch(ff_list)
qc_summary <- sw_qc_summary(qc_results, threshold = 30)
ff_clean <- qc_results$cleaned

# 4. Batch correction
uncorrected <- sw_prepare_for_correction(
  ff_clean, sample_meta, markers, cofactor = 6000
)
corrected <- sw_batch_correct(uncorrected, markers, covar = "condition")

# 5. Clustering
result <- sw_cluster(corrected, lineage_markers, n_metaclusters = 20)
assignments <- sw_get_cluster_assignments(result)
mfis <- sw_cluster_mfis(result)
sw_plot_clusters(result, "clusters.pdf")
```

### Composable Pipeline

```r
library(SpectraWeaveR)

# Build a custom pipeline with S7 steps
pip <- sw_pipeline("my_analysis", steps = list(
  sw_step("qc", sw_signal_qc, list(IT_limit = 0.55, MAD = 6)),
  sw_step("transform", function(qc) qc$FinalFF)
))

# Run it
result <- sw_pipeline_run(pip, input = my_flowframe)

# Visualize the pipeline
sw_plot_pipeline(pip)
```

### LLM-Powered Pipeline Builder

SpectraWeaveR includes an optional LLM assistant (powered by
[ellmer](https://ellmer.tidyverse.org/)) that guides you through pipeline
configuration via natural language conversation:

```r
library(SpectraWeaveR)

# Interactive session — the assistant asks questions and inspects your files
sw_assistant()

# Use a specific provider
sw_assistant(provider = "anthropic")

# Local model via Ollama (fully private, no data sent externally)
sw_assistant(provider = "ollama", model = "llama3")

# One-shot code generation from a description
code <- sw_quick_pipeline(
  "I have 30 Aurora 5L FCS files in /data/fcs/, metadata in /data/meta.csv
   with columns filename, patient_id, batch, treatment. Markers: CD3, CD4,
   CD8, CD19, CD56. Cluster on CD3, CD4, CD8."
)

# Generate code from a config list (no LLM needed)
sw_generate_pipeline_code(config, output = "file", path = "my_pipeline.R")
```

Install the optional dependency:

```r
install.packages("ellmer")
```

### Spectral Unmixing (from raw data)

```r
library(SpectraWeaveR)

# One-call unmixing pipeline
result <- sw_unmix_pipeline(
  control_dir  = "path/to/controls/",
  sample_input = "path/to/samples/",
  unstained_fcs = "path/to/unstained.fcs",
  cytometer    = "aurora",
  method       = "AutoSpectral"
)

# Or step by step:
setup    <- sw_autospectral_setup("controls/", cytometer = "aurora")
controls <- sw_prepare_controls(setup)
af       <- sw_extract_af_spectra("unstained.fcs", setup, controls$spectra)
variants <- sw_extract_spectral_variants(setup, controls$spectra, af)
unmixed  <- sw_unmix("samples/", controls$spectra, setup,
                     controls$flow_control, af_spectra = af,
                     spectra_variants = variants)
```

## Key Design Decisions

- **`truncate_max_range = FALSE`** is enforced when reading FCS files — Aurora fluorescence intensities (~4×10⁶) exceed flowCore's default range cutoffs.
- **arcsinh cofactor = 6000** is used for spectral flow cytometry (not 5 as in CyTOF).
- **RemoveMargins before transformation** — margin removal operates on raw signal intensities.
- **kohonen SOM** (not FlowSOM) for clustering — simpler dependency chain with hierarchical metaclustering via ward.D2.
- **AutoSpectral in Suggests** — unmixing is optional; users can start from pre-unmixed files.
- Format conversion utilities (`sw_flowframe_to_tibble()`, `sw_tibble_to_flowframe()`) bridge the incompatible data structures required by each tool.

## Documentation

Full documentation is available at the [Quarto website](https://mdmanurung.github.io/SpectraWeaveR/), including:
- Installation guide
- Vignettes for spectral unmixing, batch correction, and LLM assistant
- Complete API reference for all exported functions

## License

MIT
