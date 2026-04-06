# SpectraWeaveR

<!-- badges: start -->
[![R-CMD-check](https://github.com/mdmanurung/SpectraWeaveR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdmanurung/SpectraWeaveR/actions/workflows/R-CMD-check.yaml)
[![R-universe status](https://mdmanurung.r-universe.dev/badges/SpectraWeaveR)](https://mdmanurung.r-universe.dev/SpectraWeaveR)
<!-- badges: end -->

End-to-end R pipeline for 40-color spectral flow cytometry analysis.

## Overview

SpectraWeaveR orchestrates five established Bioconductor/GitHub packages into a single reproducible pipeline for spectral flow cytometry data — from spectral unmixing through clustering and visualization.

| Step | Function | Package |
|------|----------|---------|
| 1. Unmixing | `sw_unmix_pipeline()`, `sw_unmix()`, `sw_load_unmixed()` | AutoSpectral / flowCore |
| 2. Gating | `sw_gate()`, `sw_build_gating_template()` | openCyto / flowWorkspace |
| 3. Signal QC | `sw_signal_qc()`, `sw_signal_qc_batch()` | PeacoQC |
| 4. Batch correction | `sw_batch_correct()`, `sw_prepare_for_correction()` | cyCombine |
| 5. Clustering | `sw_cluster()`, `sw_plot_clusters()` | kohonen / FastPG |
| End-to-end | `run_pipeline()` | all of the above |

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
BiocManager::install(c("flowCore", "flowWorkspace", "openCyto",
                       "PeacoQC"))

# Install kohonen from CRAN
install.packages("kohonen")

# Install cyCombine from GitHub
remotes::install_github("biosurf/cyCombine")

# Install SpectraWeaveR
remotes::install_github("mdmanurung/SpectraWeaveR")
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

## Key Design Decisions

- **`truncate_max_range = FALSE`** is enforced when reading FCS files — Aurora fluorescence intensities (~4×10⁶) exceed flowCore's default range cutoffs.
- **arcsinh cofactor = 6000** is used for spectral flow cytometry (not 5 as in CyTOF).
- Format conversion utilities (`sw_flowframe_to_tibble()`, `sw_tibble_to_flowframe()`) bridge the incompatible data structures required by each tool.

## R-universe

SpectraWeaveR is available on [R-universe](https://mdmanurung.r-universe.dev/SpectraWeaveR).

To register this package on R-universe, go to <https://r-universe.dev/add> and follow the prompts. This creates a `universe` registry repo with a `packages.json` pointing to this repository.

## License

MIT
