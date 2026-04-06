#!/usr/bin/env Rscript
#
# install_dependencies.R
# ──────────────────────────────────────────────────────────────────────
# Install all SpectraWeaveR dependencies (required + optional) using pak.
#
# Usage:
#   Rscript inst/scripts/install_dependencies.R            # install everything
#   Rscript inst/scripts/install_dependencies.R --required  # required only
#
# pak resolves CRAN, Bioconductor, and GitHub packages in a single call,
# so no separate BiocManager setup is needed.
# ──────────────────────────────────────────────────────────────────────

# --- Install pak itself if not available --------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/")
}

# --- Parse command-line flags -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
install_optional <- !("--required" %in% args)

# --- Required dependencies (Imports) -----------------------------------
required_pkgs <- c(
  # Bioconductor
  "bioc::flowCore",        # FCS I/O, flowFrame/flowSet classes
  "bioc::flowWorkspace",   # GatingSet, cytoset
  "bioc::openCyto",        # CSV-driven hierarchical gating
  "bioc::PeacoQC",         # Isolation Tree + MAD signal QC

  # CRAN

  "kohonen",               # Self-organizing map clustering
  "dplyr",                 # Data manipulation
  "tibble"                 # Tibble data frames
)

# cyCombine is only on GitHub
required_github <- c(
  "biosurf/cyCombine"      # ComBat batch correction
)

# --- Optional dependencies (Suggests) ----------------------------------
optional_pkgs <- c(
  # GitHub
  "DrCytometer/AutoSpectral",     # Spectral unmixing engine
  "DrCytometer/AutoSpectralRcpp", # C++ acceleration for unmixing
  "SamGG/FastPG",                 # Graph-based clustering

  # Bioconductor
  "bioc::ggcyto",                 # Gate visualization

  # CRAN
  "S7",                           # Composable pipeline (R7 OOP)
  "ggplot2",                      # Plotting
  "uwot",                         # UMAP dimensionality reduction
  "pheatmap",                     # Cluster annotation heatmaps
  "testthat",                     # Unit testing
  "knitr",                        # Vignette support
  "rmarkdown",                    # Vignette rendering
  "quarto"                        # Quarto vignette builder
)

# --- Install ------------------------------------------------------------
all_pkgs <- c(required_pkgs, required_github)

if (install_optional) {
  all_pkgs <- c(all_pkgs, optional_pkgs)
  message("Installing all dependencies (required + optional)...")
} else {
  message("Installing required dependencies only...")
}

pak::pkg_install(all_pkgs, ask = FALSE)

message("Done! All requested dependencies are installed.")
