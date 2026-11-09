#!/usr/bin/env Rscript
# Purpose: entry point for benchmark 03 with the cells spectra set; the work is in _common/run_unmix.R.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/03_omip069v2-40color/run_cells_set.R experiments/runs/<date>_omip069v2-cells
CONFIG_PATH <- "experiments/benchmarks/03_omip069v2-40color/config_cells_set.R"
source("experiments/benchmarks/_common/run_unmix.R")
