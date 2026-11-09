#!/usr/bin/env Rscript
# Purpose: entry point for benchmark 03 with cells controls; the work is in _common/run_unmix.R.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/03_omip069v2-40color/run_cells.R experiments/runs/<date>_omip069v2-cells
CONFIG_PATH <- "experiments/benchmarks/03_omip069v2-40color/config_cells.R"
source("experiments/benchmarks/_common/run_unmix.R")
