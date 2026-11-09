#!/usr/bin/env Rscript
# Purpose: entry point for this benchmark; the work is in _common/run_unmix.R.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/01_mendeley-7color/run.R experiments/runs/<date>_<slug>
CONFIG_PATH <- "experiments/benchmarks/01_mendeley-7color/config.R"
source("experiments/benchmarks/_common/run_unmix.R")
