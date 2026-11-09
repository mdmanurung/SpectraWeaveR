#!/usr/bin/env Rscript
# Purpose: entry point for this benchmark; the work is in _common/run_unmix.R.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/02_omip109-45color/run.R experiments/runs/<date>_<slug>
CONFIG_PATH <- "experiments/benchmarks/02_omip109-45color/config.R"
source("experiments/benchmarks/_common/run_unmix.R")
