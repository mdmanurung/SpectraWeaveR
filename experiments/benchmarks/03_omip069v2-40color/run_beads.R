#!/usr/bin/env Rscript
# Purpose: entry point for benchmark 03 with beads controls; the work is in _common/run_unmix.R.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/03_omip069v2-40color/run_beads.R experiments/runs/<date>_omip069v2-beads
CONFIG_PATH <- "experiments/benchmarks/03_omip069v2-40color/config_beads.R"
source("experiments/benchmarks/_common/run_unmix.R")
