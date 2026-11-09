#!/usr/bin/env bash
# Purpose: run one method over several seeds of the simulated stages.
# Run:     bash experiments/benchmarks/05_cytether-synthetic/run_seeds.sh <sim_dir> <out_dir> <method> <seed> [<seed> ...]
# Inputs:  <sim_dir>/seed_<n>/ with observed_X.csv.gz (from export_observed.py)
# Outputs: <out_dir>/seed_<n>/<method>/<stage>/ (see run_method.R)
set -euo pipefail
sim_dir=$1
out_dir=$2
method=$3
shift 3
for seed in "$@"; do
  R_LIBS=/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/Rlib_competitors \
    /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript \
    experiments/benchmarks/05_cytether-synthetic/run_method.R \
    "$sim_dir/seed_$seed" "$out_dir/seed_$seed" "$method"
done
