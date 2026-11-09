#!/usr/bin/env bash
# Purpose: run the whole synthetic benchmark for a set of seeds and methods: generate, run, score, summarise.
# Run:     bash experiments/benchmarks/05_cytether-synthetic/run_all.sh <sim_dir> <run_dir> <reference_method> "<methods>" <seed> [<seed> ...]
#          e.g. ... cytether_sim_dev experiments/runs/20261109_synthetic_dev paired_shift "none naive_batch_center paired_shift heb" 441 442 443
# Inputs:  bundle at $SW_DATA_DIR/cytether_bundle_v1, venv at $SW_DATA_DIR/venv_cytether, pixi R env, competitor library
# Outputs (in <run_dir>):
#   results/seed_<n>/<method>/<stage>/   corrected matrices and status (run_method.R)
#   scores.csv                           evaluate.py
#   summary/                             summarize.py: primary_by_seed.csv, primary_summary.csv, guardrails.csv, primary.png
# The run dir is never overwritten. Sanity candidates: add "oracle" (must score 0) and "reversed_shift" (must score worse than none).
set -euo pipefail
sim_dir=$1
run_dir=$2
reference=$3
methods=$4
shift 4
data=${SW_DATA_DIR:-/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata}
py=$data/venv_cytether/bin/python
bench=experiments/benchmarks/05_cytether-synthetic

[ ! -e "$run_dir" ] || { echo "$run_dir exists; run dirs are never overwritten" >&2; exit 1; }
mkdir -p "$run_dir/results"

for method in $methods; do
  case $method in oracle) continue ;; esac
  bash "$bench/run_seeds.sh" "$sim_dir" "$run_dir/results" "$method" "$@"
done
$py "$bench/evaluate.py" "$sim_dir" "$run_dir/results" "$run_dir/scores.csv" $methods
$py "$bench/summarize.py" "$run_dir/scores.csv" "$run_dir/summary" "$reference"
