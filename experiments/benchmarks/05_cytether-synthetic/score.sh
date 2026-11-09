#!/usr/bin/env bash
# Purpose: score and summarise a results tree that already exists (for example after SLURM jobs added methods).
# Run:     bash experiments/benchmarks/05_cytether-synthetic/score.sh <sim_dir> <run_dir> <reference_method> <out_name> "<methods>"
# Inputs:  <run_dir>/results/seed_<n>/<method>/<stage>/ for every method listed
# Outputs: <run_dir>/<out_name>/scores.csv and <run_dir>/<out_name>/summary/ (see run_all.sh); an existing <out_name> is never overwritten
set -euo pipefail
sim_dir=$1
run_dir=$2
reference=$3
out_name=$4
methods=$5
data=${SW_DATA_DIR:-/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata}
py=$data/venv_cytether/bin/python
bench=experiments/benchmarks/05_cytether-synthetic

[ ! -e "$run_dir/$out_name" ] || { echo "$run_dir/$out_name exists; never overwritten" >&2; exit 1; }
mkdir -p "$run_dir/$out_name"
$py "$bench/evaluate.py" "$sim_dir" "$run_dir/results" "$run_dir/$out_name/scores.csv" $methods
$py "$bench/summarize.py" "$run_dir/$out_name/scores.csv" "$run_dir/$out_name/summary" "$reference"
