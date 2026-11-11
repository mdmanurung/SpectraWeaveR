#!/usr/bin/env bash
# Purpose: run one benchmark script in the pixi R environment from a frozen copy that also freezes lib/
#          (figure templates), so scripts that source lib/theme.R etc. work. Copy of run_script.sbatch's freeze
#          logic; the script reads the frozen lib/ path from SW_LIB and the repo root from SW_REPO.
# Run:     sbatch experiments/slurm/run_script_lib.sbatch <script.R> [args passed to the script]
# Inputs:  the R script; experiments/benchmarks/05_bmv-unmix-compare/*.R and lib/
# Outputs: whatever the script writes; slurm-<jobid>.out in the submit directory
#SBATCH --job-name=sw_lib
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail

script=${1:?usage: sbatch run_script_lib.sbatch <script.R> [args]}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export SW_REPO="$PWD"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Frozen copy: Rscript reads lazily, so an edit during the run would break the job.
snap_dir="$(mktemp -d)"
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$snap_dir"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/lib "$snap_dir"/lib
cp "$script" "$snap_dir"/
export SW_LIB="$snap_dir/lib"
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$snap_dir/#" "$snap_dir"/*.R
pixi run Rscript "$snap_dir/$(basename "$script")" "${@:2}"
