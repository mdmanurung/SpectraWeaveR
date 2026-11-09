#!/usr/bin/env bash
# Purpose: run one experiment script inside the project's pixi R environment.
# Run:     sbatch experiments/slurm/run_script.sbatch <script.R> <run_dir> [extra args]
#          e.g. sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/01_mendeley-7color/run.R experiments/runs/2026-11-09_mendeley-7color
# Inputs:  the R script and a run directory that does not exist yet
# Outputs: <run_dir>/ (written by the script) and slurm-<jobid>.out in the submit directory
#SBATCH --job-name=sw_beta
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail

script=${1:?usage: sbatch run_script.sbatch <script.R> <run_dir>}
run_dir=${2:?usage: sbatch run_script.sbatch <script.R> <run_dir>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
# a new run needs a new folder; figure-only reruns set ALLOW_EXISTING_RUN=1 (the script itself refuses to overwrite)
[ ! -e "$run_dir" ] || [ "${ALLOW_EXISTING_RUN:-0}" = 1 ] || { echo "run dir exists: $run_dir" >&2; exit 1; }

# threads for AutoSpectral and BLAS follow the job allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
pixi run Rscript "$script" "$run_dir" "${@:3}"
