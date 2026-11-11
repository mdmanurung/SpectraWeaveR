#!/usr/bin/env bash
# Purpose: run one integration script (experiments/integration/*.R) in the pixi R environment of the main checkout,
#          loading the package from any checkout (a worktree has no .pixi and no data).
# Run:     sbatch experiments/slurm/run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir>
#          e.g. sbatch experiments/slurm/run_integration.sbatch experiments/integration/clean_run_bmv.R \
#                 "$PWD" /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR \
#                 /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/experiments/runs/2026-11-11_clean-run-smoke
# Inputs:  the R script, a package checkout, the checkout that holds experiments/bmv_test, and a new out_dir
# Outputs: <out_dir>/ (written by the script) and experiments/slurm/logs/integration-<jobid>.out in the main checkout
#SBATCH --job-name=sw_integration
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/experiments/slurm/logs/integration-%j.out
set -euo pipefail

script=${1:?usage: sbatch run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir>}
pkg_dir=${2:?usage: sbatch run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir>}
repo_root=${3:?usage: sbatch run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir>}
out_dir=${4:?usage: sbatch run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir>}

[ ! -e "$out_dir" ] || { echo "run dir exists: $out_dir" >&2; exit 1; }

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "$script")"
cp "$pkg_dir/$script" "$snap"
# A non-zero R exit fails the job (set -e); the pixi manifest of the main checkout supplies the environment.
pixi run --manifest-path "$repo_root/pixi.toml" Rscript "$snap" "$pkg_dir" "$repo_root" "$out_dir"
