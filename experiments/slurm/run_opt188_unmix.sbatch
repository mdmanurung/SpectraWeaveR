#!/usr/bin/env bash
# Purpose: unmix one opt188 arm from a route's saved objects with 53b_opt188_unmix.R in the pixi R environment.
# Run:     sbatch experiments/slurm/run_opt188_unmix.sbatch <route> <pipeline> [n_passes n_af_passes]
#          route = pregate | legacyspec; pipeline = joint (with passes) | legacy (no passes)
# Inputs:  experiments/runs/2026-11-10_bmv-opt188-smoke/opt188_<route>_neg_same_fsc4e5_p1/ (objects)
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/
#          experiments/slurm/logs/opt188u-<jobid>.out
#SBATCH --job-name=sw_opt188u
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/slurm/logs/opt188u-%j.out
set -euo pipefail

route=${1:?usage: sbatch run_opt188_unmix.sbatch <route> <pipeline> [n_passes n_af_passes]}
pipeline=${2:?usage: sbatch run_opt188_unmix.sbatch <route> <pipeline> [n_passes n_af_passes]}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
echo "route: $route  pipeline: $pipeline  passes: ${3:-none} ${4:-none}"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/53b_opt188_unmix.R"
cp experiments/benchmarks/05_bmv-unmix-compare/53b_opt188_unmix.R "$snap"
pixi run Rscript "$snap" "$@"
