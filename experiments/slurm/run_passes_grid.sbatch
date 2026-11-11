#!/usr/bin/env bash
# Purpose: unmix one n.passes x n.af.passes grid arm (52_passes_grid.R) in the pixi R environment.
# Run:     sbatch --array=1-9 experiments/slurm/run_passes_grid.sbatch
#          task i runs line i of ARMS below.
# Inputs:  the 1.8.8 zip folder under experiments/bmv_test/compare_3D/data
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/p<np>_af<naf>/
#          experiments/slurm/logs/ppg-<arrayjob>_<task>.out
#SBATCH --job-name=sw_ppg
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/ppg-%A_%a.out
set -euo pipefail

ARMS=(p1_af1 p1_af5 p1_af10 p5_af1 p5_af5 p5_af10 p10_af1 p10_af5 p10_af10)

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
arm=${ARMS[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "arm: $arm"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/52_passes_grid.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/52_passes_grid.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$arm"
