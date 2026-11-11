#!/usr/bin/env bash
# Purpose: unmix one Alexa Fluor 405 mechanism arm (55_af405_subset_unmix.R); one array task per arm.
# Run:     sbatch --array=1-N experiments/slurm/run_af405_arm.sbatch <arm_list>
#          arm_list holds one arm per line; task i runs line i.
# Inputs:  <arm_list>; experiments/runs/2026-11-10_bmv-af405-mechanism/variant_subsets.rds [54]
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/
#          experiments/slurm/logs/af405-<arrayjob>_<task>.out
#SBATCH --job-name=sw_af405
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=experiments/slurm/logs/af405-%A_%a.out
set -euo pipefail

arm_list=${1:?usage: sbatch --array=1-N run_af405_arm.sbatch <arm_list>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
arm=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$arm_list")
[ -n "$arm" ] || { echo "no arm on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "arm: $arm"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/55_af405_subset_unmix.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/55_af405_subset_unmix.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$arm"
