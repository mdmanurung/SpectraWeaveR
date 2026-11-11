#!/usr/bin/env bash
# Purpose: run one ablation chain arm (49_ablation_chain.R) in the pixi R environment; one array task per arm.
# Run:     sbatch --array=1-6 experiments/slurm/run_ablation_chain.sbatch <arm_list>
#          arm_list holds one arm per line; task i runs line i.
# Inputs:  <arm_list>, the 1.8.8 zip folder under experiments/bmv_test/compare_3D/data
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/
#          experiments/slurm/logs/abl-<arrayjob>_<task>.out
#SBATCH --job-name=sw_abl
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/abl-%A_%a.out
set -euo pipefail

arm_list=${1:?usage: sbatch --array=1-N run_ablation_chain.sbatch <arm_list>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
arm=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$arm_list")
[ -n "$arm" ] || { echo "no arm on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "arm: $arm"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/49_ablation_chain.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/49_ablation_chain.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$arm"
