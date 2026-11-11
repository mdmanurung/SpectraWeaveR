#!/usr/bin/env bash
# Purpose: run one ablation arm script (48_ablation_swap.R, 48b_ablation_bead.R, 49_ablation_chain.R) in the pixi R
#          environment; one array task per line of an arm list "<script> <arm>".
# Run:     sbatch --array=1-N experiments/slurm/run_ablation_arm.sbatch <arm_list>
# Inputs:  <arm_list>, lines "<script basename> <arm>"; the zip folders under experiments/bmv_test/compare_3D/data
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/
#          experiments/slurm/logs/abl-<arrayjob>_<task>.out
#SBATCH --job-name=sw_abl
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/abl-%A_%a.out
set -euo pipefail

arm_list=${1:?usage: sbatch --array=1-N run_ablation_arm.sbatch <arm_list>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$arm_list")
[ -n "$line" ] || { echo "no arm on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
read -r script arm <<< "$line"
echo "script: $script  arm: $arm"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "experiments/benchmarks/05_bmv-unmix-compare/$script")"
cp "experiments/benchmarks/05_bmv-unmix-compare/$script" "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
# freeze the helper folders too (lib/theme.R and the style/ and lucid/ pointers to it)
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid} "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid}/*.R
pixi run Rscript "$snap" "$arm"
