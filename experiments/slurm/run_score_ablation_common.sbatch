#!/usr/bin/env bash
# Purpose: score ablation arms on the common SpectroFlo gates (51c_score_ablation_common.R), reading only.
# Run:     sbatch experiments/slurm/run_score_ablation_common.sbatch <run_dir> <arm1,arm2,...>
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/ (batch 0)
# Outputs: <run_dir>/ (scores, gap_closed, labels, figure); experiments/slurm/logs/abl-score-common-<jobid>.out
#SBATCH --job-name=sw_ablc
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/abl-score-common-%j.out
set -euo pipefail

run_dir=${1:?usage: sbatch run_score_ablation_common.sbatch <run_dir> <arms>}
arms=${2:?usage: sbatch run_score_ablation_common.sbatch <run_dir> <arms>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap_dir="$(mktemp -d)"
cp experiments/benchmarks/05_bmv-unmix-compare/51c_score_ablation_common.R "$snap_dir"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/lib "$snap_dir"/lib
sed -i "s#LIB <- \"experiments/benchmarks/05_bmv-unmix-compare/lib\"#LIB <- \"$snap_dir/lib\"#" "$snap_dir"/51c_score_ablation_common.R
pixi run Rscript "$snap_dir/51c_score_ablation_common.R" "$run_dir" "$arms"
