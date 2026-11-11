#!/usr/bin/env bash
# Purpose: gate and score ablation arms (51_score_ablation.R) in the pixi R environment.
# Run:     sbatch [--dependency=afterok:<jobid>] experiments/slurm/run_score_ablation.sbatch <run_dir> <arm1,arm2,...>
# Inputs:  experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/ for each arm; the base arms' GatingSets
# Outputs: <run_dir>/ (scores, gap_closed, labels, figure); experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>
#          experiments/slurm/logs/abl-score-<jobid>.out
#SBATCH --job-name=sw_ablsc
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/abl-score-%j.out
set -euo pipefail

run_dir=${1:?usage: sbatch run_score_ablation.sbatch <run_dir> <arms>}
arms=${2:?usage: sbatch run_score_ablation.sbatch <run_dir> <arms>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/51_score_ablation.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/51_score_ablation.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$run_dir" "$arms"
