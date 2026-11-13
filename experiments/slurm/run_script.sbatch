#!/usr/bin/env bash
# Purpose: run one experiment script inside the project's pixi R environment.
# Run:     sbatch experiments/slurm/run_script.sbatch <script.R> <run_dir> [extra args]
#          e.g. sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/01_mendeley-7color/run.R experiments/runs/2026-11-09_mendeley-7color
# Inputs:  the R script and a run directory that does not exist yet
# Outputs: <run_dir>/ (written by the script) and slurm-<jobid>.out in the submit directory
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes experiments/benchmarks/05_bmv-unmix-compare/datasets/ and
#          exports SW_DATASET_DIR=<frozen copy>/datasets, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_beta
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_script.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

script=${1:?usage: sbatch run_script.sbatch <script.R> <run_dir>}
run_dir=${2:?usage: sbatch run_script.sbatch <script.R> <run_dir>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
# a new run needs a new folder; figure-only reruns set ALLOW_EXISTING_RUN=1 (the script itself refuses to overwrite)
[ ! -e "$run_dir" ] || [ "${ALLOW_EXISTING_RUN:-0}" = 1 ] || { echo "run dir exists: $run_dir" >&2; exit 1; }

# threads for AutoSpectral and BLAS follow the job allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "$script")"
cp "$script" "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
# freeze the helper folders too (lib/theme.R and the style/ and lucid/ pointers to it)
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid,datasets} "$(dirname "$snap")"/
export SW_DATASET_DIR="$(dirname "$snap")/datasets" # dataset.R reads the frozen config and registry; its source() line is rewritten below
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid,datasets}/*.R
pixi run Rscript "$snap" "$run_dir" "${@:3}"
