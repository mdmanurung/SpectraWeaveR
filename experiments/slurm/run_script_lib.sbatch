#!/usr/bin/env bash
# Purpose: run one benchmark script in the pixi R environment from a frozen copy that also freezes lib/
#          (figure templates), so scripts that source lib/theme.R etc. work. Copy of run_script.sbatch's freeze
#          logic; the script reads the frozen lib/ path from SW_LIB and the repo root from SW_REPO.
# Run:     sbatch experiments/slurm/run_script_lib.sbatch <script.R> [args passed to the script]
# Inputs:  the R script; experiments/benchmarks/05_bmv-unmix-compare/*.R, lib/ and datasets/
# Outputs: whatever the script writes; slurm-<jobid>.out in the submit directory
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes experiments/benchmarks/05_bmv-unmix-compare/datasets/ and
#          exports SW_DATASET_DIR=<frozen copy>/datasets, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_lib
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_script_lib.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

script=${1:?usage: sbatch run_script_lib.sbatch <script.R> [args]}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export SW_REPO="$PWD"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Frozen copy: Rscript reads lazily, so an edit during the run would break the job.
snap_dir="$(mktemp -d)"
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$snap_dir"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/lib "$snap_dir"/lib
cp -r experiments/benchmarks/05_bmv-unmix-compare/datasets "$snap_dir"/datasets
export SW_DATASET_DIR="$snap_dir/datasets" # dataset.R reads the frozen config and registry
cp "$script" "$snap_dir"/
export SW_LIB="$snap_dir/lib"
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$snap_dir/#" "$snap_dir"/*.R "$snap_dir"/datasets/*.R
pixi run Rscript "$snap_dir/$(basename "$script")" "${@:2}"
