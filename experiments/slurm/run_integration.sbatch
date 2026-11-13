#!/usr/bin/env bash
# Purpose: run one integration script (experiments/integration/*.R) in the pixi R environment of the main checkout,
#          loading the package from any checkout (a worktree has no .pixi and no data).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_integration.sbatch \
#            <script.R> <pkg_dir> <repo_root> <out_dir> <dataset_json>
#          e.g. SW_DATASET=3d SW_RUN_DATE=2026-11-12 sbatch experiments/slurm/run_integration.sbatch \
#                 experiments/integration/clean_run_bmv.R "$PWD" /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR \
#                 /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/experiments/runs/2026-11-12_clean-run-3d \
#                 "$PWD"/experiments/benchmarks/05_bmv-unmix-compare/datasets/3d.json
#          <dataset_json> is required and has no default; its file name must be $SW_DATASET.json and its folder must
#          hold dataset.R, <ds>_runs.csv and <ds>_samples.csv. The job runs with repo_root as working directory, so
#          the registry paths experiments/runs/... resolve there.
# Inputs:  the R script, a package checkout, the checkout that holds experiments/bmv_test, a new out_dir, the dataset json
# Outputs: <out_dir>/ (written by the script) and experiments/slurm/logs/integration-<jobid>.out in the main checkout
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes the dataset_json folder and passes the frozen json path to
#          the script, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_integration
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/experiments/slurm/logs/integration-%j.out
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_integration.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

usage="usage: sbatch run_integration.sbatch <script.R> <pkg_dir> <repo_root> <out_dir> <dataset_json>"
script=${1:?$usage}
pkg_dir=${2:?$usage}
repo_root=${3:?$usage}
out_dir=${4:?$usage}
dataset_json=${5:?$usage}

out_dir="$(realpath -m "$out_dir")"
dataset_json="$(realpath "$dataset_json")"
[ ! -e "$out_dir" ] || { echo "run dir exists: $out_dir" >&2; exit 1; }
[ "$(basename "$dataset_json")" = "$SW_DATASET.json" ] || { echo "run_integration.sbatch: $dataset_json is not $SW_DATASET.json" >&2; exit 1; }
[ -f "$(dirname "$dataset_json")/dataset.R" ] || { echo "run_integration.sbatch: no dataset.R beside $dataset_json" >&2; exit 1; }

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap_dir="$(mktemp -d)"
snap="$snap_dir/$(basename "$script")"
cp "$pkg_dir/$script" "$snap"
# freeze the dataset folder and the benchmark helpers the script sources, and point its source() lines at the copies
cp -r "$(dirname "$dataset_json")" "$snap_dir/datasets"
cp "$pkg_dir"/experiments/benchmarks/05_bmv-unmix-compare/*.R "$snap_dir"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$snap_dir/#" "$snap_dir"/*.R "$snap_dir"/datasets/*.R
# registry paths are relative to the checkout that holds the data and the runs
cd "$repo_root"
# A non-zero R exit fails the job (set -e); the pixi manifest of the main checkout supplies the environment.
pixi run --manifest-path "$repo_root/pixi.toml" Rscript "$snap" "$pkg_dir" "$repo_root" "$out_dir" "$snap_dir/datasets/$(basename "$dataset_json")"
