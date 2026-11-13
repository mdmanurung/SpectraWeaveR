#!/usr/bin/env bash
# Purpose: score ablation arms on the common SpectroFlo gates (51c_score_ablation_common.R), reading only.
# Run:     sbatch experiments/slurm/run_score_ablation_common.sbatch <run_dir> <arm1,arm2,...>
# Inputs:  registry key common-gated (datasets/<ds>_runs.csv): batch 0 of the common gates
# Outputs: <run_dir>/ (scores, gap_closed, labels, figure); experiments/slurm/logs/abl-score-common-<jobid>.out
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes experiments/benchmarks/05_bmv-unmix-compare/datasets/ and
#          exports SW_DATASET_DIR=<frozen copy>/datasets, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_ablc
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/abl-score-common-%j.out
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_score_ablation_common.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

run_dir=${1:?usage: sbatch run_score_ablation_common.sbatch <run_dir> <arms>}
arms=${2:?usage: sbatch run_score_ablation_common.sbatch <run_dir> <arms>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap_dir="$(mktemp -d)"
cp experiments/benchmarks/05_bmv-unmix-compare/51c_score_ablation_common.R "$snap_dir"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/lib "$snap_dir"/lib
cp -r experiments/benchmarks/05_bmv-unmix-compare/datasets "$snap_dir"/datasets
export SW_DATASET_DIR="$snap_dir/datasets" # dataset.R reads the frozen config and registry
sed -i "s#LIB <- \"experiments/benchmarks/05_bmv-unmix-compare/lib\"#LIB <- \"$snap_dir/lib\"#" "$snap_dir"/51c_score_ablation_common.R
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R\")#source(\"$snap_dir/datasets/dataset.R\")#" "$snap_dir"/51c_score_ablation_common.R
pixi run Rscript "$snap_dir/51c_score_ablation_common.R" "$run_dir" "$arms"
