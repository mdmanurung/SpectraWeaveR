#!/usr/bin/env bash
# Purpose: run one of the clean3 scripts (97c_reclean_pooled.R, 97d_spectroflo_labels_clean3.R) in the pixi R environment.
# Run:     sbatch experiments/slurm/run_clean3.sbatch <script file name in experiments/benchmarks/05_bmv-unmix-compare>
#          e.g. sbatch experiments/slurm/run_clean3.sbatch 97c_reclean_pooled.R
# Inputs:  see the header of the script
# Outputs: see the header of the script; experiments/slurm/logs/clean3-<jobid>.out
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes experiments/benchmarks/05_bmv-unmix-compare/datasets/ and
#          exports SW_DATASET_DIR=<frozen copy>/datasets, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_clean3
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/slurm/logs/clean3-%j.out
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_clean3.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

script=${1:?usage: sbatch run_clean3.sbatch <script file name>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
# one BLAS thread per R worker; the script forks SLURM_CPUS_PER_TASK workers itself
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/${script}"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid,datasets} "$(dirname "$snap")"/
export SW_DATASET_DIR="$(dirname "$snap")/datasets" # dataset.R reads the frozen config and registry; its source() line is rewritten below
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid,datasets}/*.R
pixi run Rscript "$snap"
