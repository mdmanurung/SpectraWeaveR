#!/usr/bin/env bash
# Purpose: run one line of a job list as an array task: "<script.R> [args...]" per line, from a frozen copy.
#          Copy of run_spike_cg.sbatch for arrays; the old 88 jobs were submitted one by one with
#          run_af405_proof.sbatch.
# Run:     sbatch --array=1-N[%cap] [--time=HH:MM:SS] experiments/slurm/run_spike_cg_array.sbatch <jobs.txt>
#          task i runs line i of <jobs.txt>; any failing task exits non-zero.
#          SW_SMOKE_DONOR=<donor> in the submit environment runs on that one whole donor (smoke dirs).
# Inputs:  <jobs.txt> and the scripts it names
# Outputs: see each script's header; log in experiments/slurm/logs/spikecg-<arrayjob>_<task>.out
# Dataset: SW_DATASET (3d, 24h, exv) and SW_RUN_DATE (YYYY-MM-DD) must be set in the submit shell; the job stops when
#          either is empty. Checked: sbatch exports the submit environment by default (--export=ALL, stated below on
#          purpose), so both reach the job. The job freezes experiments/benchmarks/05_bmv-unmix-compare/datasets/ and
#          exports SW_DATASET_DIR=<frozen copy>/datasets, so dataset.R reads the frozen config and registry.
#SBATCH --export=ALL
#SBATCH --job-name=sw_spikecg
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/spikecg-%A_%a.out
set -euo pipefail

for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_spike_cg_array.sbatch: $v is empty; set SW_DATASET (3d, 24h or exv) and SW_RUN_DATE (YYYY-MM-DD) before sbatch" >&2; exit 1; }
done
echo "SW_DATASET=$SW_DATASET SW_RUN_DATE=$SW_RUN_DATE"

jobs=${1:?usage: sbatch --array=1-N run_spike_cg_array.sbatch <jobs.txt>}
cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$jobs")"
[ -n "$line" ] || { echo "empty job line ${SLURM_ARRAY_TASK_ID} in $jobs" >&2; exit 1; }
echo "task ${SLURM_ARRAY_TASK_ID}: $line"
read -r -a words <<< "$line"
script=${words[0]}
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "experiments/benchmarks/05_bmv-unmix-compare/${script}")"
cp "experiments/benchmarks/05_bmv-unmix-compare/${script}" "$snap"
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid,datasets} "$(dirname "$snap")"/
export SW_DATASET_DIR="$(dirname "$snap")/datasets" # dataset.R reads the frozen config and registry; its source() line is rewritten below
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid,datasets}/*.R
pixi run Rscript "$snap" "${words[@]:1}"
