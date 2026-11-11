#!/usr/bin/env bash
# Purpose: run one line of a job list as an array task: "<script.R> [args...]" per line, from a frozen copy.
#          Copy of run_spike_cg.sbatch for arrays; the old 88 jobs were submitted one by one with
#          run_af405_proof.sbatch.
# Run:     sbatch --array=1-N[%cap] [--time=HH:MM:SS] experiments/slurm/run_spike_cg_array.sbatch <jobs.txt>
#          task i runs line i of <jobs.txt>; any failing task exits non-zero.
#          SW_SMOKE_DONOR=<donor> in the submit environment runs on that one whole donor (smoke dirs).
# Inputs:  <jobs.txt> and the scripts it names
# Outputs: see each script's header; log in experiments/slurm/logs/spikecg-<arrayjob>_<task>.out
#SBATCH --job-name=sw_spikecg
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/spikecg-%A_%a.out
set -euo pipefail

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
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid} "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid}/*.R
pixi run Rscript "$snap" "${words[@]:1}"
