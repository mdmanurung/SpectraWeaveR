#!/usr/bin/env bash
# Purpose: run one common-gates spike script (88cg ... 88hcg, 62cg, 64cg, 61) in the pixi R environment from a frozen copy.
#          Copy of run_af405_proof.sbatch with a time limit sized for the 5-pass unmixes and a log name of its own.
# Run:     sbatch [--dependency=afterok:<id>] [--time=HH:MM:SS] experiments/slurm/run_spike_cg.sbatch <script.R> [args...]
#          script.R is a file name in experiments/benchmarks/05_bmv-unmix-compare/.
#          SW_SMOKE_DONOR=<donor> in the submit environment runs the script on that one whole donor (smoke dirs).
# Inputs:  see each script's header
# Outputs: see each script's header; log in experiments/slurm/logs/spikecg-<jobid>.out
#SBATCH --job-name=sw_spikecg
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/spikecg-%j.out
set -euo pipefail

script=${1:?usage: sbatch run_spike_cg.sbatch <script.R> [args...]}
shift

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "experiments/benchmarks/05_bmv-unmix-compare/${script}")"
cp "experiments/benchmarks/05_bmv-unmix-compare/${script}" "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
# freeze the helper folders too (lib/theme.R and the style/ and lucid/ pointers to it)
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid} "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid}/*.R
# a non-zero Rscript exit fails the job (set -e); no output is swallowed
pixi run Rscript "$snap" "$@"
