#!/usr/bin/env bash
# Purpose: run one proof-chain script (60-64, AF405 vs BV510) in the pixi R environment.
# Run:     sbatch [--dependency=afterok:<id>] experiments/slurm/run_af405_proof.sbatch <script.R> [args...]
#          script.R is a file name in experiments/benchmarks/05_bmv-unmix-compare/.
# Inputs:  see each script's header
# Outputs: see each script's header; log in experiments/slurm/logs/af405proof-<jobid>.out
#SBATCH --job-name=sw_proof
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/slurm/logs/af405proof-%j.out
set -euo pipefail

script=${1:?usage: sbatch run_af405_proof.sbatch <script.R> [args...]}
shift

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename "experiments/benchmarks/05_bmv-unmix-compare/${script}")"
cp "experiments/benchmarks/05_bmv-unmix-compare/${script}" "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$@"
