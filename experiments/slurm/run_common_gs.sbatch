#!/usr/bin/env bash
# Purpose: run one of the common-gate scripts (97e_common_gs.R, 01c, 02c, 97f figure) in the pixi R environment.
# Run:     sbatch experiments/slurm/run_common_gs.sbatch <script file name> [script args]
# Inputs:  see the header of the script
# Outputs: see the header of the script; experiments/slurm/logs/common-<jobid>.out
#SBATCH --job-name=sw_common
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/common-%j.out
set -euo pipefail

script=${1:?usage: sbatch run_common_gs.sbatch <script file name> [script args]}
shift

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
# one BLAS thread per R worker; the script forks SLURM_CPUS_PER_TASK workers itself
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/${script}"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid} "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid}/*.R
pixi run Rscript "$snap" "$@"
