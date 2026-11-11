#!/usr/bin/env bash
# Purpose: unmix one donor of one noise-model GLS arm (42_unmix_gls_noise.R) in the pixi R environment.
# Run:     sbatch experiments/slurm/run_gls_donor.sbatch <arm> <donor>
#          e.g. sbatch experiments/slurm/run_gls_donor.sbatch gls188_noise A3
# Inputs:  experiments/runs/2026-11-10_bmv-gls-noise/noise_model/noise_model.rds (42 ... noise)
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/<donor> ... AutoSpectral.fcs
#          experiments/slurm/logs/gls-<jobid>.out
#SBATCH --job-name=sw_gls
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/slurm/logs/gls-%j.out
set -euo pipefail

arm=${1:?usage: sbatch run_gls_donor.sbatch <arm> <donor>}
donor=${2:?usage: sbatch run_gls_donor.sbatch <arm> <donor>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
# one BLAS thread per R worker; the script runs SLURM_CPUS_PER_TASK workers itself
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/42_unmix_gls_noise.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/42_unmix_gls_noise.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$arm" "$donor"
