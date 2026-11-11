#!/usr/bin/env bash
# Purpose: unmix one per-dye variant arm (45_variant_per_dye.R) in the pixi R environment; one array task per dye.
# Run:     sbatch --array=1-32 experiments/slurm/run_variant_per_dye.sbatch <dye_list>
#          dye_list holds one dye per line; task i runs line i.
# Inputs:  <dye_list>, the 1.8.8 and 1.0.2 zip folders under experiments/bmv_test/compare_3D/data
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_<dye_id>/
#          experiments/slurm/logs/vpd-<arrayjob>_<task>.out
#SBATCH --job-name=sw_vpd
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/slurm/logs/vpd-%A_%a.out
set -euo pipefail

dye_list=${1:?usage: sbatch --array=1-N run_variant_per_dye.sbatch <dye_list>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
dye=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$dye_list")
[ -n "$dye" ] || { echo "no dye on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "dye: $dye"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/45_variant_per_dye.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/45_variant_per_dye.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
pixi run Rscript "$snap" "$dye"
