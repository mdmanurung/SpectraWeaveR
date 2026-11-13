#!/usr/bin/env bash
# Purpose: extract and unmix one AS 1.8.8 tuning arm of 103_af_tuning.R on the 24H smoke donor in the pixi R environment.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_af_tuning.sbatch <refine TRUE|FALSE> <n_passes> <threshold>
# Inputs:  see 103_af_tuning.R
# Outputs: experiments/runs/<SW_RUN_DATE>_bmv24h-af-tuning/<arm>/; log experiments/slurm/logs/aftune-<jobid>.out
#SBATCH --export=ALL
#SBATCH --job-name=sw_aftune
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/slurm/logs/aftune-%j.out
set -euo pipefail
for v in SW_DATASET SW_RUN_DATE; do
  [ -n "${!v:-}" ] || { echo "run_af_tuning.sbatch: $v is empty" >&2; exit 1; }
done
# first argument "unstained" runs 103c_unmix_unstained.R over every finished arm instead
refine=${1:?usage: sbatch run_af_tuning.sbatch <refine> <n_passes> <threshold> | unstained}
script=103_af_tuning.R
if [ "$refine" = unstained ]; then script=103c_unmix_unstained.R; fi
if [ "$refine" = heldout ]; then script=103f_heldout.R; fi
if [ "$refine" = unstained102 ]; then script=103h_unmix_unstained_102.R; fi
cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
export OMP_NUM_THREADS=1
snap="$(mktemp -d)/$script"
cp experiments/benchmarks/05_bmv-unmix-compare/$script "$snap"
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
cp -r experiments/benchmarks/05_bmv-unmix-compare/{datasets,lib,style,lucid} "$(dirname "$snap")"/ 2>/dev/null || true
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R
if [ "$script" = 103c_unmix_unstained.R ] || [ "$script" = 103h_unmix_unstained_102.R ]; then pixi run Rscript "$snap"; elif [ "$script" = 103f_heldout.R ]; then pixi run Rscript "$snap" "${2:?refine}" "${3:?n_passes}"; else pixi run Rscript "$snap" "$refine" "${2:?n_passes}" "${3:?threshold}"; fi
