#!/usr/bin/env bash
# Purpose: build (and for stage full, unmix) one optimised 1.8.8 arm with 53_opt188.R in the pixi R environment.
# Run:     sbatch experiments/slurm/run_opt188.sbatch <route> <neg> <n_passes> <stage>
#          route = pregate | legacy; neg = neg_same_fsc4e5 | neg_h12; stage = extract | full
# Inputs:  the 1.8.8 zip folders under experiments/bmv_test/compare_3D/data, the landmark control-file copy
# Outputs: experiments/runs/2026-11-10_bmv-opt188-smoke/<arm>/ (extract) or
#          experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/ (full)
#          experiments/slurm/logs/opt188-<jobid>.out
#SBATCH --job-name=sw_opt188
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/slurm/logs/opt188-%j.out
set -euo pipefail

route=${1:?usage: sbatch run_opt188.sbatch <route> <neg> <n_passes> <stage>}
neg=${2:?usage: sbatch run_opt188.sbatch <route> <neg> <n_passes> <stage>}
n_passes=${3:?usage: sbatch run_opt188.sbatch <route> <neg> <n_passes> <stage>}
stage=${4:?usage: sbatch run_opt188.sbatch <route> <neg> <n_passes> <stage>}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR
echo "route: $route  neg: $neg  n_passes: $n_passes  stage: $stage"
export OMP_NUM_THREADS=1
# Run a frozen copy: Rscript reads the script lazily, so an edit during the run breaks the job.
snap="$(mktemp -d)/$(basename experiments/benchmarks/05_bmv-unmix-compare/53_opt188.R)"
cp experiments/benchmarks/05_bmv-unmix-compare/53_opt188.R "$snap"
# freeze sourced helpers too: copy the benchmark scripts beside it and point source() at the copies
cp experiments/benchmarks/05_bmv-unmix-compare/*.R "$(dirname "$snap")"/
# freeze the helper folders too (lib/theme.R and the style/ and lucid/ pointers to it)
cp -r experiments/benchmarks/05_bmv-unmix-compare/{lib,style,lucid} "$(dirname "$snap")"/
sed -i "s#source(\"experiments/benchmarks/05_bmv-unmix-compare/#source(\"$(dirname "$snap")/#" "$(dirname "$snap")"/*.R "$(dirname "$snap")"/{lib,style,lucid}/*.R
pixi run Rscript "$snap" "$route" "$neg" "$n_passes" "$stage"
