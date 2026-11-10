#!/usr/bin/env bash
# Purpose: run 02_correct.R for one (method, fold) of T6; use as a job array over the 24 pairs.
# Run:     sbatch --array=0-23%5 experiments/slurm/t6_correct.sbatch <prep_dir> <out_dir>
#          index i -> method METHODS[i / 4], fold FOLDS[i % 4]; %5 caps concurrent jobs (user-approved F4).
# Inputs:  <prep_dir> from 01_prep.R
# Outputs: <out_dir>/<method>/fold_<held_out>/ and slurm-<jobid>_<index>.out in the submit directory
#SBATCH --job-name=sw_t6
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
set -euo pipefail

prep_dir=${1:?usage: sbatch t6_correct.sbatch <prep_dir> <out_dir>}
out_dir=${2:?}
METHODS=(none paired_shift heb cycombine cytonorm cytofin)
FOLDS=(tcu1 tcu2 tcs3 tcs4)
method=${METHODS[$((SLURM_ARRAY_TASK_ID / 4))]}
fold=${FOLDS[$((SLURM_ARRAY_TASK_ID % 4))]}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.claude/worktrees/cytether-benchmark
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
R_LIBS=/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/Rlib_competitors \
  /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript \
  experiments/benchmarks/06_cytether-cytobatchflagr/02_correct.R "$prep_dir" "$out_dir" "$method" "$fold"
