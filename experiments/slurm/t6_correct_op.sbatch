#!/usr/bin/env bash
# Purpose: run 02b_correct_op.R for one (method, fold, split) of the T6 OpenProblems track; a job array over 56 (cytofin_pop: 48-55).
# Run:     sbatch --array=0-55%5 experiments/slurm/t6_correct_op.sbatch <prep_dir> <out_dir>
#          index i -> method METHODS[i / 8], fold FOLDS[(i % 8) / 2], split (i % 2) + 1; %5 caps concurrent jobs.
# Inputs:  <prep_dir> from 01_prep.R and 01b_prep_op.R
# Outputs: <out_dir>/<method>/<fold>/split_<k>/ and slurm-<jobid>_<index>.out in the submit directory
#SBATCH --job-name=sw_t6op
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
set -euo pipefail

prep_dir=${1:?usage: sbatch t6_correct_op.sbatch <prep_dir> <out_dir>}
out_dir=${2:?}
METHODS=(none paired_shift heb cycombine cytonorm cytofin cytofin_pop)
FOLDS=(op_tcu1_tcs3 op_tcu1_tcs4 op_tcu2_tcs3 op_tcu2_tcs4)
method=${METHODS[$((SLURM_ARRAY_TASK_ID / 8))]}
fold=${FOLDS[$(((SLURM_ARRAY_TASK_ID % 8) / 2))]}
split=$((SLURM_ARRAY_TASK_ID % 2 + 1))

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.claude/worktrees/cytether-benchmark
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
R_LIBS=/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/Rlib_competitors \
  /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript \
  experiments/benchmarks/06_cytether-cytobatchflagr/02b_correct_op.R "$prep_dir" "$out_dir" "$method" "$fold" "$split"
