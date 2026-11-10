#!/usr/bin/env bash
# Purpose: run 04_score_op.py for one (fold, set) of the T6 OpenProblems track; a job array over 40 (10 sets per fold; jobs before 2026-11-10 21h used 9 sets).
# Run:     sbatch --array=0-39%5 experiments/slurm/t6_score_op.sbatch <prep_dir> <results_dir> <score_dir>
#          index i -> fold FOLDS[i / 10], set SETS[i % 10]; %5 caps concurrent jobs. Then run 05_scale_op.py.
# Inputs:  <prep_dir> from 01b_prep_op.R; <results_dir> from t6_correct_op.sbatch
# Outputs: <score_dir>/<fold>/<set>/scores_raw.csv and slurm-<jobid>_<index>.out in the submit directory
#SBATCH --job-name=sw_t6sc
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=8:00:00
set -euo pipefail

prep_dir=${1:?usage: sbatch t6_score_op.sbatch <prep_dir> <results_dir> <score_dir>}
results_dir=${2:?}
score_dir=${3:?}
FOLDS=(op_tcu1_tcs3 op_tcu1_tcs4 op_tcu2_tcs3 op_tcu2_tcs4)
SETS=(no_integration perfect_integration shuffle_integration_globally none paired_shift heb cycombine cytonorm cytofin cytofin_pop)
fold=${FOLDS[$((SLURM_ARRAY_TASK_ID / 10))]}
set=${SETS[$((SLURM_ARRAY_TASK_ID % 10))]}

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.claude/worktrees/cytether-benchmark
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/venv_cytether/bin/python -I \
  experiments/benchmarks/06_cytether-cytobatchflagr/04_score_op.py "$prep_dir" "$results_dir" "$score_dir" "$fold" "$set"
