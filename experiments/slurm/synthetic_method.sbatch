#!/usr/bin/env bash
# Purpose: run one method over one seed of the synthetic benchmark; use as a job array over seeds.
# Run:     sbatch --array=0-9%5 experiments/slurm/synthetic_method.sbatch <sim_dir> <results_dir> <method> <first_seed>
#          array index i runs seed <first_seed>+i; %5 caps concurrent jobs (user-approved F4).
# Inputs:  <sim_dir>/seed_<n>/ with observed_X.csv.gz
# Outputs: <results_dir>/seed_<n>/<method>/<stage>/ and slurm-<jobid>_<index>.out in the submit directory
#SBATCH --job-name=sw_synth
#SBATCH --account=div2-para
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -euo pipefail

sim_dir=${1:?usage: sbatch synthetic_method.sbatch <sim_dir> <results_dir> <method> <first_seed>}
results_dir=${2:?}
method=${3:?}
first_seed=${4:?}
seed=$((first_seed + SLURM_ARRAY_TASK_ID))

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.claude/worktrees/cytether-benchmark
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
bash experiments/benchmarks/05_cytether-synthetic/run_seeds.sh "$sim_dir" "$results_dir" "$method" "$seed"
