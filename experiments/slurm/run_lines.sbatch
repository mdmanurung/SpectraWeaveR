#!/usr/bin/env bash
# Purpose: run a chunk of a command list (one shell command per line) as one array task, in parallel inside it.
# Run:     sbatch --array=0-<n_chunks-1>%<cap> experiments/slurm/run_lines.sbatch <commands.txt> <lines_per_task>
#          task i runs lines i*k+1 .. (i+1)*k with xargs -P $SLURM_CPUS_PER_TASK; any failing line fails the task.
# Inputs:  <commands.txt>, commands written to run from the worktree root
# Outputs: whatever the commands write; slurm-<jobid>_<index>.out in the submit directory
#SBATCH --job-name=sw_lines
#SBATCH --account=div2-para
#SBATCH --partition=all
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=8:00:00
set -euo pipefail

commands=${1:?usage: sbatch run_lines.sbatch <commands.txt> <lines_per_task>}
k=${2:?}
first=$((SLURM_ARRAY_TASK_ID * k + 1))
last=$(((SLURM_ARRAY_TASK_ID + 1) * k))

cd /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.claude/worktrees/cytether-benchmark
export OMP_NUM_THREADS=1
# -d '\n': one line is one command, quotes intact (plain xargs strips them and splits quoted lists).
sed -n "${first},${last}p" "$commands" | xargs -d '\n' -P "${SLURM_CPUS_PER_TASK:-1}" -I{} bash -c '{}'
