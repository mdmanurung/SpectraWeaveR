#!/usr/bin/env bash
# Purpose: run every scenario in scenarios.tsv in parallel and append one JSON line each.
# Run:     bash experiments/benchmarks/06_anomaly-features/run_batch.sh <out.jsonl> <jobs>
# Inputs:  scenarios.tsv (tag, candidate, env overrides); R and isotree on the paths below
# Outputs: <out.jsonl> with one line per scenario: {"tag":..., <run.R output> or "error":...}
set -u
OUT=$1
JOBS=$2
RSCRIPT=${RSCRIPT:-/exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript}
DIR=experiments/benchmarks/06_anomaly-features

run_one() {
  tag=$1; cand=$2; envs=$3
  res=$(env $envs "$RSCRIPT" "$DIR/run.R" "$cand" 2>&1 | tail -1)
  case "$res" in
    "{"*) echo "{\"tag\":\"$tag\",\"result\":$res}" ;;
    *) echo "{\"tag\":\"$tag\",\"error\":\"$(echo "$res" | tr '"' "'" | cut -c1-200)\"}" ;;
  esac
}
export -f run_one
export RSCRIPT DIR

tail -n +2 "$DIR/scenarios.tsv" |
  xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r t c e <<< "{}"; run_one "$t" "$c" "$e"' >> "$OUT"
