#!/usr/bin/env bash
# Purpose: download and verify the public test datasets listed in registry.csv.
# Run:     bash experiments/datasets/download.sh [name ...]      (no name = all)
#          data go to $SW_DATA_DIR (default: /exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata)
# Inputs:  experiments/datasets/registry.csv
# Outputs: $SW_DATA_DIR/<name>/{zip,data}/ and a DONE or FAILED marker per dataset
# Data never go inside the repository.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
data_dir=${SW_DATA_DIR:-/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata}
want=("$@")

tail -n +2 "$here/registry.csv" | python3 -c '
import csv, sys
for r in csv.reader(sys.stdin):
    print("\t".join([r[0], r[2], r[3], r[4]]))
' | while IFS=$'\t' read -r name url algo sum; do
  if [ ${#want[@]} -gt 0 ] && [[ ! " ${want[*]} " =~ " $name " ]]; then continue; fi
  [ -e "$data_dir/$name/DONE" ] && { echo "skip $name (DONE)"; continue; }
  mkdir -p "$data_dir/$name/zip"
  curl -sSL --retry 3 -o "$data_dir/$name/zip/$name.zip" "$url"
  got=$("${algo}sum" "$data_dir/$name/zip/$name.zip" | cut -d' ' -f1)
  if [ "$got" != "$sum" ]; then
    echo "CHECKSUM MISMATCH $name: $got != $sum" | tee "$data_dir/$name/FAILED"; continue
  fi
  mkdir -p "$data_dir/$name/data"
  unzip -q -o "$data_dir/$name/zip/$name.zip" -d "$data_dir/$name/data"
  echo "OK $name" | tee "$data_dir/$name/DONE"
done
