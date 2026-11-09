#!/usr/bin/env bash
# Purpose: download named files of one Zenodo record, verify each md5 from the Zenodo API, and
#          extract .7z archives. Used for datasets with many files (registry rows of type zenodo).
# Run:     SW_DATA_DIR=<dir> bash experiments/datasets/download_zenodo.sh <name> <record_id> <file> [file ...]
#          use the file name "ALL" to fetch every file of the record
# Outputs: $SW_DATA_DIR/<name>/raw/<file>      downloaded files
#          $SW_DATA_DIR/<name>/data/           extracted 7z archives (other files are copied)
#          $SW_DATA_DIR/<name>/DONE or FAILED  marker
# Data never go inside the repository.
set -euo pipefail
data_dir=${SW_DATA_DIR:?set SW_DATA_DIR to a folder outside the repository}
name=$1
record=$2
shift 2
dest="$data_dir/$name"
mkdir -p "$dest/raw" "$dest/data"
listing=$(curl -sSL --retry 3 "https://zenodo.org/api/records/$record")

echo "$listing" | python3 -c '
import json, sys
want = sys.argv[1:]
for f in json.load(sys.stdin)["files"]:
    if "ALL" in want or f["key"] in want:
        print(f["key"], f["checksum"].split(":")[1], f["links"]["self"])
' "$@" | while read -r key md5 url; do
  out="$dest/raw/$key"
  if [ ! -e "$out" ] || [ "$(md5sum "$out" | cut -d' ' -f1)" != "$md5" ]; then
    curl -sSL --retry 3 -o "$out" "$url"
  fi
  got=$(md5sum "$out" | cut -d' ' -f1)
  if [ "$got" != "$md5" ]; then
    echo "CHECKSUM MISMATCH $key: $got != $md5" | tee -a "$dest/FAILED"
    continue
  fi
  case "$key" in
    *.7z) 7za x -y -bd -o"$dest/data" "$out" > /dev/null ;;
    *) cp "$out" "$dest/data/" ;;
  esac
  echo "OK $key"
done
[ -e "$dest/FAILED" ] || echo "OK $name record $record" > "$dest/DONE"
