# Purpose: per-dataset config for the Python dashboard builders of benchmark 05. Same semantics as
#          datasets/dataset.R: the dataset comes from the required env var SW_DATASET (3d, 24h, exv), with no
#          default. Standard library only; runs on the system python 3.6 (pixi python is newer).
# Run:     import dataset as ds   (the script folder dashboard/ must be on sys.path; repo root is the cwd)
#          SW_DATASET_DIR (set by the sbatch wrappers) points at a frozen copy of datasets/; unset, it is the repo copy.
# Inputs:  <dataset dir>/<ds>.json         config values
#          <dataset dir>/<ds>_runs.csv     run registry: key, path, producer, approved_on
#          <dataset dir>/<ds>_samples.csv  sample, donor, condition, well, acq_batch, stem
# Outputs: none; defines DS, DATASET, cfg(), data_path(), run_in(), run_out(), samples().

import csv
import json
import os
import re
import sys
from pathlib import Path

DATASET = os.environ.get("SW_DATASET", "")
if DATASET not in ("3d", "24h", "exv"):
    sys.exit("dataset.py: SW_DATASET must be one of 3d, 24h, exv; got '%s'" % DATASET)

_DIR = Path(os.environ.get("SW_DATASET_DIR", "experiments/benchmarks/05_bmv-unmix-compare/datasets"))
_JSON = _DIR / (DATASET + ".json")
DS = json.loads(_JSON.read_text())
if DS.get("dataset") != DATASET:
    sys.exit("dataset.py: %s does not name itself '%s'" % (_JSON, DATASET))

with (_DIR / (DATASET + "_runs.csv")).open(newline="") as _f:
    _RUNS = list(csv.DictReader(_f))


def cfg(name):
    """Config value; stops when the field is missing or null (never a silent default)."""
    v = DS.get(name)
    if v is None:
        sys.exit(
            "dataset.py: config '%s' is missing or null for dataset %s (%s); an impactful value needs the user's decision"
            % (name, DATASET, _JSON)
        )
    return v


def data_path(*parts):
    return Path(cfg("data_dir"), *parts)


def run_in(key):
    """Path of a registered run folder; stops when the key is not registered once or the folder is missing."""
    hit = [r for r in _RUNS if r["key"] == key]
    if len(hit) != 1:
        sys.exit("dataset.py: run key '%s' is not registered once in %s_runs.csv (found %d)" % (key, DATASET, len(hit)))
    path = Path(hit[0]["path"])
    if not path.is_dir():
        sys.exit("dataset.py: registered run folder is missing: %s" % path)
    return path


def run_out(key):
    """New run folder experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-<key>; stops if it exists. Does not register it."""
    day = os.environ.get("SW_RUN_DATE", "")
    if not re.match(r"^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", day):
        sys.exit("dataset.py: SW_RUN_DATE must be YYYY-MM-DD")
    out = Path("experiments/runs") / ("%s_bmv%s-%s" % (day, cfg("run_tag"), key))
    if out.exists():
        sys.exit("dataset.py: run dir exists: %s" % out)
    return out


def samples():
    with (_DIR / (DATASET + "_samples.csv")).open(newline="") as f:
        return list(csv.DictReader(f))
