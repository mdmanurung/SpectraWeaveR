#!/usr/bin/env python3
"""Purpose: make a fresh run folder that reuses finished method results by symlink, so a new arm scores beside them.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/link_methods.py <old_run_dir> <new_run_dir> <method> [<method> ...]
Inputs:  <old_run_dir>/results/seed_<n>/<method>/ for each method
Outputs: <new_run_dir>/results/seed_<n>/<method> -> absolute symlink to the old folder; <new_run_dir> must not exist
"""

import sys
from pathlib import Path


def main():
    old, new = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    methods = sys.argv[3:]
    new.mkdir(parents=True, exist_ok=False)
    for seed_dir in sorted((old / "results").glob("seed_*")):
        out = new / "results" / seed_dir.name
        out.mkdir(parents=True)
        for m in methods:
            (out / m).symlink_to(seed_dir / m, target_is_directory=True)


if __name__ == "__main__":
    main()
