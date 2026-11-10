#!/usr/bin/env python3
"""Purpose: write observed.npz as observed_X.csv.gz so R reads it without reticulate.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/export_observed.py <sim_dir>
Inputs:  <sim_dir>/seed_<n>/<stage>/observed.npz (never the truth file)
Outputs: <sim_dir>/seed_<n>/<stage>/observed_X.csv.gz, one column per marker, rows in cell_metadata order
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

for obs_path in sorted(Path(sys.argv[1]).glob("seed_*/S*/observed.npz")):
    obs = np.load(obs_path)
    out = obs_path.with_name("observed_X.csv.gz")
    if out.exists():
        continue  # written earlier; never overwritten
    pd.DataFrame(obs["X"], columns=[str(m) for m in obs["markers"]]).to_csv(out, index=False)
