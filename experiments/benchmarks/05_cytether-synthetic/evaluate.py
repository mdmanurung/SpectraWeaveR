#!/usr/bin/env python3
"""Purpose: score corrected matrices against the simulation truth. The only script that reads the truth file.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/evaluate.py <sim_dir> <results_dir> <out_csv> <method> [<method> ...]
Inputs:  <sim_dir>/seed_<n>/<stage>/{truth_EVALUATION_ONLY.npz,cell_metadata.csv.gz}
         <results_dir>/seed_<n>/<method>/<stage>/{corrected_X.csv.gz,status.csv}
Outputs: <out_csv>, one row per seed, stage and method:
  heldout_ref_rmse   RMSE of corrected held-out reference cells to the true expression (F2 primary, lower is better)
  plus the bundle's participant_rmse, donor_contrast_mae, protected_contrast_mae, future_study_rmse and
  heldout_reference_across_batch_mad, computed with the bundle's own metric().
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

bundle = Path(os.environ.get("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")) / "cytether_bundle_v1"
sys.path.insert(0, str(bundle / "src"))
from cytether_bench.evaluate import metric  # noqa: E402

sim_dir, results_dir, out_csv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
methods = sys.argv[4:]

rows = []
for seed_dir in sorted(sim_dir.glob("seed_*")):
    seed = int(seed_dir.name.split("_")[1])
    if not (results_dir / f"seed_{seed}").exists():
        continue  # only seeds that were run
    for stage_dir in sorted(seed_dir.glob("S*")):
        truth = np.load(stage_dir / "truth_EVALUATION_ONLY.npz")["X_true"]
        meta = pd.read_csv(stage_dir / "cell_metadata.csv.gz")
        holdout = meta.is_holdout_reference.to_numpy()
        for method in methods:
            if method == "oracle":
                # Sanity candidate: the truth itself must score zero error.
                row = {"seed": seed, "stage": stage_dir.name, "method": method, "supported": True, "abstained_pairs": 0,
                       "heldout_ref_rmse": 0.0 if holdout.any() else np.nan}
                row.update(metric(truth, truth, meta))
                rows.append(row)
                continue
            res = results_dir / f"seed_{seed}" / method / stage_dir.name
            status = pd.read_csv(res / "status.csv").iloc[0]
            row = {"seed": seed, "stage": stage_dir.name, "method": method,
                   "supported": bool(status.supported), "abstained_pairs": status.abstained_pairs,
                   "historical_max_delta": status.get("historical_max_delta", np.nan)}
            if status.supported:
                x = pd.read_csv(res / "corrected_X.csv.gz").to_numpy(dtype="float64")
                assert x.shape == truth.shape, f"{res}: shape {x.shape} vs truth {truth.shape}"
                row["heldout_ref_rmse"] = float(np.sqrt(np.mean((x[holdout] - truth[holdout]) ** 2))) if holdout.any() else np.nan
                row.update(metric(x, truth, meta))
            rows.append(row)

pd.DataFrame(rows).to_csv(out_csv, index=False)
