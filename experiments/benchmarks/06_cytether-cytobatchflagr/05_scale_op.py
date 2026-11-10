#!/usr/bin/env python3
"""Purpose: scale the T6 OpenProblems raw scores on the three controls (O2) and stack all folds.

Run:     <venv python> -I experiments/benchmarks/06_cytether-cytobatchflagr/05_scale_op.py <score_dir>
Inputs:  <score_dir>/<fold>/<set>/scores_raw.csv (04_score_op.py)
Outputs (in <score_dir>):
  scores_raw_all.csv     fold, set, metric, value
  scores_scaled_all.csv  fold, method, metric, scaled value, plus metric "mean"
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
from openproblems_cyto_metrics import scale_scores  # noqa: E402

CONTROLS = ["no_integration", "perfect_integration", "shuffle_integration_globally"]  # O2

root = Path(sys.argv[1])
raw = pd.concat(
    [pd.read_csv(f).assign(fold=f.parents[1].name) for f in sorted(root.glob("*/*/scores_raw.csv"))]
)
raw.to_csv(root / "scores_raw_all.csv", index=False)

rows = []
for fold, d in raw.groupby("fold"):
    by_set = {s: dict(zip(g.metric, g.value)) for s, g in d.groupby("set")}
    missing = [c for c in CONTROLS if c not in by_set]
    if missing:
        raise ValueError(f"{fold}: controls missing: {missing}")
    controls = {c: by_set[c] for c in CONTROLS}
    for method, vals in by_set.items():
        scaled, mean = scale_scores(vals, controls)
        rows += [{"fold": fold, "method": method, "metric": k, "value": v} for k, v in scaled.items()]
        rows.append({"fold": fold, "method": method, "metric": "mean", "value": mean})
pd.DataFrame(rows).to_csv(root / "scores_scaled_all.csv", index=False)
