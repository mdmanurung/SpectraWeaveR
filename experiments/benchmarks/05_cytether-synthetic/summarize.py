#!/usr/bin/env python3
"""Purpose: summarise evaluate.py scores into the benchmark card's table and figure.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/summarize.py <scores_csv> <out_dir> <reference_method>
Inputs:  <scores_csv> from evaluate.py
Outputs (in <out_dir>):
  primary_by_seed.csv     heldout_ref_rmse per seed, scenario and method, plus the mean over scenarios per seed
  primary_summary.csv     per method: mean over seeds, seed SD, paired difference to <reference_method> with a seed-level bootstrap interval
  guardrails.csv          abstention on S3, S5, S12; protected-contrast ratio to <reference_method>
  invariance.csv          largest historical_max_delta per method (must be 0)
  primary.png             one point per seed over the mean, per method
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Parameters --------------------------------------------------------------------
PRIMARY_STAGES = ["S4", "S9", "S10", "S11"]  # benchmark card: stages with reference donors
ABSTAIN_STAGES = ["S3", "S5", "S12"]  # guardrail: must abstain or leave values unchanged
N_BOOT = 10000  # resamples of seeds for the paired interval
BOOT_SEED = 20261109  # fixed so the interval reproduces

scores = pd.read_csv(sys.argv[1])
out = Path(sys.argv[2])
ref = sys.argv[3]
out.mkdir(parents=True, exist_ok=True)

# Primary: mean over the primary stages per seed -------------------------------------
prim = scores[scores.stage.isin(PRIMARY_STAGES)]
by_seed = prim.pivot_table(index=["seed", "method"], columns="stage", values="heldout_ref_rmse").reset_index()
by_seed["mean_primary"] = by_seed[PRIMARY_STAGES].mean(axis=1)
by_seed.to_csv(out / "primary_by_seed.csv", index=False)

wide = by_seed.pivot(index="seed", columns="method", values="mean_primary")
rng = np.random.default_rng(BOOT_SEED)
rows = []
for method in wide.columns:
    diff = (wide[method] - wide[ref]).dropna().to_numpy()
    boots = [rng.choice(diff, size=len(diff), replace=True).mean() for _ in range(N_BOOT)] if len(diff) > 1 else [np.nan]
    rows.append({
        "method": method, "n_seeds": int(wide[method].notna().sum()),
        "mean_primary": wide[method].mean(), "sd_over_seeds": wide[method].std(),
        f"mean_diff_vs_{ref}": diff.mean() if len(diff) else np.nan,
        "diff_ci_low": np.percentile(boots, 2.5), "diff_ci_high": np.percentile(boots, 97.5),
    })
pd.DataFrame(rows).sort_values("mean_primary").to_csv(out / "primary_summary.csv", index=False)

# Guardrails -------------------------------------------------------------------------
g = scores[scores.stage.isin(ABSTAIN_STAGES)].groupby(["method", "stage"]).agg(
    supported_share=("supported", "mean"), abstained_pairs=("abstained_pairs", "mean")).reset_index()
pc = scores[scores.stage.isin(["S4", "S9"])].groupby(["method"]).protected_contrast_mae.mean()
g["protected_contrast_ratio_to_ref"] = g.method.map(pc / pc.get(ref, np.nan))
g.to_csv(out / "guardrails.csv", index=False)
inv = scores[scores.supported].groupby("method").historical_max_delta.max().reset_index()
inv.to_csv(out / "invariance.csv", index=False)

# Figure: one point per seed over the mean ---------------------------------------------
order = list(pd.read_csv(out / "primary_summary.csv").method)
fig, ax = plt.subplots(figsize=(3.5, 2.8))
for i, method in enumerate(order):
    y = wide[method].dropna()
    ax.scatter(np.full(len(y), i) + np.linspace(-0.15, 0.15, len(y)), y, s=8, alpha=0.7)
    ax.hlines(y.mean(), i - 0.25, i + 0.25, color="black", lw=1.5)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=45, ha="right", fontsize=7)
ax.set_ylim(bottom=0)
ax.set_ylabel("held-out reference RMSE to truth\n(mean over S4, S9, S10, S11)", fontsize=7)
fig.tight_layout()
fig.savefig(out / "primary.png", dpi=200)
