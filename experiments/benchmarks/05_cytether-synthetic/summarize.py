#!/usr/bin/env python3
"""Purpose: summarise evaluate.py scores into the benchmark card's table and figure.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/summarize.py <scores_csv> <out_dir> <reference_method>
Inputs:  <scores_csv> from evaluate.py
Outputs (in <out_dir>):
  primary_by_seed.csv     heldout_ref_rmse per seed, scenario and method, plus the mean over scenarios per seed
  primary_summary.csv     per method: mean over seeds, seed SD, paired difference to <reference_method> with a seed-level bootstrap interval
  guardrails.csv          abstention on S3, S5, S12 with a content check (unsupported cells unchanged)
  protected_contrast.csv  per-seed ratio to <reference_method> with seed interval (proposed rule: upper bound <= 1.10)
  coverage.csv            90% and 95% interval coverage of each engine's own effects (proposed rule: within 5 points)
  invariance.csv          largest historical_max_delta per method (must be 0)
  primary.png             one point per seed over the mean, per method
"""
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Parameters --------------------------------------------------------------------
PRIMARY_STAGES = os.environ.get("SW_PRIMARY_STAGES", "S4,S9,S10,S11").split(",")  # benchmark card default; the shrinkage scenarios set their own
ABSTAIN_STAGES = ["S3", "S5", "S12"]  # guardrail: must abstain or leave values unchanged
N_BOOT = 10000  # resamples of seeds for the paired interval
BOOT_SEED = 20261109  # fixed so the interval reproduces

scores = pd.read_csv(sys.argv[1])
# Scenario sets differ in which columns exist (no S5 or S12, no effects file): keep the summary running.
for col in ["unsupported_unchanged_share", "heldout_share_changed", "heldout_ref_rmse_common_with_ps", "n_effects",
            "cov90", "cov95", "effect_rmse", "fits_on_heldout", "ref_mode_ok", "protected_contrast_mae", "historical_max_delta"]:
    if col not in scores.columns:
        scores[col] = np.nan
out = Path(sys.argv[2])
ref = sys.argv[3]
out.mkdir(parents=True, exist_ok=True)

# Primary: mean over the primary stages per seed -------------------------------------
prim = scores[scores.stage.isin(PRIMARY_STAGES)]
by_seed = prim.pivot_table(index=["seed", "method"], columns="stage", values="heldout_ref_rmse").reset_index()
by_seed["mean_primary"] = by_seed[PRIMARY_STAGES].mean(axis=1, skipna=False)  # a missing stage makes the seed NaN
by_seed.to_csv(out / "primary_by_seed.csv", index=False)

def seed_boot(diff):
    rng = np.random.default_rng(BOOT_SEED)
    boots = [rng.choice(diff, size=len(diff), replace=True).mean() for _ in range(N_BOOT)] if len(diff) > 1 else [np.nan]
    return np.percentile(boots, 2.5), np.percentile(boots, 97.5)

wide = by_seed.pivot(index="seed", columns="method", values="mean_primary")
# Per method: share of held-out cells changed and RMSE on cells that the method and the reference both corrected.
extra = prim.groupby(["seed", "method"]).agg(
    share_changed=("heldout_share_changed", "mean"), common=("heldout_ref_rmse_common_with_ps", "mean")).reset_index()
extra = extra.groupby("method").agg(share_changed=("share_changed", "mean"), rmse_common_with_ref=("common", "mean"))
flags = scores.groupby("method").agg(
    fits_on_heldout=("fits_on_heldout", "max"), ref_mode_ok_share=("ref_mode_ok", "mean"))
rows = []
for method in wide.columns:
    diff = (wide[method] - wide[ref]).dropna().to_numpy()
    lo, hi = seed_boot(diff)
    rows.append({
        "method": method, "n_seeds": int(wide[method].notna().sum()),
        "mean_primary": wide[method].mean(), "sd_over_seeds": wide[method].std(),
        f"mean_diff_vs_{ref}": diff.mean() if len(diff) else np.nan, "diff_ci_low": lo, "diff_ci_high": hi,
    })
summary = pd.DataFrame(rows).merge(extra, on="method", how="left").merge(flags, on="method", how="left")
summary.sort_values("mean_primary").to_csv(out / "primary_summary.csv", index=False)

# Per stage: mean over seeds with the seed-level paired difference to the reference ----------------
stage_rows = []
for stage, d in scores[scores.supported].groupby("stage"):
    w = d.pivot(index="seed", columns="method", values="heldout_ref_rmse")
    c = d.pivot(index="seed", columns="method", values="heldout_ref_rmse_common_with_ps")
    for method in w.columns:
        diff = (w[method] - w[ref]).dropna().to_numpy() if ref in w.columns else np.array([])
        lo, hi = seed_boot(diff) if len(diff) > 1 else (np.nan, np.nan)
        stage_rows.append({"stage": stage, "method": method, "rmse": w[method].mean(), "rmse_common_with_ref": c[method].mean() if method in c else np.nan,
                           f"diff_vs_{ref}": diff.mean() if len(diff) else np.nan, "diff_ci_low": lo, "diff_ci_high": hi})
pd.DataFrame(stage_rows).to_csv(out / "by_stage.csv", index=False)

# Guardrails -------------------------------------------------------------------------
# Abstention content: unsupported states must stay untouched (S5 Adaptive_NK, S12 future-study cells).
g = scores[scores.stage.isin(ABSTAIN_STAGES)].groupby(["method", "stage"]).agg(
    supported_share=("supported", "mean"), abstained_pairs=("abstained_pairs", "mean"),
    unsupported_unchanged_share=("unsupported_unchanged_share", "mean")).reset_index()
g["content_pass"] = (g.unsupported_unchanged_share == 1.0) | g.unsupported_unchanged_share.isna() & (g.supported_share == 0)
g.to_csv(out / "guardrails.csv", index=False)

# Protected contrast: per-seed ratio to the reference over S4 and S9, seed interval, proposed rule upper bound <= 1.10 (Q7 margin).
pcs = scores[scores.stage.isin(["S4", "S9"])].pivot_table(index="seed", columns="method", values="protected_contrast_mae")
pc_rows = []
for method in pcs.columns:
    ratio = (pcs[method] / pcs[ref]).dropna().to_numpy()
    lo, hi = seed_boot(ratio)
    pc_rows.append({"method": method, "mean_ratio": ratio.mean(), "ci_low": lo, "ci_high": hi, "pass_proposed_upper_le_1.10": bool(hi <= 1.10)})
pd.DataFrame(pc_rows).to_csv(out / "protected_contrast.csv", index=False)

# Interval coverage of the engine's own effects on S4 and S9 (proposed rule: within 5 points of nominal; needs confirm).
cov = scores[scores.stage.isin(PRIMARY_STAGES) & scores.n_effects.notna() & (scores.n_effects > 0)].groupby("method").agg(
    n_effects=("n_effects", "sum"), cov90=("cov90", "mean"), cov95=("cov95", "mean"), effect_rmse=("effect_rmse", "mean")).reset_index()
cov["pass_proposed_within_0.05"] = ((cov.cov90 - 0.90).abs() <= 0.05) & ((cov.cov95 - 0.95).abs() <= 0.05)
cov.to_csv(out / "coverage.csv", index=False)

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
