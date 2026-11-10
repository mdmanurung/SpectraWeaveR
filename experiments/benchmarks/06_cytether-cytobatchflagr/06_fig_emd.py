#!/usr/bin/env python3
"""Purpose: figure for the T6 EMD track: median EMD per method and fold, and HEB minus paired_shift per cell.

Run:     <venv python> -I experiments/benchmarks/06_cytether-cytobatchflagr/06_fig_emd.py <score_dir>
Inputs:  <score_dir>/emd_summary.csv, emd_long.csv, status_all.csv (03_score_emd.R)
Outputs (in <score_dir>):
  fig_emd.png, fig_emd.pdf   A: median EMD per method, one point per held-out donor, lines join folds
                             B: per (fold, population, marker) EMD difference, heb - paired_shift
  fig_emd_source_A.csv, fig_emd_source_B.csv   plotted values
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ORDER = ["none", "cytofin", "paired_shift", "heb", "cytonorm", "cycombine"]
FOLD_COLOURS = {"tcu1": "#1b9e77", "tcu2": "#d95f02", "tcs3": "#7570b3", "tcs4": "#e7298a"}

root = Path(sys.argv[1])
summ = pd.read_csv(root / "emd_summary.csv")
status = pd.read_csv(root / "status_all.csv")
long = pd.read_csv(root / "emd_long.csv")

# Panel frames -------------------------------------------------------------------------
a = summ.merge(status[["method", "fold", "fits_on_heldout"]], on=["method", "fold"], validate="1:1")
a["x"] = a.method.map({m: i for i, m in enumerate(ORDER)})
b = long.pivot_table(index=["fold", "population", "marker"], columns="method", values="emd").reset_index()
b["diff"] = b["heb"] - b["paired_shift"]
b = b.dropna(subset=["diff"])

# Plot ----------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.6, 1]})
for fold, d in a.sort_values("x").groupby("fold"):
    ax1.plot(d.x, d.median_emd, "-o", ms=3, lw=0.6, color=FOLD_COLOURS[fold], label=f"held out {fold}")
ax1.set_xticks(range(len(ORDER)))
ax1.set_xticklabels([m + (" *" if a[a.method == m].fits_on_heldout.any() else "") for m in ORDER],
                    rotation=40, ha="right", fontsize=7)
ax1.set_ylim(bottom=0)
ax1.set_ylabel("median EMD over population x marker\n(bin units, bin 0.1; lower is better)", fontsize=7)
ax1.legend(fontsize=6, frameon=False)
ax1.set_title("A  held-out donor, 18 batches", fontsize=8, loc="left")
for fold, d in b.groupby("fold"):
    ax2.hist(d["diff"], bins=60, range=(-6, 6), histtype="step", color=FOLD_COLOURS[fold], lw=0.8)
ax2.axvline(0, color="k", lw=0.5)
ax2.set_xlabel("EMD heb - paired_shift\n(per population x marker)", fontsize=7)
ax2.set_ylabel("count", fontsize=7)
ax2.set_title("B  heb vs paired_shift", fontsize=8, loc="left")
fig.text(0.01, 0.01, "* fits on the held-out donor's cells", fontsize=6)
fig.tight_layout()
fig.savefig(root / "fig_emd.png", dpi=200)
fig.savefig(root / "fig_emd.pdf")

a.to_csv(root / "fig_emd_source_A.csv", index=False)
b[["fold", "population", "marker", "heb", "paired_shift", "diff"]].to_csv(root / "fig_emd_source_B.csv", index=False)
