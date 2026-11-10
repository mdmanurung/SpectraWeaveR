#!/usr/bin/env python3
"""Purpose: 90% coverage per stage for effects backed by one vs 2+ reference donors, one point per seed.

Run:     <venv python> -I experiments/benchmarks/05_cytether-synthetic/fig_coverage_class.py <out_png> <scores.csv> [<scores.csv> ...] -- <method> ...
Inputs:  scores.csv from evaluate.py (with the cov90_nd1 / cov90_nd2p diagnostic columns)
Outputs: <out_png> and <out_png stem>_source.csv (plotted values)
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

NOMINAL = 0.90
CLASSES = {"cov90_nd1": "1 reference donor", "cov90_nd2p": "2+ reference donors"}

args = sys.argv[1:]
cut = args.index("--")
out = Path(args[0])
methods = args[cut + 1:]
scores = pd.concat([pd.read_csv(f) for f in args[1:cut]])
d = scores[scores.method.isin(methods)].melt(
    id_vars=["seed", "stage", "method"], value_vars=list(CLASSES), var_name="cls", value_name="cov90_class"
).dropna(subset=["cov90_class"])
stages = sorted(d.stage.unique(), key=lambda s: (int("".join(c for c in s if c.isdigit())), s))
colours = dict(zip(methods, ["#999999", "#1b9e77", "#d95f02", "#7570b3", "#e7298a"]))

fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True, sharey=True)
for ax, (cls, label) in zip(axes, CLASSES.items()):
    for j, m in enumerate(methods):
        dm = d[(d.cls == cls) & (d.method == m)]
        x = dm.stage.map({s: i for i, s in enumerate(stages)}) + (j - (len(methods) - 1) / 2) * 0.18
        ax.scatter(x, dm.cov90_class, s=4, alpha=0.5, color=colours[m], label=m if cls == "cov90_nd1" else None)
        mean = dm.groupby("stage").cov90_class.mean().reindex(stages)
        ax.scatter([i + (j - (len(methods) - 1) / 2) * 0.18 for i in range(len(stages))], mean,
                   marker="_", s=120, color=colours[m])
    ax.axhline(NOMINAL, color="k", lw=0.6, ls="--")
    ax.set_ylabel(f"cov90, {label}", fontsize=7)
    ax.set_ylim(0, 1.02)
axes[1].set_xticks(range(len(stages)))
axes[1].set_xticklabels(stages, fontsize=7)
axes[0].legend(fontsize=6, frameon=False, ncol=len(methods), loc="lower left")
axes[0].set_title("90% coverage per seed (points) and mean (bar); dashed = nominal", fontsize=8, loc="left")
fig.tight_layout()
fig.savefig(out, dpi=200)
d.to_csv(out.with_name(out.stem + "_source.csv"), index=False)
