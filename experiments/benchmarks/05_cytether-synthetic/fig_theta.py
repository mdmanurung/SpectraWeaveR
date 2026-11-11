#!/usr/bin/env python3
"""Purpose: plot heb's fitted variance components per arm and stage (plan A13 step 3); never reads truth.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/fig_theta.py <theta_dir>
Inputs:  <theta_dir>/<arm>_<stage>_theta.csv from diag_a9_modes.R (arm lead, anc, anc_d; stage S13, S17, S18)
Outputs (in <theta_dir>):
  fig_theta.png/.pdf       one panel per component; x = stage x arm; pooled sets (seed x marker) as points, median bar
  fig_theta_source.csv     the plotted points, long format, named by panel
  theta_medians.csv        median per arm, stage, component
  stop_rule.csv            S17 th3: lead median, anc and anc_d medians, ratio, and the A13 step 4 verdict
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parameters ---------------------------------------------------------------------
ARMS = ["lead", "anc", "anc_d"]
STAGES = ["S13", "S17", "S18"]
COMPONENTS = ["study", "run", "resid_donor", "resid_anchor_extra", "tau_d"]
STOP_RATIO = 2.0  # A13 step 4 (user Q2): no fix when lead th3 on S17 is within 2x of the anchor arm's
FLOOR = (
    1e-6  # display floor for the log axis; values at the 1e-8 variance floor plot here
)
JITTER_SEED = 1


def main():
    d = Path(sys.argv[1])
    theta = pd.concat(
        [
            pd.read_csv(f).assign(
                arm=f.name.split("_S")[0],
                stage=f"S{f.name.split('_S')[1].split('_')[0]}",
            )
            for f in sorted(d.glob("*_theta.csv"))
        ],
        ignore_index=True,
    )
    long = theta.melt(
        id_vars=["arm", "stage", "seed", "marker"],
        value_vars=COMPONENTS,
        var_name="component",
    ).dropna(subset=["value"])

    med = long.groupby(["arm", "stage", "component"], as_index=False)["value"].median()
    med.to_csv(d / "theta_medians.csv", index=False)

    # Stop rule on S17 resid_donor (th3): lead vs each anchor arm.
    s17 = med[(med.stage == "S17") & (med.component == "resid_donor")].set_index("arm")[
        "value"
    ]
    rule = pd.DataFrame(
        {
            "arm": ["anc", "anc_d"],
            "lead_th3": s17.get("lead", np.nan),
            "arm_th3": [s17.get("anc", np.nan), s17.get("anc_d", np.nan)],
        }
    )
    rule["ratio"] = rule.arm_th3 / rule.lead_th3
    rule["no_fix_stop"] = rule.ratio <= STOP_RATIO
    rule.to_csv(d / "stop_rule.csv", index=False)

    # Figure: one panel per component, groups ordered stage then arm -----------------
    groups = [(s, a) for s in STAGES for a in ARMS]
    xpos = {g: i for i, g in enumerate(groups)}
    long["x"] = [xpos[(s, a)] for s, a in zip(long.stage, long.arm)]
    long["panel"] = long.component
    long.to_csv(d / "fig_theta_source.csv", index=False)
    rng = np.random.default_rng(JITTER_SEED)
    colours = {"lead": "#4C72B0", "anc": "#DD8452", "anc_d": "#55A868"}
    fig, axes = plt.subplots(1, len(COMPONENTS), figsize=(13, 3), sharey=True)
    for ax, comp in zip(axes, COMPONENTS):
        sub = long[long.component == comp]
        for (s, a), g in sub.groupby(["stage", "arm"]):
            x = xpos[(s, a)]
            y = np.maximum(g.value.to_numpy(), FLOOR)
            ax.scatter(
                x + rng.uniform(-0.25, 0.25, len(y)),
                y,
                s=4,
                alpha=0.4,
                color=colours[a],
            )
            ax.hlines(np.median(y), x - 0.35, x + 0.35, color="black", lw=1.2)
        ax.set_yscale("log")
        ax.set_title(comp, fontsize=8)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([f"{s}\n{a}" for s, a in groups], fontsize=5)
        for b in (2.5, 5.5):
            ax.axvline(b, color="grey", lw=0.5)
    axes[0].set_ylabel(f"variance (asinh^2), log; floor shown at {FLOOR:g}", fontsize=7)
    fig.suptitle(
        "heb variance components per pooled set (seed x marker), seeds 441-443; bar = median",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(d / "fig_theta.png", dpi=200)
    fig.savefig(d / "fig_theta.pdf")


if __name__ == "__main__":
    main()
