#!/usr/bin/env python3
"""Purpose: figure and stop rule for the T6 EMD decomposition (plan T6 S3-S4); reads 08_emd_decompose.R outputs.

Run:     <python> experiments/benchmarks/06_cytether-cytobatchflagr/09_fig_emd_decompose.py <decompose_dir> <emd_v5_dir>
Inputs:  <decompose_dir>/fold_<f>.csv (4 folds); <emd_v5_dir>/emd_long.csv (check: totals must match)
Outputs (in <decompose_dir>):
  summary.csv            median total and shape EMD per method and fold
  stop_rule.csv          per fold and pooled: HEB - cycombine_ref gap in median total and shape EMD, shape share, verdict
  fig_emd_decompose.png/.pdf, fig_emd_decompose_source_A.csv, fig_emd_decompose_source_B.csv
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parameters ---------------------------------------------------------------------
SHAPE_SHARE_MIN = (
    0.5  # Q3: scale is a candidate cause only if shape explains >= half the gap
)
TOL = 1e-9  # totals must reproduce emd_scores_v5
ORDER = ["none", "cytonorm", "heb_kv_bap", "heb_kv_bap_scale", "cycombine_ref"]


def main():
    d, v5 = Path(sys.argv[1]), Path(sys.argv[2])
    long = pd.concat(
        [pd.read_csv(f) for f in sorted(d.glob("fold_*.csv"))], ignore_index=True
    )

    # Check: totals of the shared methods equal the v5 scores.
    ref = pd.read_csv(v5 / "emd_long.csv")
    chk = long.merge(
        ref, on=["method", "fold", "population", "marker"], how="inner", validate="1:1"
    )
    diff = (chk.emd_total - chk.emd).abs().max()
    if not (len(chk) > 0 and diff < TOL):
        raise ValueError(
            f"totals do not reproduce emd_scores_v5 (max diff {diff}, n {len(chk)})"
        )

    summ = long.groupby(["method", "fold"], as_index=False)[
        ["emd_total", "emd_shape"]
    ].median()
    summ.to_csv(d / "summary.csv", index=False)

    # Stop rule: gap in median EMD, HEB minus cycombine_ref, total and shape.
    w = summ.pivot(index="fold", columns="method")
    rule = pd.DataFrame(
        {
            "gap_total": w.emd_total.heb_kv_bap - w.emd_total.cycombine_ref,
            "gap_shape": w.emd_shape.heb_kv_bap - w.emd_shape.cycombine_ref,
            "scale_arm_gain": w.emd_total.heb_kv_bap - w.emd_total.heb_kv_bap_scale,
        }
    )
    rule.loc["pooled"] = rule.mean()
    rule["shape_share"] = rule.gap_shape / rule.gap_total
    rule["scale_candidate"] = rule.shape_share >= SHAPE_SHARE_MIN
    rule.to_csv(d / "stop_rule.csv")

    # Figure -----------------------------------------------------------------------
    src_a = summ.melt(
        id_vars=["method", "fold"], var_name="part", value_name="median_emd"
    )
    src_a.to_csv(d / "fig_emd_decompose_source_A.csv", index=False)
    pm = long.pivot_table(
        index=["fold", "population", "marker"],
        columns="method",
        values=["emd_total", "emd_shape"],
    )
    src_b = (
        pd.DataFrame(
            {
                "total_gap": pm.emd_total.heb_kv_bap - pm.emd_total.cycombine_ref,
                "shape_gap": pm.emd_shape.heb_kv_bap - pm.emd_shape.cycombine_ref,
            }
        )
        .dropna()
        .reset_index()
    )
    src_b.to_csv(d / "fig_emd_decompose_source_B.csv", index=False)

    folds = sorted(summ.fold.unique())
    colours = dict(zip(folds, ["#7570b3", "#e7298a", "#1b9e77", "#d95f02"]))
    fig, axes = plt.subplots(
        1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1, 1, 1.1]}
    )
    for ax, part, title in [
        (axes[0], "emd_total", "A  total EMD (as scored)"),
        (axes[1], "emd_shape", "B  shape EMD (files centred on median)"),
    ]:
        for f in folds:
            s = summ[summ.fold == f].set_index("method").reindex(ORDER)
            ax.plot(
                range(len(ORDER)),
                s[part],
                marker="o",
                ms=3,
                lw=1,
                color=colours[f],
                label=f"held out {f}",
            )
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels(ORDER, rotation=40, ha="right", fontsize=6)
        ax.set_ylim(0, None)
        ax.set_title(title, fontsize=8, loc="left")
    axes[0].set_ylabel(
        "median EMD over population x marker\n(bin units, bin 0.1; lower is better)",
        fontsize=7,
    )
    axes[1].sharey(axes[0])
    axes[0].legend(fontsize=6, frameon=False)
    ax = axes[2]
    lim = np.nanpercentile(np.abs(src_b[["total_gap", "shape_gap"]].to_numpy()), 99)
    for f in folds:
        g = src_b[src_b.fold == f]
        ax.scatter(
            g.total_gap, g.shape_gap, s=3, alpha=0.4, color=colours[f], rasterized=True
        )
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.6)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("total EMD, heb_kv_bap - cycombine_ref", fontsize=7)
    ax.set_ylabel("shape EMD, heb_kv_bap - cycombine_ref", fontsize=7)
    p = rule.loc["pooled"]
    ax.set_title(
        f"C  per population x marker; pooled shape share {p.shape_share:.2f}",
        fontsize=8,
        loc="left",
    )
    fig.tight_layout()
    fig.savefig(d / "fig_emd_decompose.png", dpi=200)
    fig.savefig(d / "fig_emd_decompose.pdf")


if __name__ == "__main__":
    main()
