#!/usr/bin/env python3
"""Purpose: figure and reading rule for the T6 OpenProblems scores (plan O8); reads 05_scale_op.py output.

Run:     <venv python> -I experiments/benchmarks/06_cytether-cytobatchflagr/07b_fig_op.py <score_dir>
Inputs:  <score_dir>/scores_scaled_all.csv (fold, method, metric, value; metric "mean" = the task's mean score)
Outputs (in <score_dir>):
  fig_op.png/.pdf, fig_op_source_A.csv (mean scaled score per fold), fig_op_source_B.csv (fold-mean per metric)
  reading_rule.csv  per fold: sc_cl and cycombine_ref mean; pooled biology-metric gaps; verdict (O8 Q3)
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parameters ---------------------------------------------------------------------
ORDER = [
    "none",
    "cytofin",
    "cytofin_pop",
    "paired_shift",
    "heb",
    "heb_kv_bap",
    "heb_kv_bap_sc",
    "heb_kv_bap_sc_cl",
    "cytonorm",
    "cycombine_ref",
    "cycombine",
]
BIO = [
    "cLisi",
    "functional_marker_preservation_cohens_d",
    "functional_marker_preservation_wilcoxon",
    "ratio_consistent_peaks",
]  # O8 Q3
BIO_MARGIN = 0.10  # O8 Q3: no biology metric more than this below cycombine_ref, pooled
MIN_FOLDS = 3  # O8 Q3: sc_cl >= cycombine_ref on the mean score in at least this many of 4 folds
CAND, REF = "heb_kv_bap_sc_cl", "cycombine_ref"
FOLD_COLOURS = {
    "op_tcu1_tcs3": "#7570b3",
    "op_tcu1_tcs4": "#e7298a",
    "op_tcu2_tcs3": "#1b9e77",
    "op_tcu2_tcs4": "#d95f02",
}


def main():
    d = Path(sys.argv[1])
    s = pd.read_csv(d / "scores_scaled_all.csv")
    s = s[s.method.isin(ORDER)]
    order = [m for m in ORDER if m in set(s.method)]

    # Reading rule -----------------------------------------------------------------
    mean = s[s.metric == "mean"].pivot(index="fold", columns="method", values="value")
    per_metric = (
        s[s.metric != "mean"].groupby(["method", "metric"]).value.mean().unstack()
    )
    bio_gap = per_metric.loc[CAND, BIO] - per_metric.loc[REF, BIO]
    rule = pd.DataFrame({"cand_mean": mean[CAND], "ref_mean": mean[REF]})
    rule["cand_ge_ref"] = rule.cand_mean >= rule.ref_mean
    for m in BIO:
        rule[f"pooled_gap_{m}"] = bio_gap[m]
    rule["passes"] = rule.cand_ge_ref.sum() >= MIN_FOLDS and bool(
        (bio_gap >= -BIO_MARGIN).all()
    )
    rule.to_csv(d / "reading_rule.csv")

    # Figure -----------------------------------------------------------------------
    src_a = s[s.metric == "mean"][["fold", "method", "value"]]
    src_a.to_csv(d / "fig_op_source_A.csv", index=False)
    metrics = [m for m in per_metric.columns if m not in BIO] + BIO
    heat = per_metric.loc[order, metrics]
    heat.to_csv(d / "fig_op_source_B.csv")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1, 1.5]}
    )
    rng = np.random.default_rng(1)
    for f, g in src_a.groupby("fold"):
        x = [order.index(m) for m in g.method]
        ax1.scatter(
            np.array(x) + rng.uniform(-0.15, 0.15, len(x)),
            g.value,
            s=12,
            color=FOLD_COLOURS[f],
            label=f,
        )
    med = src_a.groupby("method").value.median().reindex(order)
    ax1.hlines(
        med,
        np.arange(len(order)) - 0.3,
        np.arange(len(order)) + 0.3,
        color="black",
        lw=1,
    )
    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels(order, rotation=40, ha="right", fontsize=6)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel("mean scaled score (0 to 1, higher is better)", fontsize=7)
    ax1.set_title(
        "A  OpenProblems mean score per fold; bar = median", fontsize=8, loc="left"
    )
    ax1.legend(fontsize=6, frameon=False)
    im = ax2.imshow(heat.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels(order, fontsize=6)
    ax2.set_xticks(range(len(metrics)))
    ax2.set_xticklabels(metrics, rotation=40, ha="right", fontsize=6)
    ax2.axvline(len(metrics) - len(BIO) - 0.5, color="white", lw=2)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax2.text(
                j,
                i,
                f"{heat.iat[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=5,
                color="white",
            )
    ax2.set_title(
        "B  scaled metric, mean over 4 folds (right of line: biology, O8)",
        fontsize=8,
        loc="left",
    )
    fig.colorbar(im, ax=ax2, fraction=0.03)
    fig.tight_layout()
    fig.savefig(d / "fig_op.png", dpi=200)
    fig.savefig(d / "fig_op.pdf")


if __name__ == "__main__":
    main()
