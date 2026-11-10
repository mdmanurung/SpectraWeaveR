#!/usr/bin/env python3
"""Purpose: test whether A9 anchor-mode mismatch explains the anchor arm's effect errors (evaluation only: reads truth).

Run:     <python> experiments/benchmarks/05_cytether-synthetic/diag_mismatch.py <sim_dir> <stage> <results_dir> <modes_csv> <out_dir>
Inputs:  <sim_dir>/seed_<n>/<stage>/{observed_X.csv.gz, cell_metadata.csv.gz, truth_EVALUATION_ONLY.npz}
         <results_dir>/seed_<n>/<method>/<stage>/effects.csv.gz for the lead and the anchor arm
         <modes_csv> from diag_a9_modes.R (A9 mode per seed, population, marker)
Outputs (in <out_dir>, which must not exist):
  effects_mismatch.csv   one row per scored effect: seed, arm, class (nd1, nd0), batch, population, marker, error, mismatch
  slopes.csv             per arm and class: pooled OLS slope, intercept, R2, mean per-seed slope with seed-bootstrap 95% CI, verdict
  counts.csv             effects per seed, arm and class (check against n_nd1 / n_nd0 in scores.csv)
  fig_mismatch.png/.pdf  error against mismatch per panel, fitted pooled line
  fig_mismatch_source.csv  the plotted points, named by panel
Mismatch (C1, user 2026-11-10): the population's true effect minus the cell-weighted mean true effect of all populations
sharing its A9 mode, per batch x marker. Error = estimate - true effect. Truth as in evaluate.py.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parameters ---------------------------------------------------------------------
ARMS = {
    "lead": "heb_pool_pop_kv_bap",
    "anchor": "heb_pool_pop_kv_bap_anc",
    "anchor_d": "heb_pool_pop_kv_bap_anc_d",  # plan A11; skipped where not run
}
N_BOOT = 10000  # seed resamples for the mean-slope interval
BOOT_SEED = 20261109  # same as summarize.py
R2_MIN = 0.5  # C6: mechanism supported when the anchor nd1 slope CI is below 0, pooled R2 >= this, lead slope CI covers 0


def main():
    sim_dir, stage, results_dir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    modes = pd.read_csv(sys.argv[4])
    out = Path(sys.argv[5])
    out.mkdir(parents=True, exist_ok=False)

    rows = []
    for seed, m_seed in modes.groupby("seed"):
        stage_dir = sim_dir / f"seed_{seed}" / stage
        truth = np.load(stage_dir / "truth_EVALUATION_ONLY.npz")["X_true"]
        obs_df = pd.read_csv(stage_dir / "observed_X.csv.gz")
        markers = list(obs_df.columns)
        meta = pd.read_csv(stage_dir / "cell_metadata.csv.gz")
        anchor = m_seed.anchor.iloc[0]

        # True technical shift per batch x population (evaluate.py), then relative to the anchor batch.
        tech = pd.DataFrame(obs_df.to_numpy(dtype="float64") - truth, columns=markers)
        tech["batch"], tech["population"] = (
            meta.batch.to_numpy(),
            meta.cell_type.to_numpy(),
        )
        shift = tech.groupby(["batch", "population"]).mean()
        true_eff = (
            shift
            - shift.xs(anchor, level="batch").reindex(shift.index, level="population")
        ).stack()
        true_eff.index.names = ["batch", "population", "marker"]
        true_eff = true_eff.rename("true").reset_index()

        # Cell weights: the cells anchors are built from (initial, not held out), as in run_method.R.
        fit = meta[(meta.arrival == "initial") & ~meta.is_holdout_reference]
        n_cells = fit.groupby(["batch", "cell_type"]).size().rename("n").reset_index()
        n_cells = n_cells.rename(columns={"cell_type": "population"})

        # Mismatch: own true effect minus the cell-weighted mean over the populations in the same A9 mode.
        mode = m_seed.drop_duplicates(["population", "marker"])[
            ["population", "marker", "mode"]
        ].dropna()
        te = true_eff.merge(
            mode, on=["population", "marker"], how="inner", validate="m:1"
        )
        te = te.merge(
            n_cells, on=["batch", "population"], how="left", validate="m:1"
        ).fillna({"n": 0})
        te["wt"] = te.n * te.true
        grp = te.groupby(["batch", "marker", "mode"])[["wt", "n"]].transform("sum")
        te["mismatch"] = te.true - grp.wt / grp.n

        for arm, method in ARMS.items():
            eff_file = results_dir / f"seed_{seed}" / method / stage / "effects.csv.gz"
            if not eff_file.exists():
                continue  # the arm abstained on this stage (lead on S18)
            eff = pd.read_csv(eff_file)
            # Same effect filter as evaluate.py; the anchor batch has no effect to score.
            ok = (
                (eff.support != "anchor")
                & (eff.se_donor > 0)
                & (eff.df > 0)
                & (eff.batch != anchor)
            )
            e = eff[ok].merge(
                te[["batch", "population", "marker", "true", "mismatch"]],
                on=["batch", "population", "marker"],
                how="inner",
                validate="1:1",
            )
            e["cls"] = np.select(
                [e.n_donors == 1, e.n_donors == 0], ["nd1", "nd0"], default="nd2p"
            )
            e["error"] = e.estimate - e.true
            e["seed"], e["arm"] = seed, arm
            rows.append(
                e[
                    [
                        "seed",
                        "arm",
                        "cls",
                        "batch",
                        "population",
                        "marker",
                        "error",
                        "mismatch",
                    ]
                ]
            )

    em = pd.concat(rows, ignore_index=True)
    em.to_csv(out / "effects_mismatch.csv", index=False)
    em.groupby(["seed", "arm", "cls"]).size().rename("n_effects").reset_index().to_csv(
        out / "counts.csv", index=False
    )

    # Slopes: pooled OLS, then per-seed slopes with a seed-bootstrap interval of their mean.
    rng = np.random.default_rng(BOOT_SEED)
    res = []
    for (arm, cls), d in em.groupby(["arm", "cls"]):
        slope, icpt = np.polyfit(d.mismatch, d.error, 1)
        r2 = np.corrcoef(d.mismatch, d.error)[0, 1] ** 2
        per_seed = d.groupby("seed")[["mismatch", "error"]].apply(
            lambda g: np.polyfit(g.mismatch, g.error, 1)[0]
        )
        boots = [
            rng.choice(per_seed.to_numpy(), size=len(per_seed)).mean()
            for _ in range(N_BOOT)
        ]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        res.append(
            {
                "arm": arm,
                "cls": cls,
                "n_effects": len(d),
                "n_seeds": len(per_seed),
                "pooled_slope": slope,
                "pooled_intercept": icpt,
                "pooled_r2": r2,
                "mean_seed_slope": per_seed.mean(),
                "ci_low": lo,
                "ci_high": hi,
                "mismatch_sd": d.mismatch.std(),
                "error_rmse": float(np.sqrt(np.mean(d.error**2))),
            }
        )
    slopes = pd.DataFrame(res)
    a = slopes[(slopes.arm == "anchor") & (slopes.cls == "nd1")]
    lead = slopes[(slopes.arm == "lead") & (slopes.cls == "nd1")]
    verdict = (
        not a.empty
        and not lead.empty
        and bool(a.ci_high.iloc[0] < 0)
        and bool(a.pooled_r2.iloc[0] >= R2_MIN)
        and bool(lead.ci_low.iloc[0] <= 0 <= lead.ci_high.iloc[0])
    )
    slopes["mechanism_supported_c6"] = verdict
    slopes.to_csv(out / "slopes.csv", index=False)

    # Figure: one panel per arm x class present, shared axes, pooled fitted line -------------
    panels = [
        (arm, cls)
        for arm, cls in [
            ("lead", "nd1"),
            ("anchor", "nd1"),
            ("anchor", "nd0"),
            ("anchor_d", "nd1"),
            ("anchor_d", "nd0"),
        ]
        if ((em.arm == arm) & (em.cls == cls)).any()
    ]
    src = pd.concat(
        [
            em[(em.arm == a_) & (em.cls == c_)].assign(panel=f"{a_}_{c_}")
            for a_, c_ in panels
        ]
    )
    src.to_csv(out / "fig_mismatch_source.csv", index=False)
    fig, axes = plt.subplots(
        1,
        len(panels),
        figsize=(2.4 * len(panels), 2.6),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for ax, (arm, cls) in zip(axes[0], panels):
        d = em[(em.arm == arm) & (em.cls == cls)]
        s = slopes[(slopes.arm == arm) & (slopes.cls == cls)].iloc[0]
        ax.scatter(d.mismatch, d.error, s=2, alpha=0.25, rasterized=True)
        xs = np.linspace(em.mismatch.min(), em.mismatch.max(), 2)
        ax.plot(xs, s.pooled_intercept + s.pooled_slope * xs, color="black", lw=1)
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="grey", lw=0.5)
        ax.set_title(
            f"{ARMS[arm]}\n{cls}: slope {s.mean_seed_slope:.2f} [{s.ci_low:.2f}, {s.ci_high:.2f}], R2 {s.pooled_r2:.2f}",
            fontsize=6,
        )
        ax.set_xlabel("A9 mode mismatch of true effect (asinh)", fontsize=7)
    axes[0][0].set_ylabel("effect error, estimate - truth (asinh)", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fig_mismatch.png", dpi=200)
    fig.savefig(out / "fig_mismatch.pdf")


if __name__ == "__main__":
    main()
