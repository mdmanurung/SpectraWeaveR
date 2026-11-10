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
from scipy.stats import t as student_t

bundle = Path(os.environ.get("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")) / "cytether_bundle_v1"
sys.path.insert(0, str(bundle / "src"))
from cytether_bench.evaluate import metric  # noqa: E402

sim_dir, results_dir, out_csv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
methods = sys.argv[4:]

Z90, Z95 = 1.6449, 1.9600  # normal quantiles, kept only for the cov90_normal comparison column


def changed_cells(x, obs):
    """Cells whose values the method changed."""
    return np.abs(x - obs).max(axis=1) > 1e-9


def effect_coverage(eff, obs, truth, meta, markers):
    """Coverage of the engine's intervals for the true effect of each batch relative to the anchor.

    True effect of batch b, population k, marker m: mean(obs - truth) over cells of (b, k) minus the same
    over the anchor's cells of k. Cell noise averages out; this is the technical shift the engine estimates.
    """
    anchor = eff.loc[eff.support == "anchor", "batch"].iloc[0]
    tech = pd.DataFrame(obs - truth, columns=markers)
    tech["batch"], tech["population"] = meta.batch.to_numpy(), meta.cell_type.to_numpy()
    shift = tech.groupby(["batch", "population"]).mean()
    out = []
    # t quantiles with the engine's own df (user decision 2026-11-09); an engine without df gets no coverage.
    if "df" not in eff.columns:
        return {"n_effects": 0}
    ok = (eff.support != "anchor") & eff.se_donor.notna() & (eff.se_donor > 0) & eff.df.notna() & (eff.df > 0)
    out = []
    for r in eff[ok].itertuples():
        try:
            true = shift.loc[(r.batch, r.population), r.marker] - shift.loc[(anchor, r.population), r.marker]
        except KeyError:
            continue
        out.append((r.estimate - true, abs(r.estimate - true) / r.se_donor, r.df, r.n_donors, true))
    if not out:
        return {"n_effects": 0}
    err, z, df, nd, true = np.array(out).T
    hit = z <= student_t.ppf(0.95, df)
    res = {"n_effects": len(z), "cov90": float(np.mean(hit)),
           "cov95": float(np.mean(z <= student_t.ppf(0.975, df))),
           "cov90_normal": float(np.mean(z <= Z90)), "median_df": float(np.median(df)),
           "effect_rmse": float(np.sqrt(np.mean(err ** 2)))}
    # Diagnostic only (amendment 2026-11-10, not scored): split by reference donors behind the effect.
    # nd0: internal-anchor-only effects (plan A0); before anchors every effect had >= 1 donor, so nd1 is unchanged.
    for name, m in {"nd0": nd == 0, "nd1": nd == 1, "nd2p": nd >= 2}.items():
        res[f"n_{name}"] = int(m.sum())
        res[f"cov90_{name}"] = float(np.mean(hit[m])) if m.any() else np.nan
        res[f"rmse_{name}"] = float(np.sqrt(np.mean(err[m] ** 2))) if m.any() else np.nan
        res[f"medz_{name}"] = float(np.median(z[m])) if m.any() else np.nan
        # Abstaining leaves the true effect in place: its RMSE is the root mean square of the true effects.
        res[f"rmse_abstain_{name}"] = float(np.sqrt(np.mean(true[m] ** 2))) if m.any() else np.nan
    return res


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
            if method == "anchored_oracle":
                # Floor for any reference-anchored method: the truth plus the anchor batch's own technical offset,
                # per (population, marker). Anchor by the Q8 rule on initial non-holdout reference cells.
                calib = meta[meta.is_reference & ~meta.is_holdout_reference & (meta.arrival == "initial")]
                if calib.empty:
                    rows.append({"seed": seed, "stage": stage_dir.name, "method": method, "supported": False})
                    continue
                donors = calib.groupby("batch").donor_id.nunique().sort_index()
                anchor = donors[donors == donors.max()].index[0]
                obs = pd.read_csv(stage_dir / "observed_X.csv.gz").to_numpy(dtype="float64")
                x = truth.astype("float64").copy()
                for pop in meta.cell_type.unique():
                    rows_a = ((meta.batch == anchor) & (meta.cell_type == pop)).to_numpy()
                    if rows_a.sum() >= 3:
                        x[(meta.cell_type == pop).to_numpy()] += (obs[rows_a] - truth[rows_a]).mean(axis=0)
                row = {"seed": seed, "stage": stage_dir.name, "method": method, "supported": True, "abstained_pairs": 0,
                       "heldout_ref_rmse": float(np.sqrt(np.mean((x[holdout] - truth[holdout]) ** 2))) if holdout.any() else np.nan}
                row.update(metric(x, truth, meta))
                rows.append(row)
                continue
            res = results_dir / f"seed_{seed}" / method / stage_dir.name
            status = pd.read_csv(res / "status.csv").iloc[0]
            row = {"seed": seed, "stage": stage_dir.name, "method": method,
                   "supported": bool(status.supported), "abstained_pairs": status.abstained_pairs,
                   "historical_max_delta": status.get("historical_max_delta", np.nan),
                   "min_cells": status.get("min_cells", np.nan),
                   "fits_on_heldout": status.get("fits_on_heldout", False),
                   "ref_mode_ok": status.get("ref_mode_ok", np.nan)}
            if status.supported:
                x = pd.read_csv(res / "corrected_X.csv.gz").to_numpy(dtype="float64")
                assert x.shape == truth.shape, f"{res}: shape {x.shape} vs truth {truth.shape}"
                row["heldout_ref_rmse"] = float(np.sqrt(np.mean((x[holdout] - truth[holdout]) ** 2))) if holdout.any() else np.nan
                row.update(metric(x, truth, meta))
                obs_df = pd.read_csv(stage_dir / "observed_X.csv.gz")
                obs = obs_df.to_numpy(dtype="float64")
                moved = changed_cells(x, obs)
                row["heldout_share_changed"] = float(moved[holdout].mean()) if holdout.any() else np.nan
                # RMSE on held-out cells that this method and paired_shift both corrected: separates coverage from accuracy.
                ps_file = results_dir / f"seed_{seed}" / "paired_shift" / stage_dir.name / "corrected_X.csv.gz"
                if holdout.any() and ps_file.exists():
                    both = moved & changed_cells(pd.read_csv(ps_file).to_numpy(dtype="float64"), obs) & holdout
                    row["heldout_ref_rmse_common_with_ps"] = float(np.sqrt(np.mean((x[both] - truth[both]) ** 2))) if both.any() else np.nan
                    row["n_common_cells"] = int(both.sum())
                # Unsupported states must stay untouched: Adaptive_NK in S5, future-study cells in S12.
                if stage_dir.name == "S5":
                    cells = (meta.cell_type == "Adaptive_NK").to_numpy()
                    row["unsupported_unchanged_share"] = float(1 - moved[cells].mean()) if cells.any() else np.nan
                if stage_dir.name == "S12":
                    cells = (meta.arrival == "future").to_numpy()
                    row["unsupported_unchanged_share"] = float(1 - moved[cells].mean()) if cells.any() else np.nan
                if (res / "effects.csv.gz").exists():
                    row.update(effect_coverage(pd.read_csv(res / "effects.csv.gz"), obs, truth, meta, list(obs_df.columns)))
            rows.append(row)

pd.DataFrame(rows).to_csv(out_csv, index=False)
