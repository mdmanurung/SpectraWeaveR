#!/usr/bin/env python3
"""Purpose: score one control or method of the T6 OpenProblems track (O1 to O7) with the task's 9 scores.

Run:     <venv python> -I experiments/benchmarks/06_cytether-cytobatchflagr/04_score_op.py <prep_dir> <results_dir> <out_dir> <fold> <set>
         set: no_integration, perfect_integration, shuffle_integration_globally, or a method name
         One set per job: one set on ~1.8M cells takes about 2 h.
Inputs:  <prep_dir>/cells_raw.f64, obs.csv, markers.csv, populations_op.i32, split_op.i32 (01b_prep_op.R)
         <results_dir>/<method>/<fold>/split_<k>/{corrected.f64,rows.i32} (02b_correct_op.R), methods only
Outputs: <out_dir>/<fold>/<set>/scores_raw.csv  set, metric, value (unscaled); the folder must not exist
         05_scale_op.py scales them on the controls.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
from openproblems_cyto_metrics import (  # noqa: E402
    control_no_integration,
    control_shuffle_integration_globally,
    score_all,
)

# Parameters (approved, plan T6) -------------------------------------------------
PERFECT_BATCH = "batch01"  # O5
FLOWSOM_XDIM, FLOWSOM_YDIM, FLOWSOM_NCLUS = 10, 10, 20  # O5 (as N5)
SHUFFLE_SEED = 1  # N6
LINEAGE = {"CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CD103", "CCR6"}  # O4; the rest are functional
FOLDS = ["op_tcu1_tcs3", "op_tcu1_tcs4", "op_tcu2_tcs3", "op_tcu2_tcs4"]
RSCRIPT = "/exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript"


def perfect_integration_t6(X_raw, obs, perfect_batch):
    """O7: each evaluated (donor, batch) sample in split k becomes the donor's perfect_batch cells of split k.

    The task's control needs one batch per donor per split; here every donor spans all batches.
    Only evaluated cells are scored, so reference rows are left out of the output.
    """
    Xs, parts = [], []
    ev = obs[obs.is_control == 0]
    for (donor, batch, k), _ in ev.groupby(["donor", "batch", "split"], observed=True):
        take = ((obs.donor == donor) & (obs.batch == perfect_batch) & (obs.split == k)).to_numpy()
        o = obs.loc[take].assign(batch=batch)
        o.index = [f"{i}|{batch}" for i in o.index]
        parts.append(o)
        Xs.append(X_raw[take])
    return np.vstack(Xs), pd.concat(parts)


def main():
    prep, results, out_root, fold, which = (Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]),
                                            sys.argv[4], sys.argv[5])
    out = out_root / fold / which
    if out.exists():
        raise ValueError(f"{out} exists; use a new out_dir")

    markers = pd.read_csv(prep / "markers.csv")
    n_m = len(markers)
    obs = pd.read_csv(prep / "obs.csv")
    X_raw = np.fromfile(prep / "cells_raw.f64", dtype=np.float64).reshape((n_m, -1)).T
    pops = np.fromfile(prep / "populations_op.i32", dtype=np.int32).reshape((len(FOLDS), -1))[FOLDS.index(fold)]
    split = np.fromfile(prep / "split_op.i32", dtype=np.int32)
    if not (len(obs) == X_raw.shape[0] == len(pops) == len(split)):
        raise ValueError("prep exports differ in rows")

    # obs in the task's schema: evaluated = the held-out pair; controls = the two references.
    evaluated = fold.removeprefix("op_").split("_")
    is_eval = obs.control.isin(evaluated).to_numpy()
    obs = obs.assign(
        donor=obs.control, population=pops.astype(str), is_control=(~is_eval).astype(int),
        split=np.where(is_eval, split, 0),
    )
    var = pd.DataFrame({
        "to_correct": True,
        "marker_type": ["lineage" if m in LINEAGE else "functional" for m in markers.label],
    }, index=markers.channel)

    # Corrected values for the chosen set ----------------------------------------------
    obs_corr = None
    if which == "no_integration":
        X_corr, _ = control_no_integration(X_raw, obs)
    elif which == "perfect_integration":
        X_corr, obs_corr = perfect_integration_t6(X_raw, obs, PERFECT_BATCH)
    elif which == "shuffle_integration_globally":
        X_corr, _ = control_shuffle_integration_globally(X_raw, obs, SHUFFLE_SEED)
    else:
        # A method: each evaluated cell takes its value from the run of its own split.
        X_corr = X_raw.copy()
        for k in (1, 2):
            d = results / which / fold / f"split_{k}"
            rows = np.fromfile(d / "rows.i32", dtype=np.int32) - 1
            xk = np.fromfile(d / "corrected.f64", dtype=np.float64).reshape((n_m, -1)).T
            take = is_eval[rows] & (split[rows] == k)
            X_corr[rows[take]] = xk[take]

    vals = score_all(X_raw, X_corr, obs, var, flowsom_xdim=FLOWSOM_XDIM, flowsom_ydim=FLOWSOM_YDIM,
                     flowsom_n_clusters=FLOWSOM_NCLUS, obs_corrected=obs_corr, rscript=RSCRIPT)
    out.mkdir(parents=True)
    pd.DataFrame({"set": which, "metric": list(vals), "value": list(vals.values())}).to_csv(
        out / "scores_raw.csv", index=False)


if __name__ == "__main__":
    main()
