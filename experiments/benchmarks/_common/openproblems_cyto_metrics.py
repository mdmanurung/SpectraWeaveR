"""OpenProblems task_cyto_batch_integration metrics (9 scores) in numpy/scipy.

Purpose: score a batch-corrected cell table with the 9 scores of the OpenProblems
    task, and scale them between the control methods the way the task results do.
Run: python openproblems_cyto_metrics.py        (runs the synthetic self-check)
    from openproblems_cyto_metrics import score_all, scale_scores
Inputs:
    X_raw        cells x markers, the task's layer "preprocessed" (unintegrated)
    X_corrected  cells x markers, the task's layer "integrated" (same rows as X_raw,
                 unless obs_corrected is given, see score_all)
    obs          DataFrame, one row per row of X_raw, columns:
                 batch, donor, sample_id, population, is_control, split, group
                 split: 0 control sample, 1 or 2 = replicate split of the task.
                 group: biological group of the donor (2 groups for the
                 functional-marker metric, iLISI needs every group in >= 2 batches).
    var          DataFrame indexed by marker name, columns to_correct (bool) and
                 marker_type ("lineage" or "functional")
Outputs: dict of 9 floats (score_all), dict of scaled floats (scale_scores).

Source of every definition below: openproblems-bio/task_cyto_batch_integration,
tag v0.0.1-rc3, commit TASK_SHA, files cited per metric. The task runs each method
twice (split 1 and split 2, each with the control samples) and scores the two
outputs together, so X_corrected must come from two separate method runs.
Population = the task's cell_type; sample_id = the task's sample.
FlowSOM mapping similarity runs the task's R code through Rscript (FlowSOM, flowCore).
iLISI and cLISI use scib-metrics' formula with an exact kNN (cKDTree) instead of
pynndescent, because scib_metrics and pynndescent are not in the venv.
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde, mannwhitneyu, wasserstein_distance

# Settings ------------------------------------------------------------------------
# Every value below is fixed by the task code (source file named beside it).
# They are the task's choices, not ours: change one and the score no longer matches.
TASK_SHA = "6cdd9e8f6d7e74826af62eaa913e82a723e21769"  # tag v0.0.1-rc3
# src/metrics/emd/helper.py: bins -100..100 step 0.1, skip pairs with < 50 cells
EMD_BIN_EDGES = np.arange(-100, 100.1, 0.1) + 0.0000001
EMD_MIN_CELLS = 50
# src/metrics/average_batch_r2/script.py: skip donor x population with < 20 cells
R2_MIN_CELLS = 20
# src/metrics/ratio_consistent_peaks/{script,helper}.py
RCP_MIN_CELLS = 100
RCP_HEIGHT = 0.1
RCP_PROMINENCE = 0.01
RCP_GRID_POINTS = 100
# src/metrics/lisi/helper.py (kNN) and scib_metrics 0.5.7 utils/_lisi.py (LISI)
LISI_N_NEIGHBORS = 100
LISI_TOL = 1e-5
LISI_MAX_TRIES = 50
# src/metrics/functional_marker_preservation/{script,helper}.py
FMP_MIN_CELLS = 5
FMP_PVAL_THRESHOLD = 0.1
# src/metrics/flowsom_mapping_similarity/script.R
FLOWSOM_SEED = 42
# metric names and directions: config.vsh.yaml of each metric folder
METRIC_MAXIMIZE = {
    "average_batch_r2_ct": False,
    "emd_mean_ct_horiz": False,
    "emd_mean_ct_vert": False,
    "flowsom_mean_mapping_similarity": True,
    "iLisi": True,
    "cLisi": True,
    "functional_marker_preservation_wilcoxon": True,
    "functional_marker_preservation_cohens_d": False,
    "ratio_consistent_peaks": True,
}

# Shared preparation ---------------------------------------------------------------
# src/utils/helper_functions.py: subset_nocontrols, remove_unlabelled,
# subset_markers_tocorrect, get_obs_var_for_integrated


def _rename_obs(obs):
    return obs.rename(columns={"population": "cell_type", "sample_id": "sample"})


def _prep(X, obs, var):
    """Non-control, labelled cells and markers to correct. Returns (X, obs, markers)."""
    obs = _rename_obs(obs)
    if obs["cell_type"].isna().any():
        raise ValueError("population has NA, the task has no rule for it")
    unlabelled = obs["cell_type"].str.lower().isin(["unlabelled", "unlabeled"])
    keep = ((obs["is_control"].astype(int) == 0) & ~unlabelled).to_numpy()
    cols = var["to_correct"].to_numpy(dtype=bool)
    return np.asarray(X)[keep][:, cols], obs.loc[keep], list(var.index[cols])


def _by_split(X, obs, var):
    """Task: integrated split k = control samples + cells with split == k."""
    Xp, op, markers = _prep(X, obs, var)
    out = {}
    for k in (1, 2):
        m = (op["split"] == k).to_numpy()
        out[k] = (Xp[m], op.loc[m])
    return out, markers


# EMD ------------------------------------------------------------------------------
# src/metrics/emd/script.py + helper.py (SHA TASK_SHA); metrics emd_mean_ct_horiz and
# emd_mean_ct_vert, lower is better.
_BIN_CENTRES = (EMD_BIN_EDGES[:-1] + EMD_BIN_EDGES[1:]) / 2


def _bin_probs(X):
    """helper.bin_array per marker: histogram on the fixed grid, normalised."""
    out = np.empty((X.shape[1], len(_BIN_CENTRES)))
    for j in range(X.shape[1]):
        counts, _ = np.histogram(X[:, j], bins=EMD_BIN_EDGES)
        out[j] = counts / counts.sum()
    return out


def _emd_pair(pa, pb):
    return np.array(
        [wasserstein_distance(_BIN_CENTRES, _BIN_CENTRES, a, b) for a, b in zip(pa, pb)]
    )


def _emd_horizontal(s1, s2):
    (X1, o1), (X2, o2) = s1, s2
    vals = []
    for donor in o1["donor"].unique():
        d1, d2 = (o1["donor"] == donor).to_numpy(), (o2["donor"] == donor).to_numpy()
        both = np.intersect1d(o1.loc[d1, "cell_type"], o2.loc[d2, "cell_type"])
        for ct in both:
            m1 = d1 & (o1["cell_type"] == ct).to_numpy()
            m2 = d2 & (o2["cell_type"] == ct).to_numpy()
            if m1.sum() < EMD_MIN_CELLS or m2.sum() < EMD_MIN_CELLS:
                continue
            vals.append(_emd_pair(_bin_probs(X1[m1]), _bin_probs(X2[m2])))
    return np.nanmean(np.concatenate(vals))  # empty -> ValueError, as in the task


def _emd_vertical_one_split(X, obs):
    """All sample pairs within a biological group, per population."""
    combos = []
    for samples in obs.groupby("group", observed=True)["sample"].unique():
        samples = list(samples)
        combos += [(a, b) for i, a in enumerate(samples) for b in samples[i + 1 :]]
    if not combos:
        return None
    cache = {}
    for s in {s for c in combos for s in c}:
        for ct in obs["cell_type"].unique():
            m = ((obs["sample"] == s) & (obs["cell_type"] == ct)).to_numpy()
            if m.sum() >= EMD_MIN_CELLS:
                cache[(s, ct)] = _bin_probs(X[m])
    vals = [
        _emd_pair(cache[(a, ct)], cache[(b, ct)])
        for a, b in combos
        for ct in obs["cell_type"].unique()
        if (a, ct) in cache and (b, ct) in cache
    ]
    return np.concatenate(vals)


def emd_scores(s1, s2):
    horiz = _emd_horizontal(s1, s2)
    vert = [_emd_vertical_one_split(*s) for s in (s1, s2)]
    vert = [v for v in vert if v is not None]
    return horiz, (np.nanmean(np.concatenate(vert)) if vert else np.nan)


# Average batch R2 -----------------------------------------------------------------
# src/metrics/average_batch_r2/script.py + helper.py; metric average_batch_r2_ct,
# lower is better. The task fits LinearRegression(marker ~ batch dummies) on the
# donor x population cells of both splits; R2 of that fit = between-batch SS / total SS.


def r2_score_batch(s1, s2):
    (X1, o1), (X2, o2) = s1, s2
    r2 = []
    for donor in o1["donor"].unique():
        d1, d2 = (o1["donor"] == donor).to_numpy(), (o2["donor"] == donor).to_numpy()
        for ct in o1.loc[d1, "cell_type"].unique():
            m1 = d1 & (o1["cell_type"] == ct).to_numpy()
            m2 = d2 & (o2["cell_type"] == ct).to_numpy()
            if m1.sum() < R2_MIN_CELLS or m2.sum() < R2_MIN_CELLS:
                continue
            Y = np.vstack([X1[m1], X2[m2]])
            b = np.concatenate([o1.loc[m1, "batch"], o2.loc[m2, "batch"]])
            sst = ((Y - Y.mean(0)) ** 2).sum(0)
            sse = np.zeros(Y.shape[1])
            for lev in np.unique(b):
                Yb = Y[b == lev]
                sse += ((Yb - Yb.mean(0)) ** 2).sum(0)
            # sklearn r2_score: constant target and exact fit gives 1.0
            r2.append(np.where(sst == 0, 1.0, 1 - sse / np.where(sst == 0, 1, sst)))
    return np.mean(np.concatenate(r2)) if r2 else np.nan


# LISI -----------------------------------------------------------------------------
# src/metrics/lisi/script.py + helper.py; metrics iLisi and cLisi, higher is better.
# LISI per cell = scib_metrics 0.5.7 lisi_knn (src/scib_metrics/metrics/_lisi.py,
# utils/_lisi.py): perplexity floor(k/3), self edge masked, beta search tol 1e-5,
# 50 tries, simpson index 1 if H == 0 returns -1. kNN is exact here (cKDTree,
# Euclidean, k=100 incl. self); the task uses pynndescent. Arithmetic is float64,
# scib_metrics uses jax float32.


def _lisi_knn(X, labels):
    labels = np.asarray(pd.Categorical(labels).codes)
    n = X.shape[0]
    k = min(LISI_N_NEIGHBORS, n)
    dist, idx = cKDTree(X).query(X, k=k, workers=-1)
    self_free = idx != np.arange(n)[:, None]
    log_perp = np.log(np.floor(k / 3))

    def hbeta(beta):
        with np.errstate(all="ignore"):
            P = np.where(self_free, np.exp(-dist * beta[:, None]), 0.0)
            sumP = np.nansum(P, 1)
            H = np.where(
                sumP == 0, 0.0, np.log(sumP) + beta * np.nansum(dist * P, 1) / sumP
            )
            P = np.where(sumP[:, None] == 0, 0.0, P / sumP[:, None])
        return H, P

    beta = np.ones(n)
    bmin, bmax = np.full(n, -np.inf), np.full(n, np.inf)
    H, P = hbeta(beta)
    hdiff = H - log_perp
    tries = np.zeros(n, dtype=int)
    while True:
        act = (np.abs(hdiff) >= LISI_TOL) & (tries < LISI_MAX_TRIES)
        if not act.any():
            break
        up = hdiff > 0
        nbeta = np.where(
            up,
            np.where(bmax == np.inf, beta * 2, (beta + bmax) / 2),
            np.where(bmin == -np.inf, beta / 2, (beta + bmin) / 2),
        )
        nbmin, nbmax = np.where(up, beta, bmin), np.where(up, bmax, beta)
        nH, nP = hbeta(nbeta)
        beta, bmin, bmax = (
            np.where(act, a, b)
            for a, b in ((nbeta, beta), (nbmin, bmin), (nbmax, bmax))
        )
        H, P = np.where(act, nH, H), np.where(act[:, None], nP, P)
        hdiff = np.where(act, nH - log_perp, hdiff)
        tries = tries + act
    knn_lab = labels[idx]
    simpson = np.zeros(n)
    for lab in range(len(np.unique(labels))):
        simpson += np.where(knn_lab == lab, P, 0.0).sum(1) ** 2
    simpson = np.where(H == 0, -1.0, simpson)
    return 1 / simpson


def _ilisi(X, obs):
    # helper._check_batch_group_confounding: NaN if any group sits in one batch
    if (obs.groupby("group", observed=True)["batch"].nunique() < 2).any():
        return np.nan
    n_batches = obs["batch"].nunique()
    return (np.nanmedian(_lisi_knn(X, obs["batch"].to_numpy())) - 1) / (n_batches - 1)


def _clisi(X, obs):
    n_ct = obs["cell_type"].nunique()
    lisi = _lisi_knn(X, obs["cell_type"].to_numpy())
    return (n_ct - np.nanmedian(lisi)) / (n_ct - 1)


def lisi_scores(s1, s2):
    ilisi = np.mean([_ilisi(*s1), _ilisi(*s2)])
    clisi = np.mean([_clisi(*s1), _clisi(*s2)])
    return ilisi, clisi


# Ratio of consistent peaks ----------------------------------------------------------
# src/metrics/ratio_consistent_peaks/script.py + helper.py; metric
# ratio_consistent_peaks, higher is better. Uses X_raw of the split's cells of the
# same donor and population as the unintegrated reference.


def _standardise(a, b):
    pooled = np.concatenate([a, b])
    mu, sd = pooled.mean(), pooled.std()
    return (a - mu) / sd, (b - mu) / sd


def _kde_density(v):
    lo, hi = v.min(), v.max()
    kde = gaussian_kde(np.reshape(v, (1, -1)), bw_method="scott")
    density = kde(np.linspace(lo, hi, RCP_GRID_POINTS))
    if np.argmax(density) == 0 and density.size > 1:  # helper.get_kde_density shift
        grid = np.linspace(lo, hi, RCP_GRID_POINTS - 1)
        density = np.concatenate([[0.0], kde(grid)])
    return density


def _n_peaks(v):
    peaks, _ = find_peaks(_kde_density(v), prominence=RCP_PROMINENCE, height=RCP_HEIGHT)
    return len(peaks)


def ratio_consistent_peaks(raw, s1, s2):
    (Xu, ou) = raw  # unintegrated, non-control, labelled, markers to correct
    (X1, o1), (X2, o2) = s1, s2
    n1 = n2 = 0
    for donor in o1["donor"].unique():
        du = (ou["donor"] == donor).to_numpy()
        d1, d2 = (o1["donor"] == donor).to_numpy(), (o2["donor"] == donor).to_numpy()
        for ct in o1.loc[d1, "cell_type"].unique():
            m1 = d1 & (o1["cell_type"] == ct).to_numpy()
            m2 = d2 & (o2["cell_type"] == ct).to_numpy()
            if m1.sum() < RCP_MIN_CELLS or m2.sum() < RCP_MIN_CELLS:
                continue
            mu = du & (ou["cell_type"] == ct).to_numpy()
            mu1, mu2 = (
                mu & (ou["split"] == 1).to_numpy(),
                mu & (ou["split"] == 2).to_numpy(),
            )
            for j in range(Xu.shape[1]):
                u1, u2 = _standardise(Xu[mu1, j], Xu[mu2, j])
                pu1, pu2 = _n_peaks(u1), _n_peaks(u2)
                c1, c2 = X1[m1, j], X2[m2, j]
                # flat after integration but not before: automatic case 2
                if (np.std(c1) == 0 and np.std(Xu[mu1, j]) != 0) or (
                    np.std(c2) == 0 and np.std(Xu[mu2, j]) != 0
                ):
                    n2 += 1
                    continue
                i1, i2 = _standardise(c1, c2)
                p1, p2 = _n_peaks(i1), _n_peaks(i2)
                if pu1 == pu2:  # unintegrated consistent: case 1 or 2; else case NGT
                    if p1 == pu1 and pu2 == p2:
                        n1 += 1
                    else:
                        n2 += 1
    return n1 / (n1 + n2) if n1 + n2 > 0 else np.nan


# Functional marker preservation ----------------------------------------------------
# src/metrics/functional_marker_preservation/script.py + helper.py; metrics
# ..._wilcoxon (higher is better) and ..._cohens_d (lower is better). Needs exactly
# two groups and the first two batches (order of appearance) of the unintegrated data.


def _mean_expr(X, obs, markers):
    df = pd.DataFrame(X, index=obs.index, columns=markers).join(
        obs[["cell_type", "sample", "group", "batch"]]
    )
    n = df.groupby(["cell_type", "sample"], observed=True)[markers[0]].transform("size")
    df = df[n >= FMP_MIN_CELLS]
    mean = df.groupby(["cell_type", "sample", "group", "batch"], observed=True)[
        markers
    ].mean()
    return mean.reset_index().melt(
        id_vars=["cell_type", "sample", "group", "batch"],
        var_name="marker",
        value_name="mean_expr",
    )


def _group_tests(mean_df, ga, gb, kind):
    rows = []
    for (marker, ct), g in mean_df.groupby(["marker", "cell_type"], observed=True):
        a = g.loc[g["group"] == ga, "mean_expr"].to_numpy()
        b = g.loc[g["group"] == gb, "mean_expr"].to_numpy()
        if len(a) < 2 or len(b) < 2:
            continue
        if kind == "pval":
            val = mannwhitneyu(a, b, alternative="two-sided").pvalue
        else:  # Cohen's d, pooled SD with Bessel correction
            sd = np.sqrt(
                ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                / (len(a) + len(b) - 2)
            )
            val = np.nan if sd == 0 else (a.mean() - b.mean()) / sd
        rows.append({"marker": marker, "cell_type": ct, kind: val})
    return pd.DataFrame(rows, columns=["marker", "cell_type", kind])


def _sig(res):
    s = res[res["pval"] <= FMP_PVAL_THRESHOLD]
    return set(zip(s["marker"], s["cell_type"]))


def _only(df, pairs):
    return df[[(m, c) in pairs for m, c in zip(df["marker"], df["cell_type"])]]


def functional_marker_preservation(raw, s1, s2, var_u):
    """raw/s1/s2: (X, obs) restricted to functional markers among to_correct."""
    (Xu, ou), (X1, o1), (X2, o2) = raw, s1, s2
    markers = list(var_u.index)
    groups = sorted(ou["group"].unique())
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 biological groups, found: {groups}")
    ga, gb = groups
    batches = ou["batch"].unique()
    mu = {
        b: _mean_expr(
            Xu[(ou["batch"] == b).to_numpy()], ou.loc[ou["batch"] == b], markers
        )
        for b in batches
    }
    m1, m2 = _mean_expr(X1, o1, markers), _mean_expr(X2, o2, markers)
    sig_b = {b: _sig(_group_tests(mu[b], ga, gb, "pval")) for b in batches}
    sig_u = sig_b[batches[0]] & sig_b[batches[1]]
    if not sig_u:
        return np.nan, np.nan
    sig_1 = _sig(_group_tests(m1, ga, gb, "pval")) & sig_u
    sig_2 = _sig(_group_tests(m2, ga, gb, "pval")) & sig_u
    wilcoxon = len(sig_1 & sig_2) / len(sig_u)
    d_b = [
        _group_tests(_only(mu[b], sig_u), ga, gb, "d").rename(columns={"d": f"d_{i}"})
        for i, b in enumerate(batches[:2])
    ]
    d_u = d_b[0].merge(d_b[1], on=["marker", "cell_type"])
    d_u["d_u"] = d_u[["d_0", "d_1"]].mean(axis=1)
    d_1 = _group_tests(_only(m1, sig_u), ga, gb, "d").rename(columns={"d": "d_s1"})
    d_2 = _group_tests(_only(m2, sig_u), ga, gb, "d").rename(columns={"d": "d_s2"})
    comp = d_u.merge(
        d_1.merge(d_2, on=["marker", "cell_type"]), on=["marker", "cell_type"]
    )
    abs_mean = pd.concat(
        [(comp["d_s1"] - comp["d_u"]).abs(), (comp["d_s2"] - comp["d_u"]).abs()], axis=1
    ).mean(axis=1)
    return wilcoxon, abs_mean.mean()


# FlowSOM mapping similarity ---------------------------------------------------------
# src/metrics/flowsom_mapping_similarity/script.R + helper.R; metric
# flowsom_mean_mapping_similarity, higher is better (0 to 100). The R functions are
# copied below and run with Rscript, per donor and direction, on lineage markers.
_FLOWSOM_R = r"""
suppressPackageStartupMessages({library(flowCore); library(FlowSOM)})
a <- commandArgs(TRUE)
s1 <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE)
s2 <- read.csv(a[2], check.names = FALSE, stringsAsFactors = FALSE)
lineage_markers <- strsplit(a[3], "\t", fixed = TRUE)[[1]]
n_clusters <- as.integer(a[4]); xdim <- as.integer(a[5]); ydim <- as.integer(a[6])
seed <- as.integer(a[7])
marker_cols <- setdiff(names(s1), c("donor", "cell_type"))
get_ff_annotations <- function(df, donor_name) {
  m <- as.matrix(df[df$donor == donor_name, marker_cols])
  list(flowframe = flowFrame(m), ct_annotations = df$cell_type[df$donor == donor_name])
}
compute_clust_pct <- function(vector_clust, vector_celltype) {
  tab <- table(vector_clust, vector_celltype)
  pct_tab <- prop.table(tab, margin = 1) * 100
  pct_tab <- round(pct_tab, 2)
  pct_tab[is.na(pct_tab)] <- 0
  pct_tab
}
compute_clust_weights <- function(vector_clust_s1, vector_clust_s2) {
  proportions(table(c(vector_clust_s1, vector_clust_s2), exclude = NULL))
}
compute_fs_mapping_similarity <- function(fs_tree_s1, ct_ann_s1, fs_tree_s2, ct_ann_s2) {
  clust_levels <- levels(as.factor(c(GetClusters(fs_tree_s1), GetClusters(fs_tree_s2))))
  clust_s1 <- factor(paste0("clust ", GetClusters(fs_tree_s1)),
                     levels = paste0("clust ", clust_levels))
  clust_s2 <- factor(paste0("clust ", GetClusters(fs_tree_s2)),
                     levels = paste0("clust ", clust_levels))
  clust_weights <- compute_clust_weights(clust_s1, clust_s2)
  # Deviation from the task code: shared cell-type levels, so a type missing from one split gives a
  # 0 column instead of "non-conformable arrays". Identical result whenever both splits have every type.
  ct_levels <- sort(unique(c(ct_ann_s1, ct_ann_s2)))
  ct_ann_s1 <- factor(ct_ann_s1, levels = ct_levels)
  ct_ann_s2 <- factor(ct_ann_s2, levels = ct_levels)
  t1 <- compute_clust_pct(clust_s1, ct_ann_s1)
  t2 <- compute_clust_pct(clust_s2, ct_ann_s2)
  100 - sum(clust_weights * rowSums(abs(t2 - t1)))
}
res <- list()
for (donor in unique(s1$donor)) {
  d1 <- get_ff_annotations(s1, donor); d2 <- get_ff_annotations(s2, donor)
  tree1 <- FlowSOM(d1$flowframe, colsToUse = lineage_markers, nClus = n_clusters,
                   xdim = xdim, ydim = ydim, seed = seed)
  map2 <- NewData(tree1, d2$flowframe)
  tree2 <- FlowSOM(d2$flowframe, colsToUse = lineage_markers, nClus = n_clusters,
                   xdim = xdim, ydim = ydim, seed = seed)
  map1 <- NewData(tree2, d1$flowframe)
  res[[donor]] <- c(
    compute_fs_mapping_similarity(tree1, d1$ct_annotations, map2, d2$ct_annotations),
    compute_fs_mapping_similarity(tree2, d2$ct_annotations, map1, d1$ct_annotations))
}
write.csv(data.frame(similarity = unlist(res)), a[8], row.names = FALSE)
"""


def flowsom_mapping_similarity(
    s1, s2, markers, lineage, xdim, ydim, n_clusters, rscript
):
    rscript = rscript or os.environ.get("SW_RSCRIPT") or shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError(
            "Rscript with FlowSOM and flowCore is needed (set SW_RSCRIPT)"
        )
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for k, (X, o) in zip((1, 2), (s1, s2)):
            df = pd.DataFrame(X, columns=markers)
            df.insert(0, "cell_type", o["cell_type"].to_numpy())
            df.insert(0, "donor", o["donor"].to_numpy())
            df.to_csv(tmp / f"s{k}.csv", index=False)
        (tmp / "run.R").write_text(_FLOWSOM_R)
        run = subprocess.run(
            [
                rscript,
                str(tmp / "run.R"),
                str(tmp / "s1.csv"),
                str(tmp / "s2.csv"),
                "\t".join(lineage),
                str(n_clusters),
                str(xdim),
                str(ydim),
                str(FLOWSOM_SEED),
                str(tmp / "out.csv"),
            ],
            capture_output=True,
            text=True,
        )
        if run.returncode != 0:
            raise RuntimeError(f"FlowSOM mapping step failed in R:\n{run.stderr[-3000:]}")
        return pd.read_csv(tmp / "out.csv")["similarity"].mean()


# Public API ----------------------------------------------------------------------


def score_all(
    X_raw,
    X_corrected,
    obs,
    var,
    *,
    flowsom_xdim,
    flowsom_ydim,
    flowsom_n_clusters,
    obs_corrected=None,
    rscript=None,
    metrics=None,
):
    """Return the 9 task scores as a dict keyed by the task's metric ids.

    obs_corrected: only for the perfect_integration control, whose output has its own
        rows (see control_perfect_integration). Default: same rows as X_raw / obs.
    flowsom_*: the dataset's uns parameter_som_xdim, parameter_som_ydim and
        parameter_num_clusters. No default: they are the dataset's choice.
    metrics: optional list of metric ids to compute (e.g. to skip FlowSOM without R).
    """
    metrics = list(METRIC_MAXIMIZE) if metrics is None else list(metrics)
    obs_c = obs if obs_corrected is None else obs_corrected
    if obs_corrected is None and len(obs) != len(X_corrected):
        raise ValueError("X_corrected and obs rows differ; pass obs_corrected")
    splits, markers = _by_split(X_corrected, obs_c, var)
    s1, s2 = splits[1], splits[2]
    Xu, ou, _ = _prep(X_raw, obs, var)  # unintegrated reference
    out = {}
    if {"emd_mean_ct_horiz", "emd_mean_ct_vert"} & set(metrics):
        out["emd_mean_ct_horiz"], out["emd_mean_ct_vert"] = emd_scores(s1, s2)
    if "average_batch_r2_ct" in metrics:
        out["average_batch_r2_ct"] = r2_score_batch(s1, s2)
    if {"iLisi", "cLisi"} & set(metrics):
        out["iLisi"], out["cLisi"] = lisi_scores(s1, s2)
    if "ratio_consistent_peaks" in metrics:
        out["ratio_consistent_peaks"] = ratio_consistent_peaks((Xu, ou), s1, s2)
    if {
        "functional_marker_preservation_wilcoxon",
        "functional_marker_preservation_cohens_d",
    } & set(metrics):
        func = (var.loc[var["to_correct"], "marker_type"] == "functional").to_numpy()
        var_f = var.loc[var["to_correct"]].loc[func]
        w, d = functional_marker_preservation(
            (Xu[:, func], ou), (s1[0][:, func], s1[1]), (s2[0][:, func], s2[1]), var_f
        )
        out["functional_marker_preservation_wilcoxon"] = w
        out["functional_marker_preservation_cohens_d"] = d
    if "flowsom_mean_mapping_similarity" in metrics:
        lin = (var.loc[var["to_correct"], "marker_type"] == "lineage").to_numpy()
        out["flowsom_mean_mapping_similarity"] = flowsom_mapping_similarity(
            s1,
            s2,
            markers,
            [m for m, f in zip(markers, lin) if f],
            flowsom_xdim,
            flowsom_ydim,
            flowsom_n_clusters,
            rscript,
        )
    return {m: out[m] for m in metrics}


# Scaling -------------------------------------------------------------------------
# Not in the task repo. The task's results are scaled in openproblems-bio/openproblems
# commit 673fa8f63ac620fd296b506ae2a741b027f979f2, src/reporting/render_report/
# report-template.qmd (control_ranges, scaled_scores), and in the website's
# openproblems-bio/openproblems-v2 src/common/process_task_results/get_results/
# script.R (scale_scores, aggregate_scores, commit f0ef558f16a94526f16ce888f246d3a3d3986e9d
# is the last commit touching that file). Rule: per metric, min and max over the
# control methods' values; scaled = (v - min) / (max - min); 1 - scaled when the
# metric is lower-is-better. The website then sets missing/NaN to 0, needs >= 2
# control values and min != max, and the method mean clips each score to [0, 1].
# The task's control methods are no_integration, perfect_integration and
# shuffle_integration_globally; pass the ones you want as `controls`.


def scale_scores(scores, controls, clip=True):
    """scores: dict metric -> value. controls: dict control name -> dict metric -> value.

    Returns (scaled dict, mean score). clip=True is the website's aggregation.
    """
    scaled = {}
    for metric, v in scores.items():
        cv = np.array([c[metric] for c in controls.values() if not pd.isna(c[metric])])
        if len(cv) < 2 or cv.min() == cv.max() or pd.isna(v):
            scaled[metric] = 0.0
            continue
        s = (v - cv.min()) / (cv.max() - cv.min())
        scaled[metric] = s if METRIC_MAXIMIZE[metric] else 1 - s
    vals = np.array(list(scaled.values()))
    return scaled, float(np.mean(np.clip(vals, 0, 1) if clip else vals))


# Control methods ---------------------------------------------------------------------
# src/control_methods/{no_integration,perfect_integration,shuffle_integration_globally}
# /script.py. Each returns (X_corrected, obs_corrected) for score_all.


def control_no_integration(X_raw, obs):
    """Integrated = preprocessed, unchanged."""
    return np.asarray(X_raw), obs


def control_perfect_integration(X_raw, obs, perfect_batch):
    """Both splits return the cells of one batch (the task hard-codes batch "1").

    Batch and split labels of those cells are rewritten per split to the donor's
    nominal batch (src/utils/helper_functions.py get_obs_var_for_integrated).
    """
    X_raw = np.asarray(X_raw)
    obs = obs.copy()
    obs["batch"] = obs["batch"].astype(str)
    take = (
        (obs["is_control"].astype(int) == 0) & (obs["batch"] == str(perfect_batch))
    ).to_numpy()
    parts, Xs = [], []
    for k in (1, 2):
        nominal = obs[(obs["is_control"].astype(int) == 0) & (obs["split"] == k)]
        mapping = nominal.groupby("donor", observed=True)["batch"].unique()
        if (mapping.map(len) != 1).any():
            raise ValueError(f"a donor has more than one batch in split {k}")
        mapping = mapping.map(lambda b: b[0])
        o = obs.loc[take].copy()
        o["batch"] = o["donor"].map(mapping)
        o["split"] = k
        o.index = [f"{i}|split{k}" for i in o.index]
        parts.append(o)
        Xs.append(X_raw[take])
    return np.vstack(Xs), pd.concat(parts)


def control_shuffle_integration_globally(X_raw, obs, seed):
    """Shuffle all cells of each split (control cells are in the pool, as in the task)."""
    X_raw = np.asarray(X_raw)
    rng = np.random.default_rng(seed)
    ic = obs["is_control"].astype(int).to_numpy()
    X_out = X_raw.copy()
    for k in (1, 2):
        pool = np.flatnonzero((ic > 0) | (obs["split"].to_numpy() == k))
        X_out[pool] = X_raw[rng.permutation(pool)]
    return X_out, obs


# Self-check ------------------------------------------------------------------------
# Synthetic data: 3 populations, 6 donors (3 per group), 2 batches, each donor has one
# sample per split in different batches. Batch 2 adds +0.8 to every marker. Group B
# lowers functional marker F1 in populations P1 and P2.


def _synthetic(seed=1, n_per=150):
    rng = np.random.default_rng(seed)
    markers = ["L1", "L2", "L3", "L4", "F1", "F2", "Time"]
    var = pd.DataFrame(
        {
            "to_correct": [True] * 6 + [False],
            "marker_type": ["lineage"] * 4 + ["functional"] * 2 + ["lineage"],
        },
        index=markers,
    )
    centre = {"P1": [3, 0, 0, 0], "P2": [0, 3, 0, 0], "P3": [0, 0, 3, 3]}
    split1_batch = {"d1": "1", "d2": "1", "d3": "2", "d4": "1", "d5": "2", "d6": "2"}
    group = {d: "A" if d in ("d1", "d2", "d3") else "B" for d in split1_batch}
    samples = [
        (d, k, split1_batch[d] if k == 1 else str(3 - int(split1_batch[d])))
        for d in split1_batch
        for k in (1, 2)
    ]
    samples += [("ctrl", 0, "1"), ("ctrl", 0, "2")]
    X, rows = [], []
    for donor, split, batch in samples:
        s_eff = rng.normal(0, 0.1)
        for pop, c in centre.items():
            m = np.concatenate([c, [1.0, 0.5], [0.0]])
            if group.get(donor) == "B" and pop in ("P1", "P2"):
                m[4] -= 2.0
            x = m + s_eff + rng.normal(0, 0.5, (n_per, 7))
            if batch == "2":
                x[:, :6] += 0.8
            X.append(x)
            rows += [
                dict(
                    batch=batch,
                    donor=donor,
                    sample_id=f"{donor}_s{split}_b{batch}",
                    population=pop,
                    is_control=int(donor == "ctrl"),
                    split=split,
                    group=group.get(donor, "A"),
                )
            ] * n_per
    obs = pd.DataFrame(rows)
    obs.index = [f"c{i}" for i in range(len(obs))]
    return np.vstack(X), obs, var


def _selfcheck():
    X, obs, var = _synthetic()
    kw = dict(flowsom_xdim=3, flowsom_ydim=3, flowsom_n_clusters=3)
    run_flowsom = (os.environ.get("SW_RSCRIPT") or shutil.which("Rscript")) is not None
    metrics = [m for m in METRIC_MAXIMIZE if run_flowsom or "flowsom" not in m]
    if not run_flowsom:
        print("Rscript not found: FlowSOM metric skipped (set SW_RSCRIPT)")
    Xn, on = control_no_integration(X, obs)
    Xp, op = control_perfect_integration(X, obs, perfect_batch="1")
    Xs, os_ = control_shuffle_integration_globally(X, obs, seed=7)
    res = {
        "no_integration": score_all(X, Xn, obs, var, metrics=metrics, **kw),
        "perfect_integration": score_all(
            X, Xp, obs, var, obs_corrected=op, metrics=metrics, **kw
        ),
        "shuffle_integration_globally": score_all(
            X, Xs, obs, var, metrics=metrics, **kw
        ),
    }
    print(pd.DataFrame(res).round(4).to_string())
    perf, none = res["perfect_integration"], res["no_integration"]
    assert perf["emd_mean_ct_horiz"] < 1e-12, "perfect: horizontal EMD must be 0"
    assert perf["average_batch_r2_ct"] < 1e-9, "perfect: batch R2 must be 0"
    assert perf["ratio_consistent_peaks"] == 1.0, "perfect: peak ratio must be 1"
    assert perf["functional_marker_preservation_wilcoxon"] == 1.0
    if run_flowsom:
        assert abs(perf["flowsom_mean_mapping_similarity"] - 100) < 1e-9
    assert none["emd_mean_ct_horiz"] > 0.5 and none["average_batch_r2_ct"] > 0.1
    scaled, mean = scale_scores(
        none, {k: res[k] for k in ("no_integration", "perfect_integration")}
    )
    assert scaled["emd_mean_ct_horiz"] == 0.0 and scaled["average_batch_r2_ct"] == 0.0
    scaled_p, mean_p = scale_scores(
        perf, {k: res[k] for k in ("no_integration", "perfect_integration")}
    )
    assert (
        scaled_p["emd_mean_ct_horiz"] == 1.0 and scaled_p["average_batch_r2_ct"] == 1.0
    )
    print(
        "scaled (controls: no, perfect), no_integration:",
        {k: round(v, 3) for k, v in scaled.items()},
        round(mean, 3),
    )
    print(
        "scaled (controls: no, perfect), perfect_integration:",
        {k: round(v, 3) for k, v in scaled_p.items()},
        round(mean_p, 3),
    )
    print("self-check passed")


if __name__ == "__main__":
    _selfcheck()
