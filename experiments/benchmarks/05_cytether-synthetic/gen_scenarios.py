#!/usr/bin/env python3
"""Purpose: write the development scenarios S13, S13b, S14, S15, S16, S17 (shrinkage, thin reference design) for dev seeds.

Run:     <python> experiments/benchmarks/05_cytether-synthetic/gen_scenarios.py <out_dir> <seed> [<seed> ...]
Inputs:  cyTether bundle v1 (MIT) at $SW_DATA_DIR/cytether_bundle_v1, read only, never modified
Outputs (in <out_dir>/seed_<seed>/<stage>/), same layout as gen_data.py:
  observed.npz, cell_metadata.csv.gz, truth_EVALUATION_ONLY.npz, manifest.json
  (run export_observed.py afterwards to get observed_X.csv.gz)

The bundle generator draws the biology and cells. A post-hoc layer, in the style of simulate_extended.py,
then (1) swaps each batch's technical offset for a new draw and (2) drops reference rows to make the design thin.
Shift-only stages have x = x_true + offset, so a batch's old offset is exactly mean(x - x_true) over its cells.

S17 is S13 with population-specific shifts: every cell type gets its own independent study and run offset,
so an effect shared across populations is wrong. Truth offsets then have a cell-type axis (study, type, marker).
Scenario indices 13 to 17 are registered in LEVELS here; S13b reuses index 13, so S13 and S13b share truth and offsets.

S18 (plan A4, A5, A10) has no additive offsets. Each batch x channel gets an affine instrument effect on the linear
scale, x = asinh(exp(g) * sinh(truth) + o), g ~ N(0, 0.20), o ~ N(0, 0.10), zero in the anchor batch, so the
asinh-scale effect depends on intensity. 3 batches per seed (seeded) keep no training reference donor; REF_2
(held out) stays everywhere. S18w is S18 plus one participant sample per batch whose T cells get a CD3 truth
shift N(0, 0.3) (biology, not batch): the internal CD3 anchor wobbles there. Run with S18_ONLY=1.
Seeds 441 to 460 are the development seeds (451 to 460 added by F5, 2026-11-10). Test seeds 431 to 440 are never written or read here.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True  # keep the read-only bundle free of new files

# Parameters (user-approved design, 2026-11-09) -----------------------------------
CELLS_PER_PARTICIPANT = 70  # approved E2
CELLS_PER_REFERENCE_ALIQUOT = 90  # approved E2
DEV_SEEDS = range(441, 461)  # F5 adds 451-460; test seeds 431 to 440 are never touched
SD_RUN, SD_STUDY = 0.10, 0.47  # S13, S13b, S15, S16, S17 (S16 overrides two runs)
SD_FLAT_RUN = 0.47  # S14: all offset variance sits in the runs
SD_DISCREPANT_RUN = 0.60  # S16, the two discrepant runs
DISCREPANT_RUNS = ["ST2_RUN1", "ST4_RUN2"]
SD_AFFINE_OFFSET, SD_AFFINE_LOGGAIN = 0.10, 0.20  # S18: A10
N_REFLESS_S18 = 3  # S18: A5
SD_ANCHOR_WOBBLE = 0.30  # S18w: A4
WOBBLE_TYPES = ["CD4_T", "CD8_T"]  # T cells carry the CD3+ anchor
TRAINING_DONORS = ["REF_0", "REF_1"]  # REF_2 is the held-out donor and is never dropped
ANCHOR = "ST0_RUN0"  # keeps both training donors in the thin design
REFLESS_KNOWN_RUN = "ST4_RUN3"  # S15: known study, new run, no non-holdout reference rows

bundle = Path(os.environ.get("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")) / "cytether_bundle_v1"
sys.path.insert(0, str(bundle / "src"))
from cytether_bench.simulate import LEVELS, MARKERS, simulate  # noqa: E402

# Scenario registry -----------------------------------------------------------------
# Same keys as the bundle's S4 (shift-only, 3 reference donors, REF_2 held out); only studies and runs differ.
LEVELS["S13"] = dict(title="shrink_thin_nested", studies=6, runs=3, shift=1, scale=0, refs=3)
LEVELS["S14"] = dict(title="flat_large_runs", studies=6, runs=3, shift=1, scale=0, refs=3)
LEVELS["S15"] = dict(title="new_run_vs_new_study", studies=6, runs=4, shift=1, scale=0, refs=3)
LEVELS["S16"] = dict(title="discrepant_runs", studies=6, runs=3, shift=1, scale=0, refs=3)
LEVELS["S17"] = dict(title="population_specific_shift", studies=6, runs=3, shift=1, scale=0, refs=3)
LEVELS["S18"] = dict(title="affine_instrument_anchor_only_batches", studies=6, runs=3, shift=0, scale=0, refs=3)
SCENARIO_INDEX = {"S13": 13, "S14": 14, "S15": 15, "S16": 16, "S17": 17, "S18": 18}


# Post-hoc layer --------------------------------------------------------------------
def redraw_offsets(x, true, meta, n_studies, n_runs, sd_study, sd_run, index, seed, by_population=False):
    """Replace each batch's offset: x = x - old_offset + new_offset. Study 0 and run (0, 0) stay at zero.

    by_population=True (S17) draws an independent offset per cell type, so offsets have shape (..., type, marker).
    """
    rng = np.random.default_rng(seed + 1000 + index)
    n_markers = x.shape[1]
    types = sorted(meta.cell_type.unique())
    pop_shape = (len(types),) if by_population else ()
    study_offset = rng.normal(0, 1, (n_studies, *pop_shape, n_markers)) * sd_study
    sd_run = np.asarray(sd_run).reshape(n_studies, n_runs, *([1] * len(pop_shape)), 1)
    run_offset = rng.normal(0, 1, (n_studies, n_runs, *pop_shape, n_markers)) * sd_run
    study_offset[0] = 0
    run_offset[0, 0] = 0
    x = x.astype("float64")
    for s in range(n_studies):
        for r in range(n_runs):
            rows = (meta.batch == f"ST{s}_RUN{r}").to_numpy()
            old_offset = (x[rows] - true[rows]).mean(axis=0)
            new_offset = study_offset[s] + run_offset[s, r]  # (marker,) or (type, marker)
            if by_population:
                type_index = pd.Categorical(meta.cell_type[rows], categories=types).codes
                new_offset = new_offset[type_index]
            x[rows] = x[rows] - old_offset + new_offset
    return x.astype("float32"), study_offset, run_offset


def thin_reference_rows(meta, refless=()):
    """Boolean keep mask: anchor keeps both training donors, other batches alternate one, refless batches none."""
    ref_batches = list(dict.fromkeys(meta.batch[meta.is_reference]))  # generation order
    keep_donors = {ANCHOR: set(TRAINING_DONORS)}
    others = [b for b in ref_batches if b != ANCHOR and b not in refless]
    for i, b in enumerate(others):
        keep_donors[b] = {TRAINING_DONORS[i % 2]}
    for b in refless:
        keep_donors[b] = set()
    is_training_ref = (meta.is_reference & meta.donor_id.isin(TRAINING_DONORS)).to_numpy()
    wanted = np.array([d in keep_donors[b] for b, d in zip(meta.batch, meta.donor_id)]) if is_training_ref.any() else None
    keep = np.ones(len(meta), dtype=bool)
    keep[is_training_ref] = wanted[is_training_ref]
    return keep


def affine_s18(true, meta, seed, wobble):
    """S18 / S18w: affine per batch x channel on the linear scale; S18w adds the CD3 wobble to truth first."""
    rng = np.random.default_rng(seed + 1000 + SCENARIO_INDEX["S18"])
    batches = sorted(meta.batch.unique())
    n_markers = true.shape[1]
    gain = rng.normal(0, SD_AFFINE_LOGGAIN, (len(batches), n_markers))
    offset = rng.normal(0, SD_AFFINE_OFFSET, (len(batches), n_markers))
    gain[batches.index(ANCHOR)] = 0
    offset[batches.index(ANCHOR)] = 0
    true = true.astype("float64").copy()
    if wobble:
        wrng = np.random.default_rng(seed + 3000 + SCENARIO_INDEX["S18"])
        cd3 = MARKERS.index("CD3")
        for b in batches:
            samples = sorted(meta.sample_id[(meta.batch == b) & ~meta.is_reference].unique())
            pick = samples[wrng.integers(len(samples))]
            rows = ((meta.sample_id == pick) & meta.cell_type.isin(WOBBLE_TYPES)).to_numpy()
            true[rows, cd3] += wrng.normal(0, SD_ANCHOR_WOBBLE)
    bi = pd.Categorical(meta.batch, categories=batches).codes
    x = np.arcsinh(np.exp(gain[bi]) * np.sinh(true) + offset[bi])
    return x.astype("float32"), true.astype("float32"), gain, offset


def refless_s18(meta, seed):
    """A5: N_REFLESS_S18 non-anchor batches without training reference donors, seeded."""
    rng = np.random.default_rng(seed + 2000 + SCENARIO_INDEX["S18"])
    others = sorted(b for b in meta.batch.unique() if b != ANCHOR)
    return tuple(sorted(rng.choice(others, N_REFLESS_S18, replace=False)))


def build(base_stage, seed):
    """Return x, true, meta, cfg, offsets for the full (not thinned) data of one base scenario."""
    cfg = LEVELS[base_stage]
    x, true, meta, cfg = simulate(base_stage, seed, CELLS_PER_PARTICIPANT, CELLS_PER_REFERENCE_ALIQUOT)
    ns, nr = cfg["studies"], cfg["runs"]
    sd_study = {"S14": 0.0}.get(base_stage, SD_STUDY)
    sd_run = np.full((ns, nr), {"S14": SD_FLAT_RUN}.get(base_stage, SD_RUN))
    if base_stage == "S16":
        for b in DISCREPANT_RUNS:
            sd_run[int(b[2]), int(b[-1])] = SD_DISCREPANT_RUN
    x, study_off, run_off = redraw_offsets(x, true, meta, ns, nr, sd_study, sd_run, SCENARIO_INDEX[base_stage], seed,
                                         by_population=base_stage == "S17")
    meta = meta.copy()
    extra = {}
    if base_stage == "S15":
        new_study = meta.study == "ST5"
        new_run = meta.batch.str.endswith("RUN3") & ~new_study
        meta["arrival"] = np.where(new_study | new_run, "future", "initial")
        meta["arrival_kind"] = np.select([new_study, new_run], ["new_study", "known_study_new_run"], "initial")
        extra["refless_known_study_runs"] = [REFLESS_KNOWN_RUN]
    if base_stage == "S16":
        extra["discrepant_runs"] = DISCREPANT_RUNS  # manifest only; no method reads it
    return x, true, meta, cfg, dict(study=study_off, run=run_off), extra


def write_stage(seed_dir, stage, x, true, meta, cfg, offsets, extra, seed, design, base_stage):
    path = seed_dir / stage
    path.mkdir(parents=True)
    np.savez_compressed(path / "observed.npz", X=x, markers=np.asarray(MARKERS))
    np.savez_compressed(path / "truth_EVALUATION_ONLY.npz", X_true=true, markers=np.asarray(MARKERS),
                        study_offset=offsets["study"], run_offset=offsets["run"])
    meta.to_csv(path / "cell_metadata.csv.gz", index=False, compression="gzip")
    settings = dict(cfg, title=stage_titles[stage], base=base_stage, scenario_index=SCENARIO_INDEX[base_stage],
                    reference_design=design, training_donors=TRAINING_DONORS, heldout_donor="REF_2",
                    posthoc_offset_seed=seed + 1000 + SCENARIO_INDEX[base_stage])
    manifest = dict(stage=stage, description=stage_titles[stage], settings=settings, seed=seed, n_cells=len(meta),
                    n_donors=int(meta.donor_id.nunique()), n_markers=len(MARKERS), markers=MARKERS,
                    n_studies=int(meta.study.nunique()), n_batches=int(meta.batch.nunique()),
                    ground_truth_separate=True, units="simulated transformed marker expression",
                    warning="DO NOT use truth_EVALUATION_ONLY.npz for correction or selection", **extra)
    (path / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


stage_titles = {"S13": "shrink_thin_nested", "S13b": "shrink_balanced_twin", "S14": "flat_large_runs",
                "S15": "new_run_vs_new_study", "S16": "discrepant_runs", "S17": "population_specific_shift",
                "S18": "affine_instrument_anchor_only_batches", "S18w": "affine_instrument_cd3_wobble"}

BASE_STAGES = (["S17"] if os.environ.get("S17_ONLY") else ["S18"] if os.environ.get("S18_ONLY")
               else ["S13", "S14", "S15", "S16", "S17"])

# Generate --------------------------------------------------------------------------
out_dir = Path(sys.argv[1])
seeds = [int(s) for s in sys.argv[2:]]
for seed in seeds:
    if seed not in DEV_SEEDS:
        sys.exit(f"seed {seed} is not a development seed (441 to 460)")
    seed_dir = out_dir / f"seed_{seed}"
    if seed_dir.exists():
        sys.exit(f"{seed_dir} exists; run dirs are never overwritten")
    for base in BASE_STAGES:
        if base == "S18":
            _, true0, meta, cfg = simulate("S18", seed, CELLS_PER_PARTICIPANT, CELLS_PER_REFERENCE_ALIQUOT)
            refless = refless_s18(meta, seed)
            keep = thin_reference_rows(meta, refless)
            for stage, wobble in (("S18", False), ("S18w", True)):
                x, true, gain, offset = affine_s18(true0, meta, seed, wobble)
                extra = dict(refless_batches=list(refless), anchor_wobble_sd=SD_ANCHOR_WOBBLE if wobble else 0.0)
                write_stage(seed_dir, stage, x[keep], true[keep], meta[keep].reset_index(drop=True), cfg,
                            dict(study=offset, run=gain), extra, seed, "thin", base)
            continue
        x, true, meta, cfg, offsets, extra = build(base, seed)
        refless = (REFLESS_KNOWN_RUN,) if base == "S15" else ()
        keep = thin_reference_rows(meta, refless)
        thin = lambda a: a[keep]
        write_stage(seed_dir, base, thin(x), thin(true), meta[keep].reset_index(drop=True), cfg, offsets, extra,
                    seed, "thin", base)
        if base == "S13":  # twin: same truth and offsets, both training donors in every batch
            write_stage(seed_dir, "S13b", x, true, meta, cfg, offsets, extra, seed, "balanced", base)
