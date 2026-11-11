#!/usr/bin/env python3
"""Purpose: write the bundle's simulated stages S0 to S12 for a list of seeds.

Run:     pixi run python experiments/benchmarks/05_cytether-synthetic/gen_data.py <out_dir> <seed> [<seed> ...]
Inputs:  cyTether bundle v1 (MIT) at $SW_DATA_DIR/cytether_bundle_v1, sha256 checked
Outputs (in <out_dir>/seed_<seed>/<stage>/):
  observed.npz                 transformed expression X and marker names
  cell_metadata.csv.gz         sample, donor, study, batch, reference and arrival flags
  truth_EVALUATION_ONLY.npz    truth; only evaluate.R reads it
  manifest.json                stage settings

The bundle's own generator is used unchanged, so seed 431 reproduces the bundle's files.
"""
import os
import sys
from pathlib import Path

# Parameters ------------------------------------------------------------------
CELLS_PER_PARTICIPANT = 70  # approved E2: bundle configs/benchmark.yaml
CELLS_PER_REFERENCE_ALIQUOT = 90  # approved E2

bundle = Path(os.environ.get("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")) / "cytether_bundle_v1"
sys.path.insert(0, str(bundle / "src"))
from cytether_bench.simulate import LEVELS, save  # noqa: E402
from cytether_bench.simulate_extended import EXTRAS, save_extra  # noqa: E402

# Generate --------------------------------------------------------------------
out_dir = Path(sys.argv[1])
seeds = [int(s) for s in sys.argv[2:]]
for seed in seeds:
    seed_dir = out_dir / f"seed_{seed}"
    if seed_dir.exists():
        sys.exit(f"{seed_dir} exists; run dirs are never overwritten")
    for stage in list(LEVELS) + list(EXTRAS):
        writer = save if stage in LEVELS else save_extra
        writer(stage, seed_dir, seed, CELLS_PER_PARTICIPANT, CELLS_PER_REFERENCE_ALIQUOT)
