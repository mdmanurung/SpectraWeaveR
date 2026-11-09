# 05_cytether-synthetic: synthetic benchmark for reference-anchored correctors

Plan: `plans/cytether_benchmark_plan.md` and `plans/cytether_benchmark_strategy.md`. Frozen benchmark card: `$SW_DATA_DIR/tusoskill_runs/cytether_heb_20261109/benchmark/benchmark_card.md`.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `gen_data.py` | Bundle simulator for S0 to S12 per seed (70 cells per participant, 90 per reference aliquot). Seed 431 reproduces the bundle files exactly. |
| 2 | `export_observed.py` | Writes `observed_X.csv.gz` so R reads the data without the truth file. |
| 3 | `run_method.R` | One method over all stages of one seed. Fit uses initial, non-holdout reference cells only. Held-out reference cells never calibrate. |
| 4 | `evaluate.py` | The only script that reads `truth_EVALUATION_ONLY.npz`. |
| 5 | `summarize.py` | Primary table with seed-level paired bootstrap, guardrail table, figure. |
| all | `run_all.sh` | Steps 3 to 5 for a seed set. |

Seeds: development 441 to 450 (`cytether_sim_dev/`), test 431 to 440 (`cytether_sim/`, opened once at the final audit).

## Methods

| Method | Role |
|---|---|
| `none` | Negative control |
| `naive_batch_center` | Unsafe stress control (bundle baseline): uses all cells |
| `paired_shift` | Baseline to beat (C2) |
| `heb` | HEB engine under development (C6) |
| `oracle` | Sanity: truth itself; must score 0 |
| `reversed_shift` | Sanity: `paired_shift` applied with the wrong sign; must score worse than `none` |

## Settings that are the user's (approved)

`min_cells = 10`, oracle `cell_type` labels, anchor rule (most reference donors, ties to the first sorted batch), 10 seeds, cells per unit 70 and 90. Nominal values set in code: `TRANSFORM` (label only), `MIN_BATCHES = 2`, `VARIANCE_FLOOR = 1e-8` (heb).

## Metrics

Primary: RMSE of corrected held-out reference cells to the true expression, mean over S4, S9, S10, S11. Lower is better. Guardrails: abstention on S3, S5 and S12; protected-contrast ratio to `paired_shift`; frozen invariance (to add); interval coverage (HEB, to add).

## Baselines on the development seeds (441 to 450, approved settings)

Primary metric: held-out reference RMSE to truth, mean over S4, S9, S10, S11, 10 seeds. Source: `experiments/runs/20261109_synthetic_dev_v0/scored_all/summary/` (figure inspected: `primary.png`).

| Method | Mean | SD over seeds | Diff vs `paired_shift` (95% seed bootstrap) | historical_max_delta |
|---|---|---|---|---|
| oracle | 0.000 | 0.000 | -0.253 | n/a |
| naive_batch_center | 0.148 | 0.015 | -0.105 [-0.126, -0.087] | n/a (uses all cells; protected-contrast ratio 1.70) |
| heb (v0) | 0.241 | 0.032 | -0.012 [-0.018, -0.006] | 0 |
| paired_shift | 0.253 | 0.035 | 0 | 0 |
| cyCombine 0.3.0 | 0.424 | 0.020 | +0.172 [+0.152, +0.191] | 5.32 (refits) |
| none | 0.734 | 0.034 | +0.481 | n/a |
| CytoNorm 2.0.12 | 0.820 | 0.289 | +0.568 [+0.464, +0.746] | 0 |
| reversed_shift (sanity) | 1.411 | 0.059 | +1.158 | 0 |

Reading the table:
- Sanity passes: the oracle scores 0 and `reversed_shift` is worse than `none`.
- `naive_batch_center` wins this metric but has the worst protected-contrast ratio among the corrective methods that use references only on all-cell grounds (1.70 vs 1.00 for `paired_shift`). It is the unsafe stress control, as in the bundle.
- HEB-v0 beats `paired_shift` by 0.012 at `min_cells = 10`. In a labelled diagnostic at `min_cells = 3` (not the approved setting) both score about 0.185 and the difference is -0.0007 [-0.0022, +0.0011]. So HEB-v0's gain comes from using summaries `paired_shift` drops, not from better shrinkage.
- CytoNorm with the approved package defaults (15x15 SOM, 10 clusters) computes no quantiles where a batch has under 50 reference cells in a cluster. At 70 and 90 cells per unit it is data-limited here, and one seed scores 1.63. This is not a verdict on CytoNorm; the Z247 track is its fair test.
- cyCombine ran in reference mode first; where ComBat failed because a cluster lacked anchor cells it ran with no reference batch (recorded in each `status.csv`).

## Reproduce

```bash
# 1. data (once): python gen_data.py <sim_dir> <seeds...>; python export_observed.py <sim_dir>
# 2. local methods
bash experiments/benchmarks/05_cytether-synthetic/run_all.sh <sim_dir> experiments/runs/<date>_synthetic_dev paired_shift \
  "none naive_batch_center paired_shift heb reversed_shift oracle" 441 442 443 444 445 446 447 448 449 450
# 3. competitors on SLURM (5 concurrent jobs at most)
sbatch --mem=4G --cpus-per-task=2 --array=0-9%5 experiments/slurm/synthetic_method.sbatch <sim_dir> <run_dir>/results cytonorm 441
sbatch --mem=4G --cpus-per-task=2 --array=0-9%5 experiments/slurm/synthetic_method.sbatch <sim_dir> <run_dir>/results cycombine 441
# 4. score everything
bash experiments/benchmarks/05_cytether-synthetic/score.sh <sim_dir> <run_dir> paired_shift scored_all \
  "none naive_batch_center paired_shift heb cycombine cytonorm reversed_shift oracle"
```
