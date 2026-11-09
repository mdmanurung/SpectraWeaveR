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
