# Handoff — cyTether/HEB benchmark: HEB beats fair cyCombine on T6 EMD (F12 + F13); O8 OpenProblems check running

**Written:** 2026-11-11 · **Branch:** cytether-benchmark @ 72d648d · **Status:** F12 scale term in the package; F13 sub-cluster arm (benchmark-side prototype) beats cycombine_ref on held-out EMD in 4 of 4 T6 folds; O8 (OpenProblems scores for sc and sc_cl) scoring on SLURM, figure script written, not yet run

## Check first
Run `git log --oneline 72d648d..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Show HEB (reference-anchored hierarchical empirical Bayes corrector, `R/corrector_heb.R`) is not inferior to cyCombine, CytoNorm and CytofIn on synthetic, T6 (CytoBatchFlagR) and OpenProblems tasks; tusoskill loop (50-iteration floor, ledger at 22/50).

## Next action
Check O8 jobs 25994245 (corrections), 25994248 (scoring), 25994249 (scaling) with `sacct -j 25994245,25994248,25994249 -X`. When all are COMPLETED, run `<venv python> -I experiments/benchmarks/06_cytether-cytobatchflagr/07b_fig_op.py experiments/runs/20261111_t6_op_sc/scores`, open `fig_op.png`, read `reading_rule.csv`, and apply the O8 rule (plan row O8). Label the result "same data, new metrics and held-out pairs", not out-of-sample.

## State
- Uncommitted: O8 plan row in `plans/cytether_t6_cytobatchflagr_plan.md`; new `experiments/benchmarks/06_cytether-cytobatchflagr/07b_fig_op.py`; this file. Pre-existing, not this work: `experiments/datasets/registry.csv`, untracked `man/*.Rd`, `__pycache__/`.
- Tests: pass — corrector_heb 51, corrector_heb_shared 18, corrector 66 (`pkgload::load_all` with `R_LIBS=.../Rlib_competitors`; devtools needs newer rlang).
- Run: all compute on SLURM via `experiments/slurm/run_lines.sbatch`; generators in `experiments/runs/20261110_slurm_lists/` (git-ignored), O8 lists from `make_op_sc_lists.py`.
- Package: F12 `scale_model = "none" | "log_sd"` (required, no default) — REML mixed model on per-aliquot log SD, cells x' = c + exp(−δ_b)(x − shift − c); A13 `anchor_batches = "all" | "unreferenced"` (required when anchors given).
- Results: T6 held-out EMD per fold (tcs3, tcs4, tcu1, tcu2) — heb_kv_bap 5.40/4.37/4.99/4.99; sc 5.24/4.21/4.66/4.71; sc_cl 4.87/4.11/4.29/4.50; cycombine_ref 5.11/4.51/4.50/4.64 — `experiments/runs/20261111_t6_emd_cluster/scores/fig_emd.png`
- Results: shape (spread) explains 51% of the heb_kv_bap − cycombine_ref EMD gap; MAD-ratio scale prototype made total EMD worse — `experiments/runs/20261111_t6_emd_decompose/fig_emd_decompose.png`
- Results: F12 on synthetic costs +0.0007 to +0.0011 held-out RMSE vs lead (passes) — `experiments/runs/20261111_sc_s17/scored/summary/primary.png`
- Results: A13 anc_u passes every anchor rule on dev, shrink, S17, S18, S18w — `experiments/runs/20261111_anc_u_s17/scored/summary/primary.png`, `experiments/runs/20261111_s17_mismatch_u/fig_mismatch.png`
- Flag: sc_cl tcs4 fell back to stage 1 for 136,200 cell-markers (sub-clusters without an effect); unsplit populations differ from the plain sc run by ~9e-7 (refit noise).
- Flag: T6 methods were tuned on these 18 batches (F9–F13); every T6 win risks overfitting until Z247 or another dataset confirms it. EMD is the only T6 metric scored for sc/sc_cl so far.
- Flag: S17 anc_u anchor-only (nd0) cov90 0.43 — table fallback: `experiments/runs/20261111_anc_u_s17/scored/scores.csv`, no class-coverage figure.

## Locked decisions
- Plans: `plans/cytether_t6_cytobatchflagr_plan.md` (N1–N7, O1–O8, S1–S4), `plans/cytether_fix_plan.md` (F1–F13), `plans/cytether_anchor_plan.md` (A0–A13, leak rule) — user, 2026-11-09/11.
- F12 log-SD scale term; F13 sub-clusters k = max(1, min(4, floor(median calibrating cells per aliquot / 100))), SOM seed 1, rlen 10, fallback to stage 1 — user, 2026-11-11.
- O8 rule: sc_cl ≥ cycombine_ref mean scaled score in ≥ 3 of 4 folds, and no biology metric more than 0.10 below cycombine_ref pooled — user, 2026-11-11.
- All compute on SLURM; max 5 concurrent jobs. Test seeds 431–440 and Z247 sealed; only `evaluate.py` and `diag_mismatch.py` read truth; impactful choices need an approved plan row; plot-first.

## Dead ends — do not redo
- MAD-ratio scale prototype (S1–S4) — worse total EMD; replaced by F12 log-SD model.
- Boundary prior alone (F6, c022) — residual collapses instead.
- Anchors in every batch (A0–A10 anc) and anchor discrepancy tau_d (A11 anc_d) on S17 — mode mismatch, slope −0.95.
- Exact CD3-relativisation; pre-normalising channels by CD3; cyCombine `correct_data` matrix path.

## Read first
- `plans/cytether_fix_plan.md` rows F12, F13; `plans/cytether_t6_cytobatchflagr_plan.md` row O8.
- `R/corrector_heb.R` (`.sw_heb_fit_designs`, `.sw_heb_scale`), `R/corrector.R` (`.sw_apply_effects` scale step).
- `experiments/benchmarks/06_cytether-cytobatchflagr/correct_lib.R` (`run_heb_subclusters`), `07b_fig_op.py`.
- Ledger: `/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/tusoskill_runs/cytether_heb_20261109/experiments/results.jsonl` (c020–c024; c025+ not logged).

## Open items
- O8 result — see Next action.
- Move sub-clustering (F13) into the package — needs a user decision and plan row.
- T6 strip array 25992174 failed; dependants 25992175, 25992197 at DependencyNeverSatisfied — scancel, read logs, fix, resubmit.
- N1–N3 HEB candidates — rescore (3 FAILED tasks in 25991343).
- Ledger records c025+ (anc, anc_d, anc_u, sc, sc_cl) — not logged.
- EMD non-inferiority margin vs CytoNorm/cyCombine — needs user plan row.
- Independent test on Z247 — sealed until the final audit.
