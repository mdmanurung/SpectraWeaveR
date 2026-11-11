# Handoff — cyTether/HEB benchmark: HEB beats fair cyCombine on T6 (EMD and OpenProblems); anchors fail on real data

**Written:** 2026-11-11 · **Branch:** main (cytether-benchmark merged at d03f422, branch deleted) · **Status:** F12 + F13 HEB (sc_cl) beats cycombine_ref on T6 held-out EMD 4/4 folds and OpenProblems mean 4/4 folds (O8 pass); internal anchors fail A8 on T6 (A14), so HEB now abstains on unbridged batches (A15, `unbridged` argument)

## Check first
Run `git log --oneline d03f422..HEAD -- R/corrector* tests/testthat/test-corrector* experiments/benchmarks/0[56]_cytether* plans/cytether*` and `git status --short`. This file is the cyTether handoff; `HANDOFF.md` belongs to the BMV unmix work. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Show HEB (reference-anchored hierarchical empirical Bayes corrector, `R/corrector_heb.R`) is not inferior to cyCombine, CytoNorm and CytofIn on synthetic, T6 (CytoBatchFlagR) and OpenProblems tasks; tusoskill loop (50-iteration floor, ledger at 22/50).

## Next action
N3 min_cells 50 diagnostic: scoring 25991343 tasks 6–8 (`scores.txt` lines 7–9, runs `experiments/runs/20261110_n3_mc50_{dev,shrink,s17}`) crashed in `summarize.py:44` with `KeyError: 'heldout_ref_rmse'` — `scored/scores.csv` lacks that column for this run; find why (likely every method abstains at min_cells 50 so evaluate.py writes no held-out rows), then rescore. N1, N2 and N3 mc20 scored fine (`scored_n12`).

## State
- Runs now on main: the worktree `.claude/worktrees/cytether-benchmark` was removed; SLURM lists in `experiments/runs/20261110_slurm_lists/` were git-ignored and went with it unless copied (generators are small; rewrite as needed).
- Man pages `man/sw_corrector_heb*.Rd` not generated (devtools needs newer rlang); run `devtools::document()` once rlang is updated.
- Tests: pass — corrector_heb 58, corrector_heb_shared 18, corrector 66, corrector_paired_shift 63 (`pkgload::load_all` with `R_LIBS=.../Rlib_competitors`).
- Package: A15 `unbridged = "abstain" | "stop"` (required for heb, no default); anchors documented as experimental.
- Results: O8 OpenProblems mean scaled score, sc_cl vs cycombine_ref 0.542/0.533, 0.536/0.532, 0.541/0.535, 0.539/0.533 (4/4 folds; same data, new metrics and held-out pairs, not out-of-sample). Caveats: gaps 0.004–0.009 with no uncertainty; functional_marker_preservation scaled 0 and cLisi 1 for every method, so only ratio_consistent_peaks discriminates (sc_cl 0.93 vs 0.85); some scaled values fall outside [0, 1] despite O2 clipping (unchecked) — `experiments/runs/20261111_t6_op_sc/scores/fig_op.png`
- Results: A14 strip check fails — anc lowers stripped-batch EMD 4/4 folds but leak share 0.66; anc_u worse EMD 3/4, leak 0.61; anchor effects compressed toward 0 — `experiments/runs/20261111_t6_strip/emd_strip/fig_emd_strip.png`, `experiments/runs/20261111_t6_strip/leak_heb_kv_bap_anc_u/fig_leak.png`
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
- Plans: `plans/cytether_t6_cytobatchflagr_plan.md` (N1–N7, O1–O8, S1–S4), `plans/cytether_fix_plan.md` (F1–F13), `plans/cytether_anchor_plan.md` (A0–A15, leak rule) — user, 2026-11-09/11.
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
- N3 min_cells 50 rescore — see Next action.
- Move sub-clustering (F13) into the package — needs a user decision and plan row.
- OpenProblems metric check: why functional_marker_preservation is 0 and cLisi 1 for every method, and why scaled values leave [0, 1].
- Ledger records c025+ (anc, anc_d, anc_u, sc, sc_cl, A14, A15) — not logged.
- EMD non-inferiority margin vs CytoNorm/cyCombine — needs user plan row.
- Independent test on Z247 — sealed until the final audit.
