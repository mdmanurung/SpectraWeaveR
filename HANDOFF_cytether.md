# Handoff — cyTether/HEB benchmark: HEB beats fair cyCombine on T6; anchors fail on real data

**Written:** 2026-11-11 15:01 CET · **Branch:** main @ f51f17c · **Status:** cytether-benchmark merged into main (f51f17c); N3 min_cells 50 rescore paused mid-diagnosis

## Check first
Run `git log --oneline f51f17c..HEAD` and `git status --short`, with `f51f17c` from the line above. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it. `HANDOFF.md` is the BMV unmix session's file, not this one.

## Goal
Show HEB (`R/corrector_heb.R`) is not inferior to cyCombine, CytoNorm and CytofIn on synthetic, T6 (CytoBatchFlagR) and OpenProblems tasks; tusoskill loop (50-iteration floor, ledger at 22/50).

## Next action
Find why `experiments/runs/20261110_n3_mc50_dev/scored/scores.csv` has no `heldout_ref_rmse` column (crash at `experiments/benchmarks/05_cytether-synthetic/summarize.py:44`); start at `evaluate.py:90-125`, where the column is written only when holdout rows exist. Then rescore lines 7–9 of the old `scores.txt` (N3 min_cells 50: dev, shrink, s17).

## State
- Uncommitted: none of this work. `R/gate.R`, `R/utils.R`, `tests/testthat/test-utils.R`, `man/sw_channel_is_*.Rd` belong to another session; do not commit them.
- Tests: pass on f51f17c — corrector_heb 58, corrector_heb_shared 18, corrector 66, corrector_paired_shift 63 (`R_LIBS=/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/Rlib_competitors .pixi/envs/default/bin/Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file(...)'`).
- Run: all outputs in `experiments/runs/` (moved from the worktree; 891 symlinks repointed, 0 broken). SLURM via `experiments/slurm/run_lines.sbatch`; the old git-ignored lists in `20261110_slurm_lists/` moved too.
- Worktree `.claude/worktrees/cytether-benchmark` removed and branch `cytether-benchmark` deleted (was f51f17c).
- Results: T6 held-out EMD, sc_cl beats cycombine_ref 4/4 folds (4.87/4.11/4.29/4.50 vs 5.11/4.51/4.50/4.64) — `experiments/runs/20261111_t6_emd_cluster/scores/fig_emd.png`
- Results: O8 OpenProblems mean, sc_cl ≥ cycombine_ref 4/4 folds, gaps 0.004–0.009, no uncertainty; same data, not out-of-sample — `experiments/runs/20261111_t6_op_sc/scores/fig_op.png`
- Results: A14 anchors fail A8 (leak share anc 0.66, anc_u 0.61; need 0.85); anc_u worse EMD 3/4 folds — `experiments/runs/20261111_t6_strip/emd_strip/fig_emd_strip.png`, `experiments/runs/20261111_t6_strip/leak_heb_kv_bap_anc_u/fig_leak.png`
- Results: F12 log-SD scale costs +0.0007 to +0.0011 synthetic RMSE (passes) — `experiments/runs/20261111_sc_s17/scored/summary/primary.png`
- Assumed, not run: the N3 mc50 cause (every method abstains at min_cells 50) is a guess; why OpenProblems functional_marker_preservation is 0 and cLisi 1 for all methods; whether `05_scale_op.py` clips to [0, 1] (values 1.09 and −0.07 seen).
- Man pages `man/sw_corrector_heb*.Rd` not generated (devtools needs newer rlang).

## Locked decisions
- Plans: `plans/cytether_t6_cytobatchflagr_plan.md` (N1–N7, O1–O8, S1–S4), `plans/cytether_fix_plan.md` (F1–F13), `plans/cytether_anchor_plan.md` (A0–A15) — user, 2026-11-09/11.
- A15: heb `unbridged = "abstain" | "stop"` required; anchors experimental; heb only — because anchor rows are per sample × marker, so population effects in anchor-only batches are not identified (user, 2026-11-11).
- All compute on SLURM, at most 5 concurrent jobs; test seeds 431–440 and Z247 sealed; only `evaluate.py` and `diag_mismatch.py` read truth; impactful choices need a plan row; plot-first (user).

## Dead ends — do not redo
- Internal anchors standing in for references (A0–A14: anc, anc_d, anc_u) — S17 nd0 cov90 0.43; T6 leak 0.61–0.66.
- MAD-ratio scale prototype (S1–S4) — worse total EMD; replaced by F12.
- Boundary prior alone (F6) — residual collapses instead.
- CD3 relativisation and pre-normalising channels by CD3 — no single-stain controls; adds donor×batch wobble.

## Read first
- `plans/cytether_anchor_plan.md` rows A14, A15; `plans/cytether_fix_plan.md` rows F12, F13.
- `R/corrector.R` `.sw_check_audit` (unbridged gate); `R/corrector_heb.R` roxygen header.
- `experiments/benchmarks/06_cytether-cytobatchflagr/correct_lib.R` `run_heb_subclusters` (F13 prototype).

## Open items
- Move F13 sub-clustering into the package — needs a plan row; lean (not a decision): yes, it carries the T6 win.
- OpenProblems metric sanity (zeros, saturation, clipping) — before any OP claim; lean: check `05_scale_op.py` first.
- Ledger records c025+ (anc, anc_d, anc_u, sc, sc_cl, A14, A15) — not logged; lean: log before more iterations.
- EMD non-inferiority margin vs CytoNorm/cyCombine — needs a plan row; lean: reuse Q7's 10%.
- Independent test on Z247 — sealed until the final audit; lean: the only out-of-sample check for the T6 tuning.
- Push main to origin — not done; lean: ask the user.
