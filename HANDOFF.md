# Handoff — batch-effect and reference diagnostics build (workflow A done except P4b)

**Written:** 2026-11-12 15:20 CET · **Branch:** main @ fce4e30 · **Status:** 16 of 17 PRs merged into main; P4b not built; Workflow B (benchmarks, docs) not started. Nothing pushed.

## Check first
Run `git log --oneline fce4e30..HEAD` and `git status --short`, with `fce4e30` from the line above. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Close the batch-effect and technical-diagnostic gaps found in the 2026-11-12 status review. Done = every export in the approved plan exists with tests and a CLAUDE.md row, benchmarks 12 and 13 and the scorecard rows are run and their figures inspected, and ROADMAP and docs are synced.

## Next action
Build P4b alone: no branch `feat/p4b-holdout-influence` exists. Its agent died on "Could not read the repository git config to neutralize filter drivers" in wave 4. Spec: `PRS.P4b` in the workflow script (path under Read first). Write a one-PR Workflow from that spec, or ask the user the open questions below first.

## State
- Uncommitted: untracked only: `.mycelium-extra/`, `plans/2026-11-12_unmix-zip-script-mapping.md`, two bundle zips under `plans/`. `tests/testthat/_problems/test-batch_correct-643.R` is a leftover failed-test file (not shown by git status); delete it.
- Tests: integrators report no failures against the baseline (0 failures, 42 skips; skips are missing isotree, mclust, ggiraph and others). The reporter prints no pass total.
- Merged (`--no-ff`): P0a to P0d, P1, P2a to P2c, P3a to P3c, P4a, P5a to P5c, P6. Wave 5 left 68 warnings "Unknown or uninitialised column: source_batch" in `test-ref_composition.R` and `test-ref_vs_study.R`.
- Results: 18 of 32 real-file fault classes meet the D1 bar — `reference/figures/scorecard.png` (opened 2026-11-12). Nothing from the new functions is benchmarked yet: `unsupported: no inspected figure`.
- Assumed, not run: I did not rerun any test after the workflow. The pass and skip numbers are the integrator agents' reports. The critic checked exports against the plan by reading, not by running.
- 33 git worktrees exist, including the `feat/p*` ones. Remove them when done.

## Locked decisions
- Unfreeze A2b, A2c, A3, A5 and A6; report-only on the six weak classes (user, 2026-11-12).
- DR24 four shift families; DR26 and DR27 built; DR28 unit donor or aliquot, never cell; DR37 caller transforms; DR38 rank within check × marker; DR30 plot built; DR33 five diagnose modules (user, 2026-11-12).
- A5 and A6 are ranked score rows, no flags; kBET and LISI built from the paper definitions only (user, 2026-11-12).
- DR12 re-extract CBF at 20,000 events; DR15 `min_events` 500; `n_splits` and `n_draws` 20; probs c(0.01, 0.99); `min_batches_per_reference` 3; k and perplexity 30 (user, 2026-11-12; benchmark values only).
- Heavy runs go to SLURM; code merges one branch per PR, `--no-ff` (user, 2026-11-12).
- DR39 `spread_ratio` ranks by abs(log) — my addition, approved with the plan.

## Dead ends — do not redo
- `pixi run` inside a worktree — worktrees have no `.pixi/envs`. Use the root `Rscript` with `R_LIBS_USER=<root>/.pixi/Rlib` and `R_ENVIRON_USER=<root>/.pixi/Renviron.none`.
- Resuming `wf_76616695-3d1` to rerun P4b — replay would also rerun cached integrators. Use a fresh one-PR workflow.

## Read first
- `/home/mdmanurung/.claude/plans/atomic-baking-journal.md` — approved plan (outside the repo; copy into `plans/` if it should be kept).
- `plans/2026-11-11_reference-and-batch-diagnostics.md` — inventory and DR1 to DR38 (status header is out of date).
- Workflow script with `PRS` specs: `/home/mdmanurung/.claude/projects/-exports-para-lipg-hpc-mdmanurung-SpectraWeaveR/9d388268-c044-43a6-b8c7-9f60404f05fc/workflows/scripts/batch-diagnostics-build-wf_76616695-3d1.js`.
- Full workflow result with all implementer open questions and the critic: `/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-SpectraWeaveR/9d388268-c044-43a6-b8c7-9f60404f05fc/tasks/wbd17h8d3.output` (temp file; may be gone).
- `plans/2026-11-12_unmix-zip-script-mapping.md` — zip 1.8.8 script to `sw_unmix_*` mapping.

## Open items
- P4b (`sw_ref_split`, `sw_ref_influence`, `sw_corrector_compare(exclude=, n_quantiles)` required) — rerun; migrate `test-corrector.R` lines 339, 381, 393, 409 and `test-corrector_paired_shift.R` line 489; grep `experiments/benchmarks/05_cytether-synthetic` and `06_cytether-cytobatchflagr`. Lean (not a decision): fresh workflow, one PR.
- Implementer questions (first ten seen): `keep_order` required in `sw_io_subsample`? `source_batch` must name a manifest batch? `sw_corrector_fit` warn on synthetic references? `effect` beyond the fence or from the median? Count `below_min_effect` in `sw_flag_iqr_summary`? `sw_synref_as_manifest` wrap `sw_sample_manifest` now? Does "equal events per file" mean `strata = "file"` only? Lean (not a decision): yes to wrapping, otherwise keep as built.
- DR9 `min_effect` values — need the CBF null arm (a: pseudo-aliquot floor, b: shuffled batch labels; 90th, 95th, 99th percentiles); the user picks null and percentile.
- Workflow B: CBF re-extract at 20,000 events, install mclust and isotree (ask first), benchmarks 12 and 13 on SLURM, material swap, `scorecard_map.csv` rows, vignettes, ROADMAP and PROGRESS sync. Population label source for benchmarks (gate or k-means k) not yet chosen.
- Housekeeping: P0c merge commit lacks the Co-Authored-By line; CLAUDE.md export count and `plans/2026-11-11...` header need a sync.
- Previous handoff (BMV unmix comparison, closed 2026-11-11) is replaced; it is in git at `6845d1d:HANDOFF.md`. `HANDOFF-diagnostics.md` is untouched.
