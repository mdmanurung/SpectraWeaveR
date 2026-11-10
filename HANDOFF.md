# Handoff — diagnostics for technical problems and batch effects (cytoscan-diagnostics worktree)

**Written:** 2026-11-10 15:57 CET · **Branch:** worktree-cytoscan-diagnostics @ a3711c0 · **Status:** gap round in progress; two-way residual (uncommitted, agent run) and FAUST full run (hours) still running; ledger 112/150

## Check first
Run `git log --oneline a3711c0..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Active /goal: "work autonomously until we can consistently diagnose for any technical problems including but not limited to batch effects extensively". Bar per fault class: real-file AUROC mean >= 0.90, lower 95% CI >= 0.80 (D1). tusoskill ledger target 150 (G4); stop gate: `/usr/bin/python3.12 /home/mdmanurung/.claude/skills/tusoskill/scripts/biomni_method_state.py stop-check --root experiments/runs/2026-11-09_tusoskill`.

## Next action
Wait for the FAUST full run (`tail full.log` in the scratchpad below; ends `exit N`), then do the FAUST open item. The two-way residual finished, failed the bar and is committed.

## State
- Committed 2026-11-10: two-way residual (G10), G11 FAUST forest depth, G13 row. Untracked: `experiments/integration/aurora_exvivo/15_faust_before_after.R` (commit after the FAUST run).
- Tests after the two-way edits: anomaly_features 43, batch_scores 33, cluster_scores 41, fcs_signal_audit 68, faust_wrap 49, 0 fail (via `pkgload::load_all`; `devtools` fails on an rlang 1.1.5 vs >= 1.1.7 clash in the scratchpad rlib).
- Tests: last full check before the two-way edits: 348 expectations over 11 diagnostic test files, 0 fail; negative_audit 45, faust_wrap 49.
- Running: FAUST full run, launched 15:5x via `setsid nohup bash <scratchpad>/full.sh` (scratchpad = `/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-SpectraWeaveR/99bb6407-3499-4914-a4c4-455aad92c547/scratchpad`); log `full.log` ends `exit N`; output `experiments/runs/2026-11-10_faust`; expected 4–6 h (guess from smoke test). Two-way benchmark run by subagent (pid 2503603 at 15:40), output `experiments/runs/2026-11-10_cbf-refgroup-v2`.
- Run commands: `R_LIBS=<scratchpad>/rlib /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript`; after `roxygen2::roxygenise()` revert churn: `git checkout -- DESCRIPTION man/sw_channel_is_fluor.Rd man/sw_channel_is_signal.Rd man/sw_diff_abundance.Rd man/sw_unmix_extract_af.Rd man/sw_unmix_save_figures.Rd man/sw_unmix_setup.Rd` (if DESCRIPTION has real edits, drop only the added RoxygenNote line).
- Results: doublets `area_height_z` 1.00 from 2% at q0.99 — `reference/figures/cbf_doublets.png`
- Results: spillover spread does not explain the Fortessa AF700/PB residual (all ratios ~0.9–1.5 both instruments) — `reference/figures/glehr_spread_af700_pb.png`
- Results: `spread_z` finds injected spread 0.98 (1 neg SD, q0.95), 0.76 at 0.5 SD — `reference/figures/aurora_spread_real.png`
- Results: one-way `reference_group` (donor) best 0.71; four donor/batch-limited faults stay below bar — `reference/figures/cbf_reference_group.png`
- Results: full coverage table with figures in `reference/diagnostics.qmd` (lines ~73–105).
- Results: two-way `reference_group` best 0.68, no fault at bar — same figure, `reference/diagnostics.qmd` after the one-way paragraph.
- Assumed, not run: wider Fortessa detector noise floor as AF700/PB cause; FAUST full-run outcome; that the FAUST run survives (setsid).

## Locked decisions
- G1–G13 in `plans/2026-11-10_diagnostics-next.md` (last section): spread SD ratio; reference_group; FAUST conventions; keep 150; median over groups; group EMD pool; FAUST rank norm; 37 markers; confirm refgroup details; two-way residual; forest depth 3; prop_by_batch range; correct on 37 — user, 2026-11-10.
- Earlier rows D1–D11, GI1–GI9 in `plans/2026-11-10_diagnostics-next.md`, `plans/2026-11-10_gating.md` — user.
- Standing: impactful choices are required args (no default) and the user's; GPL/AGPL wrapped never copied; plot-first; no sample IDs or data committed (rename f1.., d1..); real data read-only in `/exports/para-lipg-hpc/mdmanurung/sw_data`; never edit a running script; absolute script paths in background runs; commits end with the Co-Authored-By and Claude-Session lines; no bare git stash.

## Dead ends — do not redo
- One-way within-donor z for small file faults on CBF — 0.47–0.71 (`reference/figures/cbf_reference_group.png`); batch variation stays inside each donor.
- `reference` clusters to fix the rare-population gap — 0.86, no gain; the cause was the max over clusters.
- Spread as the AF700/PB cause — refuted above.
- More events for the ±0.5 SD two-mode split — 5000 events gave 0.69, donor variation is the limit.
- Two-way (donor × batch medpolish) residual for small file faults on CBF — best 0.68, below one-way on both faults where one-way helped (`reference/figures/cbf_reference_group.png`); code kept, committed 2026-11-10.
- Two-way residual for the batch-split fault — removes the effect under test; sw_batch_scores must stop on a data-frame reference_group.

## Read first
- `plans/2026-11-10_diagnostics-next.md` — decisions and plan.
- `reference/diagnostics.qmd` — coverage table and every real-data result with figures.
- `plans/2026-11-09_faust-wrapper.md` — FAUST plan.
- `experiments/runs/2026-11-09_tusoskill/reports/final_round2.md` — outdated; write `final_round3.md` at the end.

## Open items
- After FAUST ends: open `fig_faust_ari.png`, `fig_faust_prop_by_batch.png`, `run_info.txt`; add results to `reference/faust.qmd`; commit the 15 script; ledger entry. lean (not a decision): if `exit` is nonzero, fix and rerun into a new dir.
- Ledger 112 → 150: 38 iterations left; replenish from diagnosis (e.g. test detector-noise cause with unstained/negative SD per instrument — needs user approval of the design). lean (not a decision): ask before any new impactful choice.
- Batch split ±0.5 SD (0.68) — stays a documented limit unless the user wants another design. lean (not a decision): document.
- Old `HANDOFF.md` (unmixing beta, main @ 6c99bd6) is replaced here; its open items live in git: `git show a3711c0:HANDOFF.md`.
