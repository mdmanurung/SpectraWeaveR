# Handoff — diagnostics for technical problems and batch effects (cytoscan-diagnostics worktree)

**Written:** 2026-11-10 (updated after FAUST) · **Branch:** worktree-cytoscan-diagnostics @ c16607f · **Status:** gap round; two-way and FAUST done and committed; ledger 116/150

## Check first
Run `git log --oneline c16607f..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Active /goal: "work autonomously until we can consistently diagnose for any technical problems including but not limited to batch effects extensively". Bar per fault class: real-file AUROC mean >= 0.90, lower 95% CI >= 0.80 (D1). tusoskill ledger target 150 (G4); stop gate: `/usr/bin/python3.12 /home/mdmanurung/.claude/skills/tusoskill/scripts/biomni_method_state.py stop-check --root experiments/runs/2026-11-09_tusoskill`.

## Next action
Pick the next gap from Open items with the user. Nothing is running.

## State
- All work committed through c16607f: two-way residual (G10), FAUST forest depth (G11), G13, `match_by = "marker_state"` (G14), FAUST script and results.
- Tests: anomaly_features 43, batch_scores 33, cluster_scores 41, fcs_signal_audit 68, faust_wrap 54, 0 fail (via `pkgload::load_all`; `devtools` fails on an rlang 1.1.5 vs >= 1.1.7 clash in the scratchpad rlib).
- Runs finished: FAUST `experiments/runs/2026-11-10_faust` (exit 0, ~1 h 50 min; its compare used exact match and is superseded); compare redone from its labels in `experiments/runs/2026-11-10_faust-markerstate-v2` (`-markerstate` without v2 has wrong panel labels). Two-way benchmark `experiments/runs/2026-11-10_cbf-refgroup-v2`.
- Results: FAUST per-file ARI median 0.77 (0.19–0.89), reference files B3 0.49, B4 0.28 — `reference/figures/faust_ari.png`; 43/43 phenotypes matched, donor batch-median shares move ≤ 2.7 points — `reference/figures/faust_prop_by_batch.png`; text in `reference/faust.qmd`.
- Run commands: `R_LIBS=<scratchpad>/rlib /exports/para-lipg-hpc/mdmanurung/SpectraWeaveR/.pixi/envs/default/bin/Rscript`; after `roxygen2::roxygenise()` revert churn: `git checkout -- DESCRIPTION man/sw_channel_is_fluor.Rd man/sw_channel_is_signal.Rd man/sw_diff_abundance.Rd man/sw_unmix_extract_af.Rd man/sw_unmix_save_figures.Rd man/sw_unmix_setup.Rd` (if DESCRIPTION has real edits, drop only the added RoxygenNote line).
- Results: doublets `area_height_z` 1.00 from 2% at q0.99 — `reference/figures/cbf_doublets.png`
- Results: spillover spread does not explain the Fortessa AF700/PB residual (all ratios ~0.9–1.5 both instruments) — `reference/figures/glehr_spread_af700_pb.png`
- Results: `spread_z` finds injected spread 0.98 (1 neg SD, q0.95), 0.76 at 0.5 SD — `reference/figures/aurora_spread_real.png`
- Results: one-way `reference_group` (donor) best 0.71; four donor/batch-limited faults stay below bar — `reference/figures/cbf_reference_group.png`
- Results: full coverage table with figures in `reference/diagnostics.qmd` (lines ~73–105).
- Results: two-way `reference_group` best 0.68, no fault at bar — same figure, `reference/diagnostics.qmd` after the one-way paragraph.
- Results: detector noise floor accounts for PB (unstained/panel ratio 1.23) but only part of AF700 (0.55; unstained SD 0.09 vs panel 0.53) — `reference/figures/glehr_unstained_noise.png`; run `experiments/runs/2026-11-10_glehr-unstained-noise-v2` (the run without v2 used a 1–99% gate that kept granulocytes; discarded).
- Results: AF700 single stains — Fortessa AF700 SD unstained 0.09, AF700-stain CD3- 0.13, panel negative 0.53; no single dye spreads more on the Fortessa (APC 5.9 vs 4.3 only) — `reference/figures/glehr_af700_stain_sd.png`, `reference/figures/glehr_af700_single_spread.png`; run `experiments/runs/2026-11-10_glehr-af700-stain`.

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
- Ledger 116 → 150: 34 iterations left; replenish from diagnosis (candidates: CD8-dim cells merging into the Fortessa AF700 negative, or summed spread in the full panel). lean (not a decision): ask before any new impactful choice.
- Batch split ±0.5 SD (0.68) — stays a documented limit unless the user wants another design. lean (not a decision): document.
- Old `HANDOFF.md` (unmixing beta, main @ 6c99bd6) is replaced here; its open items live in git: `git show a3711c0:HANDOFF.md`.
