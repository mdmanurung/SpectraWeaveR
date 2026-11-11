# Handoff — diagnostics for technical problems and batch effects (cytoscan-diagnostics worktree)

**Written:** 2026-11-11 (after the four-gap round) · **Branch:** worktree-cytoscan-diagnostics @ fdccbbc + uncommitted G37 · **Status:** G36 committed; rolling and bridge reference (G37) run, uncommitted; ledger 124/150

## Check first
Run `git log --oneline fdccbbc..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Active /goal: "work autonomously until we can consistently diagnose for any technical problems including but not limited to batch effects extensively". Bar per fault class: real-file AUROC mean >= 0.90, lower 95% CI >= 0.80 (D1). tusoskill ledger target 150 (G4); stop gate: `/usr/bin/python3.12 /home/mdmanurung/.claude/skills/tusoskill/scripts/biomni_method_state.py stop-check --root experiments/runs/2026-11-09_tusoskill`.

## Next action
Commit G37 (08_rolling_bridge.R, figures, docs) when the user asks. Then pick the next gap from Open items with the user. Nothing is running.

## State
- G37 rolling and bridge reference (2026-11-11, uncommitted): no arm meets the rule. Means over 4 bridge controls, false flags at |z| >= 3: frozen 10.9%, rolling 9.7%, bridge 14.8%, in-sample 5.7%; 80% fault AUROC 0.88 / 0.59 / 0.93; only tcu1 as bridge lowers flags (8.5%, AUROC 0.90, wide CI); excess sits in batches 11, 14, 15, 17 — `reference/figures/rolling_bridge_cutoff.png`, `reference/figures/rolling_bridge_null_z.png`; run `experiments/runs/2026-11-11_rolling-bridge`.
- G36 frozen calibration (2026-11-11, committed fdccbbc): clean later batches give frozen |z| >= 3 in 9.4% of file x marker values (reference in-sample 4.1%, cohort 5.2%), no time trend (Spearman 0.12/0.03); hit marker 75% vs 1% cohort at |z| >= 3 — `reference/figures/frozen_cutoff.png`, `reference/figures/frozen_null_z.png`; run `experiments/runs/2026-11-11_frozen-calibration`. Stop rule applied: documented limit.
- Round G20–G35 (2026-11-10/11), committed 7d2fea7..74e456d: (A) `bimodality`/`bimodality_z` in `sw_fcs_signal_audit` — 07 split 0.58, dip fallback 0.51 (not packaged); real 1.00/0.65/0.85 at 0.5/1/1.5 SD, only on already-bimodal markers — `reference/figures/bimodal_sim_tools.png`, `reference/figures/bimodal_limits_real.png`; runs `experiments/runs/2026-11-10_{bimodal-07,bimodal-limits-real,dip-07,dip-limits-real}`. (B) `min_separation` gate in `sw_negative_audit` (per-channel median, `split_ok`), 10 of 40 Aurora channels scored, 06/08/14 AUROC identical — `reference/figures/negative_gate.png`; run `experiments/runs/2026-11-10_negative-gate`. (C) `rest_emd`/`rest_emd_z` in `sw_batch_scores` — CBF batch split 0.48/0.60 at 0.5/1 SD, `pair_emd` 0.50 at 20 seeds (the 0.68 was 5-seed noise) — `reference/figures/cbf_batch_rest.png`; run `experiments/runs/2026-11-10_cbf-batch-rest`. (D) `sw_signal_reference_fit/save/load` + required `reference` arg on `sw_fcs_signal_audit` — 80% faulty: frozen 0.93 ± 0.06 vs cohort 0.16, nulls 0.51 — `reference/figures/frozen_reference.png`; run `experiments/runs/2026-11-10_frozen-reference`.
- Tests after the round: signal_reference 36, fcs_signal_audit 80, diagnose 38, negative_audit 62, batch_scores 40, 0 fail.
- `min_separation = 3` was added to Glehr 08/09 and benchmark 10 so they still run; they were not rerun, and a rerun gates channels below 3 (e.g. Aurora PB-A).
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
- Results: CD8-dim refuted — Fortessa CD4+ T cells (true CD8 negatives) have AF700 SD 0.51 vs panel negative 0.53 (Aurora 0.094), shifted up — `reference/figures/glehr_af700_subset_sd.png`, `reference/figures/glehr_af700_subset_density.png`; run `experiments/runs/2026-11-10_glehr-af700-cd8dim-v2` (the run without v2 split CD3 inside CD3+ files; discarded).
- Results: APC→AF700 under-compensation on the Fortessa (+1374 vs -142 in the CD3-APC single stain) predicts ~407 of the 536-unit AF700 offset in CD4+ T (table fallback); linear fit does not narrow AF700 (0.51 → 0.54) — `reference/figures/glehr_af700_comp_shift.png`, `reference/figures/glehr_af700_comp_slope.png`; run `experiments/runs/2026-11-10_glehr-af700-comp-v2` (same numbers as the run without v2; only a label changed).

## Locked decisions
- G20–G35 (four-gap round) and G1–G19 in `plans/2026-11-10_diagnostics-next.md`; G1–G13: spread SD ratio; reference_group; FAUST conventions; keep 150; median over groups; group EMD pool; FAUST rank norm; 37 markers; confirm refgroup details; two-way residual; forest depth 3; prop_by_batch range; correct on 37 — user, 2026-11-10.
- Earlier rows D1–D11, GI1–GI9 in `plans/2026-11-10_diagnostics-next.md`, `plans/2026-11-10_gating.md` — user.
- Standing: impactful choices are required args (no default) and the user's; GPL/AGPL wrapped never copied; plot-first; no sample IDs or data committed (rename f1.., d1..); real data read-only in `/exports/para-lipg-hpc/mdmanurung/sw_data`; never edit a running script; absolute script paths in background runs; commits end with the Co-Authored-By and Claude-Session lines; no bare git stash.

## Dead ends — do not redo
- Further AF700 width tests on Glehr BAD — thread closed by the user 2026-11-10 after noise floor, single stains, CD8-dim and compensation (`reference/diagnostics.qmd`); reopen only with new data (e.g. a CD8-AF700 control).
- One-way within-donor z for small file faults on CBF — 0.47–0.71 (`reference/figures/cbf_reference_group.png`); batch variation stays inside each donor.
- `reference` clusters to fix the rare-population gap — 0.86, no gain; the cause was the max over clusters.
- Spread as the AF700/PB cause — refuted above.
- More events for the ±0.5 SD two-mode split — 5000 events gave 0.69, donor variation is the limit.
- Two-way (donor × batch medpolish) residual for small file faults on CBF — best 0.68, below one-way on both faults where one-way helped (`reference/figures/cbf_reference_group.png`); code kept, committed 2026-11-10.
- Two-way residual for the batch-split fault — removes the effect under test; sw_batch_scores must stop on a data-frame reference_group.
- Rolling reference for the frozen z — absorbs a persistent fault (AUROC 0.59) and barely cuts false flags (G37). Single-file bridge — lifts detection but adds false flags (14.8%) (G37).
- Shape scores (bimodality coefficient, Hartigan dip) for a ±0.5 split in simulation — the mix stays one peak; 0.58 and 0.51 (`reference/figures/bimodal_sim_tools.png`).
- `rest_emd_z` and more seeds for the CBF batch split — 20 seeds: at chance at 0.5 SD, 0.60 at 1 SD (`reference/figures/cbf_batch_rest.png`); closed by stop rule G31.

## Read first
- `plans/2026-11-10_diagnostics-next.md` — decisions and plan.
- `reference/diagnostics.qmd` — coverage table and every real-data result with figures.
- `plans/2026-11-09_faust-wrapper.md` — FAUST plan.
- `experiments/runs/2026-11-09_tusoskill/reports/final_round2.md` — outdated; write `final_round3.md` at the end.

## Open items
- Ledger 124 → 150: 26 iterations left; replenish from diagnosis. The AF700 thread is closed (user, 2026-11-10): shift fits APC under-compensation, width stays a documented open cause. lean (not a decision): ask before any new impactful choice.
- Frozen reference: offset doubles clean false flags (G36); rolling and single-file bridge do not fix it (G37). A bridge averaging several files per batch is untested. lean (not a decision): document, do not pursue without new bridge data.
- `spread_z` is not in the signed-check list of `sw_diagnosis_top` (ranks by raw value) — parked from plan B.
- `fig_batch_rest.png` and `fig_frozen_reference.png` have a clipped x label (cosmetic).
- Old `HANDOFF.md` (unmixing beta, main @ 6c99bd6) is replaced here; its open items live in git: `git show a3711c0:HANDOFF.md`.
