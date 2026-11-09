# Handoff — SpectraWeaveR unmixing beta, dashboard and plan cleanup

**Written:** 2026-11-09 13:34 CET · **Branch:** main @ 6c99bd6 · **Status:** beta unmixing, report, diagnostics, correction wrappers and dashboard work on real data; nothing committed; 3 protocol-audit follow-ups wait for the user

## Check first
Run `git log --oneline 6c99bd6..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
One-stop spectral-flow package. Active /goal: working unmixing wrapper, unmixing report and diagnostics, and correction of unmixing errors, tested on the Mendeley, OMIP-069v2 and OMIP-109 data. Latest AutoSpectral (1.8.8, commit 3c504c4) is the gold standard.

## Next action
Ask the user the three protocol-audit questions in Open items (control-file review, `speed` check, channel/gate columns). Do not rerun any run before they answer.

## State
- Uncommitted: 120 paths (all work this session; the user has not asked for a commit).
- Tests: full suite 1847 expectations, 11 known failing older files (annotate, assistant, batch_correct, cluster, composable_pipeline, differential, gate, gating_utils, pipeline, qc, targets_backend); the same 11 failed at the original commit by count, not compared by name. New test files pass (dashboard 112, report 152, figures 122).
- Run: finished runs in `experiments/runs/2026-11-09_*` (never edit). Current dashboards: `experiments/runs/2026-11-09_dashboard_v3/`. Run via SLURM (account div2-para, partition all) with `pixi run Rscript`.
- Results: Mendeley, OMIP-109, OMIP-069v2 cells and beads, and the correction run all completed — see `experiments/benchmarks/RESULTS.md` for inspected figure paths per conclusion.
- Assumed, not run: speed ("fast" vs "slow") takes effect in AutoSpectral's joint pipeline; the runner's new `report/dashboard.html` step end to end; dashboards in browsers other than Firefox; that the 11 failing files are the same by name as at the original commit.
- Shiny app left running in the background on this node (127.0.0.1:38211); stop with `pkill -f serve.R`.

## Locked decisions
- Settings for OMIP runs: OLS + AutoSpectral fast + slow, seed 1, 50,000 events, pseudo-log sigma 1000, 200 bins — user, 2026-11-09.
- OMIP-069v2: one control file `experiments/benchmarks/03_omip069v2-40color/control_combined.csv`, two runs (cells set, beads set), both unmixed — user, 2026-11-09.
- SIF limits: warn 4, fail 6 (AutoSpectral article 01), passed explicitly via `sw_unmix_limits_autospectral()`; other metrics stay not_judged — user, 2026-11-09.
- The seven `sw_unmix_assess_*`/`sw_unmix_compare*` wrappers are internal; roadmap B3 to B10 parked; flowstate and batch-side code frozen — user, 2026-11-09.
- `sw_pipeline_run_all` defaults (cofactor 6000, n_metaclusters 20, seed 42, gate_node "lymphocytes", unmix_method "AutoSpectral") stay as they are — user, 2026-11-09.
- Dashboard: HTML file plus Shiny app; what-if SIF sliders display-only; several runs side by side; AutoSpectral tab uses tables interactively and links the images (not embedded); no variants view — user, 2026-11-09.
- Standing: wrap AGPL/GPL code, never copy; impactful choices are required arguments with no default; plot-first (cite an inspected figure or write `unsupported`); no commits unless asked.

## Dead ends — do not redo
- Renaming bead variants in place after extraction — the variants object has many name-keyed parts; the runner now derives a per-set control file instead (`experiments/benchmarks/_common/run_unmix.R`).
- SIF as the largest entry per row of the hotspot matrix — Mage et al. define SIF as the diagonal only; code now uses the diagonal.
- `sw_unmix_sample_events` on all unmixed files unaligned — OLS output has no AF channels; align on shared `-A` channels.

## Read first
- `experiments/benchmarks/RESULTS.md` — results and open items.
- `ROADMAP.md` — tracks, gold-standard rules, register of user choices, unconfirmed defaults table.
- `R/unmix_report.R` (`max_sif`, `sw_unmix_limits_autospectral`), `R/unmix_dashboard.R`, `R/unmix_dashboard_data.R`, `inst/dashboard/dashboard.html`, `R/unmix_app.R`.
- `PROGRESS.md` Session 14 — what changed this session.

## Open items
- Review the OMIP-109 control file (4 blank markers, AF `universal.negative` empty, cFluor BYG750 fell back to the legacy path on unreviewed gates), then rerun into a new folder — needs the user's go; lean (not a decision): review and rerun.
- Check whether `speed` takes effect in the joint pipeline; if not, relabel fast vs slow — needs the user's go; lean (not a decision): check first, it is read-only.
- Drop the filled `channel` and gate columns from the OMIP-069v2 control file (article 02) and rerun — needs the user's go; lean (not a decision): drop only if AutoSpectral 1.8.8's automated path ignores them.
- Nineteen pre-existing defaults on impactful choices are listed in ROADMAP as unconfirmed — return when the user confirms each; lean (not a decision): confirm in one table.
- Run-folder gaps: the SpectraWeaveR sha in `provenance.csv` is NA; runs lack uniform `sessionInfo.txt` — return at the next runner edit.
- Benchmark 01 README may still carry draft wording; the OMIP-109 dashboard tabs were not inspected.
