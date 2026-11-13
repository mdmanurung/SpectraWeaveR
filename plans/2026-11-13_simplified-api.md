# Plan: a small user-facing API (about 35 functions instead of 241)

Date: 2026-11-13. Status: **awaiting approval**. No code written. Mapping of every export: `plans/2026-11-13_api_mapping.csv`.
Decisions already answered (user, 2026-11-13): three verbs per step (run, check, report); un-export now (no deprecation period).

## Why

`NAMESPACE` exports 241 `sw_` functions. A user cannot tell which to call first. The package mixes
steps (run a thing) with parts (a score, a plot, a reader) at the same level. Users should see one
verb per action per step; the parts stay as internal functions that the verbs call.

## Target surface

| Step | User-facing functions | Replaces (now internal) |
|---|---|---|
| Input | `sw_sample_manifest`, `sw_io_read_fcs` | `sw_manifest_from_sample_meta`; ten other `sw_io_*`, `sw_channel_*`, `sw_transform_*` |
| Import | `sw_import_check`, `sw_import_report` (new; `sw_report_import` is renamed) | `sw_manifest_bind_files`, `sw_fcs_*` audits, `sw_cytotools_header`, `sw_io_event_audit` |
| Unmix | `sw_unmix_run`, `sw_unmix_check`, `sw_unmix_report`, `sw_unmix_app` | setup, prepare, extract_*, spectracle, correct_*; checks, spike, screens; figures, dashboard |
| Clean and signal QC | `sw_clean_run`, `sw_qc_check`, `sw_clean_report`, `sw_qc_app` | filters, `sw_qc_run/batch/summary`; signal audit, anomaly, negative audit, signal reference, cytoscan; QC plots |
| Gate | `sw_gate_run`, `sw_gate_check`, `sw_gate_report`, `sw_gate_read`, `sw_gate_plot` | density cuts, polygons, transfer, flowMagic, FAUST, GatingML/FlowJo converters |
| Batch | `sw_batch_check`, `sw_batch_correct`, `sw_batch_decide`, `sw_batch_report` | scores, OTD, mixing, flags, reference checks, synthetic references; cyCombine and cyTether correctors; relativisation |
| Cluster | `sw_cluster_run`, `sw_cluster_check`, `sw_cluster_report`, `sw_dimred_run` | scores, stability, MFI, predict, plots |
| Annotate and differential | `sw_annotate_run`, `sw_diff_run`, `sw_diff_report` | manual annotation, plots, abundance and expression tests |
| Across steps | `sw_pipeline_run_all`, `sw_diagnose`, `sw_pipeline`, `sw_pipeline_run`, `sw_assistant` | 28 composable-pipeline methods and step constructors, `sw_diagnosis_top`, five assistant helpers |

Count: 19 names that exist today stay as they are; 16 are new verbs that wrap internals
(`sw_import_check`, `sw_import_report`, `sw_unmix_check`, `sw_clean_report`, `sw_qc_check`,
`sw_gate_check`, `sw_gate_report`, `sw_gate_read`, `sw_batch_check`, `sw_batch_correct`,
`sw_batch_report`, `sw_cluster_check`, `sw_cluster_report`, `sw_diff_run`, `sw_diff_report`, and the
rename `sw_report_import`). About 35 in all. Some existing names change meaning (see risks).

## How a verb is built

- `*_check`: the pattern `sw_diagnose()` already uses. One call takes the data plus one settings list
  per module (`NULL` skips a module). Every impactful value stays a required argument inside its list.
  The result is one long table with a `module` column; no pass or fail; no per-file roll-up.
  Example: `sw_import_check(paths, bind = list(...), integrity = list(...), signatures = list(...), compensation = list(...), batch_date = list(...), dye_map = list(...), header = list(...), signal_summary = list(...))`.
- `*_report(result, out_dir)`: new run folder, tables, figures, `source_data/`, `provenance.json`,
  `sessionInfo.txt` (parent plan contract). One per step.
- `*_run`: does the work (read, clean, gate, correct, cluster). Options that pick between methods
  (`method`, `engine`) are required arguments with no default.
- The internal functions keep their code, tests and roxygen. Only `@export` is removed, and
  `@keywords internal` is added, so `?sw_fcs_integrity_audit` still works for developers.

## What breaks and what does not

- Tests (84 files) use `devtools::load_all()`, which sees internal functions: no change needed.
- `experiments/` (88 files mention an affected name): 85 scripts use `load_all()`, none uses
  `library(SpectraWeaveR)`. They keep working. Scripts that become wrong are the ones calling a
  function whose meaning changes (below).
- Vignettes (11 files) call affected names. Each is rewritten to the new verbs.
- `R/` (87 files) call each other directly; no change.
- Anyone outside this repo who calls an un-exported name with `library()` breaks. This is accepted:
  no deprecation period (your answer).

## Risks and open design points

1. **`sw_unmix_run` changes meaning.** Today it is the low-level run; `sw_unmix_pipeline()` is the
   whole chain (setup, prepare, run). The new `sw_unmix_run` should be the whole chain, so the old
   `sw_unmix_run` becomes internal and `sw_unmix_pipeline` is retired. Callers: grep before the change.
2. **`sw_gate_run(method = )`** merges template gating, density cuts, transfer, flowMagic and FAUST.
   Each method has its own required values; the argument list must not become one long unreadable list.
   Proposed: `method` selects and a `settings` list holds that method's required values.
3. **`sw_batch_correct(engine = )`** covers cyCombine (`sw_correct_*`), the cyTether correctors and
   relativisation. The corrector contract (fit, integrate, diagnose, save, load) must survive as
   arguments (`model =`, `save =`). Batch work is frozen in ROADMAP A2; this is a renaming and
   wrapping of existing code, with no new method. Confirm that is allowed.
4. **Composable pipeline.** 28 step constructors and helpers go internal; `sw_pipeline_step(fun, ...)`
   stays the one way to add a step. The `targets` backend and `sw_plot_pipeline` stay reachable
   through `sw_pipeline_run(backend = )`. The assistant tool list in `R/assistant.R` names functions
   directly; each tool must call the internal function.
5. **Help and discovery.** `_quarto.yml` and `reference/*.qmd` list exports by step; they shrink to
   the 35. A `?SpectraWeaveR` page shows the step table above.
6. **Hidden power.** Advanced users lose one-call access to a single score (e.g. `sw_batch_otd`). They
   can still call it with `SpectraWeaveR:::`, or select it with a module setting in `sw_batch_check()`.
7. Work in progress: the composition-check workflow is still editing `R/composable_pipeline.R`,
   `R/assistant.R` and `R/fcs_header_audit.R`. Start only after it finishes.

## Build order (each step: tests first, then code, then docs)

1. Mapping review: you read `plans/2026-11-13_api_mapping.csv` and change any row.
2. Import step first (already built): `sw_import_check` and `sw_import_report`, un-export its parts.
   This is the template for the other steps.
3. Clean and QC, gate, batch, unmix, cluster, annotate and differential, one step at a time.
4. Composable pipeline and assistant last.
5. Rewrite the 11 vignettes, the reference pages, `CLAUDE.md`, `ROADMAP.md` register rows (argument
   names), and `NAMESPACE` via `devtools::document()`.
6. Check: whole suite, `R CMD check`-style example run, rerun the 8-file real-data composition script
   through `sw_import_check` and `sw_import_report`, compare its tables with the earlier run.

## Decisions for the user

| # | Choice | Recommended | Why |
|---|---|---|---|
| SA1 | The mapping CSV (which names stay, which verb absorbs the rest) | Review and edit | One row per export |
| SA2 | `sw_unmix_run` becomes the whole chain; `sw_unmix_pipeline` retires | Yes | One obvious call |
| SA3 | `method` / `engine` plus a `settings` list per verb | Yes | Keeps the required-argument rule readable |
| SA4 | Wrapping batch code under `sw_batch_correct` despite the ROADMAP freeze | Yes, rename and wrap only | No new method |
| SA5 | `sw_pipeline_step()` as the only step constructor | Yes | 11 constructors go internal |
| SA6 | Start after the composition check ends | Yes | Same files |
