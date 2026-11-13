# Plan: import and harmonise step (B1) and `sw_report_import` (revision 2)

Date: 2026-11-13. Status: **built 2026-11-13; arguments approved by the user the same day**.
Parent plan: `plans/2026-11-13_scenarios-and-step-reports.md`.
Sources: `experiments/repo_crawls/spectraweaver_incorporation_proposals.md` (ids M17, M21, M22, M30,
M38, M39, M40, N01, N02, N08, N09, N05); `R/fcs_header_audit.R`; two read-only audits of revision 1
(code claims; ROADMAP, decisions and source caveats).

## What the audit changed

Revision 1 said no manifest join exists, that the report copies `sw_unmix_report()`, that no row is
inherited, and that nine new functions are needed. All four were wrong or too broad.

- A manifest-to-file join exists: `.sw_ref_manifest()` (R/ref_vs_study.R:367-399) and the join in
  `sw_correct_prepare()` (R/batch_correct.R:93-175). `.sw_ref_manifest()` drops manifest rows with no
  file without saying so. New in I1: both orphan lists and `$NEXTDATA`.
- `sw_unmix_report()` (R/unmix_report.R:1742) writes CSVs and HTML; it stops on an existing folder only
  when `overwrite = FALSE`, and writes no `provenance.json`, `sessionInfo.txt`, figures or `source_data/`.
  The report contract is assembled from `sw_unmix_save_figures()` (R/unmix_figures.R:1000; PNG +
  `source_data/`, refuses to overwrite), the `sessionInfo.txt` lines in `sw_unmix_dashboard_quarto()`
  (R/unmix_dashboard_quarto.R:279-280) and the provenance helpers (R/unmix_report.R:217-345).
- IR11 (event scale) is inherited: arcsinh 6000 (CLAUDE.md design constraints; ROADMAP.md:214,
  confirmed 2026-11-09).
- `sw_fcs_header_audit()` already covers part of I3, I4 and I8 (needs a reference file; at least two
  paths; one unreadable file fails the call; no `not_assessed`). `sw_fcs_signal_audit()`
  (R/fcs_signal_audit.R:143) already gives `n_events` and `nonfinite_share`.

## Rules for every function

- Clean rewrite from the idea. Licences of the sources: M30, M38, M39, M17, M21 have no licence found;
  M40, N01, N02 and N08 sources are GPL-3. All are idea only; no code is copied (ROADMAP.md:25).
- Read headers or events only; never edit a file. Arguments are validated before dependency checks
  (CLAUDE.md). Required values use `.sw_stop_missing()` (R/corrector.R:142), so tests match
  `"no default"`.
- Event reads go through `sw_io_read_fcs()` (R/utils.R:60), which sets `truncate_max_range = FALSE`.
- No pass/fail and no per-file roll-up (Q14). One row per file and item. An unreadable or thin file is a
  `not_assessed` row with a reason. Nothing is dropped without a row.
- Every result carries a provenance block: package version, seeds, `sessionInfo()` summary (ROADMAP D6).

## Functions

| # | Function | What it does | Reuses |
|---|---|---|---|
| I1 | `sw_manifest_bind_files` (N01) | Exact-key join of a `sw_sample_manifest()` to FCS paths. Returns manifest rows with no file, files with no row, multi-dataset files (`$NEXTDATA`). Duplicate keys stop. | `.sw_ref_manifest()` logic, `sw_sample_manifest()` |
| I2 | `sw_fcs_integrity_audit` (M22) | Per file: readable, `$TOT` vs rows read, time span, channel set. Never edits `$TOT`. Reads every event unless `header_only = TRUE`, which skips the `$TOT` check and is stated in the output. Time-span part rests on a source file not opened in the survey. | `sw_io_read_fcs()` |
| I3 | `sw_fcs_header_signatures` (M40, M21) | Groups files by identical values of the named keywords and channel sets, with counts. No reference file needed. Optional `panel` list (one or more reference marker sets) gives per-file markers missing and extra (this is I8, folded in). Unreadable file is a `not_assessed` row. | `.sw_hdr_items()` (R/fcs_header_audit.R:159) |
| I4 | `sw_fcs_compensation_state` (M30) | Per file, separate columns: matrix present, matrix is identity, apply flag value. No precedence rule. | `.sw_hdr_items()` |
| I5 | `sw_fcs_batch_date_check` (N02) | Declared batch (and manifest `run`) crossed with date, begin-time and serial keywords; table for a tile map. | `sw_sample_manifest()` |
| I6 | `sw_fcs_marker_dye_map` (M38) | Marker x batch table of dyes from `$PnN`/`$PnS`; names that do not split are listed, not guessed. | `.sw_hdr_items()` |
| I9 | extend `sw_fcs_signal_audit` (N08) | Add per-channel mean, quantiles and counts; `min_events` floor gives `not_assessed`. Changes an existing exported function and its tests (approved). | existing function |

Deferred: I7 `sw_marker_names` (M39) waits for ROADMAP D5 so the panel shape is fixed once; when it is
built it must log markers that leave the shared set and from which files (ROADMAP.md:158 bug).
`sw_fcs_harmonize` (M17, renames), `sw_bead_peak_audit` (N05), `sw_plate_layout_read` (N09) stay later.

## Report: `sw_report_import(result, out_dir)`

New run folder; stops if it exists (no `overwrite`). Writes `tables/*.csv`, `figures/*.png`,
`source_data/<panel>.csv`, `provenance.json`, `sessionInfo.txt`. Returns
`list(tables, figures, source_data, dir)`. Numbers only.

| Output | Content | Figure |
|---|---|---|
| `bound_files.csv` | I1 | `table fallback: names and counts, nothing to plot` |
| `header_audit.csv` | `sw_fcs_header_audit()` output (differences against a named reference) | file x keyword tile map |
| `integrity.csv` | I2 | one dot per file for `$TOT` minus rows read, ordered by value, y in events |
| `signatures.csv` | I3 incl. panel match | file x keyword tile map grouped by signature; file x marker presence tile map |
| `compensation_state.csv` | I4 | one point per file in a state x file grid (not a bare count) |
| `batch_date.csv` | I5 | batch x date tile map, fill = file count |
| `dye_map.csv` | I6 | marker x batch tile map |
| `channel_summary.csv` | I9 | per-file median and IQR dots, one panel per channel, free y axes with units, files in acquisition or batch order |

Figures are proposals. Per the plot-first rule, no conclusion from a real run rests on a figure until
it is saved and inspected. Figure panels for cytometry events use ggcyto hexbins where events are drawn
(memory: gate-figures-ggcyto); the tile maps and dot panels here show tables, not events.

## Impactful values (plan rows; each a required argument, no default)

Status column: **new** = no earlier decision; **inherited** = earlier decision, one-line confirm.

| Row | Function | Value | Status |
|---|---|---|---|
| IR1 | I1 | match key; case rule; recursion; duplicate-file and extra-file policy; `$NEXTDATA` policy | new |
| IR2 | I2 | `tot_tolerance` (source used exact equality) | new |
| IR3 | I2 | `time_keyword` (which keyword or channel gives the span) | new |
| IR4 | I3, I4, I5 | `keywords` / `matrix_keywords` / `apply_keywords` / date keywords | new (same idea as `sw_fcs_header_audit(keywords=)`) |
| IR5 | I4 | `zero_tol` (distance from identity that counts as identity; source used exactly 0) | new |
| IR6 | I4 | `expected_state` per file or group | new |
| IR8 | I5 | `batch` labels; date format; time zone; `run_gap_hours` | new |
| IR10 | I3 | `panel` marker sets | new; shape revisited at D5 |
| IR11 | I9 | `probs`; `min_events`; median vs mean; event scale and cofactor | `cofactor` **inherited**: arcsinh 6000 (ROADMAP.md:214), confirmed here; rest new |
| IR12 | I2, I9 | `n_events` and `seed` for I9 only; I2 reads every event | new, per function |
| IR13 | I2 | `header_only` (TRUE skips the `$TOT` check) | new |

IR7 (precedence) is gone: I4 reports the facts in separate columns. IR9 (alias map) moves to the
deferred I7. Rows IR1 to IR13 are added to the ROADMAP choice register (ROADMAP.md:237 asks for a row
and a one-line confirm when a later step inherits a choice).

## Tests (written first)

`tests/testthat/`, synthetic FCS written with `flowCore::write.FCS()` into `withr::local_tempdir()` as in
`test-fcs_header_audit.R:5-27`. For each function: missing required argument stops with `"no default"`;
an unreadable file gives a `not_assessed` row, not a crash; every dropped or unmatched file appears in
the output (no silent filter). Specific cases: duplicate keys stop (I1); orphan files and manifest rows
listed both ways (I1); `$TOT` mismatch found and `$TOT` never edited (I2); `header_only = TRUE` says so
(I2); identity vs non-identity matrix and apply flag in separate columns (I4); new `sw_fcs_signal_audit`
columns and unchanged old columns (I9); report refuses an existing folder.

## Vignette

`vignettes/import-and-harmonise.qmd`: synthetic FCS files in a temp dir (varied `$DATE`, `$TOT`, matrix
and marker spellings); runs in CI; ends with `sw_report_import()`. Example values labelled
"illustrative, choose your own". The real public set (FR-FCM-Z4KT) goes in a section with
`eval: false`; its availability is unverified and FlowRepository HTTPS has had an expired certificate
(HANDOFF-diagnostics.md).

## Decisions answered (user, 2026-11-13)

| Row | Answer |
|---|---|
| Q1 panel | Fold I8 into I3; defer I7 to D5. |
| Q2 I9 | Extend `sw_fcs_signal_audit`. |
| Q3 inherited | Confirm arcsinh 6000; add IR rows to the ROADMAP register. |
| Q4 I2 full read | Full read, plus a required `header_only` argument (TRUE skips the `$TOT` check). |

## Estimate

I1 to I6 and I9, report, tests and vignette: about two working days. Build order: tests, then I1, I2,
I3, I4, I5, I6, I9, report, vignette.

## Arguments approved after the build (user, 2026-11-13: "fix functions and approve arguments")

Choices the builders made that were not plan rows. All approved as built; each is also in the ROADMAP register.

| Row | Function | Approved argument or shape |
|---|---|---|
| IA1 | `sw_fcs_compensation_state` | `apply_true` and `apply_false`: required flag-word vectors, no built-in vocabulary |
| IA2 | `sw_fcs_integrity_audit` | `time_keyword`: one channel name (needs events) or two header keys, begin and end (header only, seconds); an end before the begin is `not_assessed` |
| IA3 | `sw_fcs_compensation_state` | `expected_state`: NULL or a data frame, file plus logical `matrix_present`, `matrix_identity`, `apply_parsed`; NA means no expectation |
| IA4 | `sw_manifest_bind_files` | `key` (file_stem or full name), `case`, `recursive`, `on_duplicate_file` (stop or unbind), `on_extra_file` (list or stop), `nextdata`; a path that does not exist stops |
| IA5 | `sw_fcs_batch_date_check` | batch and run come from the manifest; no `batch` argument (IR8 changed) |
| IA6 | `sw_fcs_signal_audit` | `summary` has no `n_events`, `seed` or cofactor (IR11, IR12 changed): the function uses the events it is given |
| IA7 | `sw_fcs_integrity_audit` | The row count for a file whose `$TOT` disagrees with its data comes from flowCore's error text; if that wording changes the row becomes `not_assessed` with the error. Accepted as a known fragility |

Fixes made after the audit: `experiments/` callers of `sw_fcs_signal_audit()` pass `summary = NULL`
(16 calls); one shared header reader (`.sw_read_fcs_header()`, R/import_util.R) closes the connections
flowCore leaves open on a bad file; ROADMAP register rows match the built arguments.
