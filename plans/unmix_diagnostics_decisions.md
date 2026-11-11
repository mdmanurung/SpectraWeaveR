# Unmix diagnostics suite: grill decisions (user, 2026-11-11)

Sources: `unmix_diagnostics_grill.md` (Fable draft), `unmix_diagnostics_critique.md` (Opus critique). Every row below was answered by the user in popups on 2026-11-11, except the "default" rows, which were stated as settled unless the user objected (no objection).

| ID | Decision | Note |
|---|---|---|
| R1 | Full plan: build, cospike, score, shape, screen, figures, plus diagnose and report integration and helpers | departs from recommendation (standalone) |
| R2 | Independence check required in score: every fluorescence channel, arcsinh, level-0 rows and untouched rows reported apart, required `baseline_tol`, per-file status `independent`/`coupled` carried onto summaries; recovery values kept, never NA'd | 88h drift: `experiments/runs/2026-11-11_bmv-spike-passes/drift_summary.csv` |
| R3 | Required `spike_fraction`; realised fraction recorded per file | |
| R4 | Populations must be disjoint (error on overlap); levels drawn per population | |
| R5 | Build accepts a precomputed `assignment` (sample, raw_row, pop, level) beside `seed` + `spike_fraction`; regression feeds 88's spike keys | |
| R6 | Error unless each spike-spectrum row has max 1; store cosine between spike row and each setup's row | |
| R7 | Required `loss_cut`; add recovery quantiles 5/25/50/75/95% per donor | |
| R8 | Co-spike as 88d: rebuild from raw, spike only `level1` cells with dye 1 + dye 2, sets `cospike_<dye2>_<level2>`, level2 = 0 reported as the fraction check | |
| R9 | Negatives = each setup's own below-cut events (helper75) | |
| R10 | Required `min_negatives` on the below-cut count, for setup and reference; no default | renamed from `min_events` (collides with existing limit) |
| R11 | Flag = squeezed share rises (> 0) in every donor and donor-median rise >= `min_delta`; required `min_delta`, `min_donors`; NA when too few donors | differs from 86b NA->FALSE; list expected differences |
| R12 | Pure ggplot2 hexbin (no ggcyto); one builder, modes "scatter" and "spike"; file stem per marker x population; ggh4x strip fills only when installed | |
| Q15 | `sw_unmix_diagnose(spike = , screen = )`; check ids `spike_recovery`, `marker_screen` | |
| R2.4 | One row per setup x marker x pop; `item = "setup|marker|pop"`; value = donor median; reason = n donors and range | |
| Q16 | No built-in limits: `not_judged` unless the caller passes a limit | |
| Q17 | `sw_unmix_report()` sections "Spike-in recovery" and "Marker screen"; CSVs spike_summary, marker_screen, baseline_check | |
| Q20 | Nature theme for new figures AND migrate F1-F4 in this change; before/after figures for user approval per panel; opt-in `check = "nature"` | departs from recommendation |
| R3.3 | Theme = the user's own prism theme code (pasted 2026-11-11), included under MIT; user states they hold the copyright. Header records this and that the code matches `theme_prism()` in ggprism (github.com/csdaw/ggprism). Made standalone: hard-code the `black_and_white` palette colours, drop `prism.ticks.length`, base-R replacements for `deframe()`, `is_bool()`; no ggprism dependency | user decision after licence discussion |
| R2.2 | Ship `sw_unmix_spike_rows()` (GatingSet nodes -> raw rows; duplicate-key stop) | |
| R2.8 | Write `$SW_SPIKE_DESIGN`, `$SW_SPIKE_DYE`, `$SW_SPIKE_LEVELS` FCS keywords; test survival through write.FCS | |
| R2.9 | Co-spike `level1` required | |
| R3.4 | Ship spike summary panels F7 (% lost + recovery quantiles, donors as points) and F8 (recovery ECDF), with source CSVs; plot-first approval rows | |
| R3.5 | Add spike/screen tab to `sw_unmix_dashboard()` via `R/unmix_dashboard_data.R` | |
| R3.6 | Also ship a Quarto HTML dashboard generator `sw_unmix_dashboard_quarto()` (user, 2026-11-11: "i want the quarto HTML dashboard as well"): the 96_dashboard.qmd pattern packaged as a template in `inst/` — verdict cards, per-section Figures / Data (reactable) / AutoSpectral-notes tabs, theme: none with own CSS (dark mode), embed-resources, 300-dpi images from the panels' PDFs (user, 2026-11-11: "figures should be 300 dpi"); quarto and reactable in Suggests. Added to PR 4 | extra to plan v2 |
| default | Decoupled unmix (package never runs or submits unmixing); `sw_unmix_*` names; key-uniqueness stop; required output dir; mocks in tests (real OLS mock + coupled mock unmixer); roxygen examples; vignette with mock section always run and real-data section gated on `SW_DATA_DIR`; R2.1 no callback; R2.3 pairing helper default suffix " AutoSpectral"; R2.5 job-table helper shipped; R3.1 table-only report sections; R3.2 vignette gating | |
