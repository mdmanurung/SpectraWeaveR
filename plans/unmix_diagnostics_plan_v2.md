# Implementation plan v2: unmix diagnostics suite (spike-in, negative shape and screen, figures)

Status: **revised to the user's grill answers** (`plans/unmix_diagnostics_decisions.md`, 2026-11-11). Sources: v1 plan (`unmix_diagnostics_plan.md`, kept), Fable grill (`unmix_diagnostics_grill.md`), Opus critique (`unmix_diagnostics_critique.md`). No package code is written until the user approves this plan. Decision IDs in the last column of each step refer to the decisions file.

House rules applied: AGPL code wrapped, never copied (ROADMAP rule 1; AutoSpectral, spectracle and flowstate are untouched). Impactful choices are required arguments with register rows (rule 2). No figure ships without its approved plot-first row (rule 5). Unit tests run without flowCore or AutoSpectral (`experiments/README.md` rule 3). Unmixing is decoupled: no package function runs or submits unmixing (default row).

## Function signatures

Required arguments have no default and fail by name. Defaults appear only on display-only or structural arguments, marked `# display` or `# structural`.

```r
# Theme (R3.3, Q20) -- internal
.sw_theme_prism(base_size, base_family = "sans", base_fontface = "bold",  # display
                base_line_size = base_size / 14, base_rect_size = base_size / 14,
                axis_text_angle = 0, border = FALSE)
.sw_theme_nature(font_family = "sans")   # display; 7 pt, plain text, bold titles/strips/axis titles,
                                         # line 0.4, legend title kept, legend right, strip.clip "off"

# Module 1: spike-in (R2-R8, R2.8, R2.9)
sw_unmix_spike_build(
  raw,            # named character (sample = raw FCS path) or named list of matrices
  populations,    # list: populations[[sample]][[pop]] integer raw rows; pops must be disjoint (R4)
  spectra,        # fluorophores x detectors; spiked rows must have max 1 (R6)
  dyes,           # named character: marker = spectra rowname
  levels,         # numeric, must contain 0
  cofactor,       # k = cofactor * sinh(level)
  spike_fraction, # share of each population's cells that get a level > 0 (R3)
  seed,           # used only when assignment is NULL (R5)
  assignment = NULL,  # structural: tibble sample, raw_row, pop, level; replaces seed draw (R5)
  key_cols,       # alignment keys; must be disjoint from spectra columns; must be unique per file
  output_dir      # must not exist; NULL only with matrix input
)
# -> sw_spike_design: design_id, key (sample, raw_row, pop, marker, dye, level, k, raw_file,
#    spiked_file), files, realised_fraction (sample, pop, set, fraction), spectra_used, levels,
#    cofactor, spike_fraction, seed, key_cols, provenance. FCS keywords $SW_SPIKE_DESIGN,
#    $SW_SPIKE_DYE, $SW_SPIKE_LEVELS written (R2.8); $PnR raised where the spiked max exceeds it.

sw_unmix_spike_cospike(design, dye2, levels2, level1, output_dir)
# 88d construction (R8): rebuild from raw, spike only level1 cells with dye 1 + dye 2;
# sets "cospike_<dye2>_<level2>"; level1 required (R2.9); levels2 may contain 0 (fraction check).

sw_unmix_spike_score(
  design,
  unmixed,        # tibble: setup, sample, set, path (or matrix list column)
  loss_cut,       # recovery below this = lost (R7)
  baseline_tol,   # arcsinh shift above which an unspiked cell counts as moved (R2)
  results = NULL  # structural: named list setup -> sw_unmix_result; enables cosine (R6) and notes
)
# -> sw_spike_score: cells; summary per (setup, set, marker, pop, sample, level) with n, share_lost,
#    recovery quantiles q05 q25 q50 q75 q95 (R7), file_status; by_setup donor medians and ranges;
#    baseline_check (setup, sample, set, group = level0 | untouched, channel, n, max_abs_diff,
#    share_moved) and file_status (independent | coupled) (R2); spike_cosine (R6); provenance.

# Helpers (R2.2, R2.3, R2.5)
sw_unmix_spike_rows(gs, nodes, cleaned_files, raw_files, key_cols)  # duplicate-key stop
sw_unmix_spike_pairs(design, setups, file_suffix = " AutoSpectral")  # structural default
sw_unmix_spike_jobs(design, setups)  # one row per (setup, set): raw_dir, out_dir

# Module 2: negative shape and screen (R9-R11)
sw_unmix_negative_shape(
  x,              # x[[setup]][[sample]]: matrix or FCS path; rows paired across setups
  reference,      # setup name
  populations,    # populations[[sample]][[pop]] integer rows on the paired rows
  channels,
  cofactor, cut, squeeze_frac,
  min_negatives,  # floor on below-cut events, applied to setup and reference (R10)
  key_cols
)
# -> sw_negative_shape (S3 list): metrics tibble setup, sample, pop, channel, n_pop, n_neg, rsd_neg,
#    ref_rsd, rsd_ratio, squeezed_share, pos_share; negatives = each setup's own below-cut events (R9)

sw_unmix_marker_screen(shape, min_delta, min_donors, channel_populations = NULL)
# flag = squeezed share rises (> 0) in every donor AND donor-median rise >= min_delta;
# NA when fewer than min_donors (R11). -> sw_marker_screen: flags, rule, reference.

# Figures (R12, R3.4, Q20)
sw_unmix_figure_scatter_grid(            # panel F6
  events,         # tibble setup, sample, <x>, <y>, optional lost
  x, y,
  mode,           # "scatter" (% above cut) or "spike" (before/after, % lost, reference curves)
  cut,            # required in scatter mode
  cofactor,
  setup_fill,     # named colours, one per setup; unknown setup = error
  stem,           # file stem, e.g. "scatter-grid_TNF_CD8T"
  bins = 50, limits = NULL, font_family = "sans"  # display
)
sw_unmix_figure_spike_summary(score, readout_level, stem, font_family = "sans")  # F7
sw_unmix_figure_spike_ecdf(score, readout_level, stem, font_family = "sans")     # F8

# Changed signatures
sw_unmix_diagnose(result, limits = list(), unmixed_files = NULL, unmixed_data = NULL,
                  spike = NULL, screen = NULL, ...)                            # Q15, R2.4, Q16
sw_unmix_save_figures(figures, dir, sizes = NULL, check = c("none", "nature")) # Q20
```

## Steps

Hours are one developer, tests included, review excluded.

| # | Step | Files | Tests | h | Decisions |
|---|---|---|---|---|---|
| 1 | ROADMAP item B11 and register rows (table below); note in the gold-standard section that AutoSpectral has `sim.flow.data()` but no real-cell spike, so module 1 is not a rebuild | `ROADMAP.md` | none | 1.5 | R1, all required-arg rows |
| 2 | Standalone prism theme from the user's code. Header: "Theme code supplied by the package author (Mikhael D. Manurung), who states they hold its copyright, included under MIT (user decision 2026-11-11). The code matches `theme_prism()` in ggprism (github.com/csdaw/ggprism)." Hard-code `black_and_white` (pageBackgroundColor, plottingAreaColor `#FFFFFF`; axisColor, axisLabelColor, axisTitleColor, graphTitleColor `#000000`); drop `palette`, `prism.ticks.length`, `legend.text.align`, `legend.title.align`; base-R replacements for `deframe()` and `is_bool()`; parent `theme_grey(base_size, base_family) %+replace% t` instead of ggprism's `all_null` | `R/theme_prism.R` | `complete = TRUE`; every element in `t` set; axis angle error; border switch; parity with `ggprism::theme_prism()` on the elements `t` sets (skip if ggprism absent) | 5 | R3.3 |
| 3 | `.sw_theme_nature()` on top of step 2 (r-0017 values) | `R/theme_prism.R` | text sizes 5 to 7 pt; legend right | 1 | Q20 |
| 4 | Migrate F1 to F4: `.sw_fig_theme()` returns `.sw_theme_nature()`; facets `axes = "all", axis.labels = "margins"`; legend right; widths snapped to 89 or 183 mm in `.sw_fig_sizes`. Render before and after PNGs of each panel from the test fixtures into `experiments/runs/<date>_fig-migration/` for user approval; no merge before approval | `R/unmix_figures.R`, `tests/testthat/test-unmix_figures.R` | existing figure tests updated to new sizes | 8 | Q20 |
| 5 | Nature save check: `check = "nature"` on `sw_unmix_save_figures()`; internal width check, grob-walk pt check (cairo device), font embedding via an internal `pdffonts` wrapper (skipped with a note when absent) | `R/unmix_figures.R` | width 120 refused; 9 pt label refused; wrapper mocked | 4 | Q20 |
| 6 | Spike build: required-arg errors; disjoint populations (error on overlap, R4); per-population level draw from `spike_fraction` and `seed`, or `assignment` as given (R5); key uniqueness and `key_cols` disjoint from spectra; spectrum max 1 (R6); write FCS with keywords (R2.8) and raised `$PnR`; optional zero-spike round-tripped base set; `spike_key.csv`, `design.rds` | `R/unmix_spike.R` | below, step 9 | 10 | R2-R6, R2.8 |
| 7 | Co-spike as 88d (R8, R2.9) | `R/unmix_spike.R` | below | 4 | R8, R2.9 |
| 8 | Score: identity checks on `key_cols`; recovery per cell; `summary` with quantiles (R7); baseline check over every fluorescence channel in arcsinh, level-0 and untouched rows apart, `share_moved` above `baseline_tol`, `file_status` carried onto summaries, recovery never set to NA (R2); spike cosine when `results` given (R6); note naming `refine.af.quantile` when a result records `n.af.passes > 1` | `R/unmix_spike.R` | below | 11 | R2, R6, R7 |
| 9 | Module 1 tests: OLS mock unmixer (`qr.solve` on overlapping spectra) gives recovery 1; **coupled mock unmixer** (subtracts a file-level quantile of one channel from every cell) gives `file_status = "coupled"` and `share_moved > 0`; dropped spike gives `share_lost` 1; quantiles on a known recovery vector; per-pop level counts exact; same seed same key; `assignment` path bypasses seed; overlap error; duplicate-key error; spectrum max error; cospike recovers both dyes and keeps levels 2 and 3 unspiked; keywords survive `write.FCS` (flowCore skip); FCS tolerances not `== 0` | `tests/testthat/test-unmix_spike.R` | as listed | 8 | R2-R8, R2.8 |
| 10 | Helpers `sw_unmix_spike_rows()` (flowWorkspace skip), `sw_unmix_spike_pairs()`, `sw_unmix_spike_jobs()` | `R/unmix_spike_helpers.R` | mapping on a mock GatingSet-like list; pairs names; jobs rows | 6 | R2.2, defaults R2.3, R2.5 |
| 11 | Negative shape: pairing check; own below-cut events (R9); Nguyen rSD; `min_negatives` floor on setup and reference, `NA` below it; guard `ref_rsd = 0` | `R/unmix_screen.R` | below | 5 | R9, R10 |
| 12 | Marker screen: strict rise in all donors, median >= `min_delta`, NA under `min_donors` | `R/unmix_screen.R` | below | 3 | R11 |
| 13 | Module 2 tests: rSD of a known sample; reference ratio 1; squeezed share of a point mass 1; NA under `min_negatives` (setup and reference); flag cases: rise in 5/5 above delta TRUE, 4/5 FALSE, fall in 5/5 FALSE, below delta FALSE, too few donors NA | `tests/testthat/test-unmix_screen.R` | as listed | 5 | R9-R11 |
| 14 | F6 builder, both modes: pure ggplot2 `geom_hex`; facets `axes = "all"`; ggh4x `facet_grid2` + `strip_themed` fills when installed, plain strips otherwise; `stem` per marker x pop; source CSV of % per (setup, sample) | `R/unmix_figures_scatter.R`, `R/unmix_figures.R` (`F6` in panel list) | below | 8 | R12 |
| 15 | F7 (% lost and recovery quantiles, donors as points, setups down y) and F8 (recovery ECDF per setup, every setup drawn, linetype per setup checked against levels) | `R/unmix_figures_spike.R` | below | 6 | R3.4 |
| 16 | Figure tests F6 to F8: ggplot returned; source % equals a hand count; unknown setup error; spike mode needs `lost`; F8 errors when a setup has no linetype (the 88b/88e bug); `check = "nature"` passes on all three | `tests/testthat/test-unmix_figures_scatter.R`, `test-unmix_figures_spike.R` | as listed | 5 | R12, R3.4 |
| 17 | Diagnose: `spike =`, `screen =`; check ids `spike_recovery`, `marker_screen`; `item = "setup|marker|pop"` from one internal delimiter constant; value donor median, reason n and range; `not_judged` unless a limit is passed; limit spec rebuilt from a named direction vector (fixes the positional `direction`, `unmix_report.R:62`); new limit names `max_share_lost`, `min_median_recovery`, `max_share_moved`, `max_delta_squeezed` | `R/unmix_report.R` | below | 6 | Q15, Q16, R2.4 |
| 18 | Report sections "Spike-in recovery" and "Marker screen" (tables only), CSVs `spike_summary.csv`, `marker_screen.csv`, `baseline_check.csv`, JSON entries, `.sw_unmix_report_files` extended | `R/unmix_report.R` | below | 4 | Q17, R3.1 |
| 19 | Diagnose and report tests: rows absent without arguments; item round trip; status from supplied limit; every limit name's direction; CSVs written and replaced on overwrite | `tests/testthat/test-unmix_report.R` | as listed | 4 | Q15-Q17 |
| 20 | Dashboard tab: data-model entry reading the three CSVs from a run folder; tab with tables and F7/F8 PNGs; ignored when absent | `R/unmix_dashboard_data.R`, `R/unmix_dashboard.R`, `inst/dashboard/dashboard.html` | tab present only with CSVs; existing dashboard tests pass | 7 | R3.5 |
| 21 | `DESCRIPTION`: `ggplot2 (>= 3.5.0)`; Suggests `hexbin`, `ggh4x`; no ggprism; `inst/scripts/install_dependencies.R`; `devtools::document()` | metadata | `R CMD check` notes | 1 | R12, R3.3 |
| 22 | Real-data regression in `experiments/integration/`: feed saved spike keys through `assignment` for 88 (T cells), 88b (panel 2), 88e (CD3-), 88d (co-spike) and 88h (passes); compare `summary` to each run's `spike_summary.csv` and baseline to its CSV; one 88h `p1_af5` case must be `coupled` and `p5_af1` `independent`; screen vs 86b `flags.csv` only where 86b did not NA-coerce | `experiments/integration/spike_bmv.R`, `screen_bmv.R` | expected differences listed below | 8 | R5, R2, R11 |
| 23 | Docs: vignette (mock section always run; BMV section gated on `SW_DATA_DIR`), reference page, roxygen examples, `CLAUDE.md` pipeline row, `PROGRESS.md` | `vignettes/unmix-diagnostics.qmd`, `reference/`, `CLAUDE.md`, `PROGRESS.md` | vignette renders | 5 | defaults R3.2 |
| 24 | `devtools::test()` and `devtools::check()` in pixi; fix findings | all | green | 3 | all |

**Total: about 128 h.** The v1 estimate (86 h) grows by the F1 to F4 migration and Nature check (Q20, +13 h), the dashboard tab (+7 h), F7 and F8 (+6 h), the independence and quantile redesign (+5 h), the theme port (+6 h) and the wider regression (+2 h); the critique's deferrals do not apply under R1 = full plan.

### Expected regression differences (step 22)

| Difference | Cause | Decision |
|---|---|---|
| Screen flag `NA` where 86b wrote `FALSE` | 86b `coalesce(flag, FALSE)` (`86b:121`) | R11 |
| Screen `NA` where below-cut count < `min_negatives` | helper75 had no floor | R10 |
| Per-population level counts differ from 88 when run through `seed` | 88 drew pooled across populations (`88:128-147`) | R4; regression uses `assignment`, so no difference there |
| Baseline values differ from 88's `max_abs_diff` | new metric: arcsinh, all channels, share above `baseline_tol` | R2 |
| Summary gains quantile columns and `file_status` | new readouts | R2, R7 |

## Register rows to add to ROADMAP

Source for every row: "user, 2026-11-11, grill <ID>". No defaults.

| Argument | Function | Grill ID | Benchmark value (for reference only) |
|---|---|---|---|
| `levels` | spike_build | R8 default row / Q8 | 0, 1, 2, 3 |
| `cofactor` | spike_build, negative_shape, figures | Q8 | 6000 |
| `spike_fraction` | spike_build | R3 | 0.75 implied by levels 0-3 |
| `seed` | spike_build | R5 | 42 |
| `populations` | spike_build, negative_shape | R4 | CD4T, CD8T, CD3- |
| `loss_cut` | spike_score | R7 | 0.5 |
| `baseline_tol` | spike_score | R2 | none chosen yet |
| `level1` | spike_cospike | R2.9 | 1 |
| `levels2`, `dye2` | spike_cospike | R8 | 0, 1, 2; GLUT1 |
| `cut` | negative_shape, F6 scatter | R9 | arcsinh 1 |
| `squeeze_frac` | negative_shape | R9 | from 51 |
| `min_negatives` | negative_shape | R10 | none (new) |
| `min_delta` | marker_screen | R11 | 0.10 |
| `min_donors` | marker_screen | R11 | 5 |
| `reference` | negative_shape | R9 | as188_joint |
| verdict limits | diagnose | Q16 | none: `not_judged` |

ROADMAP item text:

| # | Item | Function(s) | Status |
|---|---|---|---|
| B11 | Spike-in recovery with an independence check, negative shape and marker screen, scatter and spike figures (F6 to F8), Nature theme for all unmix figures, as diagnostics on AutoSpectral outputs (from benchmark 05). Unmixing stays with the user. Limits unsourced: `not_judged`. | `sw_unmix_spike_*()`, `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()`, `sw_unmix_figure_scatter_grid()`, `sw_unmix_figure_spike_*()` | **next**, after plan approval |

## Plot-first approval rows

Each needs the user's approval of a rendered example before merge (ROADMAP rule 5).

| Panel | What it shows | Source table | Example to render |
|---|---|---|---|
| F1 migrated | spectra overlay, Nature theme | fixture spectra | before/after PNG |
| F2 migrated | similarity, Nature theme | fixture | before/after PNG |
| F3 migrated | distributions, Nature theme | fixture events | before/after PNG |
| F4 migrated | biplot, Nature theme | fixture events | before/after PNG |
| F6 scatter mode | marker vs culprit dye per setup x donor, % above cut | events + source CSV | BMV per-dye TNF vs GLUT1 |
| F6 spike mode | before vs after per cell, % lost, reference curves | `spike_cells` | BMV Th2 CD8T |
| F7 | % lost and recovery quantiles, donors as points | `summary` | BMV Th2 |
| F8 | recovery ECDF per setup | `cells` | BMV panel 2 |

## Risks

| Risk | Where | Mitigation |
|---|---|---|
| Coupled unmix methods mix spike response with file drift | score | required independence check and `file_status` (R2); coupled-mock test |
| Drift may grow with spiked fraction; fraction differs between single, co-spike and 88e designs | build, cospike | `spike_fraction` required, realised fraction recorded (R3); level2 = 0 fraction check (R8) |
| AutoSpectral may clip values above `$PnR` (spike adds up to ~60,000 at level 3) | build | raise `$PnR`; FCS round-trip test (reader behaviour not verified) |
| Raw and unmixed rows stop pairing | score, shape | `key_cols` identity and uniqueness stops |
| FCS keywords may not survive AutoSpectral unmixing | build | test survival through `write.FCS`; note result for AutoSpectral (unverified) |
| Positional `direction` in the limit spec | diagnose | rebuild as a named vector; test every name |
| `item` encoding parsed by report and dashboard | diagnose | one delimiter constant; round-trip test |
| Theme port differs from ggprism where `t` leaves elements unset (`theme_grey` parent instead of `all_null`) | theme | parity test on the elements `t` sets |
| Copyright of the theme code rests on the user's statement; the code matches ggprism (GPL >= 3) | theme | header records the statement and the match, as decided (R3.3); not verified here |
| F1 to F4 sizes and legends change approved figures | figures | before/after approval rows; existing tests updated |
| ggh4x and ggplot2 facet argument names differ | F6 | one internal facet helper; gtable axis test per backend |
| ggplot2 < 3.5.0 lacks `axes =` | figures | version floor |
| `pdffonts` or Nimbus Sans absent | save, theme | font check skipped with a note; `font_family` display arg, default "sans" |
| Spiked FCS sets are large | build | `output_dir` must not exist; sizes recorded; vignette notes disk use |
| `experiments/runs/` is git-ignored | regression | expected values in the script header |
| ggcyto stays in Suggests but fails under ggplot2 4.0.3 | DESCRIPTION | out of scope; ROADMAP note |

## Suggested PR chunks

| PR | Steps | Hours | Review focus |
|---|---|---|---|
| 1. Theme and figure base | 1-5, 21 | 20.5 | theme port, F1-F4 before/after approval, Nature check |
| 2. Spike-in module | 6-10 | 39 | method, independence check, mocks |
| 3. Screen and new figures | 11-16 | 32 | screen rule, F6-F8 approval |
| 4. Integration and docs | 17-20, 22-24 | 37 | diagnose/report schema, dashboard, regression, check |
