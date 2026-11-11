# Implementation plan: spike-in recovery, negative shape and marker screen, raw scatter grid

Status: **draft, conditional on the recommended answers in `plans/unmix_diagnostics_grill.md`**. Each step names the question ids it depends on. If an answer differs from the recommendation, the step changes as noted in its "if not" line. No package code is written until the user approves. Written 2026-11-11.

Conventions: AGPL code is wrapped, never copied (ROADMAP rule 1); impactful choices are required arguments (rule 2); no plot without its approved row (rule 5); unit tests pass without flowCore or AutoSpectral (`experiments/README.md` rule 3).

## Files to add or change

| File | Change | Depends on |
|---|---|---|
| `R/unmix_spike.R` | new: build, cospike, score, `print()` methods, internal helpers | Q1, Q2, Q26 |
| `R/unmix_screen.R` | new: negative shape, marker screen, `print()` | Q12, Q13, Q14 |
| `R/unmix_figures_scatter.R` | new: scatter-grid builder, Nature theme, events check | Q18, Q19, Q20 |
| `R/unmix_figures.R` | edit: add `F6 = "scatter-grid"` to `.sw_fig_panels` and `.sw_fig_sizes`; `check` argument on `sw_unmix_save_figures()` | Q19, Q20 |
| `R/unmix_report.R` | edit: limit-spec rows; `spike`, `screen` arguments on `sw_unmix_diagnose()`; two report sections and two CSVs; `.sw_unmix_report_files` | Q15, Q16, Q17 |
| `DESCRIPTION` | edit: `ggplot2 (>= 3.5.0)`; Suggests `hexbin`, `ggh4x`, `ggprism` | Q18, R2.6, Q20 |
| `NAMESPACE`, `man/` | regenerate with `devtools::document()` | all |
| `tests/testthat/test-unmix_spike.R`, `test-unmix_screen.R`, `test-unmix_figures_scatter.R` | new | Q24 |
| `tests/testthat/test-unmix_report.R` | edit: new limit names, new check rows, report sections | Q15, Q17 |
| `experiments/integration/spike_bmv.R`, `screen_bmv.R` | new: real-data regression against 88 and 86b outputs | Q24 |
| `vignettes/unmix-diagnostics.qmd`, `reference/unmix-diagnostics.qmd` | new | Q25 |
| `CLAUDE.md`, `ROADMAP.md`, `PROGRESS.md` | edit: pipeline-table row, item B11, register rows | Q25 |

## Function signatures

Required arguments have no default and error by name (`.sw_fig_require()` style). Optional arguments show their default.

```r
# Module 1 (Q1 A, Q4 A, Q5 A, Q6 A, Q8 A, Q9 A, Q22 A)
sw_unmix_spike_build(
  raw,          # named character (sample = raw FCS path) or named list of matrices (events x detectors)
  populations,  # list: populations[[sample]][[pop]] integer raw rows
  spectra,      # fluorophores x detectors; the spiked dye rows must have max 1 (Q9)
  dyes,         # named character: marker = rowname of spectra
  levels,       # numeric, must contain 0
  cofactor,     # k = cofactor * sinh(level)
  seed,
  allocation,   # "equal" (sample(rep_len(levels, n))) or list(n_per_level = <int>)
  output_dir,   # must not exist; NULL only with matrix input (returns spiked matrices)
  key_cols      # stored for the identity check in score
)
# -> sw_spike_design: list(design_id, key (tibble: sample, raw_row, pop, marker, dye, level, k,
#      raw_file, spiked_file), files (tibble: marker, sample, raw_file, spiked_file),
#      spectra_used, levels, cofactor, seed, allocation, key_cols, spiked (matrices, matrix
#      input only), provenance)

sw_unmix_spike_cospike(               # Q11 A, R2.9
  design,
  dye2,         # named character: marker = rowname of design$spectra_used
  levels2,      # numeric, second dye; 0 allowed (set with dye 1 only)
  level1,       # which level of dye 1 the co-spiked cells carry (88d: 1)
  output_dir
)
# -> sw_spike_design with cospike = list(dye2, levels2, level1, files)

sw_unmix_spike_score(                 # Q7 A, Q10 A+B, Q27 A
  design,
  unmixed,      # tibble: setup, sample, set ("base" or a marker, or a cospike set), path or matrix
  loss_cut,     # recovery below this = lost
  results = NULL  # named list setup -> sw_unmix_result, for provenance
)
# -> sw_spike_score: list(cells, summary, by_setup, baseline_check, loss_cut, design_id, provenance)

# Module 2 (Q12 A, Q13 A, Q14 A)
sw_unmix_negative_shape(
  x,            # list: x[[setup]][[sample]] matrix or FCS path; rows paired across setups
  reference,    # setup name
  populations,  # list: populations[[sample]][[pop]] integer rows on the paired rows
  channels,     # character
  cofactor, cut, squeeze_frac, min_events,
  key_cols
)
# -> tibble: setup, sample, pop, channel, n, rsd_neg, ref_rsd, rsd_ratio, squeezed_share, pos_share
#    (NA metrics where n < min_events), attr "settings"

sw_unmix_marker_screen(
  shape,        # output of sw_unmix_negative_shape()
  reference,
  min_delta,    # donor-median rise in squeezed_share
  min_donors,   # donors that must all rise; fewer -> NA
  channel_populations = NULL  # tibble channel, pop: designated parents
)
# -> sw_marker_screen: list(flags (tibble: setup, pop, channel, n_donors, median_d_squeezed,
#      min_d_squeezed, all_donors_up, median_rsd_ratio, flag), rule, reference)

# Module 3 (Q18 A, Q19 A, Q20 A, Q21 A, R2.6 A, R2.10)
sw_unmix_figure_scatter_grid(
  events,       # tibble: setup, sample, <x>, <y>, optional lost (logical)
  x, y,         # column names
  cut,          # arcsinh on y; drawn as a line; NA allowed only with pct_label = "lost"
  cofactor,
  setup_fill,   # named colours, one per setup level; unknown setup = error
  transform_cols,           # columns to arcsinh(x / cofactor) before drawing; character() if none
  pct_label = c("above_cut", "lost", "none"),
  reference_curves = NULL,  # tibble x, y, curve
  limits = NULL, bins = 50, label_pt = 6, font_family = NULL
)
# -> list(panel = "F6", plot, source (tibble: setup, sample, n, pct))

# Changed signatures
sw_unmix_diagnose(result, limits = list(), unmixed_files = NULL, unmixed_data = NULL,
                  spike = NULL, screen = NULL, ...)
sw_unmix_save_figures(figures, dir, sizes = NULL, check = c("none", "nature"))
```

## Steps

Hours are drafter estimates for one developer, tests included, excluding review.

| # | Step | Depends on | Hours | Deliverable |
|---|---|---|---|---|
| 1 | ROADMAP: add item B11 with the three modules; add the register rows from the grill | Q8, Q13, Q14, Q16 | 1 | `ROADMAP.md` |
| 2 | `R/unmix_spike.R` build: argument checks (names, integer rows inside range, `levels` contains 0, spectrum rows max 1, detector names match the raw columns), seeded equal split per sample, `k = cofactor * sinh(level)`, add `outer(k, spectrum)` to the listed detectors, write FCS (path input) or return matrices (matrix input), write `spike_key.csv` and `design.rds`, provenance via `.sw_unmix_provenance()` | Q4, Q5, Q8, Q9, Q22 | 8 | build + print |
| 3 | Cospike: extend a design on its `level1` cells, one set per `levels2`, same write path | Q11, R2.9 | 4 | cospike |
| 4 | Score: read pairs (paths via flowCore, matrices as given), identity check on `key_cols` for base vs spiked and vs the design, per-cell recovered and recovery, `baseline_check` on level 0, `summary` per donor, `by_setup` donor medians and ranges, provenance from `results` | Q6, Q7, Q10, Q27 | 8 | score + print |
| 5 | Unit tests for module 1 on matrices: identity unmix gives recovery 1 and `share_lost` 0; dropped spike gives 1; level 0 gives `max_abs_diff` 0; permuted rows stop with the key message; `levels` without 0 errors; spectrum max != 1 errors; equal split counts; cospike recovers both dyes; FCS round trip under `skip_if_not_installed("flowCore")` | step 2 to 4, Q24 | 6 | `test-unmix_spike.R` |
| 6 | `R/unmix_screen.R` negative shape: pairing check on `key_cols`, Nguyen rSD below the cut, reference rSD on the same rows, squeezed share, positive share, `NA` under `min_events` | Q12, Q13 | 5 | shape |
| 7 | Marker screen: per-donor delta vs reference joined `many-to-one` with `unmatched = "error"`, donor median, all-donors-up, `min_donors`, optional designated parents | Q14 | 3 | screen + print |
| 8 | Unit tests for module 2: synthetic 5-donor, 2-setup table where one channel is squeezed in every donor (flag TRUE), one in 4 of 5 (FALSE), one below `min_delta` (FALSE), one with `n < min_events` (NA); rSD on known quantiles; reference setup ratio is 1 | step 6, 7 | 4 | `test-unmix_screen.R` |
| 9 | Diagnose integration: limit-spec rows `max_share_lost -> share_lost`, `min_median_recovery -> median_recovery`, `max_baseline_diff -> baseline_max_abs_diff`, `max_delta_squeezed -> median_d_squeezed`; `spike` rows `check_id = "spike_recovery"`, `item = "<marker>|<pop>"`, donor median in `value`, n and range in `reason`; `screen` rows `check_id = "marker_screen"`; `not_judged` without limits; notes updated | Q15, Q16, R2.4 | 5 | `sw_unmix_diagnose()` |
| 10 | Report sections "Spike-in recovery" and "Marker screen" (tables only), `spike_summary.csv`, `marker_screen.csv`, `.sw_unmix_report_files` extended, JSON gets both tables | Q17, R3.1 | 4 | `sw_unmix_report()` |
| 11 | Tests for 9 and 10: parse of the new limit names, status by limit, rows absent when `spike = NULL`, report writes the two CSVs, `overwrite` replaces them | step 9, 10 | 3 | `test-unmix_report.R` |
| 12 | `R/unmix_figures_scatter.R`: events check (columns, one population per call, setups in `setup_fill`), arcsinh of `transform_cols`, `% above cut` or `% lost` per (setup, sample) on every event, `geom_hex`, cut line, label, `ggh4x::facet_grid2(axes = "all", remove_labels = "all")` with `strip_themed` fills (ggh4x 0.3.1 names; fallback `ggplot2::facet_grid(axes = "all", axis.labels = "margins")` with plain strips and a message), `aspect.ratio = 1`, Nature theme, reference curves; both modes of the approved row Q21b | Q18, Q19, Q20, Q21, Q21b, R2.6, R2.10 | 8 | builder |
| 13 | `sw_unmix_save_figures()`: register `F6`, size `c(183, 160)` default, `check = "nature"` (width in 89 or 183; text 5 to 7 pt via grob walk; font embedding when `pdffonts` is present) | Q19, Q20, R2.7 | 4 | save |
| 14 | Tests for 12 and 13: returns ggplot and source with correct %; unknown setup errors; two populations error; `check = "nature"` refuses width 120; font check skipped without `pdffonts` (mock `Sys.which`) | step 12, 13 | 4 | `test-unmix_figures_scatter.R` |
| 15 | `DESCRIPTION`: `ggplot2 (>= 3.5.0)`; Suggests `hexbin`, `ggh4x`, `ggprism`; `inst/scripts/install_dependencies.R` updated; `devtools::document()` | Q18, R2.6 | 1 | metadata |
| 16 | `experiments/integration/spike_bmv.R`: rebuild the 88 TNF design from the BMV raw files with the same seed and levels, score the existing `spk_c188_*` unmixes, assert `summary` equals `runs/2026-11-10_bmv-spike-in-v2/spike_summary.csv` for setup `c188` (tolerance 1e-6); `screen_bmv.R`: assert `flags` equals `runs/2026-11-10_bmv-marker-flags-v2/flags.csv` (source `screen`) | step 4, 7, Q24 | 6 | regression scripts, run once on SLURM |
| 17 | Helpers R2.2 (`sw_unmix_spike_rows()` from a GatingSet, tests skipped without flowWorkspace) and R2.3 (`sw_unmix_spike_pairs()`) | R2.2 = A, R2.3 = B | 5 | helpers |
| 18 | Vignette, reference page, `CLAUDE.md` row, `PROGRESS.md` entry | Q25 | 4 | docs |
| 19 | `devtools::test()` and `devtools::check()` in the pixi env; fix findings | all | 3 | green check |

Total: about 86 hours; steps 2 to 5, 6 to 8 and 12 to 14 are independent of each other and can run in parallel.

If not: Q3 = B adds a `sw_unmix_spike_unmix()` step (3 h). Q15 = B replaces step 9 with two standalone functions (same hours) and drops the report sections' automatic path. Q18 = B or C adds a ggcyto path that cannot be tested in the main env and a documented ggplot2 pin (6 h, plus the risk below). Q20 = B adds four plot-first rows and a restyle of F1 to F4 (8 h). Q5 = B moves the key mapping into build (needs flowWorkspace in build; tests skip). Q9 = B adds a recorded normalisation factor to the design (1 h).

## Test list

| Test file | What it checks | Needs |
|---|---|---|
| `test-unmix_spike.R` | required-argument errors by name; `levels` must contain 0; spectrum max != 1 error; detector mismatch error; equal split gives `n %/% length(levels)` per level within 1; the same seed gives the same key; identity unmix recovery 1; dropped spike `share_lost` 1; half-lost spike `share_lost` matches `loss_cut`; level-0 `max_abs_diff` 0; permuted base rows stop; cospike reports two recoveries; `output_dir` exists error; FCS round trip | none; flowCore for the round trip |
| `test-unmix_screen.R` | pairing stop; rSD of a known sample; `ref_rsd` equals `rsd_neg` in the reference; squeezed share of a point mass at 0 is 1; `NA` under `min_events`; flag rule four cases; `min_donors` larger than available gives NA; designated-parent filter | none |
| `test-unmix_report.R` (edit) | new limit names parse; unknown still errors; `spike` rows absent without the argument; status pass/warn/fail from a supplied limit; `item` encoding; report CSVs written and replaced on overwrite; JSON has the tables | none; jsonlite optional |
| `test-unmix_figures_scatter.R` | builder returns `panel = "F6"`, ggplot, source with % per (setup, sample) equal to a hand count; unknown setup errors; `pct_label = "lost"` needs `lost` column; `check = "nature"` width refusal; text-size refusal with a 9 pt label; font check skipped without `pdffonts` | ggplot2, hexbin; ggh4x optional |
| `experiments/integration/spike_bmv.R` | summary equals the 88 run for setup `c188`; `baseline_check` max 3e-8 (`INDEX.md:252`) | BMV data, flowCore, SLURM |
| `experiments/integration/screen_bmv.R` | flags equal 86b's `flags.csv` (`source == "screen"`) | BMV data, gated runs |

## Register rows to add to ROADMAP

See the grill's "Register rows" table; it is copied into `ROADMAP.md` under "Register of choices the user must set" in step 1. Plus one row for ROADMAP item B11:

| # | Item | Function(s) | Status |
|---|---|---|---|
| B11 | Spike-in recovery, negative shape and marker screen, raw scatter grid, as reusable diagnostics on AutoSpectral outputs (from benchmark 05). Verdict limits unsourced: `not_judged`. | `sw_unmix_spike_*()`, `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()`, `sw_unmix_figure_scatter_grid()` | **next** after grill approval |

Also add to the gold-standard section a line that AutoSpectral 1.8.8 has `sim.flow.data()` (synthetic cells with ground truth) and `get.spread.thresholds()` (spread-scaled positivity cuts); neither adds a known signal to real cells, so the spike-in is not a rebuild of an AutoSpectral routine (rule 3).

## Risks

| Risk | Where | Mitigation |
|---|---|---|
| AutoSpectral reorders or drops events in some configuration, so raw and unmixed rows no longer pair | score, shape | `key_cols` identity stop (Q6); document `include.raw = TRUE` as the stronger check; the pairing was verified on five BMV files only |
| Spiked FCS files are large (one raw copy per marker and per co-spike level) | build | `output_dir` must not exist; the design records sizes; document disk use in the vignette |
| `item = "<marker>|<pop>"` encoding is parsed by the report | diagnose | keep the delimiter in one internal constant; test the round trip |
| `.sw_unmix_limit_spec$direction` is positional (`rep("max", 7L), rep("min", 3L)`, `unmix_report.R:62`); a new limit row in the wrong block shifts every direction | diagnose | build the spec from a named vector `c(max_cosine = "max", ...)`; add a test asserting the direction of every limit name |
| ggh4x `facet_grid2()` and ggplot2 `facet_grid()` differ in argument names (`remove_labels` vs `axis.labels`) and may differ in which axes they draw | figure | one internal facet helper chooses the backend; a snapshot test of the built gtable's axis grobs per backend |
| ggcyto look differs from the pure ggplot2 figure (hex binning and colour ramp are the same calls; the gate box becomes a line) | figure | the restyle already made that change (`95:257-260`); the integration script redraws one BMV panel for visual comparison, cited in `PROGRESS.md` |
| `facet_grid(axes = )` needs ggplot2 >= 3.5.0; users on older ggplot2 fail at draw time | figure | version floor in `DESCRIPTION`; `.sw_fig_need_ggplot2()` checks the version |
| `ggprism` or `Nimbus Sans` absent on a user's machine | figure | fallbacks: `theme_classic(7)` and device default font; the Nature check reports the font used |
| `pdffonts` absent | save | the font check is skipped with a note in the returned table (R2.7) |
| A channel's name in the unmixed file is `<dye>-A` in AutoSpectral and `<dye>` elsewhere | score | reuse `.sw_match_channels()` (`unmix_report.R:1038`) |
| Donor count below `min_donors` after a failed unmix job | screen | `NA` flag with the reason, never FALSE (ROADMAP rule 4) |
| Regression numbers in `experiments/runs/` are git-ignored | integration | the integration script reads them from the run folder and states the expected values in its header; `INDEX.md` keeps the table |
| The schema of `sw_unmix_diagnostics$checks` is consumed by `sw_unmix_dashboard()` | dashboard | new `check_id` values are ignored by the dashboard reader; verify in step 11 by running `test-unmix_dashboard.R` |
