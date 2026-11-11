# Grill: spike-in recovery, negative shape and marker screen, raw scatter grid

Status: **draft, nothing decided**. Every question below is for the user. The recommendations are the drafter's picks with their confidence; none is a default. Written 2026-11-11.

Scope: three new SpectraWeaveR modules built on top of AutoSpectral outputs, reusing the hand-built checks of benchmark 05 (`experiments/benchmarks/05_bmv-unmix-compare/`).

| Module | Functions asked for | Benchmark source |
|---|---|---|
| 1 Spike-in recovery | `sw_unmix_spike_build()`, `sw_unmix_spike_score()`, `sw_unmix_spike_cospike()` | 62, 64, 88, 88b, 88c, 88d, 88e, 88h |
| 2 Negative shape and marker screen | `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()` | 51, 75, helper75, 86b |
| 3 Raw scatter grid | `sw_unmix_figure_scatter_grid()` | 72, 88c, 93, 94, 95_restyle_events_scatter, lib/theme.R |

How to answer: reply with the question id and the option letter, or a correction. Round 2 opens once its round-1 dependencies are answered; round 3 once round 2 is.

## Sources opened for this grill

| Source | Status | Used for |
|---|---|---|
| `CLAUDE.md`, `ROADMAP.md` | opened | rules 1 to 6, register format (`ROADMAP.md:167-190`), gold-standard rules (`ROADMAP.md:31-45`) |
| `R/unmix_report.R` | opened: 30-80, 360-780, 1160-1390, 1500-1791 | result, check rows, limit spec, diagnose, report |
| `R/unmix_figures.R` | opened: 1-135, 203-302, 930-1110 | panel keys, palette, theme, event contract, save |
| `R/unmix.R` | opened: 1040-1380 | `sw_unmix_run()` signature and output list |
| `R/unmix_dashboard.R`, `R/unmix_dashboard_data.R`, `R/unmix_app.R` | function list only (grep) | scope decision Q28 |
| `R/gate.R` | signatures only (grep, 54-66, 194-206) | population spec Q5 |
| `tests/testthat/test-unmix_report.R` 1-140, `test-unmix_figures.R` 1-60 | opened | mock style |
| `experiments/benchmarks/05_bmv-unmix-compare/README.md` | opened, all | approved choices, decision table |
| `experiments/runs/2026-11-10_bmv-summary/INDEX.md` | opened, all | results, AutoSpectral quotes (section 14) |
| Scripts 62, 88, 88c, helper75, 75, 86b, 72, lib/theme.R | opened, all | methods |
| Scripts 88d 1-130, 88e 1-110, 88h 1-120, 64 1-75, 51 1-110, 93 1-90, 94 1-60, 95_restyle_events_scatter 1-327 | opened, partial | methods |
| `pixi.toml`, `DESCRIPTION`, `experiments/README.md` 1-80 | opened | dependencies, rules |
| `10_ggcyto_strategy.R:9-10` | grep only | ggcyto failure note |
| ggplot2 4.0.3 `NEWS.md` (installed) | grep | `axes` argument arrived in 3.5.0 (`NEWS.md:531, 801`) |
| AutoSpectral 1.8.8 help (installed, `.pixi/Rlib`) | `unmix.fcs`, `sim.flow.data`, `get.spread.thresholds` pages | formals; `include.raw`; no spike-in routine exists |
| AutoSpectral articles 01, 10, 11, 19 | **not opened**; quoted via `INDEX.md:171-183` | verdict limits Q16 |
| Mage, Konecny, Mair 2026 (Cytometry A) | **not opened**; cited via `ROADMAP.md:186` | Q16 |
| Nguyen robust SD source | **not opened**; name only from `README.md:56` | Q13 |

Facts verified this session: installed AutoSpectral is 1.8.8, `RemoteSha 3c504c45...` (matches the pin). Main env ggplot2 is 4.0.3, ggcyto 1.34.0; the private lib `experiments/bmv_test/compare_3D/data/lib_ggplot2_35` holds ggplot2 3.5.2. `experiments/integration/` is empty. `ggh4x`, `ggprism`, `hexbin`, `ggpubr` are installed through `pixi.toml` but not in `DESCRIPTION`.

---

## Round 1: answerable now

Answer in the listed order. A later question is moot when its dependency picks the option that removes it (for example Q16 and Q17 fall away if Q15 = C). Each question's `Evidence:` line cites the files by path and line; their opened or not-opened status is in the "Sources opened" table above.

### Q1 API shape for module 1

Impact: high. Depends on: none.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Three functions, one per stage (build, score, cospike), plain result objects | Build writes spiked files and a key; user unmixes; score reads pairs | Matches the benchmark's `build / unmix / score / figures` modes (`88:13-17`); each stage runs on a different node and library | Three objects to pass around |
| B. One `sw_spike` object with methods (`$build()`, `$score()`) | Stateful object holds paths | Fewer arguments | No precedent in the package; hides state; hard to resume after a SLURM step |
| C. One function `sw_unmix_spike()` that does everything, unmix via callback | One call | Simplest for small data | Cannot span SLURM jobs or a second AutoSpectral library (`88:185-187`) |

Recommendation: **A**, confidence high. The unmix step must sit outside the package (Q3). Would change my mind: if the user only ever runs one in-process setup and never on SLURM.

Evidence: `88_spike_in_markers.R:13-17` (four modes, run per mode on SLURM); `88:185-187` (`.libPaths` swap to the 1.0.2 library); `R/unmix_report.R:409-571` (`sw_unmix_result()` as the precedent for a plain S3 list).

### Q2 Object system

Impact: medium. Depends on: Q1.

| Option | What happens | For | Against |
|---|---|---|---|
| A. S3 lists with a class and `print()` method (`sw_spike_design`, `sw_spike_score`, `sw_marker_screen`) | Same as `sw_unmix_result`, `sw_unmix_diagnostics` | Consistent with every unmix module; no new dependency | Less validation than S7 |
| B. S7 classes | Typed properties | Validation at construction | S7 is in Suggests and used only by the composable pipeline; unmix modules would then need it |
| C. Bare tibbles only, no class | Score returns a list of tibbles | Simplest | No `print`, no provenance slot |

Recommendation: **A**, confidence high. Would change my mind: a decision to move the whole unmix module to S7.

Evidence: `R/unmix_report.R:547-570` (S3 `structure(list(...), class = "sw_unmix_result")`); `DESCRIPTION` Suggests `S7 (>= 0.2.0)`; `CLAUDE.md` composable pipeline section (S7 only there).

### Q3 Where the unmix step plugs in

Impact: high. Depends on: Q1.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Decoupled: build writes spiked raw FCS; the user unmixes with `sw_unmix_run()`, raw AutoSpectral or any tool; score takes explicit base/spiked pairs | Package never runs or submits unmixing | Spans SLURM jobs, two AutoSpectral libraries, SpectroFlo exports; the heavy step stays outside the package | User writes the pairing table |
| B. Optional `unmix_fn` callback on score (signature `function(raw_paths, output_dir) -> unmixed paths`) | Score runs the callback per set in-process | Convenient for one setup | Runs hours of unmixing inside a scoring call; cannot span jobs |
| C. Score calls `sw_unmix_run()` itself from an `sw_unmix_result`'s setup | One object in | Reuses the package wrapper | Ties the diagnostic to one in-process setup; the benchmark setups could not be expressed this way (`88:199-329`) |

Recommendation: **A**, confidence high. A callback (B) can be added later as sugar (round 2, R2.1). Would change my mind: nothing seen so far.

Evidence: `88:168-338` (unmix mode per setup and set, three different call shapes, one SLURM job each); `experiments/README.md` rule 7 ("Compute goes through SLURM"); `R/unmix.R:1113-1310` (`sw_unmix_run()` runs one in-process AutoSpectral, returns `output_files`).

### Q4 Build input: how raw events arrive

Impact: high. Depends on: none.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Named character vector of raw FCS paths (names = sample ids) **or** a named list of numeric matrices (events x detectors) | Paths are read with `flowCore::read.FCS(truncate_max_range = FALSE)`; matrices allow mock tests | Same dual contract as `sw_unmix_sample_events()`; tests pass without flowCore | Two code paths |
| B. A `flowSet` | One object | Native flowCore | Tests need flowCore; `experiments/README.md` rule 3 says unit tests pass without it |
| C. A `GatingSet` | Carries populations | One object for Q5 | A GatingSet is built on cleaned, transformed data; the spike needs raw detector values (`62:2-3`) |

Recommendation: **A**, confidence high. Would change my mind: a decision to drop the "tests pass without flowCore" rule.

Evidence: `R/unmix_figures.R:255-258, 295-302` (paths or matrices); `experiments/README.md` rule 3; `88:130-152` (raw read, exprs edit, `write.FCS`); `CLAUDE.md` design constraint `truncate_max_range = FALSE`.

### Q5 How populations are specified

Impact: high. Depends on: Q4.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `populations`: a list `populations[[sample]][[pop]]` of integer **raw row** indices (required) | Build samples cells from these rows | Backend-free; the benchmark's own structure (`88:100-121` builds `raw_row, pop`) | The user maps gated rows to raw rows (round 2 helper, R2.2) |
| B. A `GatingSet` plus node names | Build reads `gh_pop_get_indices()` | One object | Gated rows index the cleaned file, not the raw file; mapping needs the key columns (`88:103-108`) |
| C. Logical matrix (events x pops) per sample | Dense | Easy to check | Big; no gain over A |

Recommendation: **A**, confidence medium-high. Add an optional helper `sw_unmix_spike_rows()` (GatingSet + cleaned file + raw file + key columns to raw rows) in round 2. Would change my mind: if every user gates on the raw unmixed file so no mapping is needed.

Evidence: `88:100-121` (raw rows per pop via `kept_rows[gh_pop_get_indices()]`); `88e:108-110` (CD3- as live minus CD3pos: a set operation on rows, which option A supports and B does not without a node); `R/gate.R:194-195` (`sw_gate_run()` works on FCS files, returns a GatingSet).

### Q6 Alignment of raw and unmixed rows

Impact: high. Depends on: none.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Required `key_cols` (e.g. `c("Time", "FSC-A", "FSC-H", "SSC-A")`); score stops when the key of base and spiked unmixed rows differs, and when the raw key differs from the design's | One identity check per file, as the benchmark did | Catches reordering or dropped events; the only check the rule set allows is exactly this one | A key that is not unique per event gives a false pass; user must know the columns |
| B. Row count only | Cheap | No arguments | Silent permutation |
| C. `include.raw = TRUE` on the unmix so raw detectors travel with the output | Strongest check | Needs no key choice | Doubles output size; only AutoSpectral offers it (`unmix.fcs` help: `include.raw`) |

Recommendation: **A**, confidence high; document C as the stronger check for small runs. Would change my mind: evidence that AutoSpectral can reorder events (none seen: `README.md:20` says the four arms keep identical `Time, FSC-A, FSC-H, SSC-A`, verified on five files only).

Evidence: `88:51, 97, 135-140, 386-391`; `helper75:13, 40-45`; `README.md:20`; AutoSpectral `unmix.fcs` help (`include.raw`, default `FALSE`).

### Q7 Pairing base and spiked unmixed files in score

Impact: high. Depends on: Q3.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Explicit tibble `unmixed` with columns `setup, sample, set, base, spiked` (paths or matrices) | Score reads exactly what it is told | No filename convention baked in; works for 1.0.2, 1.8.8, any tool | The user builds the table |
| B. Filename rule (`<raw>.fcs` to `<raw> AutoSpectral.fcs`) plus folder per setup | Fewer arguments | Matches AutoSpectral 1.8.8 output | `88:359` hard-codes the rule; breaks with `file_suffix` or other tools |
| C. `sw_unmix_run()` output lists per setup | Uses `output_files` | Package-native | Only for in-process runs (Q3) |

Recommendation: **A**, with a round-2 helper that fills the table from AutoSpectral's naming (R2.3). Confidence high. Would change my mind: nothing seen.

Evidence: `88:359` (`unmixed_name <- \(raw_file) sub("\\.fcs$", " AutoSpectral.fcs", raw_file)`); `R/unmix.R:1082, 1291-1295` (`file_suffix`, `output_files` by mtime); `88:362-401` (score joins by `raw_file`, `raw_row`).

### Q8 Required arguments of build (impactful, no default)

Impact: high. Depends on: none.

| Argument | Benchmark value | Why impactful |
|---|---|---|
| `spectra` (fluorophores x detectors, the dye rows used for the spike) | 1.8.8 spectra `..._143803.csv` in every setup | The spike's physical shape; `88:8-9` chose the 1.8.8 spectrum for all setups |
| `dyes` (named: marker = dye row) | `c(TNF = "BV510", ...)` | Which channels are tested |
| `levels` (arcsinh units; must include 0) | `c(0, 1, 2, 3)` | Dose; level 0 is the baseline check |
| `cofactor` | 6000 | `k = cofactor * sinh(level)` |
| `seed` | 42 | Which cells get which level |
| `allocation` | equal split `sample(rep_len(levels, n))` | Alternative: fixed `n_per_level`; changes cell counts per level |
| `output_dir` | one folder per marker | Where spiked FCS go (Q22) |
| `key_cols` | `c("Time", "FSC-A", "FSC-H", "SSC-A")` | Q6 |

| Option | What happens | For | Against |
|---|---|---|---|
| A. All eight required, no defaults; `levels` must contain 0 | Error names the missing argument | ROADMAP rule 2; matches `sw_unmix_sample_events()` style | Verbose calls |
| B. Defaults for `levels`, `allocation`, `key_cols`; others required | Shorter calls | Convenience | A default level grid is an impactful choice (`README.md:81` "user-approved grid") |

Recommendation: **A**, confidence high. Would change my mind: a written decision to grandfather the benchmark values as package defaults (then they go to the register as "pre-existing defaults").

Evidence: `ROADMAP.md:25` (rule 2); `88:45-54` (parameters with reasons); `88:144-146` (equal split); `62:29-31`; `R/unmix_figures.R:279-293` (`.sw_fig_require`).

### Q9 Spectrum row scale check

Impact: high. Depends on: Q8.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Error unless each spiked dye row has `max == 1` (within `1e-6`) | User normalises on purpose | No hidden transform; the recovery formula assumes it (`62:3-4`) | spectracle spectra may need a normalise call first |
| B. Normalise rows to max 1 inside build and record the factor | Always works | Convenient | A silent rescale changes `k`; against the impactful-choices rule |
| C. Argument `normalise = TRUE/FALSE`, required | Explicit either way | Explicit | One more argument |

Recommendation: **A**, confidence medium-high. Would change my mind: if spectracle-derived spectra are a common input; then C.

Evidence: `62:3-4` ("Spectra rows have max 1, so an exact unmix raises BV510 by k"); `88:149-150` (`outer(key$k, s188[dye, ])`); `R/unmix_report.R:82-128` (`.sw_check_named_matrix`, `.sw_align_detectors` for the detector-name check to reuse).

### Q10 Score arguments and output contract

Impact: high. Depends on: Q7.

Required: `design` (from build), `unmixed` (Q7 table), `loss_cut` (benchmark 0.5: "lost" = recovery below 0.5). Optional: `results` (named list setup -> `sw_unmix_result`, Q27). The channel name is `<dye>-A` or `<dye>` (reuse `.sw_match_channels()`).

| Output | Columns |
|---|---|
| `cells` | `setup, sample, pop, marker, dye, level, k, raw_row, base_value, recovered, recovery` |
| `summary` | `setup, marker, pop, sample, level, n_cells, median_recovery, share_lost` |
| `baseline_check` | `setup, marker, sample, n, max_abs_diff` (level-0 cells) |
| `loss_cut`, `design_id`, `provenance` | scalars and list |

| Option | What happens | For | Against |
|---|---|---|---|
| A. As above; donors stay separate rows; no donor-median column | Donor is the replicate; the figure or report takes medians | Matches `88:417-423` | The report must aggregate |
| B. Add a `by_setup` table with donor medians and ranges | Ready for the report | One less step | Hides donor spread unless both kept |

Recommendation: **A plus B** (both tables; B is derived and cheap), confidence high. Would change my mind: nothing.

Evidence: `88:361-424`; `88:50` (`LOSS_FRAC`); `INDEX.md:252-260` (readout: share losing >50%, donor median, median recovery in brackets).

### Q11 Co-spike shape

Impact: medium. Depends on: Q8, Q10.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `sw_unmix_spike_cospike(design, dye2, levels2, output_dir)` adds a second dye to the already-spiked cells of an existing design, one spiked set per `levels2` value; score reports both recoveries | Same cells as the single spike (`88d:38, 82-89`) | Direct comparison with the single spike; no resampling | Needs the first design |
| B. `sw_unmix_spike_build()` takes a list of dye grids and does the cross | One function | Fewer exports | Grid semantics get complicated; the task names a separate function |
| C. Independent resample | Simpler code | Independent | Loses the paired comparison |

Recommendation: **A**, confidence high. Would change my mind: nothing.

Evidence: `88d:1-9, 82-113` (reads `spike_key_TNF.csv`, level-1 cells get TNF + GLUT1 at 0/1/2, level-0 stay unspiked); `INDEX.md:320-335`.

### Q12 Module 2 input contract

Impact: high. Depends on: Q4, Q6.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `x`: list `x[[setup]][[sample]]` of matrices or FCS paths, rows paired across setups (same raw events); `reference`: a setup name (required); `populations[[sample]][[pop]]` integer rows on the paired rows; `key_cols` for the identity check | One function scores every channel in every population | Matches `helper75:22-47` and `75:85-110` | Memory: all setups of one sample at once |
| B. Per-setup call with the reference matrix passed each time | Lower memory | Simple | The user loops; easy to pass the wrong reference |
| C. GatingSets per setup | Native | Gates travel along | Pairing across setups is still by key, and `GatingSet` rows are the cleaned rows |

Recommendation: **A**, confidence high. Would change my mind: datasets where one sample's matrices do not fit in memory.

Evidence: `helper75:22-47` (rows paired on `as188_joint`'s cleaned rows, identity stop); `75:89-110`; `README.md:93` (test 1 settings).

### Q13 Module 2 metrics and required arguments

Impact: high. Depends on: Q12.

Metrics per `(setup, sample, pop, channel)`, on linear values, cut on arcsinh: `n`, `rsd_neg` (Nguyen rSD = q84 minus q50 of events below the cut), `ref_rsd` (same in the reference setup, same rows), `rsd_ratio`, `squeezed_share` (share with `|x| < squeeze_frac * ref_rsd`), `pos_share` (share above the cut). Required: `cofactor` (6000), `cut` (1), `squeeze_frac` (0.5), `min_events` (100; metric `NA` below it).

| Option | What happens | For | Against |
|---|---|---|---|
| A. Exactly helper75's metrics, plus `min_events` to `NA` | Reproduces the benchmark numbers | Regression oracle exists (`scores.csv`) | Nguyen rSD source not opened |
| B. Add MAD-based spread as a second rSD | More | Comparison | README dropped MAD spread as a ranking metric (`README.md:55`) |

Recommendation: **A**, confidence high. Would change my mind: nothing.

Evidence: `helper75:49-69`; `75:43-45`; `README.md:56-57` (Nguyen rSD primary; `MIN_N = 100`); `README.md:55` (MAD spread dropped).

### Q14 Screen rule and its arguments

Impact: high. Depends on: Q13.

Rule (revised, `README.md:103`): flag `(setup, pop, channel)` when the donor-median rise in squeezed share vs the reference is `>= min_delta` (0.10) **and** the rise has the same sign in every donor. The rSD ratio is reported, not used.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `sw_unmix_marker_screen(shape, reference, min_delta, min_donors, channel_populations = NULL)`: `min_donors` required (the number of donors that must all agree; rows with fewer donors are `NA`); `channel_populations` optional table `channel, pop` to restrict to designated parents | Reproduces 86b | Donor count comes from the user, not a constant `N_DONORS` | Two required numbers |
| B. Derive `min_donors` from the data (all donors present) | One argument less | Automatic | A missing donor silently relaxes the rule |
| C. Keep the dropped rSD criterion as an option | Flexibility | History | README superseded it (`README.md:102-103`) |

Recommendation: **A**, confidence high. Would change my mind: nothing.

Evidence: `86b:30-31, 83-124`; `README.md:102-103`; `75:112-135` (designated parents via `marker_parent.csv`).

### Q15 Integration into `sw_unmix_diagnose()` check rows

Impact: high. Depends on: Q10, Q14.

The check table has fixed columns `check_id, scope, item, metric, value, limit, status, reason` (`unmix_report.R:600-619`), limits are a closed spec (`:37-64`, `:638-647`), and `sw_unmix_diagnose()` takes one `result`.

| Option | What happens | For | Against |
|---|---|---|---|
| A. New optional arguments `spike = NULL`, `screen = NULL` on `sw_unmix_diagnose()`; rows `check_id = "spike_recovery"` (metrics `share_lost`, `median_recovery`, `baseline_max_abs_diff`) and `check_id = "marker_screen"` (metrics `median_d_squeezed`, `median_rsd_ratio`); new limit names `max_share_lost`, `min_median_recovery`, `max_baseline_diff`, `max_delta_squeezed`; `not_judged` unless a limit is supplied | One checks table, one report | Fits the existing plumbing; report and CSV come free | `item` must encode `marker|pop` (round 2, R2.4); donor aggregation needed per row (R2.4); the screen compares setups, so the `result` is the setup being judged and `reference` is named in `reason` |
| B. Separate `sw_unmix_diagnose_spike()` and `sw_unmix_diagnose_screen()` returning `sw_unmix_diagnostics` objects the user binds | No change to `sw_unmix_diagnose()` | Isolation | Two tables; report needs a merge path anyway |
| C. No integration; module tables stand alone | Nothing changes | Least work | The user asked for diagnostics "added on top"; the report would not show them |

Recommendation: **A**, confidence medium-high. Would change my mind: if widening the check schema (a `group` column) is refused; then B.

Evidence: `R/unmix_report.R:600-619, 37-64, 622-740, 1288-1365`; `R/unmix_report.R:1038` (`.sw_match_channels`).

### Q16 Verdict limits: sourced cutoffs or not judged

Impact: high. Depends on: Q15.

Values in use: loss cut 0.5 (`88:50`), pre-registered recovery 0.9 (`64:41`), flag rule 0.10 and all donors (`README.md:103`), `MIN_N` 100 (`README.md:57`). All are user choices from the benchmark. AutoSpectral article 11 says the joint pipeline "can understate a dimmer collinear fluorophore's true abundance" and gives no number (`INDEX.md:176`; article not opened here).

| Option | What happens | For | Against |
|---|---|---|---|
| A. No built-in limit; status `not_judged`; no `sw_unmix_limits_*()` helper for these metrics; roxygen records the benchmark values as "user choice, 2026-11-10, uncalibrated" | Same policy as condition number and cosine (`ROADMAP.md:186`) | Honest: nothing is sourced | The report shows numbers without verdicts |
| B. Helper `sw_unmix_limits_bmv()` returning the benchmark values | One call to reuse them | Convenience | Would present a one-dataset choice as a package limit |
| C. Hard-coded limits | Simple | None | Against ROADMAP rule 2 |

Recommendation: **A**, confidence high. Would change my mind: a published cutoff for spike loss or squeezed share.

Evidence: `ROADMAP.md:186`; `R/unmix_report.R:1167-1191` (the SIF helper exists only because article 01 gives numbers); `INDEX.md:171-183`.

### Q17 Report sections

Impact: medium. Depends on: Q15.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `sw_unmix_report()` gains sections "Spike-in recovery" and "Marker screen" (tables only, as the report draws no figures) plus `spike_summary.csv`, `marker_screen.csv`; sections appear only when the diagnostics carry those rows | Tables of donor medians and ranges; per-donor CSV | Consistent with the HTML-only report | More report code |
| B. Only the checks table (already written) | No new sections | Minimal | Per-donor detail lost |

Recommendation: **A**, confidence medium. Would change my mind: a decision to move reporting to the dashboard (Q28).

Evidence: `R/unmix_report.R:1530-1556, 1737-1788`; `.sw_unmix_report_files` (`:69-76`) must list the new CSVs so `overwrite` replaces them.

### Q18 Figure backend: ggcyto or ggplot2

Impact: high. Depends on: none.

Facts: ggcyto 1.34.0 with ggplot2 4.0.3 fails (`10_ggcyto_strategy.R:9-10`, "add_ggcyto error"; not reproduced this session); the benchmark used a private ggplot2 3.5.2 lib (`72:19`, `88c:13`, `93:17`, `94:16`, `95:14`). `facet_grid(axes = "all", axis.labels = "margins")` needs ggplot2 >= 3.5.0 (ggplot2 `NEWS.md:531, 801`). `geom_hex()` needs `hexbin`. The restyle uses `ggcyto()` only as a data adapter plus `ggcyto_par_set(limits)`; the geoms are plain ggplot2 (`95:258-304`).

| Option | What happens | For | Against |
|---|---|---|---|
| A. Pure ggplot2 (>= 3.5.0) on a tibble of events; `hexbin` to Suggests; ggcyto not used | Works in the main env; no private lib | One dependency floor; the package's other figures are pure ggplot2 | Loses `ggcyto` gate drawing (the restyle already replaced `geom_gate` with `geom_hline`, `95:257-260`) |
| B. ggcyto required; document the ggplot2 pin | Same look as the benchmark | Flow-native | Breaks in the main env; pins the package to an older ggplot2 |
| C. Both, `backend = c("ggplot2", "ggcyto")` | Choice | Flexibility | Two code paths to test; the ggcyto path untestable here |

Recommendation: **A**, confidence high. Would change my mind: a ggcyto release that supports ggplot2 4.x together with a feature only ggcyto gives.

Evidence: `10_ggcyto_strategy.R:9-10`; `95_restyle_events_scatter.R:14, 258-304`; `DESCRIPTION` (`ggplot2` unversioned, `ggcyto (>= 1.22.0)`); ggplot2 4.0.3 `DESCRIPTION:38` (hexbin in Suggests).

### Q19 Figure placement and panel key

Impact: medium. Depends on: Q18.

The figures module keys panels F1 to F4 (`unmix_figures.R:36-50`), `sw_unmix_save_figures()` rejects other keys (`:974-984`), the palette accepts six fixed method labels (`:53-82`), and `.sw_fig_check_events()` forbids more than one file per method (`:217-231`), which a setup x donor grid needs.

| Option | What happens | For | Against |
|---|---|---|---|
| A. New file `R/unmix_figures_scatter.R`; register panel `F6 = "scatter-grid"` and a size in `.sw_fig_panels` / `.sw_fig_sizes`; a new events contract (`setup, sample, x, y`) with its own check; `sw_unmix_save_figures()` unchanged otherwise | One save path for all panels | Reuses PNG, PDF, source CSV writing and the overwrite refusal | Touches two internal constants |
| B. Own `sw_unmix_save_scatter_grid()` | No change to existing constants | Isolation | Duplicates saving code |
| C. Put the builder inside `R/unmix_figures.R` | One file | Fewer files | The file is already 1110 lines |

Recommendation: **A**, confidence medium-high. Would change my mind: if `sw_unmix_save_figures()` adopts the Nature checks (Q20) and the two panels families diverge.

Evidence: `R/unmix_figures.R:36-50, 53-82, 203-233, 966-1095`.

### Q20 Figure theme and the Nature checks

Impact: medium. Depends on: Q18, Q19.

Package theme: `theme_pubr(8)` or `theme_classic(8)`, legend bottom (`unmix_figures.R:96-111`). Benchmark theme: `ggprism::theme_prism(7, "Nimbus Sans")`, bold titles and strips, legend right, `strip.clip = "off"` (`theme.R:14-29`); `save_panel()` stops unless width is 89 or 183 mm, text is 5 to 7 pt, and every PDF font is embedded (checked with `pdffonts`, `theme.R:66-97`).

| Option | What happens | For | Against |
|---|---|---|---|
| A. The new figure uses a new internal `.sw_fig_theme_nature(font_family)` (theme_prism when ggprism is installed, else theme_classic 7 pt; legend right); F1 to F4 untouched; a `check = c("none", "nature")` argument on `sw_unmix_save_figures()` runs the width and text-size checks, and the font check only when `pdffonts` is on the PATH | Benchmark look for the new figure; checks opt-in | No poppler dependency; no change to approved F1 to F4 | Two themes in one package until F1 to F4 migrate (round 3) |
| B. Migrate F1 to F4 to the Nature theme now | One theme | Consistent | Changes approved figures without a plan row (ROADMAP rule 5) |
| C. Keep the package theme for the new figure | No new theme | Minimal | Contradicts the user's figure rules (7 pt, legend right, strips) |

Recommendation: **A**, confidence medium. Would change my mind: a decision to restyle F1 to F4 in the same change.

Evidence: `R/unmix_figures.R:96-111`; `lib/theme.R:1-29, 42-110`; `pixi.toml` (ggprism, ggh4x present; not in `DESCRIPTION`).

### Q21 Figure arguments: required vs display-only defaults

Impact: medium. Depends on: Q18.

| Argument | Role | Proposal |
|---|---|---|
| `events` (tibble: `setup, sample, x, y`) | data | required |
| `x`, `y` (column names) | data | required |
| `cut` (arcsinh, on y; the line and the % label) | changes a printed number | required |
| `cofactor` (for the axis label and linear-to-arcsinh of listed columns) | changes numbers | required |
| `setup_fill` (named vector, one colour per setup; strips filled; error on an unknown setup) | colour key of the data | required |
| `transform_cols` (which columns are linear and get arcsinh) | changes numbers | required (may be `character()`) |
| `reference_curves` (tibble `x, y, curve`; e.g. the full-recovery curve of 88c) | display | optional, `NULL` |
| `pct_label = c("above_cut", "none")` | display | default `"above_cut"` |
| `bins = 50`, `limits = NULL`, `label_pt = 6`, `font_family = NULL` | display only | defaults allowed |

| Option | What happens | For | Against |
|---|---|---|---|
| A. As the table: anything that changes a printed number is required; pure display has defaults | Short calls for display | ROADMAP rule 2 scope ("changes a ranking or a selected setting") | Judgement calls on the border (bins) |
| B. Everything required | Uniform | No judgement | Noisy calls for display knobs |

Recommendation: **A**, confidence high. Would change my mind: nothing.

Evidence: `72:100-106` (`AXIS`, `BINS` marked "display only"; `CUT` user choice); `95:149-154`; `ROADMAP.md:25`.

### Q22 Disk and I/O for spiked files

Impact: medium. Depends on: Q4.

Spiked FCS are the size of the raw files, one copy per marker (and per co-spike level). The benchmark wrote them next to the data (`88:23-24`), outside `runs/`.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `output_dir` required; must not exist; one subfolder per marker (`<output_dir>/<marker>/<raw basename>`); `design.rds` and `spike_key.csv` written in `output_dir`; no temp dirs (files must outlive the R session for SLURM) | Same layout as the benchmark | Refuses to overwrite; resumable | Disk use = n_markers x raw size |
| B. `tempdir()` by default | No path to choose | Convenience | Files vanish before the SLURM unmix runs; against "never a hidden default" |
| C. Write spiked matrices only to RDS | Smaller | Fast | AutoSpectral needs FCS files |

Recommendation: **A**, confidence high. With matrix input (Q4) and `output_dir = NULL`, build returns spiked matrices in memory for tests. Would change my mind: nothing.

Evidence: `88:23-24, 125-127` (dir must not exist); `R/unmix_report.R:1549-1553` (same "must not exist" policy); `R/unmix_figures.R:1025-1036`.

### Q23 Runtime and HPC

Impact: medium. Depends on: Q3.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Package functions never submit jobs or call AutoSpectral for the spike; build and score are single-process; `experiments/slurm/` keeps the sbatch wrappers; an optional `sw_unmix_spike_jobs(design, setups)` returns a tibble of (setup, set, raw_dir, out_dir) the user feeds to sbatch (round 2, R2.5) | Package stays portable | Matches `experiments/README.md` rule 1 and 7 | The user runs the jobs |
| B. `future`/`callr`-based in-package runner | One-call runs | Convenience | New dependency; cannot use a second AutoSpectral library |

Recommendation: **A**, confidence high. Would change my mind: nothing.

Evidence: `experiments/README.md` rules 1, 7; `experiments/slurm/run_af405_proof.sbatch` (exists; not opened); `88:13-17`.

### Q24 Testing

Impact: medium. Depends on: Q4, Q12.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Unit tests with matrices only (no flowCore): identity "unmix" gives recovery exactly 1 and `share_lost` 0; a "dropped spike" unmix gives `share_lost` 1; level-0 cells give `max_abs_diff` 0; key mismatch stops; screen rule on a synthetic 5-donor table reproduces 86b's logic; figure builder returns a ggplot and a source table with the right % | Fast, no stack | Rule 3 of `experiments/README.md` | FCS reading untested in unit tests |
| B. Plus an FCS round trip guarded by `skip_if_not_installed("flowCore")` | Covers the path branch | Realistic | Slower |
| C. Plus `experiments/integration/spike_bmv.R` on the BMV data reproducing 88's `spike_summary.csv` (1.8.8 zip TNF CD8T: 2.1% lost, median 1.00) and 86b's `flags.csv` | Regression oracle | Real numbers | Needs the data folder and SLURM; not run in `R CMD check` |

Recommendation: **A + B + C**, confidence high. Would change my mind: nothing.

Evidence: `tests/testthat/test-unmix_report.R:1-2`; `tests/testthat/test-unmix_figures.R:1-60`; `INDEX.md:256, 283-289`; `experiments/integration/` (empty).

### Q25 Documentation

Impact: low. Depends on: Q1 to Q3.

| Option | What happens | For | Against |
|---|---|---|---|
| A. New vignette `vignettes/unmix-diagnostics.qmd` on mock matrices (always builds), with the BMV workflow as a non-evaluated section; a `reference/` page; a `CLAUDE.md` pipeline-table row; ROADMAP item B11 and register rows | Complete | Users find it | More files |
| B. A section in `spectral-unmixing.qmd` | One vignette | Less | That vignette is about the pipeline, not diagnostics |

Recommendation: **A**, confidence medium. Would change my mind: nothing.

Evidence: `vignettes/` listing; `CLAUDE.md` pipeline table; `ROADMAP.md:60-77`.

### Q26 Naming

Impact: low. Depends on: none.

| Option | What happens |
|---|---|
| A. `sw_unmix_spike_build()`, `sw_unmix_spike_cospike()`, `sw_unmix_spike_score()`, `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()`, `sw_unmix_figure_scatter_grid()`; classes `sw_spike_design`, `sw_spike_score`, `sw_marker_screen`; files `R/unmix_spike.R`, `R/unmix_screen.R`, `R/unmix_figures_scatter.R` | As asked, `sw_unmix_*` prefix, verb last as in `sw_unmix_check_controls()` |
| B. `sw_spike_*` prefix | Shorter; breaks the `sw_unmix_` grouping in the reference index |

Recommendation: **A**, confidence high. Would change my mind: a package-wide rename of the unmix module.

Evidence: `NAMESPACE` conventions as listed in `CLAUDE.md` pipeline table.

### Q27 Provenance link from score to the setups

Impact: medium. Depends on: Q10.

| Option | What happens | For | Against |
|---|---|---|---|
| A. Optional `results = list(<setup> = sw_unmix_result)`; score stores each `run_id`, AutoSpectral version, commit and backend in `provenance$setups`; setups without a result get `NA` and a note | Baseline recorded per ROADMAP gold-standard rule 5 | No new objects | Optional, so often empty |
| B. Required | Always recorded | Strict | SpectroFlo or 1.0.2 setups have no `sw_unmix_result` |

Recommendation: **A**, confidence high. Would change my mind: a rule that every scored setup must be an in-package run; then B.

Evidence: `ROADMAP.md:41-42`; `R/unmix_report.R:305-358` (`.sw_unmix_provenance`).

### Q21b Plot-first row for the scatter grid (ROADMAP rule 5)

Impact: high. Depends on: Q18 to Q21, R2.10. A plotting function is added only after the user approves its row. Format as `plans/beta_unmix_report_figure_plan.md:14-21`.

| # | Figure | Question it answers | Source data (saved by the run) | Visual check | Plot or table? |
|---|---|---|---|---|---|
| F6a | Scatter grid, mode `above_cut`: hexbin of marker y against marker x, one panel per setup (columns, strips filled by the setup key) x sample (rows), the cut drawn as a line, `% above cut` printed per panel; one population per figure | Does a setup move cells across the positivity cut, and where in x do the moved cells sit (a band at y = 0 at high x, as in `INDEX.md:232`)? | `figures/panels/source_data/scatter-grid.csv`: `setup, sample, n, pct` (the printed numbers), as `95:240-246, 326` | The printed % equals the source row; a band at y = 0 is visible only in the setups the screen flags | plot |
| F6b | Scatter grid, mode `lost`: per spiked cell, y = marker in the spiked unmix against x = marker in the unspiked unmix; dashed full-recovery curve `asinh(sinh(x) + sinh(level))`, solid total-loss diagonal; `% lost` printed per panel | Which cells lose the spike: are they the dim ones near x = 0, and do they sit on the diagonal (total loss) or between the curves (partial)? | `spike_cells.csv` rows at the shown level, plus the same `scatter-grid.csv` with `pct` = share lost | Cells of a clean setup sit on the dashed curve; `% lost` equals `summary$share_lost` for that setup, sample and level | plot |

| Option | What happens |
|---|---|
| A. Approve F6a and F6b as one builder with two modes (R2.10) | One panel key `F6`, source table shared |
| B. Approve F6a only; F6b later with its own row | Spike raw companions stay in `experiments/` for now |
| C. Approve neither | Module 3 is dropped from this change |

Recommendation: **A**, confidence medium-high. Would change my mind: a preference to keep spike figures out of the package until the summary figures (R3.4) are planned too.

Evidence: `plans/beta_unmix_report_figure_plan.md:14-21` (format); `72:196-249`, `95:248-327` (F6a method); `88c:63-162` (F6b method and curves); `INDEX.md:264` (what F6b showed).

### Q28 Scope: what stays out

Impact: medium. Depends on: none.

Proposed out of scope for this change (each could be its own later plan): running or submitting unmixing (Q3); gating and GatingSet creation; the ablation gap-closed scoring and labels (51); the spillover spread matrix (82); UMAP kNN (92); the cluster raw companions (94: cluster overlay dots on a hexbin, which is a different figure); a dashboard tab for spike or screen results (`R/unmix_dashboard.R`); handling of a second AutoSpectral library; the spike summary figures (donors as points, ECDF: `88:453-551`), which need their own plot-first rows (round 3, R3.4).

| Option | What happens |
|---|---|
| A. As listed | Three modules, one figure |
| B. Also include the spike summary and ECDF figures now | Two more builders; two more plan rows |
| C. Also include the cluster overlay variant (94) as a `overlay` argument of the scatter grid | One more argument; a different data contract (cluster rows) |

Recommendation: **A**, confidence medium. Would change my mind: if the user wants the full 88 figure set in the package in one go (then B, with rows in round 3).

Evidence: `88:432-556`; `94:1-14`; `R/unmix_dashboard.R:67-101`.

---

## Round 2: depends on round-1 answers

### R2.1 Callback sugar for in-process runs

Depends on: Q3 = A. Impact: low.

| Option | What happens |
|---|---|
| A. No callback in this change | Users call `sw_unmix_run()` themselves |
| B. `sw_unmix_spike_unmix(design, setup_fn, out_dir)` where `setup_fn(raw_dir, out_dir)` returns unmixed paths; it fills the Q7 pairing table | Convenience for one in-process setup |

Recommendation: A now, B later. Confidence medium.

### R2.2 Helper from GatingSet to raw rows

Depends on: Q5 = A, Q6 = A. Impact: medium.

`sw_unmix_spike_rows(gs, nodes, cleaned_files, raw_files, key_cols)` returns `populations[[sample]][[pop]]` by mapping cleaned-row keys to raw-row keys (`88:100-121`). Set operations (live minus CD3pos, `88e:108`) are the user's, on the returned integer vectors.

| Option | What happens |
|---|---|
| A. Ship it in this change (needs flowWorkspace; tests skip without it) | One fewer manual step |
| B. Document the mapping; ship later | Smaller change |

Recommendation: A, confidence medium.

### R2.3 Pairing helper for AutoSpectral output names

Depends on: Q7 = A. Impact: low.

`sw_unmix_spike_pairs(design, setups = c(<setup> = <unmix root>), suffix = " AutoSpectral")` builds the `unmixed` table from the benchmark's folder layout (`88:177-184, 359`). Required `suffix` (no default) or default `" AutoSpectral"` as in AutoSpectral 1.8.8?

| Option | What happens |
|---|---|
| A. `suffix` required | Explicit |
| B. Default `" AutoSpectral"` with `file_suffix` passthrough | Matches `sw_unmix_run()` output |

Recommendation: B (a filename convention is not an impactful analysis choice), confidence medium.

### R2.4 Check-row encoding and donor aggregation

Depends on: Q15 = A. Impact: high.

| Option | What happens | For | Against |
|---|---|---|---|
| A. One row per `(setup, marker, pop)`: `item = "<marker>|<pop>"`, `value` = donor median, `reason` = "n donors, range a to b" | Fits the schema | No schema change | Donor spread only in `reason` |
| B. Add columns `group` and `n` to the check table | Clean | Richer | Every existing check, test and the report table change |
| C. One row per donor | Nothing aggregated | Full detail | Many rows; limits judge each donor |

Recommendation: A, confidence medium; the per-donor CSV stays in the report (Q17). Would change my mind: a decision to widen the schema for other future checks.

### R2.5 Job table helper

Depends on: Q23 = A. Impact: low.

`sw_unmix_spike_jobs(design, setups)` returns one row per (setup, set) with `raw_dir`, `out_dir`. The user maps it to sbatch. Ship now or later? Recommendation: later, confidence medium.

### R2.6 Strip colouring implementation

Depends on: Q18 = A. Impact: low.

| Option | What happens | For | Against |
|---|---|---|---|
| A. `ggh4x::facet_grid2(axes = "all", remove_labels = "all", strip = strip_themed(background_x = elem_list_rect(fill = setup_fill)))`; `ggh4x` to Suggests; fall back to plain strips with a message when absent. Note: ggh4x 0.3.1 names the label argument `remove_labels`, not ggplot2's `axis.labels` (checked: `formals(ggh4x::facet_grid2)`) | Declarative; the user's own rule for colour-key strips | No gtable surgery | One more Suggests; the two APIs differ in argument names |
| B. gtable surgery as in `95:314-325` | No dependency | Works today | Fragile across ggplot2 versions |

Recommendation: A, confidence medium-high.

### R2.7 Nature save checks

Depends on: Q20 = A. Impact: low.

Which checks when `check = "nature"`: width in `c(89, 183)`; drawn text 5 to 7 pt via a grob walk (`theme.R:42-61`); font embedding only when `Sys.which("pdffonts")` is non-empty, else a note in the returned table. Recommendation: all three as described, confidence medium.

### R2.8 Spiked FCS metadata

Depends on: Q22 = A. Impact: low.

Write `$SW_SPIKE_DESIGN` (design id), `$SW_SPIKE_DYE`, `$SW_SPIKE_LEVELS` keywords into each spiked FCS header so a file can be traced without the key. Recommendation: yes, confidence medium (flowCore `keyword<-`; no AutoSpectral effect expected, not verified).

### R2.9 Co-spike levels and recovery of both dyes

Depends on: Q11 = A. Impact: medium.

Score reports `recovery` for dye 1 and `recovery_2` for dye 2 when `levels2 > 0` (`88d:8-9`). The co-spiked cells are the `level == level1` cells of the design (88d used level 1 only). Argument `level1` required? Recommendation: `level1` required, confidence medium.

### R2.10 Before/after mode of the scatter grid

Depends on: Q19 = A, Q21 = A. Impact: medium.

88c's figure is the same grid with `x = before`, `y = after`, two reference curves and a "lost >50%" label. Options: one builder with `pct_label = c("above_cut", "lost")` and `reference_curves`, or a second builder `sw_unmix_figure_spike_cells()`. Recommendation: one builder with the two modes; the "lost" mode takes the `loss_cut` through the events table column `lost`. Confidence medium.

### R2.11 Population column in the figure

Depends on: Q21. Impact: low.

Grid is setup x sample; one population per figure (`72:100`, `88c:75`). Keep one population per call, or add `pop` as a facet row group? Recommendation: one population per call; the caller loops. Confidence high.

---

## Round 3: depends on round-2 answers

The items below are stubs: a title, dependency and the decision at stake. They get options, for and against, a recommendation and evidence once round 2 is answered. The same holds for the shorter round-2 items (R2.5, R2.7 to R2.9, R2.11).

### R3.1 Report layout for the new sections

Depends on: Q17 = A, R2.4. Tables: donor median and range per (setup, marker, pop); link to `spike_summary.csv`. Figures stay out of the HTML report.

### R3.2 Vignette content and data gating

Depends on: Q25 = A, R2.2. Mock section always evaluated; BMV section `eval: false` unless `SW_DATA_DIR` is set.

### R3.3 Migrate F1 to F4 to the Nature theme

Depends on: Q20 = A. A separate plan row each (ROADMAP rule 5).

### R3.4 Spike summary figures (donors as points, ECDF)

Depends on: Q28. Each needs its own plot-first row: source table `spike_summary.csv`, question answered, visual check.

### R3.5 Dashboard tab

Depends on: Q28, Q17. `sw_unmix_dashboard()` reads run folders; a spike or screen tab needs a data-model entry in `R/unmix_dashboard_data.R`.

---

## Register rows this grill would add to ROADMAP

Each is an impactful choice with no package default; the starting value is the benchmark's and is marked `inherited`.

| Choice | Where | Source of the starting value |
|---|---|---|
| Spike spectra set (which dye spectra are added) | `sw_unmix_spike_build()` | inherited: `88:8-9` (1.8.8 spectra in every setup) |
| Spike levels (must include 0) | `sw_unmix_spike_build()` | inherited: `62:29`, `88:46` (0, 1, 2, 3) |
| Spike seed | `sw_unmix_spike_build()` | inherited: `62:31`, `88:49` (42) |
| Spike cofactor | `sw_unmix_spike_build()` | inherited: `CLAUDE.md` (6000), confirmed by user |
| Allocation of cells to levels | `sw_unmix_spike_build()` | inherited: `88:144` (equal split) |
| Key columns for row identity | build, score, shape | inherited: `88:51` (`Time, FSC-A, FSC-H, SSC-A`) |
| Loss cut (recovery below it = lost) | `sw_unmix_spike_score()` | inherited: `88:50` (0.5) |
| Co-spike levels of the second dye | `sw_unmix_spike_cospike()` | inherited: `88d:46` (0, 1, 2) |
| Negative cut (arcsinh) | `sw_unmix_negative_shape()` | inherited: `75:44` (1) |
| Squeeze fraction | `sw_unmix_negative_shape()` | inherited: `75:45` (0.5) |
| Minimum events per population | `sw_unmix_negative_shape()` | inherited: `README.md:57` (100) |
| Reference setup | `sw_unmix_negative_shape()`, `sw_unmix_marker_screen()` | inherited: `README.md:28` (`as188_joint`) |
| Screen minimum delta | `sw_unmix_marker_screen()` | inherited: `86b:30` (0.10) |
| Screen minimum donors that must agree | `sw_unmix_marker_screen()` | inherited: `86b:31` (5) |
| Figure cut and cofactor | `sw_unmix_figure_scatter_grid()` | inherited: `72:101-102` (1, 6000) |
| Setup colour key | `sw_unmix_figure_scatter_grid()` | user: `72:46, 61, 76, 90` |
| Verdict limits for spike and screen metrics | `sw_unmix_diagnose()` | none sourced; stay not judged (Q16) |

## Not verified in this session

AutoSpectral articles 01, 10, 11, 19 and Mage et al. 2026 (quoted through `INDEX.md` only); the Nguyen robust SD paper; the ggcyto failure under ggplot2 4.0.3 (one script comment); row-order preservation by AutoSpectral beyond the five BMV files; whether FCS keywords added to spiked files pass through `unmix.fcs()`; `pdffonts` availability outside this HPC; `experiments/slurm/run_af405_proof.sbatch` contents; whether `ggh4x::facet_grid2(axes = "all", remove_labels = "all")` draws the same axes as ggplot2's `facet_grid(axes = "all", axis.labels = "margins")` (argument names checked, output not drawn).
