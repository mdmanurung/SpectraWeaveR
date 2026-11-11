# Critique: unmix diagnostics grill and plan

Reviewed 2026-11-11, cold, against the sources. Files reviewed: `plans/unmix_diagnostics_grill.md` and `plans/unmix_diagnostics_plan.md`. Nothing else was edited.

## Verdict

The architecture is sound: decoupled build and score, an explicit pairing table, required impactful arguments, `not_judged` limits and pure ggplot2 are all right. The spike method has one unstated assumption: each cell unmixes independently of the other cells in its file. Run 88h shows this fails at `n.af.passes = 5`. The draft's baseline check (max linear difference, spiked channel, level-0 rows) cannot detect that failure well, and it is not a verdict. Three more method details are misdescribed or missing: the co-spike construction, the screen's NA and sign rule, and the `min_events` source. The scope is too wide. 45 questions and about 86 h can drop to 12 round-1 questions and about 55 to 60 h by deferring diagnose and report integration, the Nature save checks and the helpers.

## Top issues, ranked

### 1. Independence of cells is assumed, and the baseline check cannot carry it

**Issue.** Spike recovery is `spiked unmix - base unmix` of the same cell. That difference measures the spike only if unspiked cells do not move when other cells are spiked. 88 states this in a comment ("unmixing is per cell", `88:403`, `88h:214`) and never tests it as a pass/fail. Run 88h broke it: at `n.af.passes = 5` unspiked cells moved, up to about 0.26% of unspiked CD4T and CD8T cells by more than 0.05 arcsinh in any channel, and up to 0.02% by more than 0.5 in the spiked channel (`INDEX.md:360`). At `n.af.passes = 1` they did not move. Table fallback: the drift is a count, not a figure (`experiments/runs/2026-11-11_bmv-spike-passes/baseline_check.csv`; `INDEX.md:360`).

**Why it matters.** Under a method that is not per cell, recovered values mix the spike response with file-level drift. Recovery is then not a property of the spiked cell alone.

The draft's `baseline_check` has four gaps (grill Q10, plan step 4):
- It reports `max_abs_diff` in linear units. A linear max such as 5627 (`baseline_check.csv`, row `p1_af5,TNF,CD4T`) has no scale, and one cell sets it.
- It checks only the spiked marker's channel. The 0.26% figure counted any channel.
- It checks only level-0 cells inside the spiked populations, not the rest of the file. 88h's own score never computes the share figures: the CSV holds only the max (`88h:214-221`). The INDEX numbers come from "a per-cell drift summary computed from the unmixed files" that no script writes, so plan step 16 cannot reproduce them.
- It is not a verdict. Recovery tables are produced whatever the baseline shows.

**Likely mechanism (hedged).** The `unmix.folder()` 1.8.8 help defines `refine.af.quantile` (default 0.5) as the "Fraction of cells taken forward for additional AF passes". A quantile over the file's cells is a file-level selection. That fits the 88h result, but the selection code is compiled (`AutoSpectralRcpp::unmix.autospectral.rcpp` passes `refine_af_quantile` to C++) and I did not read it. `chunk.size` (default `2e6` events) is a second file-level boundary to document.

**Proposed fix.**
- `baseline_check` becomes a required output of score. Compute it in arcsinh (the design's `cofactor`) over every fluorescence channel shared by base and spiked files.
- Cover two row groups, reported apart: level-0 rows of the spiked populations, and all rows outside the design (`untouched`).
- Columns: `setup, sample, set, group, channel, n, max_abs_diff, share_moved`. `share_moved` is the share above a **required** `baseline_tol` (arcsinh; a register row, no default).
- Carry a per-file status (`independent` / `coupled`) onto `summary` and `by_setup`. Recovery numbers are kept, never set to NA, because rule 4 covers missing input, not a failed assumption. `print()` and any report row lead with the coupled status.
- With Q3 = A the package cannot know the unmix method, so it measures drift instead of inferring it. When `results[[setup]]` records `n.af.passes > 1`, add a note naming `refine.af.quantile` as the likely cause.

**Touches:** Q10, Q16, Q24, plan steps 4, 5 and 16, and the register (new row `baseline_tol`).

### 2. The spiked fraction per file is an unstated impactful choice

**Issue.** `allocation = "equal"` gives each population a 1/length(levels) share per level, so with levels `0:3`, 75% of T cells are spiked in every file (`88:144`, `62:86`). Under a coupled method, drift likely grows with the spiked fraction. The 88h data cannot show that, because only one fraction was run.

**Why it matters.** Recoveries are comparable across setups, sets and co-spikes only at the same spiked fraction. The co-spike sets spike about 25% (issue 3), and 88e spikes CD3- cells in separate files.

**Proposed fix.** Add `spike_fraction` (or explicit `n_per_level`) as a required build argument. Record the realised fraction per file in the design and in `baseline_check`. Offer an option to put one level per file set, which lowers the fraction per file at the cost of more unmix jobs. Ask the user (round-1 question 3 below).

**Touches:** Q8 (`allocation`), Q22 (disk use grows with sets), plan step 2.

### 3. The co-spike construction is misdescribed

**Issue.** Grill Q11 A says cospike "adds a second dye to the already-spiked cells of an existing design". 88d instead re-reads the **raw** file. It keeps only level-0 and level-1 rows of the TNF key (`88d:83-85`) and spikes the level-1 rows with TNF + GLUT1 (`88d:103-109`). Cells at TNF levels 2 and 3 are unspiked in every co-spike set and are absent from `cospike_key.csv`.

**Why it matters.** The co-spike files have a different spiked fraction (about 25%, not 75%), so issue 2 applies. The `levels2 = 0` set is a TNF-only spike at that lower fraction, which makes it a free positive control for the fraction effect, if the plan says so. The plan's `unmixed$set` contract does not say how a co-spike set is named.

**Proposed fix.** Plan step 3 should start from raw, take the `level1` rows only and treat the other levels as unspiked and unscored. It should name sets `cospike_<dye2>_<level2>`. Score should report the `levels2 = 0` set beside the single-spike level-1 recovery as a fraction check.

**Touches:** Q11, R2.9, plan step 3, test list.

### 4. Reproducing 88 by RNG is fragile; build needs an explicit assignment input

**Issue.** 88 calls `set.seed(SEED)` once per marker, then draws `sample(rep_len(LEVELS, n()))` sample after sample in manifest order (`88:128-147`). The draw depends on sample order, on the row order of `gh_pop_get_indices()` per population, on `distinct()` keeping the first population of an overlapping row (`88:116`), and on R's `sample.kind`. Plan step 16 asserts equality with 88's `spike_summary.csv` after rebuilding "with the same seed and levels". A package build that seeds per sample, or orders rows differently, fails that test without any method error.

**Proposed fix.**
- Build accepts either `seed` + `allocation` or a precomputed `assignment` tibble (`sample, raw_row, pop, level`). The regression feeds 88's `spike_key_<marker>.csv` directly.
- Cospike takes the same `assignment` route.
- Overlapping populations become an error, or need a stated precedence argument. The draft's silent "first wins" inherits `distinct()`.
- Per-population level counts are not equal under the pooled draw, so the plan's test "equal split gives `n %/% length(levels)` per level within 1" holds per sample only, not per population. Say which one build guarantees.

**Touches:** Q5, Q8, plan steps 2 and 16, test list.

### 5. Module 2's rule differs from both the README and the regression oracle

**Issue.** Four mismatches:
- (a) Grill Q14 calls the rule "the rise has the same sign in every donor", copying `README.md:103`. The code tests a rise: `all_donors_up = all(d_squeezed > 0)` (`86b:102`). A consistent fall has "the same sign" but is not flagged.
- (b) 86b requires `n_donors == N_DONORS` exactly and then turns NA into FALSE with `coalesce(flag, FALSE)` (`86b:107-121`). The plan says fewer donors gives NA, never FALSE. That follows ROADMAP rule 4, but it breaks the step-16 claim that the package flags equal `flags.csv`.
- (c) `min_events` (100) is cited as inherited from `README.md:57`. That row approves `MIN_N = 100` for scripts 30 and 31, not for the screen, and helper75 applies no minimum (`helper75:52-68`). In the screen it is a new choice, not an inherited one.
- (d) helper75 sets `n = length(rows)`, the population size. The rSD runs on the events **below the cut** (`helper75:58, 61`), so a large population with 3 negatives passes any population-level minimum.

**Why it matters.** Each mismatch changes which `(setup, pop, channel)` cells are flagged or NA.

**Proposed fix.**
- Ask the rule as the code runs it: rise in all donors, strict `> 0`.
- Put the floor on the below-cut count (`min_negatives`) and on the reference's below-cut count, which also guards `ref_rsd = 0` giving `Inf` ratios.
- Mark `min_negatives` as a new choice.
- Make step 16 compare only where 86b's flag was not NA-coerced, and list the NA-versus-FALSE rows as expected differences.
- `README.md` also has a "Shared negatives" row: "below the reference arm's cut only". helper75 instead defines each arm's negatives by that arm's own values (`helper75:61`). Which definition the package uses is a missing question.

**Touches:** Q13, Q14, register rows, plan steps 6 to 8 and 16.

### Further issues (lower rank)

6. **Over-recovery is invisible.** `share_lost` is one-sided (`recovery < loss_cut`). The 1.0.2 code with 1.8.8 spectra overshoots Th2 in CD8T at a median of 1.16 (`INDEX.md:263`). Add a two-sided readout (share outside `[loss_cut, 1/loss_cut]`, or quantiles of recovery per donor) to `summary`. This is a user question. Touches Q10.
7. **Spike spectrum versus unmix spectrum is not recorded.** For a setup that unmixes with the same spectrum as the spike, an OLS-like solve recovers 1 by construction. For a setup with a different spectrum (s102_c102), loss mixes spectrum mismatch with AF and variant choice. When `results` is supplied, store the cosine between the spike row and each setup's row for that dye. Touches Q8, Q27.
8. **The spike is noise-free.** The added signal has no shot noise, so recovery is an optimistic test for weighted or Poisson methods. Document it in roxygen; it needs no argument.
9. **Key uniqueness is not checked.** The identity stops in score use `identical()`, which is positional and fine. But the raw-to-cleaned mapping uses `%in%` (`88:108`, `helper75:38`), and a duplicated `Time|FSC-A|FSC-H|SSC-A` key maps silently. Add a uniqueness stop where keys are built. Touches Q6, R2.2, plan step 6.
10. **The diagnose integration does not fit one-run diagnose.** `sw_unmix_diagnose(result, ...)` judges one run. Spike and screen tables span setups and need a mapping from the user's setup label to `result` that the plan never defines. With Q16 = A every row is `not_judged`, so integration adds schema risk (`item = "marker|pop"`, the positional `direction` at `unmix_report.R:62`) for no verdict. Prefer Q15 B or C now (standalone objects with `print()` and CSV writers), and integrate later when a sourced limit exists. Flag for the user, do not decide.
11. **`min_events` collides with an existing limit name.** `.sw_unmix_limit_spec` already has `min_events` judging `n_events` (`unmix_report.R:46, 58`). Rename the module 2 argument (`min_negatives`, issue 5).
12. **One panel key, many files.** Q19 A registers `F6 = "scatter-grid"`, and R2.11 says "one population per call; the caller loops". `sw_unmix_save_figures()` names files from the panel stem and refuses existing files (`unmix_figures.R:1016-1033`), so a second population or marker in the same `dir` is refused. Add a `stem` suffix (for example `scatter-grid_<marker>_<pop>`) or one `dir` per call.
13. **Rule 3 ("call, don't rebuild") needs a written reason for module 3.** AutoSpectral 1.8.8 exports `compare.unmixed.data()`, a grid of biplots over already-unmixed data sets, one row per set. Module 3 rebuilds a similar layout. The reason exists: the hexbin with a % label against a cut, a donor x setup grid, the `lost` mode with reference curves, a source CSV, and no cowplot or JPEG side effects. The plan should state it next to Q18.
14. **API vocabulary drifts.** The figure module's events contract uses `file, method, event_id` with a fixed method palette. The new code uses `setup` with a free `setup_fill`. `sw_unmix_negative_shape()` returns a tibble with an attribute, while every other new function returns an S3 list. Pick one name (`setup` is better for multi-arm work; say that `method` is the single-run term) and one return style.

## Per-question notes (only where I disagree)

| Q | Note |
|---|---|
| Q5 | Agree with A. Add: populations must be disjoint or precedence stated (issue 4). The level draw is pooled across populations in 88, so per-population counts are random; ask whether the draw is per population. |
| Q6 | Agree with A. Add a key-uniqueness stop (issue 9). The draft cites "verified on five files only" under `README.md:20`, but that line says "all four arms ... same event count per sample" and does not say five files. |
| Q8 | Missing: `spike_fraction` (issue 2) and `assignment` (issue 4). Keep `key_cols` required, and require it to be disjoint from the spectrum's detector columns. A spiked detector in the key would break the identity check. |
| Q10 | `baseline_check` needs a redesign (issue 1). Add the two-sided readout (issue 6). Agree on keeping donors as rows. |
| Q11 | Mis-framed (issue 3). |
| Q13 | Option A does not "exactly" reproduce helper75: helper75 has no minimum. Floor the below-cut count (issue 5). |
| Q14 | Option text says "same sign"; the code says "rise" (issue 5). The NA rule differs from 86b. |
| Q15 | Prefer B or C now (issue 10). The draft is right that `direction` is positional and fragile (`unmix_report.R:62`). Keep that fix if A is chosen. |
| Q18 | Agree with A, with more confidence than the draft had. I reproduced the failure: ggcyto 1.34.0 under ggplot2 4.0.3 stops at plot construction with `object '+.gg' not found`. `grep ggcyto:: R/` is empty, so every package figure is already pure ggplot2. Add the rule-3 reason (issue 13). |
| Q20 | The Nature width, pt and `pdffonts` checks are the user's journal spec for benchmark 05 (`lib/theme.R:1-3`), not a package contract. Defer `check = "nature"` and ggprism to a later plan. For now, the new builder uses the existing `.sw_fig_theme()` with a `legend.position` override. |
| Q21 | Agree. `bins` is display-only. |
| Q21b | Agree that F6a and F6b are one builder with two modes. F6b depends on module 1, so approve it with Q10. |
| Q22 | Agree. Add: build writes a zero-spike round-tripped copy of each raw file as the `base` set (optional). Base and spiked then pass through the same `write.FCS` path, so a rewrite artefact cannot appear as recovery. The float32 FCS round trip is visible as `2.98e-08` in `baseline_check.csv`. |
| Q24 | Identity "unmix" gives recovery 1 only when dye spectra are unit vectors on distinct detectors. Use a real OLS mock (`qr.solve(t(S), t(x))`) on overlapping spectra. Add a coupled mock unmixer: subtract a file-level quantile of one channel from every cell. Its test asserts `baseline_check` flags `coupled` and `share_moved > 0`, which is the test that matters for issue 1. Use tolerances, not `== 0`, on FCS round trips. Mock `pdffonts` through an internal wrapper, not `Sys.which` itself (moot if Q20 is deferred). |
| Q28 | Agree with A. Also move Q15 integration, Q17 report sections, Q20 Nature checks and R2.2, R2.3 and R2.5 helpers out of this change. |

## Questions: add, merge, drop

### Proposed round 1 (12 questions, one sitting)

1. **Scope of this change** (merges Q15, Q17, Q20, Q28, R2.1 to R2.5). Ship build, cospike, score, shape, screen and one figure builder as standalone S3 objects with `print()` and CSV writers. Defer diagnose and report integration, Nature checks and helpers.
2. **Baseline independence check** (new; issue 1). Make it a required output with a per-file `independent`/`coupled` status? Required `baseline_tol` in arcsinh? Score all channels and the untouched rows?
3. **Spiked fraction per file** (new; issue 2). A required `spike_fraction` or `n_per_level`? One level per file set as an option?
4. **Populations and allocation** (Q5 + issue 4). Raw row indices, disjoint required, draw pooled or per population.
5. **Assignment input** (new; issue 4). Accept a precomputed `assignment` beside `seed`, used by the regression and by cospike.
6. **Spectrum checks** (Q9 + issue 7). Error unless max 1; store the spike-to-setup cosine when `results` is given.
7. **Score readouts** (Q10 + issue 6). Required `loss_cut`; add a two-sided readout or recovery quantiles?
8. **Co-spike semantics** (Q11 + R2.9 + issue 3). Rebuild from raw on `level1` rows, other levels unspiked, `levels2 = 0` as the fraction control.
9. **Negative definition** (new; issue 5). Each setup's own below-cut events (helper75), or the reference's below-cut rows (README "Shared negatives")?
10. **Screen metrics and floor** (Q13 + issue 5). helper75 metrics; floor on below-cut count (`min_negatives`), new choice, no inherited value.
11. **Screen rule** (Q14 + issue 5). Rise (`> 0`) in all donors and median rise `>= min_delta`; fewer than `min_donors` gives NA, not FALSE. Accept the stated differences from 86b's `flags.csv`.
12. **Figure** (Q18 + Q19 + Q21 + Q21b + R2.6 + R2.10 + R2.11). Pure ggplot2; one builder, two modes; a stem per marker and population; strip fills via ggh4x with a plain fallback, or no fills.

### Decided unless the user objects (not impactful)

Q1 A, Q2 A, Q3 A, Q4 A, Q6 A (+ uniqueness), Q7 A, Q22 A, Q23 A, Q25 (roxygen examples now; vignette deferred), Q26 A, Q27 A. Each had "would change my mind: nothing", or it is a structural or naming choice that moves no number.

### Deferred

R2.1 to R2.5, R2.7, R2.8, all of round 3 and Q17 move to a follow-up plan opened after the first real-data run of the package functions.

### Dropped

Q21 as a separate question (folded into round-1 question 12) and Q2 (no live alternative).

## Plan changes

| Step | Change | Hours |
|---|---|---|
| 1 | Add register rows `baseline_tol`, `spike_fraction`, `min_negatives`. Remove the "inherited" `min_events` row. Mark the screen rule source as `86b` code, not the README wording. | 1 |
| 2 | Add `assignment`, `spike_fraction`, disjoint-population and key-uniqueness stops, `key_cols` disjoint from spectrum columns, optional round-tripped base set. | 9 |
| 3 | Rewrite to the 88d construction (issue 3). | 4 |
| 4 | Redesign `baseline_check` (issue 1); two-sided readout; per-file status on `summary`. | 10 |
| 5 | Add the OLS mock and the coupled-unmixer test; tolerances on FCS round trips. | 7 |
| 6 to 8 | `min_negatives`, strict-rise rule, NA not FALSE; tests for a ref rSD of 0. | 12 |
| 9 to 11 | **Defer** (diagnose and report integration). | -12 |
| 12 | Keep both modes. Use `.sw_fig_theme()` with an override; stem suffix; no ggprism. | 7 |
| 13 | Register F6 and the stem suffix only; **defer** the Nature check. | 1 |
| 14 | Drop the Nature tests. | 3 |
| 15 | `ggplot2 (>= 3.5.0)`; Suggests `hexbin` (and `ggh4x` only if fills are kept). | 1 |
| 16 | Feed 88's saved keys through `assignment`. Compare `summary` to `spike_summary.csv` and `baseline_check` to its CSV. Compare flags only where 86b did not NA-coerce. Add one 88h `p1_af5` case that must come out `coupled`. | 6 |
| 17 | **Defer** (helpers). | -5 |
| 18 | Roxygen examples and the `CLAUDE.md` row only; vignette deferred. | 1 |
| 19 | Unchanged. | 3 |

The plan's per-step hours sum to 86, as stated. The revised core is about 65 h before review, or about 55 h if the figure builder ships one mode first.

### Risks to add

- A coupled unmix method makes recovery mix spike response with file drift (issue 1).
- The spiked fraction differs between single-spike, co-spike and 88e designs (issues 2 and 3).
- AutoSpectral may clip spiked values above `$PnR`. The spike adds up to `6000 * sinh(3)` ≈ 60,000 per detector. Not checked: whether AutoSpectral's reader honours `$PnR`. Build should raise `$PnR` where the spiked maximum exceeds it, and the FCS round-trip test should check it.
- `ggcyto` stays in `DESCRIPTION` Suggests (`>= 1.22.0`) but no `R/` file calls it, and it fails under the installed ggplot2. This is out of scope; say in `ROADMAP.md` whether it stays in Suggests.

## Evidence spot-check

| Citation | What it says | Status |
|---|---|---|
| `88:13-17` | four modes, run per mode on SLURM | ok |
| `88:185-187` | `.libPaths` swap to `lib_as102` for non-1.8.8 setups | ok |
| `88:144-146` | `level = sample(rep_len(LEVELS, n()))`, `k = COFACTOR * sinh(level)` | ok |
| `88:359` | `sub("\\.fcs$", " AutoSpectral.fcs", raw_file)` | ok |
| `88:50`, `88:51` | `LOSS_FRAC <- 0.5`; `KEY_COLS` | ok |
| `62:3-4` | "Spectra rows have max 1, so an exact unmix raises BV510 by k" | ok as quoted; it is a comment, and no code checks the max |
| `64:41` | `RECOVERY_MIN <- 0.9` pre-registered | ok |
| `88d:46` | `SETS <- c(cospike_g0 = 0, cospike_g1 = 1, cospike_g2 = 2)` | ok |
| `88d:82-89` (cited for Q11 "already-spiked cells") | keeps levels 0 and 1, rebuilds from raw | ok as a line ref, but the grill's reading is wrong (issue 3) |
| `88e:108-110` | live minus CD3pos set operation | content ok, lines off: the operation is at `88e:119-121` |
| `helper75:49-69` | Nguyen rSD of below-cut events; `n = length(rows)` | ok; the `n` meaning is not stated in the grill (issue 5) |
| `86b:30-31` | `D_SQUEEZE_MIN <- 0.10`, `N_DONORS <- 5` | ok |
| `README.md:57` | `MIN_N = 100` for scripts 30, 31 | wrong source for the screen (issue 5) |
| `README.md:103` | "same sign in all 5 donors" | quote ok; differs from code `86b:102` |
| `README.md:20` | four arms, same event count, key columns identical | ok; "five files only" is not in that line |
| `unmix_report.R:37-64`, `:62` | limit spec; positional `direction` | ok |
| `unmix_report.R:600-619` | `.sw_check_rows()` fixed columns | ok |
| `unmix_report.R:1038` | `.sw_match_channels()` exact name, else `-A` | ok |
| `unmix_report.R:1167-1191` | SIF limits helper from article 01 | ok |
| `unmix_report.R:1549-1553` | report path "must not exist unless `overwrite = TRUE`" | ok (roxygen; it allows overwrite, unlike build) |
| `unmix_figures.R:36-50` | panels F1 to F4, sizes | ok |
| `unmix_figures.R:53-82` | "palette accepts six fixed method labels" | wrong: 10 labels, including 4 correction labels |
| `unmix_figures.R:217-231` | one file per method, else error | ok |
| `unmix.R:1291-1295` | `output_files` by mtime | ok |
| `10_ggcyto_strategy.R:9-10` | ggcyto 1.34 fails with ggplot2 4.0.3 | ok, and reproduced today: `object '+.gg' not found` |
| `95:258-304`, `95:314-325` | ggcyto as data adapter, plain geoms; gtable strip surgery | ok |
| `theme.R:42-61`, `:66-97` | grob walk for pt sizes; `save_panel` checks | ok |
| `ROADMAP.md:25`, `:186` | rule 2; verdict-limit register row | ok |
| `INDEX.md:252` | unspiked max difference 3e-8 | ok (`spike-in-v2/baseline_check.csv` agrees: `2.98e-08` for c188 TNF) |
| grill line 38 | "`hexbin` ... installed through `pixi.toml`" | wrong: `hexbin` 1.28.6 is installed but not listed in `pixi.toml` (transitive) |
| R2.6 | ggh4x 0.3.1 `facet_grid2` has `remove_labels`, not `axis.labels` | ok (`formals()` checked) |
| Not checked | `INDEX.md:232`, `:320-335`; ggplot2 `NEWS.md:531, 801`; `72:46, 61, 76, 90` | not checked |
