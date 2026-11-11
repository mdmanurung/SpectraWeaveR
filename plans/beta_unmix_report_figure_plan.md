# Figure plan: beta unmixing report (approved)

Superseded 2026-11-11 by `plans/unmix_diagnostics_plan_v2.md`: F1 to F5 are built (`R/unmix_figures.R`, `R/unmix_report.R`) and extended by the Nature theme, F6 to F8 and the dashboard tabs. Kept as the record of the approved wording.

Status: **F1, F2, F3, F4 approved by the user on 2026-11-08. F5 stays a table.** The user stated the question and claim in their own words on 2026-11-09. The F3/F4 event sample is 50,000 events, seed 1, as confirmed on 2026-11-09.

Question the report answers (user, 2026-11-09): "Are the reference spectra and unmixed data good enough, or is correction needed?"

Claim (user, 2026-11-09): "I can evaluate unmixing results and determine whether correction is needed."

Every figure below cites its source table, so the figure and the number can be checked against each other. A table replaces a figure only where a reason is stated.

| # | Figure | Question it answers | Source data (saved by the run) | Visual check | Plot or table? |
|---|---|---|---|---|---|
| F1 | Spectra overlay, one panel per fluorophore, detectors on x, peak-normalised, ours vs spectracle vs SpectroFlo where available | Do our reference spectra match the comparators in shape? | `figures/panels/source_data/spectra-overlay.csv` (long form) | Peaks sit at the same detector; angle in the panel title matches `spectra_compare.csv` | plot |
| F2 | Cosine-similarity matrix of reference spectra (fluorophores × fluorophores), diverging-free sequential scale 0 to 1 | Which pairs are close enough to cause spillover trouble? | `checks.csv` rows `max_cosine`, plus the full matrix `figures/panels/source_data/similarity-matrix.csv` | Diagonal is 1; known close pairs (for example BV421 vs a violet neighbour) stand out | plot |
| F3 | Per-channel distribution of unmixed values for one stained file: ours vs SpectroFlo, one density per channel, symmetric log axis | Is the negative population centred near zero, and is the spread similar? | `figures/panels/source_data/unmixed-distributions.csv` (the sampled events behind each density) | Negative peak near 0 in both; spread compared on the same axis | plot |
| F4 | Biplot of two markers (the vignette's pair, for example BUV395 vs BV421) ours vs SpectroFlo, same axes | Do populations separate the same way? | `figures/panels/source_data/biplot.csv` | Same gate regions occupied; no diagonal smear present in one and not the other | plot |
| F5 | Condition number and AF library size | Is the unmixing matrix well conditioned? | `checks.csv` | One number each | **table** (a single value per run does not need a plot) |

Notes for approval:
- F3 and F4 need a sample of events written by the run (`figures/panels/source_data/unmixed-distributions.csv` and `biplot.csv`, fixed seed, size set by you). The sample size is an impactful choice and is a required argument.
- F1 to F4 use a colour-blind-safe discrete palette keyed by method (ours, spectracle, SpectroFlo); an unknown method name is an error.
- Each panel gets the number it shows in its title or caption (angle, n events); no number is typed by hand.
- Nothing is drawn for a dataset until its source table exists.

Decisions: rows F1 to F4 are plots, F5 is a table. Question and claim are the user's own words.
