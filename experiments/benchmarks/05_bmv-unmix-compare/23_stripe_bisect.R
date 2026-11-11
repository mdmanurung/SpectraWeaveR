# Purpose: find which unmixing stage creates two artefacts in as188_joint donor D1: the stripe (very bright
#          Alexa Fluor 594, BV605, PE-Cy5 events pushed to IFNg ~ -4.4) and the diagonal TNF / CD8 negative streak.
#          Donor D1 raw file is unmixed in 3 stages with the as188_joint spectra, and the same event rows are compared.
#          stage 1 = OLS, no AF
#          stage 2 = AutoSpectral per-cell AF (the documented route to AF with af.spectra), no variants
#          stage 3 = AutoSpectral joint with AF and variants (= as188_joint)
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/23_stripe_bisect.R <run_dir>
#          from the repo root. <run_dir> is e.g. experiments/runs/2026-11-10_bmv-diagnostics (created if absent;
#          the subfolder 23_stripe_bisect inside it must not exist).
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw/D1 *.fcs                 raw sample [user zip]
#          .../Reference Group Clean (flow.control), table_spectra (spectra, AF), figure_spectral_variants [as 06_unmix_c188.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/repro_s188_c188_fast/D1 *.fcs         only to check stage 3
# Outputs: <run_dir>/23_stripe_bisect/stage<k>/*.fcs   unmixed D1 per stage
#          <run_dir>/23_stripe_bisect/stage_events.rds  arcsinh values of the plotted channels, one table, all stages
#          <run_dir>/23_stripe_bisect/stage_counts.csv  counts per stage
#          <run_dir>/23_stripe_bisect/stripe_events.csv stripe events per stage (IFNg and the three bright channels)
#          <run_dir>/23_stripe_bisect/check_vs_repro.csv  stage 3 against the repro file
#          <run_dir>/23_stripe_bisect/hexbin_ifng_tnf.png, hexbin_cd8_tnf.png   one figure per marker pair, one panel per stage
#          <run_dir>/23_stripe_bisect/args.txt, sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
DONOR <- "D1" # user, 2026-11-10
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
IFNG_BAND <- c(-4.8, -4.0) # IFNg arcsinh window of the stripe (user, 2026-11-10)
STREAK_CUT <- -2 # TNF and CD8 both below this arcsinh value = streak (user, 2026-11-10)
STRIPE_BRIGHT_CUT <- 6.5 # stripe events: Alexa Fluor 594, BV605 and PE-Cy5 all >= this in stage 3. Lower end of the user's "~6.5-6.9"
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
PIPELINE <- "joint" # same source
BINS <- 150 # hexagons per panel, display only
CH_IFNG <- "BV570-A"
CH_TNF <- "BV510-A"
CH_CD8 <- "SBUV445-A"
CH_BRIGHT <- c("Alexa Fluor 594-A", "BV605-A", "PE-Cy5-A")

# Setup -----------------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_dir <- file.path(run_dir, "23_stripe_bisect")
stopifnot("23 output folder must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
out_dir <- normalizePath(out_dir)
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
spectra_dir <- file.path(d188, "table_spectra")
spectra_file <- "Automated_autospectral_spectra_20261105_143803.csv"
af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
variants_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
raw_file <- list.files(
  file.path(d188, "Raw"),
  pattern = paste0("^", DONOR, " .*\\.fcs$"),
  full.names = TRUE
)
repro_file <- list.files(
  file.path(DATA, "unmix_cells", "repro_s188_c188_fast"),
  pattern = paste0("^", DONOR, " .*\\.fcs$"),
  full.names = TRUE
)
stopifnot(length(raw_file) == 1, length(repro_file) == 1)

asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in 06_unmix_c188.R
af.spectra <- read.spectra(af_file, spectra.dir = spectra_dir)
spectra.variants <- readRDS(variants_rds)

# Stages ----------------------------------------------------------------------------------
STAGES <- list(
  list(
    id = "1_OLS_noAF",
    method = "OLS",
    af = NULL,
    variants = NULL
  ),
  list(
    id = "2_AS_AF_novariants",
    method = "AutoSpectral",
    af = af.spectra,
    variants = NULL
  ),
  list(
    id = "3_AS_AF_variants",
    method = "AutoSpectral",
    af = af.spectra,
    variants = spectra.variants
  )
)
writeLines(
  c(
    paste("raw:", raw_file),
    paste("spectra:", file.path(spectra_dir, spectra_file)),
    paste("af:", file.path(spectra_dir, af_file)),
    paste("variants:", variants_rds),
    paste("control dir:", control_dir),
    paste("speed:", SPEED, " pipeline:", PIPELINE, " (AutoSpectral stages)"),
    paste(
      "stage 1: method OLS, no af.spectra, no variants (pipeline and speed are not used by OLS)"
    ),
    "stage 2: method AutoSpectral, af.spectra given, spectra.variants NULL",
    "stage 3: method AutoSpectral, af.spectra and spectra.variants given",
    paste(
      "stripe rule: IFNg in",
      paste(IFNG_BAND, collapse = " to "),
      "and all of",
      paste(CH_BRIGHT, collapse = ", "),
      ">=",
      STRIPE_BRIGHT_CUT,
      "in stage 3"
    )
  ),
  file.path(out_dir, "args.txt")
)

n_raw <- as.integer(flowCore::read.FCSheader(raw_file)[[1]][["$TOT"]])
read_stage <- function(f) {
  ff <- flowCore::read.FCS(f, truncate_max_range = FALSE)
  stopifnot("event count equals raw" = nrow(ff) == n_raw)
  exprs(ff)
}
tick <- function(x) asinh(x / COFACTOR)

# Unmix each stage and keep the plotted channels --------------------------------------------
ev <- list()
for (s in STAGES) {
  sdir <- file.path(out_dir, paste0("stage", substr(s$id, 1, 1)))
  dir.create(sdir)
  args <- list(
    fcs.file = raw_file,
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = s$method,
    output.dir = sdir,
    parallel = TRUE,
    verbose = TRUE
  )
  if (s$method == "AutoSpectral") {
    args <- c(
      args,
      list(
        af.spectra = s$af,
        spectra.variants = s$variants,
        speed = SPEED,
        pipeline = PIPELINE
      )
    )
  }
  t0 <- Sys.time()
  do.call(unmix.fcs, args)
  cat(
    s$id,
    "elapsed_min:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "\n"
  )
  f <- list.files(sdir, pattern = "\\.fcs$", full.names = TRUE)
  stopifnot("one output file per stage" = length(f) == 1)
  x <- read_stage(f)
  need <- c(CH_IFNG, CH_TNF, CH_CD8, CH_BRIGHT)
  stopifnot("plotted channels exist" = all(need %in% colnames(x)))
  ev[[s$id]] <- as_tibble(tick(x[, need, drop = FALSE])) |>
    mutate(row = row_number(), stage = s$id)
  if (s$id == "3_AS_AF_variants") {
    xr <- read_stage(repro_file)
    common <- intersect(colnames(x), colnames(xr))
    d <- abs(x[, common, drop = FALSE] - xr[, common, drop = FALSE])
    write.csv(
      tibble(
        n_events = nrow(x),
        n_common_cols = length(common),
        same_colnames = identical(colnames(x), colnames(xr)),
        max_abs_diff = max(d),
        n_cols_with_diff = sum(apply(d, 2, max) > 0)
      ),
      file.path(out_dir, "check_vs_repro.csv"),
      row.names = FALSE
    )
  }
}
events <- bind_rows(ev) |>
  mutate(stage = factor(stage, levels = vapply(STAGES, `[[`, "", "id")))
saveRDS(events, file.path(out_dir, "stage_events.rds"))

# Counts ----------------------------------------------------------------------------------
stripe_rows <- events |>
  filter(
    stage == "3_AS_AF_variants",
    .data[[CH_IFNG]] > IFNG_BAND[[1]],
    .data[[CH_IFNG]] < IFNG_BAND[[2]],
    if_all(all_of(CH_BRIGHT), \(v) v >= STRIPE_BRIGHT_CUT)
  ) |>
  pull(row)
events <- events |> mutate(is_stripe = row %in% stripe_rows)
counts <- events |>
  group_by(stage) |>
  summarise(
    n_events = n(),
    n_ifng_band = sum(
      .data[[CH_IFNG]] > IFNG_BAND[[1]] & .data[[CH_IFNG]] < IFNG_BAND[[2]]
    ),
    n_stripe_events = sum(is_stripe),
    stripe_ifng_median = median(.data[[CH_IFNG]][is_stripe]),
    stripe_ifng_min = min(.data[[CH_IFNG]][is_stripe]),
    stripe_ifng_max = max(.data[[CH_IFNG]][is_stripe]),
    n_stripe_in_band = sum(
      is_stripe &
        .data[[CH_IFNG]] > IFNG_BAND[[1]] &
        .data[[CH_IFNG]] < IFNG_BAND[[2]]
    ),
    n_streak = sum(.data[[CH_TNF]] < STREAK_CUT & .data[[CH_CD8]] < STREAK_CUT),
    .groups = "drop"
  )
write.csv(counts, file.path(out_dir, "stage_counts.csv"), row.names = FALSE)
write.csv(
  events |>
    filter(is_stripe) |>
    select(stage, row, all_of(c(CH_IFNG, CH_BRIGHT))),
  file.path(out_dir, "stripe_events.csv"),
  row.names = FALSE
)
print(as.data.frame(counts))

# Figures: one per marker pair, one panel per stage -----------------------------------------
lims <- range(unlist(events[c(CH_IFNG, CH_TNF, CH_CD8)]))
hex_fig <- function(ycol, ylab, file) {
  p <- ggplot(events, aes(x = .data[[CH_TNF]], y = .data[[ycol]])) +
    geom_hex(bins = BINS) +
    geom_point(
      data = \(d) filter(d, is_stripe),
      shape = 1,
      colour = "red",
      size = 1.2
    ) +
    scale_fill_viridis_c(trans = "log10", name = "events") +
    facet_wrap(vars(stage), nrow = 1) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(
      x = paste0("TNF (", CH_TNF, ") arcsinh(x/", COFACTOR, ")"),
      y = ylab,
      title = paste0(
        DONOR,
        ": all events, red circles = stripe events from stage 3"
      )
    ) +
    theme_bw()
  ggsave(file, p, width = 15, height = 5.5, dpi = 150)
}
hex_fig(
  CH_IFNG,
  paste0("IFNg (", CH_IFNG, ") arcsinh"),
  file.path(out_dir, "hexbin_ifng_tnf.png")
)
hex_fig(
  CH_CD8,
  paste0("CD8 (", CH_CD8, ") arcsinh"),
  file.path(out_dir, "hexbin_cd8_tnf.png")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
