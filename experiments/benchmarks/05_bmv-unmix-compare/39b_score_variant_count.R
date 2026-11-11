# Purpose: score the arms of 39_variant_count.R against as188_joint, swap_var102 (37) and s102_c188.
#          Same method as 37b_score_swap.R. Also tabulates variant counts per dye and the angle spread of
#          the BV510 variants around the 1.8.8 BV510 signature.
#          On CD4T and CD8T cells defined by as188_joint's gates (rows paired across arms), per donor, for
#          BV510-A (TNF) and BV570-A (IFNg): Nguyen rSD of events below arcsinh 1, share of events with
#          |x| < 0.5 x as188_joint's rSD, and share of events above arcsinh 1.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/39b_score_variant_count.R
#            from the repo root, after the four arms of 39_variant_count.R finished.
# Inputs:  experiments/runs/2026-11-10_bmv-gated/gs/as188_joint        GatingSet [19a_gate_arm.R]
#          experiments/runs/2026-11-10_bmv-gated/manifest.csv          cleaned FCS per arm [12_cleanup.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/swap_var102 [37_swap_af_variants.R], the four arms [39_variant_count.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed   raw as188_joint
# Outputs: experiments/runs/2026-11-10_bmv-variant-count/variant_counts.csv  variants per dye, 1.0.2 and 1.8.8
#          experiments/runs/2026-11-10_bmv-variant-count/bv510_angles.csv    angle of each BV510 variant to the 1.8.8 signature
#          experiments/runs/2026-11-10_bmv-variant-count/scores.csv           one row per arm, donor, population, channel
#          experiments/runs/2026-11-10_bmv-variant-count/bv510_events.csv.gz  source data of the density figure
#          experiments/runs/2026-11-10_bmv-variant-count/fig_bv510_density.pdf, .png
#          experiments/runs/2026-11-10_bmv-variant-count/fig_scores.pdf, .png
#          experiments/runs/2026-11-10_bmv-variant-count/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
RUN_DIR <- "experiments/runs/2026-11-10_bmv-variant-count"
GATED <- "experiments/runs/2026-11-10_bmv-gated"
DATA <- "experiments/bmv_test/compare_3D/data"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1 # negatives are events below arcsinh 1 (user, 2026-11-10)
POS_CUT <- 1 # positives are events above arcsinh 1 (user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (user, 2026-11-10)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A")
POPS <- c("CD4T", "CD8T") # gates of as188_joint
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
ARM_COLORS <- c(
  as188_joint = "grey30",
  s102_c188 = "#D55E00",
  swap_var102 = "#009E73",
  var102_n26 = "#CC79A7",
  var188_n100 = "#56B4E9",
  var102_bv510only = "#F0E442",
  var102_except_bv510 = "#999933"
)
NEW_ARMS <- c(
  "swap_var102",
  "var102_n26",
  "var188_n100",
  "var102_bv510only",
  "var102_except_bv510"
)
VARIANT_DYE <- "BV510" # dye whose variant angles are reported
stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
theme_set(ggpubr::theme_pubr(base_size = 11))
source("experiments/benchmarks/05_bmv-unmix-compare/helper41_plot-scores.R") # plot_bv510_density(), plot_scores()

# Variant counts per dye, and BV510 variant angles to the 1.8.8 signature --------------------
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
v188 <- readRDS(file.path(
  d188,
  "figure_spectral_variants/Spectral_variants.rds"
))
v102 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
))
spectra <- as.matrix(read.csv(
  file.path(
    d188,
    "table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  row.names = 1,
  check.names = FALSE
))
tibble(
  dye = names(v188$variants),
  n_188 = map_int(v188$variants, nrow),
  n_102 = map_int(v102$variants[dye], nrow)
) |>
  write_csv(file.path(RUN_DIR, "variant_counts.csv"))
angle_deg <- \(m, ref) {
  acos(pmin(1, as.vector(m %*% ref) / sqrt(rowSums(m^2) * sum(ref^2)))) *
    180 /
    pi
}
bind_rows(
  tibble(
    version = "1.8.8",
    angle = angle_deg(v188$variants[[VARIANT_DYE]], spectra[VARIANT_DYE, ])
  ),
  tibble(
    version = "1.0.2",
    angle = angle_deg(v102$variants[[VARIANT_DYE]], spectra[VARIANT_DYE, ])
  )
) |>
  write_csv(file.path(RUN_DIR, "bv510_angles.csv"))

# Rows kept by the cleaning, from raw as188_joint vs its cleaned file ---------------------------
manifest <- read_csv(
  file.path(GATED, "manifest.csv"),
  show_col_types = FALSE
) |>
  filter(arm %in% names(ARM_COLORS)) |>
  mutate(arm = factor(arm, levels = names(ARM_COLORS)))
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
raw_ref_dir <- file.path(DATA, "autoSpectral_v1.8.8/AutoSpectral_unmixed")
ref <- manifest |> filter(arm == "as188_joint")
keep <- map2(ref$file, ref$sample, \(clean_file, sample) {
  raw_key <- key_of(read_x(file.path(raw_ref_dir, basename(clean_file))))
  stopifnot("keys unique in raw file" = !anyDuplicated(raw_key))
  raw_key %in% key_of(read_x(clean_file))
}) |>
  set_names(ref$sample)

# Gates of as188_joint ------------------------------------------------------------------
gs <- load_gs(file.path(GATED, "gs/as188_joint"))
gate_rows <- imap(keep, \(k, sample) {
  gh <- gs[[sample]]
  map(POPS, \(p) gh_pop_get_indices(gh, p)) |> set_names(POPS)
})

# Values per arm, donor, population ---------------------------------------------------------
arm_file <- \(arm, sample, clean_file) {
  if (arm %in% NEW_ARMS) {
    file.path(DATA, "unmix_cells", arm, basename(clean_file))
  } else {
    clean_file
  }
}
events <- manifest |>
  bind_rows(
    expand_grid(arm = NEW_ARMS, sample = ref$sample) |>
      left_join(
        ref |> select(sample, file),
        by = "sample",
        relationship = "many-to-one"
      ) |>
      mutate(arm = factor(arm, levels = names(ARM_COLORS)))
  ) |>
  mutate(
    path = pmap_chr(list(as.character(arm), sample, file), arm_file),
    values = pmap(list(as.character(arm), path, sample), \(a, p, s) {
      e <- read_x(p)
      # swap arms are raw unmixed files (mask the cleaned-out rows); the other arms are already cleaned
      e <- if (a %in% NEW_ARMS) e[keep[[s]], , drop = FALSE] else e
      if (nrow(e) != sum(keep[[s]])) {
        stop("unexpected row count in ", p, call. = FALSE)
      }
      e[, CHANNELS, drop = FALSE]
    })
  )
# identity check: masking a raw crossed-arm file reproduces its cleaned file (tests the row mapping)
chk <- events |> filter(arm == "s102_c188", sample == "A3")
raw_s102 <- read_x(file.path(
  DATA,
  "unmix_cells/s102_c188",
  basename(chk$file)
))[keep[["A3"]], CHANNELS]
stopifnot(
  "row mapping reproduces cleaned file" = isTRUE(all.equal(
    unname(raw_s102),
    unname(chk$values[[1]])
  ))
)
stopifnot(
  "gate rows match cleaned rows" = all(map2_lgl(
    keep,
    gate_rows,
    \(k, g) all(lengths(g) == sum(k))
  ))
)

long <- events |>
  select(arm, sample, values) |>
  mutate(
    values = map2(values, sample, \(v, s) {
      map(POPS, \(p) {
        as_tibble(v[gate_rows[[s]][[p]], , drop = FALSE]) |>
          mutate(pop = p)
      }) |>
        list_rbind()
    })
  ) |>
  unnest(values) |>
  pivot_longer(
    all_of(unname(CHANNELS)),
    names_to = "channel",
    values_to = "x"
  ) |>
  mutate(
    pop = factor(pop, levels = POPS),
    channel = factor(channel, levels = CHANNELS)
  )

# Scores ------------------------------------------------------------------------------------
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84)))) # q84 - q50
scores <- long |>
  summarise(
    n = n(),
    rsd_neg = nguyen(x[asinh(x / COFACTOR) < NEG_CUT]),
    pos_share = mean(asinh(x / COFACTOR) > POS_CUT),
    .by = c(arm, sample, pop, channel)
  )
ref_rsd <- scores |>
  filter(arm == "as188_joint") |>
  select(sample, pop, channel, ref_rsd = rsd_neg)
squeezed <- long |>
  left_join(
    ref_rsd,
    by = c("sample", "pop", "channel"),
    relationship = "many-to-one"
  ) |>
  summarise(
    squeezed_share = mean(abs(x) < SQUEEZE_FRAC * ref_rsd[1]),
    .by = c(arm, sample, pop, channel)
  )
scores <- scores |>
  left_join(
    squeezed,
    by = c("arm", "sample", "pop", "channel"),
    relationship = "one-to-one"
  ) |>
  left_join(
    ref_rsd,
    by = c("sample", "pop", "channel"),
    relationship = "many-to-one"
  ) |>
  mutate(rsd_ratio = rsd_neg / ref_rsd)
write_csv(scores, file.path(RUN_DIR, "scores.csv"))
long |>
  filter(channel == CHANNELS[["TNF"]]) |>
  mutate(arcsinh = round(asinh(x / COFACTOR), 4)) |>
  select(arm, sample, pop, arcsinh) |>
  write_csv(file.path(RUN_DIR, "bv510_events.csv.gz"))

# Figures: BV510 density and the three scores per arm (helper41_plot-scores.R) ---------------
dens <- long |>
  filter(channel == CHANNELS[["TNF"]]) |>
  mutate(arcsinh = asinh(x / COFACTOR))
p_dens <- plot_bv510_density(dens, ARM_COLORS, NEG_CUT)
p_scores <- plot_scores(scores, ARM_COLORS)
ggsave(
  file.path(RUN_DIR, "fig_bv510_density.pdf"),
  p_dens,
  width = 220,
  height = 120,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  file.path(RUN_DIR, "fig_bv510_density.png"),
  p_dens,
  width = 220,
  height = 120,
  units = "mm",
  dpi = 200,
  bg = "white"
)
ggsave(
  file.path(RUN_DIR, "fig_scores.pdf"),
  p_scores,
  width = 200,
  height = 200,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  file.path(RUN_DIR, "fig_scores.png"),
  p_scores,
  width = 200,
  height = 200,
  units = "mm",
  dpi = 200,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
