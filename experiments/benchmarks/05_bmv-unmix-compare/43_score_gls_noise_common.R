# Purpose: score the noise-model GLS arms of 42_unmix_gls_noise.R against as188_joint and as102_var, with
#          Common-gates copy (plans/common_gates_plan.md, batch 3): the gates are the common SpectroFlo labels (clean3
#          rows) and the above-cut readout uses arcsinh 1.5 (G7a). Negatives for the rSD are below arcsinh 1.5 (G13).
#          the method of 39b_score_variant_count.R: on CD4T and CD8T cells defined by as188_joint's gates
#          (rows paired across arms), per donor, for BV510-A (TNF) and BV570-A (IFNg): Nguyen rSD of events
#          below arcsinh 1.5, share of events with |x| < 0.5 x as188_joint's rSD, share above arcsinh 1.5 (G7a).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/43_score_gls_noise_common.R
#            from the repo root, after both arms of 42_unmix_gls_noise.R finished.
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/gs/as188_joint        common-gate GatingSet [batch 0, 97d]
#          experiments/runs/2026-11-11_bmv-common-gated/manifest.csv          clean3 FCS per arm [97c]
#          experiments/bmv_test/compare_3D/data/unmix_cells/gls188_noise, gls188_noise_spill [42_unmix_gls_noise.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed   raw as188_joint
# Outputs: experiments/runs/2026-11-11_bmv-gls-noise-common/scores/scores.csv           one row per arm, donor, population, channel
#          experiments/runs/2026-11-11_bmv-gls-noise-common/scores/bv510_events.csv.gz  source data of the density figure
#          experiments/runs/2026-11-11_bmv-gls-noise-common/scores/fig_bv510_density.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-gls-noise-common/scores/fig_scores.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-gls-noise-common/scores/sessionInfo.txt

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
RUN_DIR <- "experiments/runs/2026-11-11_bmv-gls-noise-common/scores"
GATED <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0: gs/<arm> on SpectroFlo labels, manifest.csv
DATA <- "experiments/bmv_test/compare_3D/data"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1.5 # negatives (rSD, squeeze cut) are events below arcsinh 1.5, as the positivity cut (G13, user 2026-11-11; was 1)
POS_CUT <- 1.5 # positives are events above arcsinh 1.5 (G7a, user 2026-11-11; was 1)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (user, 2026-11-10)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A")
POPS <- c("CD4T", "CD8T") # gates of as188_joint
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
ARM_COLORS <- c(
  as188_joint = "grey30",
  as102_var = "#D55E00",
  gls188_noise = "#0072B2",
  gls188_noise_spill = "#CC79A7"
)
NEW_ARMS <- c("gls188_noise", "gls188_noise_spill")
stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_common_scores.R") # density_curves(), plot_*()

d188 <- file.path(DATA, "autoSpectral_v1.8.8")

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
      # GLS arms are raw unmixed files (mask the cleaned-out rows); the other arms are already cleaned
      e <- if (a %in% NEW_ARMS) e[keep[[s]], , drop = FALSE] else e
      if (nrow(e) != sum(keep[[s]])) {
        stop("unexpected row count in ", p, call. = FALSE)
      }
      e[, CHANNELS, drop = FALSE]
    })
  )
# identity check: masking the raw as188_joint file reproduces its cleaned file (tests the row mapping)
chk <- events |> filter(arm == "as188_joint", sample == "A3")
raw_ref <- read_x(file.path(raw_ref_dir, basename(chk$file)))[
  keep[["A3"]],
  CHANNELS
]
stopifnot(
  "row mapping reproduces cleaned file" = isTRUE(all.equal(
    unname(raw_ref),
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

# Figures: BV510 density and the three scores per arm (helper98_common_scores.R) --------------
curves <- long |>
  filter(channel == CHANNELS[["TNF"]]) |>
  transmute(arm, sample, pop, arcsinh = asinh(x / COFACTOR)) |>
  density_curves(from = -1.5, to = 4)
p_dens <- plot_density_curves(curves, ARM_COLORS, NEG_CUT)
save_panel(p_dens, curves, RUN_DIR, "fig_bv510_density", 183, 90)
scores_plot <- plot_arm_scores(scores, ARM_COLORS, POS_CUT)
save_panel(scores_plot$plot, scores_plot$frame, RUN_DIR, "fig_scores", 183, 160)
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
