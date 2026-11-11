# Copy of 64_af405_spike_anatomy.R for the common-gates spike (plans/common_gates_plan.md G5): reads the new 62cg key
#          and the new spiked unmixes (af405proof_cgspike188, af405proof_cgspike102af405), writes a new run dir.
#          Unspiked baselines af405proof_base188 and var102_only_Alexa_Fluor_405 are reused. Method unchanged.
# Purpose: proof chain steps 5 and 6 (user, 2026-11-10).
#          Step 5: recovery of a known BV510 (TNF) spike. Per spiked CD4T / CD8T cell, recovered = BV510 in
#          the spiked unmix minus BV510 of the same cell in the unspiked unmix of the same setup; recovery =
#          recovered / added. Setups: 1.8.8 (af405proof_spike188 vs af405proof_base188) and 1.8.8 with the
#          1.0.2 AF405 variants (af405proof_spike102af405 vs var102_only_Alexa_Fluor_405). Unspiked cells
#          must match their baseline exactly; GLUT1 = Alexa Fluor 405-A of the unspiked 1.8.8 unmix.
#          Step 6: for 10 CD8T cells (set.seed(42), user-approved) with TNF above arcsinh 1 under 1.8.8 and
#          squeezed under the 1.0.2 AF405 variants, rebuild the raw spectrum from each unmix (reference
#          spectra for the other dyes, AF row by AF Index x AF-A) and find which AF405 variant fits best.
# Run:     sbatch experiments/slurm/run_spike_cg.sbatch 64cg_af405_spike_anatomy.R
# Inputs:  experiments/runs/2026-11-11_bmv-af405-proof-cg/spike_key.csv   [62_af405_spike_build.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/af405proof_{base188,spike188,spike102af405} [61]
#          experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_Alexa_Fluor_405 [45_variant_per_dye.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw, table_spectra, Spectral_variants.rds
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/.../Spectral_variants.rds
# Outputs: experiments/runs/2026-11-11_bmv-af405-proof-cg/spike_cells.csv       per spiked cell: added, recovered
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/spike_summary.csv     recovery and squeezed share by level
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/fig5_spike_recovery.png (+ _source.csv)
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/fig5b_recovery_by_glut1.png (+ _source.csv)
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/predictions_spike.csv  pass / fail
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/anatomy_cells.csv     10 cells: abundances, best variant
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/fig6_anatomy.png (+ _source.csv)
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/sessionInfo_64.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # one whole donor for a smoke run; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else ""
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
RUN_DIR <- paste0(
  "experiments/runs/2026-11-11_bmv-af405-proof-cg",
  SMOKE_SUFFIX
) # 62cg wrote spike_key.csv here
UNMIX <- file.path(DATA, "unmix_cells")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1.5 # negatives for the rSD: BV510 below this arcsinh (G13, user 2026-11-11; was 1)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below 0.5 x the 1.8.8 rSD of BV510 negatives (39b; user, 2026-11-10)
RECOVERY_MIN <- 0.9 # pre-registered: 1.8.8 recovers >= 0.9 of the spike (user-approved plan)
N_DECILES <- 10
N_ANATOMY <- 10 # cells shown in step 6 (user-approved plan)
SEED <- 42 # user-approved
DYE <- "Alexa Fluor 405"
SETUPS <- tribble(
  ~setup                         , ~spiked                                             , ~baseline                     ,
  "1.8.8"                        , paste0("af405proof_", SMOKE_TAG, "cgspike188")      , "af405proof_base188"          ,
  "1.8.8 + 1.0.2 AF405 variants" , paste0("af405proof_", SMOKE_TAG, "cgspike102af405") , "var102_only_Alexa_Fluor_405"
)

# Inputs ----------------------------------------------------------------------------------
stopifnot(
  "64cg outputs must not exist" = !file.exists(file.path(
    RUN_DIR,
    "spike_cells.csv"
  ))
)
key <- read_csv(file.path(RUN_DIR, "spike_key.csv"), show_col_types = FALSE)
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
unmixed_name <- \(raw_file) sub("\\.fcs$", " AutoSpectral.fcs", raw_file)
files <- distinct(key, sample, raw_file)

# Step 5: per spiked cell, recovered BV510 and GLUT1 -----------------------------------------
cells <- pmap(files, \(sample, raw_file) {
  key_s <- filter(key, .data$sample == .env$sample)
  base188 <- read_x(file.path(
    UNMIX,
    "af405proof_base188",
    unmixed_name(raw_file)
  ))
  per_setup <- pmap(SETUPS, \(setup, spiked, baseline) {
    sp <- read_x(file.path(UNMIX, spiked, unmixed_name(raw_file)))[
      key_s$raw_row,
      "BV510-A"
    ]
    bl <- read_x(file.path(UNMIX, baseline, unmixed_name(raw_file)))[
      key_s$raw_row,
      "BV510-A"
    ]
    tibble(setup = setup, bv510_base = bl, bv510_spiked = sp)
  }) |>
    list_rbind()
  bind_cols(key_s[rep(seq_len(nrow(key_s)), nrow(SETUPS)), ], per_setup) |>
    mutate(
      glut1 = rep(
        asinh(base188[key_s$raw_row, "Alexa Fluor 405-A"] / COFACTOR),
        nrow(SETUPS)
      )
    )
}) |>
  list_rbind() |>
  mutate(recovered = bv510_spiked - bv510_base, recovery = recovered / k)

# unspiked cells must be unchanged by the spike, else the baseline is not the same unmix
unspiked_diff <- cells |>
  filter(level == 0) |>
  summarise(max_abs_diff = max(abs(recovered)), .by = setup)
write_csv(unspiked_diff, file.path(RUN_DIR, "spike_unspiked_check.csv"))
stopifnot(
  "unspiked cells equal their baseline unmix" = all(
    unspiked_diff$max_abs_diff < 1
  )
)

# 1.8.8 rSD of BV510 negatives per donor and population, for the squeezed share
ref_rsd <- cells |>
  filter(setup == "1.8.8", asinh(bv510_base / COFACTOR) < NEG_CUT) |>
  summarise(
    ref_rsd = unname(diff(quantile(bv510_base, c(0.5, 0.84)))),
    .by = c(sample, pop)
  )
cells <- cells |>
  left_join(ref_rsd, by = c("sample", "pop"), relationship = "many-to-one") |>
  mutate(
    squeezed = abs(bv510_spiked) < SQUEEZE_FRAC * ref_rsd,
    glut1_decile = ntile(glut1, N_DECILES),
    .by = c(sample, setup)
  )
write_csv(cells, file.path(RUN_DIR, "spike_cells.csv"))

spike_summary <- cells |>
  filter(level > 0) |>
  summarise(
    n = n(),
    median_recovery = median(recovery),
    q25 = quantile(recovery, 0.25),
    q75 = quantile(recovery, 0.75),
    squeezed_share = mean(squeezed),
    .by = c(setup, sample, pop, level)
  )
write_csv(spike_summary, file.path(RUN_DIR, "spike_summary.csv"))

# Figure 5: recovery per level and setup, donors as points --------------------------------------
fig5 <- spike_summary |>
  mutate(level = factor(level, labels = paste0("+", 1:3, " arcsinh")))
write_csv(fig5, file.path(RUN_DIR, "fig5_spike_recovery_source.csv"))
p5 <- ggplot(fig5, aes(setup, median_recovery)) +
  geom_hline(yintercept = c(0, 1), colour = "grey60") +
  geom_line(aes(group = sample, colour = sample), alpha = 0.5) +
  geom_point(aes(colour = sample), size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_grid(pop ~ level) +
  labs(x = NULL, y = "Median share of added BV510 recovered") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(RUN_DIR, "fig5_spike_recovery.png"),
  p5,
  width = 11,
  height = 6,
  dpi = 130,
  bg = "white"
)

# Figure 5b: recovery by GLUT1 decile, pooled over donors ---------------------------------------
fig5b <- cells |>
  filter(level > 0) |>
  summarise(
    median_recovery = median(recovery),
    .by = c(setup, pop, level, glut1_decile, sample)
  )
write_csv(fig5b, file.path(RUN_DIR, "fig5b_recovery_by_glut1_source.csv"))
p5b <- ggplot(fig5b, aes(glut1_decile, median_recovery, colour = setup)) +
  geom_hline(yintercept = c(0, 1), colour = "grey60") +
  geom_point(size = 1, alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = setup),
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c("1.8.8" = "#1b9e77", "1.8.8 + 1.0.2 AF405 variants" = "#d95f02"),
    name = NULL
  ) +
  scale_x_continuous(breaks = 1:N_DECILES) +
  facet_grid(
    pop ~ level,
    labeller = labeller(level = \(l) paste0("+", l, " arcsinh"))
  ) +
  labs(
    x = "GLUT1 (Alexa Fluor 405) decile within donor",
    y = "Median share recovered (points = donors)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(
  file.path(RUN_DIR, "fig5b_recovery_by_glut1.png"),
  p5b,
  width = 11,
  height = 6,
  dpi = 130,
  bg = "white"
)

# Pre-registered spike predictions --------------------------------------------------------------
rec_by_setup <- cells |>
  filter(level > 0) |>
  summarise(median_recovery = median(recovery), .by = setup)
glut_trend <- cells |>
  filter(level > 0, setup == "1.8.8 + 1.0.2 AF405 variants") |>
  summarise(r = median(recovery), .by = glut1_decile) |>
  summarise(rho = cor(glut1_decile, r, method = "spearman"))
rec_of <- \(s) rec_by_setup$median_recovery[rec_by_setup$setup == s]
predictions <- tibble(
  prediction = c(
    "1.8.8 median recovery >= 0.9",
    "1.0.2 AF405 variants recover less than 1.8.8",
    "under 1.0.2 AF405 variants, recovery falls with GLUT1 decile (Spearman rho < 0)"
  ),
  value = c(
    sprintf("%.3f", rec_of("1.8.8")),
    sprintf(
      "%.3f vs %.3f",
      rec_of("1.8.8 + 1.0.2 AF405 variants"),
      rec_of("1.8.8")
    ),
    sprintf("rho = %.2f", glut_trend$rho)
  ),
  pass = c(
    rec_of("1.8.8") >= RECOVERY_MIN,
    rec_of("1.8.8 + 1.0.2 AF405 variants") < rec_of("1.8.8"),
    glut_trend$rho < 0
  )
)
write_csv(predictions, file.path(RUN_DIR, "predictions_spike.csv"))

# Step 6: anatomy of 10 squeezed CD8T cells ------------------------------------------------------
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
s188 <- read_raw(file.path(
  d188,
  "table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
))
af <- read_raw(file.path(
  d188,
  "table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
))
v102 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
))$variants[[DYE]]
v188 <- readRDS(file.path(
  d188,
  "figure_spectral_variants/Spectral_variants.rds"
))$variants[[DYE]]
fluor_cols <- paste0(rownames(s188), "-A")
b510 <- s188["BV510", ]
cand <- cells |>
  filter(
    setup == "1.8.8 + 1.0.2 AF405 variants",
    pop == "CD8T",
    level == 0
  ) |>
  select(sample, raw_file, raw_row)
# TNF above 1 under 1.8.8 and squeezed under the 1.0.2 AF405 variants, unspiked cells only
anatomy_pool <- pmap(distinct(cand, sample, raw_file), \(sample, raw_file) {
  rows <- cand$raw_row[cand$sample == sample]
  x188 <- read_x(file.path(
    UNMIX,
    "af405proof_base188",
    unmixed_name(raw_file)
  ))[rows, , drop = FALSE]
  x102 <- read_x(file.path(
    UNMIX,
    "var102_only_Alexa_Fluor_405",
    unmixed_name(raw_file)
  ))[rows, , drop = FALSE]
  rsd <- ref_rsd$ref_rsd[ref_rsd$sample == sample & ref_rsd$pop == "CD8T"]
  ok <- asinh(x188[, "BV510-A"] / COFACTOR) > 1 &
    abs(x102[, "BV510-A"]) < SQUEEZE_FRAC * rsd
  tibble(sample = sample, raw_file = raw_file, raw_row = rows[ok])
}) |>
  list_rbind()
set.seed(SEED)
picked <- slice_sample(anatomy_pool, n = N_ANATOMY)
fit_one <- \(raw_vec, unmixed, variants) {
  a <- unmixed[fluor_cols]
  af_idx <- as.integer(round(unmixed[["AF Index"]])) # AF Index = row of the AF matrix (0 = no AF row)
  af_part <- if (af_idx >= 1) {
    unmixed[["AF-A"]] * af[af_idx, colnames(s188)]
  } else {
    0
  }
  rest <- colSums(
    a[names(a) != paste0(DYE, "-A")] *
      s188[rownames(s188) != DYE, , drop = FALSE]
  )
  resid_by_variant <- apply(variants, 1, \(v) {
    sum((raw_vec - rest - af_part - a[[paste0(DYE, "-A")]] * v)^2)
  })
  best <- which.min(resid_by_variant)
  list(
    recon = rest + af_part + a[[paste0(DYE, "-A")]] * variants[best, ],
    best = names(best),
    best_coef = sum((variants[best, ] - s188[DYE, ]) * b510) / sum(b510^2),
    bv510 = a[["BV510-A"]],
    af405 = a[[paste0(DYE, "-A")]]
  )
}
anatomy <- pmap(picked, \(sample, raw_file, raw_row) {
  raw_vec <- read_x(file.path(d188, "Raw", raw_file))[raw_row, colnames(s188)]
  u188 <- read_x(file.path(
    UNMIX,
    "af405proof_base188",
    unmixed_name(raw_file)
  ))[raw_row, ]
  u102 <- read_x(file.path(
    UNMIX,
    "var102_only_Alexa_Fluor_405",
    unmixed_name(raw_file)
  ))[raw_row, ]
  f188 <- fit_one(raw_vec, u188, v188)
  f102 <- fit_one(raw_vec, u102, v102)
  cell_id <- paste(sample, raw_row)
  list(
    table = tibble(
      cell = cell_id,
      bv510_188 = f188$bv510,
      bv510_102af405 = f102$bv510,
      af405_188 = f188$af405,
      af405_102af405 = f102$af405,
      best_variant_102 = f102$best,
      best_variant_bv510_coef = f102$best_coef,
      bv510_lost = f188$bv510 - f102$bv510
    ),
    spectra = tibble(
      cell = cell_id,
      detector = colnames(s188),
      raw = raw_vec,
      `fit 1.8.8` = f188$recon,
      `fit 1.0.2 AF405 variants` = f102$recon
    )
  )
})
anatomy_table <- map(anatomy, "table") |> list_rbind()
write_csv(anatomy_table, file.path(RUN_DIR, "anatomy_cells.csv"))
fig6 <- map(anatomy, "spectra") |>
  list_rbind() |>
  pivot_longer(
    c(raw, `fit 1.8.8`, `fit 1.0.2 AF405 variants`),
    names_to = "series",
    values_to = "value"
  ) |>
  mutate(detector = factor(detector, levels = colnames(s188)))
write_csv(fig6, file.path(RUN_DIR, "fig6_anatomy_source.csv"))
p6 <- ggplot(fig6, aes(detector, value, colour = series, group = series)) +
  geom_line(linewidth = 0.4) +
  scale_colour_manual(
    values = c(
      raw = "black",
      `fit 1.8.8` = "#1b9e77",
      `fit 1.0.2 AF405 variants` = "#d95f02"
    ),
    name = NULL
  ) +
  facet_wrap(~cell, ncol = 2, scales = "free_y") +
  labs(x = "Detector", y = "Raw signal") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
ggsave(
  file.path(RUN_DIR, "fig6_anatomy.png"),
  p6,
  width = 11,
  height = 12,
  dpi = 120,
  bg = "white"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_64.txt")
)
