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
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_spike_cg.sbatch 64cg_af405_spike_anatomy.R
#          Both variables are required. The dataset sets the data dir, the stamps, the BV510-slot channel, the AF405
#          dye, the negative cut (cut_of(..., "neg")) and the run folder.
#          Here "BV510" is the TNF-slot dye (cfg("targets")$tnf_dye) and "AF405" the GLUT1 dye (cfg("targets")$glut1_dye);
#          the column names bv510_* are kept for every dataset.
# Inputs:  run_in("af405-proof-cg")/spike_key.csv   [62cg_af405_spike_build.R; registry key of datasets/<ds>_runs.csv]
#          data_path("unmix_cells")/af405proof_{base188,cgspike188,cgspike102af405} [61cg]
#          data_path("unmix_cells")/var102_only_<GLUT1 dye> [45_variant_per_dye.R]
#            unmixed folders are found recursively and matched to samples_tbl() on the file stem (24H keeps one
#            folder per condition)
#          data_path("autoSpectral_v1.8.8")/Raw, table_spectra (cfg("stamps") names), Spectral_variants.rds
#          data_path("autoSpectral_v1.0.2")/.../Spectral_variants.rds
# Outputs: all into run_in("af405-proof-cg") (= experiments/runs/<date>_bmv<run_tag>-af405-proof-cg, written by 62cg;
#          "-smoke" is added to the key in a smoke run):
#          <run dir>/spike_cells.csv       per spiked cell: added, recovered
#          <run dir>/spike_summary.csv     recovery and squeezed share by level
#          <run dir>/fig5_spike_recovery.png (+ _source.csv)
#          <run dir>/fig5b_recovery_by_glut1.png (+ _source.csv)
#          <run dir>/predictions_spike.csv  pass / fail
#          <run dir>/anatomy_cells.csv     10 cells: abundances, best variant
#          <run dir>/fig6_anatomy.png (+ _source.csv)
#          <run dir>/sessionInfo_64.txt
#          Datasets with conditions (24H) add donor and condition columns to spike_cells.csv, spike_summary.csv and
#          the fig5 and fig5b source files, and facet figures 5 and 5b by condition.

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ------------------------------------------------------------------------------
DATA <- data_path()
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # one whole donor for a smoke run; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else ""
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
RUN_DIR <- run_in(paste0("af405-proof-cg", SMOKE_SUFFIX)) # 62cg wrote spike_key.csv here; 64cg adds its tables
UNMIX <- file.path(DATA, "unmix_cells")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_ARM <- "as188_joint" # arm whose BV510 negatives set the rSD; the cut per sample comes from cut_of(..., "neg")
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below 0.5 x the 1.8.8 rSD of BV510 negatives (39b; user, 2026-11-10)
RECOVERY_MIN <- 0.9 # pre-registered: 1.8.8 recovers >= 0.9 of the spike (user-approved plan)
N_DECILES <- 10
N_ANATOMY <- 10 # cells shown in step 6 (user-approved plan)
SEED <- 42 # user-approved
DYE <- cfg("targets")$glut1_dye
GLUT1_CH <- paste0(DYE, "-A")
TNF_DYE <- cfg("targets")$tnf_dye # the "BV510" of the column names below
TNF_CH <- cfg("channels")$TNF_slot
VAR102_DIR <- paste0("var102_only_", gsub(" ", "_", DYE))
UNMIXED_SUFFIX <- " AutoSpectral\\.fcs$"
SETUPS <- tribble(
  ~setup                         , ~spiked                                             , ~baseline            ,
  "1.8.8"                        , paste0("af405proof_", SMOKE_TAG, "cgspike188")      , "af405proof_base188" ,
  "1.8.8 + 1.0.2 AF405 variants" , paste0("af405proof_", SMOKE_TAG, "cgspike102af405") , VAR102_DIR
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
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
sample_info <- samples_tbl() |>
  filter(sample %in% key$sample)
# adds the file of each raw file found by its stem under dir (recursive: 24H has a folder per condition)
add_path <- \(files, dir, suffix, col) {
  found <- tibble(
    path = list.files(
      dir,
      pattern = suffix,
      recursive = TRUE,
      full.names = TRUE
    )
  ) |>
    mutate(stem = sub(suffix, "", basename(path))) |>
    filter(stem %in% files$stem)
  stopifnot(
    "64cg: one file per raw file under dir" = setequal(found$stem, files$stem)
  )
  files |>
    inner_join(
      found,
      by = "stem",
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
    rename(!!col := path)
}
# raw file -> registered sample (join on the full file stem) -> its raw file and every unmixed file used below
files <- distinct(key, sample, raw_file) |>
  mutate(stem = sub("\\.fcs$", "", basename(raw_file)))
stopifnot(
  "64cg: every key raw file is a registered sample (samples_tbl stem)" = all(
    files$stem %in% sample_info$stem
  )
)
files <- files |>
  inner_join(
    select(sample_info, sample, stem),
    by = c("sample", "stem"),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  add_path(file.path(d188, "Raw"), "\\.fcs$", "raw_path")
files <- reduce(
  unique(c("af405proof_base188", SETUPS$baseline, SETUPS$spiked)),
  \(f, d) add_path(f, file.path(UNMIX, d), UNMIXED_SUFFIX, d),
  .init = files
)
# path of one kind of file for one raw file; files holds exactly one (one-to-one joins above)
path_of <- \(raw, col) files[[col]][files$raw_file == raw]

# Step 5: per spiked cell, recovered BV510 and GLUT1 -----------------------------------------
cells <- pmap(distinct(files, sample, raw_file), \(sample, raw_file) {
  key_s <- filter(key, .data$sample == .env$sample)
  base188 <- read_x(path_of(raw_file, "af405proof_base188"))
  per_setup <- pmap(SETUPS, \(setup, spiked, baseline) {
    sp <- read_x(path_of(raw_file, spiked))[key_s$raw_row, TNF_CH]
    bl <- read_x(path_of(raw_file, baseline))[key_s$raw_row, TNF_CH]
    tibble(setup = setup, bv510_base = bl, bv510_spiked = sp)
  }) |>
    list_rbind()
  bind_cols(key_s[rep(seq_len(nrow(key_s)), nrow(SETUPS)), ], per_setup) |>
    mutate(
      glut1 = rep(
        asinh(base188[key_s$raw_row, GLUT1_CH] / COFACTOR),
        nrow(SETUPS)
      )
    )
}) |>
  list_rbind() |>
  mutate(recovered = bv510_spiked - bv510_base, recovery = recovered / k)
if (!is.null(DS$conditions)) {
  cells <- cells |>
    inner_join(
      select(sample_info, sample, donor, condition),
      by = "sample",
      relationship = "many-to-one",
      unmatched = "error"
    )
}

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
neg_cuts <- tibble(sample = unique(cells$sample)) |>
  mutate(neg_cut = map_dbl(sample, \(s) cut_of(NEG_ARM, s, TNF_CH, "neg")))
ref_rsd <- cells |>
  inner_join(
    neg_cuts,
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  filter(setup == "1.8.8", asinh(bv510_base / COFACTOR) < neg_cut) |>
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
    .by = c(setup, sample, any_of(c("donor", "condition")), pop, level)
  )
write_csv(spike_summary, file.path(RUN_DIR, "spike_summary.csv"))

# Figure 5: recovery per level and setup, donors as points --------------------------------------
fig5 <- spike_summary |>
  mutate(level = factor(level, labels = paste0("+", 1:3, " arcsinh")))
write_csv(fig5, file.path(RUN_DIR, "fig5_spike_recovery_source.csv"))
# colour by sample (one per donor) in 3D and EXV, by donor in 24H (3 samples per donor); Dark2 holds 8 colours
colour_by <- if (is.null(DS$conditions)) "sample" else "donor"
n_colours <- n_distinct(fig5[[colour_by]])
fig_rows <- if (is.null(DS$conditions)) vars(pop) else vars(pop, condition)
fig_height <- 6 * max(1, length(DS$conditions)) # one block of rows per condition
p5 <- ggplot(fig5, aes(setup, median_recovery)) +
  geom_hline(yintercept = c(0, 1), colour = "grey60") +
  geom_line(aes(group = sample, colour = .data[[colour_by]]), alpha = 0.5) +
  geom_point(aes(colour = .data[[colour_by]]), size = 1.8) +
  (if (n_colours <= 8) {
    scale_colour_brewer(palette = "Dark2", name = "Donor")
  } else {
    scale_colour_manual(
      values = hcl.colors(n_colours, "Dark 3"),
      name = "Donor"
    )
  }) +
  facet_grid(rows = fig_rows, cols = vars(level)) +
  labs(
    x = NULL,
    y = sprintf("Median share of added %s recovered", TNF_DYE)
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(RUN_DIR, "fig5_spike_recovery.png"),
  p5,
  width = 11,
  height = fig_height,
  dpi = 130,
  bg = "white"
)

# Figure 5b: recovery by GLUT1 decile, pooled over donors ---------------------------------------
fig5b <- cells |>
  filter(level > 0) |>
  summarise(
    median_recovery = median(recovery),
    .by = c(setup, pop, level, glut1_decile, sample, any_of("condition"))
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
    rows = fig_rows,
    cols = vars(level),
    labeller = labeller(level = \(l) paste0("+", l, " arcsinh"))
  ) +
  labs(
    x = sprintf("GLUT1 (%s) decile within donor", DYE),
    y = "Median share recovered (points = donors)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(
  file.path(RUN_DIR, "fig5b_recovery_by_glut1.png"),
  p5b,
  width = 11,
  height = fig_height,
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
s188 <- read_raw(file.path(d188, "table_spectra", cfg("stamps")$as188_spectra))
af <- read_raw(file.path(d188, "table_spectra", cfg("stamps")$as188_af))
v102 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
))$variants[[DYE]]
v188 <- readRDS(file.path(
  d188,
  "figure_spectral_variants/Spectral_variants.rds"
))$variants[[DYE]]
fluor_cols <- paste0(rownames(s188), "-A")
b510 <- s188[TNF_DYE, ]
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
  x188 <- read_x(path_of(raw_file, "af405proof_base188"))[rows, , drop = FALSE]
  x102 <- read_x(path_of(raw_file, VAR102_DIR))[rows, , drop = FALSE]
  rsd <- ref_rsd$ref_rsd[ref_rsd$sample == sample & ref_rsd$pop == "CD8T"]
  ok <- asinh(x188[, TNF_CH] / COFACTOR) > 1 &
    abs(x102[, TNF_CH]) < SQUEEZE_FRAC * rsd
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
    bv510 = a[[TNF_CH]],
    af405 = a[[paste0(DYE, "-A")]]
  )
}
anatomy <- pmap(picked, \(sample, raw_file, raw_row) {
  raw_vec <- read_x(path_of(raw_file, "raw_path"))[raw_row, colnames(s188)]
  u188 <- read_x(path_of(raw_file, "af405proof_base188"))[raw_row, ]
  u102 <- read_x(path_of(raw_file, VAR102_DIR))[raw_row, ]
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
