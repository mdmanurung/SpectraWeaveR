# Purpose: find which cells get squeezed to 0 in BV510 (TNF) when Alexa Fluor 405 (GLUT1) uses its
#          AutoSpectral 1.0.2 variants. Pairs the same events in as188_joint and var102_only_Alexa_Fluor_405
#          (45_variant_per_dye.R) inside as188_joint's CD4T and CD8T gates, and relates the BV510 change to
#          the cell's GLUT1 (Alexa Fluor 405) abundance and autofluorescence (AF-A), and to whether GLUT1
#          exceeds the 1.8.8 AF405 positivity threshold (per-cell variant optimisation active).
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/56_af405_cells.R \
#            experiments/runs/2026-11-10_bmv-af405-mechanism/cells
#          (from the repo root; run dir must not exist)
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv, gs/as188_joint [19a_gate_arm.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed  raw as188_joint
#          experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_Alexa_Fluor_405/*.fcs [45]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds
# Outputs: <run_dir>/cell_summary.csv          per donor, pop, GLUT1 decile: share newly squeezed, median
#                                              BV510 change, median AF-A
#          <run_dir>/squeezed_vs_not.csv       per donor and pop: medians of GLUT1, AF-A, BV510 for cells that
#                                              become squeezed vs the rest, and share above the AF405 threshold
#          <run_dir>/fig_bv510_change_vs_glut1.png  hexbin, CD8T, all donors (+ _source.csv)
#          <run_dir>/fig_squeezed_by_glut1_decile.png  donors as points (+ _source.csv)
#          <run_dir>/sessionInfo.txt

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
DATA <- "experiments/bmv_test/compare_3D/data"
POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled" # cleaned files and as188_joint GatingSet
SWAP_ARM <- "var102_only_Alexa_Fluor_405"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1 # negatives below arcsinh 1 define the reference rSD (39b; user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
POPS <- c("CD4T", "CD8T")
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
CHANNELS <- c(TNF = "BV510-A", GLUT1 = "Alexa Fluor 405-A", AF = "AF-A")

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
theme_set(ggpubr::theme_pubr(base_size = 11))

# Rows: raw as188_joint file -> cleaned rows -> CD4T / CD8T gate indices ----------------------------
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
ref <- read_csv(file.path(POOLED, "manifest.csv"), show_col_types = FALSE) |>
  filter(arm == "as188_joint")
raw_ref_dir <- file.path(DATA, "autoSpectral_v1.8.8/AutoSpectral_unmixed")
gs <- load_gs(file.path(POOLED, "gs/as188_joint"))
af405_threshold <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds"
))$thresholds[["Alexa Fluor 405"]]

# Paired events per donor and population ----------------------------------------------------------
cells <- map2(ref$file, ref$sample, \(clean_file, donor) {
  raw_ref <- read_x(file.path(raw_ref_dir, basename(clean_file)))
  raw_swap <- read_x(file.path(
    DATA,
    "unmix_cells",
    SWAP_ARM,
    basename(clean_file)
  ))
  stopifnot(
    "paired files share event order" = identical(
      key_of(raw_ref),
      key_of(raw_swap)
    )
  )
  keep <- key_of(raw_ref) %in% key_of(read_x(clean_file))
  gh <- gs[[donor]]
  map(POPS, \(p) {
    idx <- gh_pop_get_indices(gh, p)
    stopifnot("gate rows match cleaned rows" = length(idx) == sum(keep))
    r <- raw_ref[keep, CHANNELS][idx, , drop = FALSE]
    s <- raw_swap[keep, CHANNELS][idx, , drop = FALSE]
    tibble(
      sample = donor,
      pop = p,
      tnf_ref = r[, "BV510-A"],
      tnf_swap = s[, "BV510-A"],
      glut1 = r[, "Alexa Fluor 405-A"],
      glut1_swap = s[, "Alexa Fluor 405-A"],
      af = r[, "AF-A"]
    )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(pop = factor(pop, levels = POPS))

# Newly squeezed cells: outside the squeeze band in as188_joint, inside it with the 1.0.2 variants -----
ref_rsd <- cells |>
  filter(asinh(tnf_ref / COFACTOR) < NEG_CUT) |>
  summarise(
    ref_rsd = unname(diff(quantile(tnf_ref, c(0.5, 0.84)))),
    .by = c(sample, pop)
  )
cells <- cells |>
  left_join(ref_rsd, by = c("sample", "pop"), relationship = "many-to-one") |>
  mutate(
    squeezed_ref = abs(tnf_ref) < SQUEEZE_FRAC * ref_rsd,
    squeezed_swap = abs(tnf_swap) < SQUEEZE_FRAC * ref_rsd,
    newly_squeezed = squeezed_swap & !squeezed_ref,
    tnf_change = asinh(tnf_swap / COFACTOR) - asinh(tnf_ref / COFACTOR),
    glut1_above_thr = glut1 > af405_threshold,
    glut1_decile = ntile(glut1, 10),
    .by = c(sample, pop)
  )
cell_summary <- cells |>
  summarise(
    n = n(),
    share_newly_squeezed = mean(newly_squeezed),
    median_tnf_change = median(tnf_change),
    median_glut1 = median(glut1),
    median_af = median(af),
    share_glut1_above_thr = mean(glut1_above_thr),
    .by = c(sample, pop, glut1_decile)
  ) |>
  arrange(sample, pop, glut1_decile)
write_csv(cell_summary, file.path(run_dir, "cell_summary.csv"))
cells |>
  summarise(
    n = n(),
    median_glut1 = median(glut1),
    median_glut1_swap = median(glut1_swap),
    median_af = median(af),
    median_tnf_ref = median(tnf_ref),
    share_glut1_above_thr = mean(glut1_above_thr),
    .by = c(sample, pop, newly_squeezed)
  ) |>
  arrange(sample, pop, newly_squeezed) |>
  write_csv(file.path(run_dir, "squeezed_vs_not.csv"))

# Figure: BV510 change vs GLUT1 abundance, CD8T ---------------------------------------------
hex_frame <- cells |>
  filter(pop == "CD8T") |>
  transmute(sample, glut1_arcsinh = asinh(glut1 / COFACTOR), tnf_change)
write_csv(
  hex_frame,
  file.path(run_dir, "fig_bv510_change_vs_glut1_source.csv.gz")
)
p_hex <- ggplot(hex_frame, aes(glut1_arcsinh, tnf_change)) +
  geom_hex(bins = 80) +
  geom_vline(
    xintercept = asinh(af405_threshold / COFACTOR),
    linetype = "dashed"
  ) +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_fill_viridis_c(trans = "log10", name = "Cells") +
  facet_wrap(~sample, nrow = 1) +
  labs(
    x = "GLUT1 (Alexa Fluor 405), arcsinh, as188_joint",
    y = "TNF change, arcsinh (1.0.2 AF405 variants - 1.8.8)"
  )
ggsave(
  file.path(run_dir, "fig_bv510_change_vs_glut1.png"),
  p_hex,
  width = 280,
  height = 90,
  units = "mm",
  dpi = 150,
  bg = "white"
)

# Figure: share newly squeezed by GLUT1 decile, donors as points --------------------------------
write_csv(
  cell_summary,
  file.path(run_dir, "fig_squeezed_by_glut1_decile_source.csv")
)
p_dec <- ggplot(
  cell_summary,
  aes(factor(glut1_decile), share_newly_squeezed, colour = sample)
) +
  geom_line(aes(group = sample), alpha = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~pop) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  labs(
    x = "GLUT1 decile within donor and population (1 = lowest)",
    y = "Share newly squeezed in TNF"
  )
ggsave(
  file.path(run_dir, "fig_squeezed_by_glut1_decile.png"),
  p_dec,
  width = 220,
  height = 90,
  units = "mm",
  dpi = 150,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
