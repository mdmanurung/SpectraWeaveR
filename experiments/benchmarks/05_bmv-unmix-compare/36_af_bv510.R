# Purpose: test whether the autofluorescence (AF) spectra explain why BV510 (TNF) looks different between
#          unmixing arms. Arms on the AutoSpectral 1.0.2 spectra set push many CD4 T / CD8 T cells to ~0 in
#          BV510-A; arms on the 1.8.8 set do not. Hypothesis: the AF spectra travelling with the 1.0.2 set
#          absorb BV510 signal when each cell is assigned an AF spectrum.
#          1. cosine of every AF spectrum to the BV510 signature, and its 3 most similar fluor signatures
#          2. cells by AF Index per arm and donor: n, median BV510, fraction near zero
#          3. paired view against the reference arm, and AF-A abundance against BV510
#          4. BV510 content of each AF spectrum (least squares on the fluor signatures) times AF-A
#          Nothing is re-unmixed. Evidence only; no method ranking.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/36_af_bv510.R <run_dir>
#          from the repo root. run_dir must not exist; use experiments/runs/2026-11-10_bmv-af-bv510.
# Inputs:  experiments/runs/2026-11-10_bmv-gated/manifest.csv      unmixed FCS per arm and donor [01/19a]
#          experiments/runs/2026-11-10_bmv-gated/gs/as188_joint    GatingSet with CD4T, CD8T [19a_gate_arm.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/*.csv
#            spectra (Matching Negative, AF clean up) and the zip's AF set [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/*.csv
#            run 143803 spectra and its AF set 143933 [user zip]
#          experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/af_spectra.rds
#            refined AF sets of best102_neg, best102_af, best188_p1 [07_best_c102.R, 07_best_c188.R]
# Outputs: <run_dir>/af_spectra_cosines.csv        one row per AF spectrum: cosine to BV510, top 3 fluors, BV510 content
#          <run_dir>/bv510_signature_angle.csv     angle between the BV510 signatures of the spectra sets
#          <run_dir>/reference_cut.csv             reference rSD and the near-zero cut (linear BV510)
#          <run_dir>/key_af_indices.csv            AF Indices holding a large share of an arm's CD4T cells
#          <run_dir>/arm_pop_summary.csv           per arm, population, donor: near-zero and exact-zero fractions, rSD
#          <run_dir>/af_index_by_donor.csv         per arm, population, donor, AF Index: n, median BV510, near-zero fraction
#          <run_dir>/af_index_pooled.csv           same, donors pooled, with the AF spectrum's cosine and share of cells
#          <run_dir>/near_zero_by_afa_decile.csv   fig5 top panel source
#          <run_dir>/near_zero_by_cos_bin.csv      fig5 bottom panel source
#          <run_dir>/absorption_summary.csv        fig6 source: predicted absorbed BV510 vs observed change
#          <run_dir>/fig1_af_spectra.png           AF spectra with the BV510 signature, one panel per AF set
#          <run_dir>/fig1b_cosine_by_set.png       cosine of each AF spectrum to BV510; heavily used indices labelled
#          <run_dir>/fig2_bv510_by_af_index.png    BV510 distribution by AF Index, one panel per arm
#          <run_dir>/fig3_paired_bv510.png         BV510 in the reference arm vs three arms, split by AF Index
#          <run_dir>/fig4_afa_vs_bv510.png         AF-A abundance vs BV510 per arm
#          <run_dir>/fig5_near_zero_by_afa_cos.png near-zero fraction by AF-A decile and by AF cosine to BV510
#          <run_dir>/fig6_absorption.png           BV510 content of AF spectra; predicted absorption vs observed change
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggridges)
})
theme_set(ggpubr::theme_pubr(base_size = 11))

# Parameters ------------------------------------------------------------------------------
REF_ARM <- "as188_joint" # reference arm (README, approved choices); defines CD4T/CD8T membership and the rSD
POPS <- c("CD4T", "CD8T") # populations scored in the tables
FIG_POP <- "CD4T" # population shown in the figures (task brief); CD8T numbers are in the tables
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md)
RSD_CEIL <- 1 # Nguyen rSD uses reference CD4T events below this arcsinh value (task brief)
NEAR_ZERO_FRACTION <- 0.5 # near zero = |BV510| < this x reference rSD (task brief)
EXACT_ZERO_TOL <- 1e-6 # |BV510| below this counts as exactly zero (clamped value), linear scale
PAIRED_ARMS <- c("as102_var", "s102_c188", "s188_c102") # arms in the paired view (task brief)
TOP_INDEX <- 4 # AF Indices drawn per arm in figs 2 and 3 (others pooled); keeps panels readable
KEY_SHARE <- 0.05 # an AF Index is "heavily used" when it holds at least this share of an arm's CD4T cells
HEX_BINS <- 50 # hexagon bins per axis; display only
N_DECILES <- 10 # AF-A groups per arm for fig 5; deciles of the pooled CD4T cells
COS_BIN_WIDTH <- 0.1 # cosine bin width for fig 5
MIN_CELLS_BIN <- 30 # a donor-bin point is drawn only with at least this many cells; fewer is noise
XLIM <- c(-1.5, 3.5) # arcsinh display range for BV510
NEG_CLIP <- -0.3 # fig 1 y-axis floor; AF spectra dipping below it (row max = 1) are clipped (min_value in af_spectra_cosines.csv)
JITTER_SEED <- 42 # jitter only moves points sideways; same seed gives the same picture
DATA_DIR <- "experiments/bmv_test/compare_3D/data"
GATED_DIR <- "experiments/runs/2026-11-10_bmv-gated"

# Spectra sets and arms (settings from 06_unmix_c188.R, 06_unmix_c102.R, 07_best_*.R args.txt) --
# spectra_set = which spectra csv the arm unmixed with; af_set = which AF matrix its AF Index points into
SPECTRA_FILES <- c(
  sp102_neg = "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv",
  sp102_afclean = "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells AF clean up_autospectral_spectra.csv",
  sp188 = "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
)
AF_SETS <- tribble(
  ~af_set                 , ~spectra_set    , ~af_path                                                                                                                                       , ~label                           ,
  "af188"                 , "sp188"         , "autoSpectral_v1.8.8/table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"                   , "1.8.8 zip AF (143933)"          ,
  "af188_refined"         , "sp188"         , "unmix_cells/best188_p1/af_spectra.rds"                                                                                                        , "1.8.8 refined AF (best188)"     ,
  "af102"                 , "sp102_neg"     , "autoSpectral_v1.0.2/Reference Group/table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv" , "1.0.2 zip AF"                   ,
  "af102_neg_refined"     , "sp102_neg"     , "unmix_cells/best102_neg/af_spectra.rds"                                                                                                       , "1.0.2 refined AF (best102_neg)" ,
  "af102_afclean_refined" , "sp102_afclean" , "unmix_cells/best102_af/af_spectra.rds"                                                                                                        , "1.0.2 refined AF (best102_af)"
)
ARMS <- tribble(
  ~arm             , ~af_set                 , ~source     ,
  "as102_basic"    , "af102"                 , "1.0.2 set" ,
  "as102_var"      , "af102"                 , "1.0.2 set" ,
  "best102_neg"    , "af102_neg_refined"     , "1.0.2 set" ,
  "best102_af"     , "af102_afclean_refined" , "1.0.2 set" ,
  "s102_c188"      , "af102"                 , "1.0.2 set" ,
  "as188_joint"    , "af188"                 , "1.8.8 set" ,
  "as188_fast_csv" , "af188"                 , "1.8.8 set" ,
  "best188_p1"     , "af188_refined"         , "1.8.8 set" ,
  "best188_p5"     , "af188_refined"         , "1.8.8 set" ,
  "s188_c102"      , "af188"                 , "1.8.8 set"
)
PALETTE <- c(`1.0.2 set` = "#C0392B", `1.8.8 set` = "#2E6DA4")

# Output folder ---------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run_dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Spectra and AF matrices (raw csv or rds as the arms read them; read.spectra() only rescales rows)
read_mat <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
SPECTRA <- SPECTRA_FILES |>
  map(\(f) read_mat(file.path(DATA_DIR, f))) |>
  map(\(m) m[rownames(m) != "AF", , drop = FALSE]) # as in ReRunUnmix.R and 06_unmix_*.R
AF <- AF_SETS$af_path |>
  set_names(AF_SETS$af_set) |>
  map(\(f) {
    if (str_detect(f, fixed(".rds"))) {
      readRDS(file.path(DATA_DIR, f))
    } else {
      read_mat(file.path(DATA_DIR, f))
    }
  })
stopifnot(
  "detector order is shared" = all(map_lgl(
    c(SPECTRA, AF),
    \(m) identical(colnames(m), colnames(SPECTRA$sp188))
  )),
  "best188_p5 reused the best188_p1 AF set" = identical(
    readRDS(file.path(DATA_DIR, "unmix_cells/best188_p5/af_spectra.rds")),
    AF$af188_refined
  )
)

# Test 1: cosine of each AF spectrum to BV510 and its most similar fluor signatures -------
# cosine ignores row scale, so the raw csv and the rescaled rows give the same values
unit_rows <- \(m) m / sqrt(rowSums(m^2))
rescale_max <- \(m) m / apply(m, 1, max) # what read.spectra() does; used for the BV510 content below
af_cos <- AF_SETS |>
  select(af_set, spectra_set) |>
  pmap(\(af_set, spectra_set) {
    af <- AF[[af_set]]
    sp <- SPECTRA[[spectra_set]]
    cos_fl <- unit_rows(af) %*% t(unit_rows(sp)) # AF spectra x fluor signatures
    # least squares of each AF spectrum on the 33 fluor signatures: BV510 coefficient (test 4)
    coefs <- qr.coef(qr(t(rescale_max(sp))), t(rescale_max(af)))
    tibble(
      af_set = af_set,
      af_index = seq_len(nrow(af)),
      cos_bv510 = cos_fl[, "BV510"],
      top3_similar = apply(cos_fl, 1, \(x) {
        o <- order(x, decreasing = TRUE)[1:3]
        str_c(colnames(cos_fl)[o], " (", round(x[o], 3), ")", collapse = "; ")
      }),
      bv510_content = coefs["BV510", ],
      min_value = apply(af, 1, min)
    )
  }) |>
  list_rbind()
write_csv(af_cos, file.path(run_dir, "af_spectra_cosines.csv"))

bv_sig <- SPECTRA |> map(\(m) m["BV510", ])
angle_deg <- \(a, b) acos(sum(a * b) / sqrt(sum(a^2) * sum(b^2))) * 180 / pi
tibble(
  pair = c(
    "sp102_neg vs sp188",
    "sp102_afclean vs sp188",
    "sp102_neg vs sp102_afclean"
  ),
  angle_deg = c(
    angle_deg(bv_sig$sp102_neg, bv_sig$sp188),
    angle_deg(bv_sig$sp102_afclean, bv_sig$sp188),
    angle_deg(bv_sig$sp102_neg, bv_sig$sp102_afclean)
  )
) |>
  write_csv(file.path(run_dir, "bv510_signature_angle.csv"))

# Cells: one row per CD4T/CD8T event of the reference arm, for every arm -------------------
gs <- load_gs(file.path(GATED_DIR, "gs", REF_ARM))
samples <- sampleNames(gs)
ref_pop <- samples |>
  set_names() |>
  map(\(s) {
    is_cd4 <- gh_pop_get_indices(gs[[s]], "CD4T")
    is_cd8 <- gh_pop_get_indices(gs[[s]], "CD8T")
    stopifnot("CD4T and CD8T are disjoint" = !any(is_cd4 & is_cd8))
    case_when(is_cd4 ~ "CD4T", is_cd8 ~ "CD8T", .default = "other")
  })
ref_time <- samples |>
  set_names() |>
  map(\(s) unname(exprs(gh_pop_get_data(gs[[s]], "root"))[, "Time"]))

manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
) |>
  filter(arm %in% ARMS$arm)
stopifnot(
  "every arm has every donor once" = nrow(manifest) ==
    nrow(ARMS) * length(samples) &&
    all(table(manifest$arm, manifest$sample) == 1)
)

cells <- manifest |>
  select(file_arm = arm, file_sample = sample, file) |>
  pmap(\(file_arm, file_sample, file) {
    ff <- read.FCS(file, truncate_max_range = FALSE)
    ex <- exprs(ff)
    stopifnot(
      "event count equals the reference gs" = nrow(ex) ==
        length(ref_pop[[file_sample]]),
      "Time identical to the reference (rows are paired)" = identical(
        unname(ex[, "Time"]),
        ref_time[[file_sample]]
      )
    )
    tibble(
      arm = file_arm,
      sample = file_sample,
      event = seq_len(nrow(ex)),
      pop = ref_pop[[file_sample]],
      bv510 = ex[, "BV510-A"],
      afa = ex[, "AF-A"],
      af_index = as.integer(ex[, "AF Index"])
    ) |>
      filter(pop %in% POPS)
  }) |>
  list_rbind() |>
  left_join(
    ARMS,
    join_by(arm),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  left_join(
    af_cos |> select(af_set, af_index, cos_bv510, bv510_content, min_value),
    join_by(af_set, af_index),
    relationship = "many-to-one"
  ) |>
  mutate(
    bv_asinh = asinh(bv510 / COFACTOR),
    afa_asinh = asinh(afa / COFACTOR),
    arm = factor(arm, levels = ARMS$arm)
  )
stopifnot(
  "no missing values in cells (NA cos_bv510 = AF Index outside its arm's AF matrix)" = !anyNA(
    cells |> select(bv510, afa, af_index, cos_bv510)
  )
)

# Reference rSD and near-zero cut ---------------------------------------------------------
ref_neg <- cells |>
  filter(arm == REF_ARM, pop == "CD4T", bv_asinh < RSD_CEIL) |>
  pull(bv510)
RSD_REF <- unname(quantile(ref_neg, 0.84) - quantile(ref_neg, 0.5))
ZERO_CUT <- NEAR_ZERO_FRACTION * RSD_REF
cells <- cells |>
  mutate(
    near_zero = abs(bv510) < ZERO_CUT,
    exact_zero = abs(bv510) < EXACT_ZERO_TOL
  )
tibble(
  ref_arm = REF_ARM,
  rsd_ref = RSD_REF,
  near_zero_cut = ZERO_CUT,
  n_ref_negatives = length(ref_neg)
) |>
  write_csv(file.path(run_dir, "reference_cut.csv"))

# Arm summary: per arm, population, donor -------------------------------------------------
rsd_below <- \(x, ax) {
  unname(quantile(x[ax < RSD_CEIL], 0.84) - quantile(x[ax < RSD_CEIL], 0.5))
}
arm_pop <- cells |>
  summarise(
    n = n(),
    frac_near_zero = mean(near_zero),
    frac_exact_zero = mean(exact_zero),
    rsd_own = rsd_below(bv510, bv_asinh),
    median_bv510 = median(bv510),
    .by = c(arm, source, pop, sample)
  ) |>
  arrange(pop, arm, sample)
write_csv(arm_pop, file.path(run_dir, "arm_pop_summary.csv"))

# Test 2: cells by AF Index per arm, population and donor ---------------------------------
af_by_donor <- cells |>
  summarise(
    n = n(),
    median_bv510 = median(bv510),
    frac_near_zero = mean(near_zero),
    .by = c(arm, pop, sample, af_index)
  ) |>
  arrange(pop, arm, sample, desc(n))
write_csv(af_by_donor, file.path(run_dir, "af_index_by_donor.csv"))

af_pooled <- cells |>
  summarise(
    n = n(),
    median_bv510 = median(bv510),
    frac_near_zero = mean(near_zero),
    frac_exact_zero = mean(exact_zero),
    median_afa = median(afa),
    cos_bv510 = first(cos_bv510),
    bv510_content = first(bv510_content),
    min_value = first(min_value),
    .by = c(arm, af_set, pop, af_index)
  ) |>
  mutate(share = n / sum(n), .by = c(arm, pop)) |>
  arrange(pop, arm, desc(n))
write_csv(af_pooled, file.path(run_dir, "af_index_pooled.csv"))

# Top AF Indices per arm (CD4T), used by figs 1b, 2 and 3
top_idx <- af_pooled |>
  filter(pop == FIG_POP) |>
  mutate(rank = row_number(), .by = arm) |>
  filter(rank <= TOP_INDEX) |>
  mutate(
    class = str_c("#", rank),
    label = sprintf(
      "AF%d  cos %.2f  %.0f%% of cells  near0 %.0f%%",
      af_index,
      cos_bv510,
      100 * share,
      100 * frac_near_zero
    )
  )
class_levels <- c(str_c("#", seq_len(TOP_INDEX)), "other")
cells_fig <- cells |>
  filter(pop == FIG_POP) |>
  left_join(
    top_idx |> select(arm, af_index, class),
    join_by(arm, af_index),
    relationship = "many-to-one"
  ) |>
  mutate(class = factor(replace_na(class, "other"), levels = class_levels))

# Figure 1: AF spectra with the BV510 signature, one panel per AF set ---------------------
set_label <- AF_SETS |>
  mutate(
    panel = sprintf("%s (n = %d)", label, map_int(af_set, \(s) nrow(AF[[s]])))
  )
af_top3 <- af_cos |>
  mutate(top3 = rank(-cos_bv510, ties.method = "first") <= 3, .by = af_set) |>
  select(af_set, af_index, top3)
fig1_lines <- AF_SETS |>
  pull(af_set) |>
  map(\(s) {
    as_tibble(AF[[s]], rownames = "spectrum") |>
      mutate(af_set = s, af_index = row_number()) |>
      pivot_longer(
        -c(spectrum, af_set, af_index),
        names_to = "detector",
        values_to = "value"
      )
  }) |>
  list_rbind() |>
  left_join(
    af_top3,
    join_by(af_set, af_index),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  left_join(
    set_label |> select(af_set, panel, spectra_set),
    join_by(af_set),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(
    detector = factor(detector, levels = colnames(SPECTRA$sp188)),
    panel = factor(panel, levels = set_label$panel)
  )
fig1_bv <- set_label |>
  mutate(
    sig = map(spectra_set, \(s) {
      tibble(detector = names(bv_sig[[s]]), value = unname(bv_sig[[s]]))
    })
  ) |>
  select(panel, sig) |>
  unnest(sig) |>
  mutate(
    detector = factor(detector, levels = colnames(SPECTRA$sp188)),
    panel = factor(panel, levels = set_label$panel)
  )
p1 <- ggplot(fig1_lines, aes(detector, value, group = af_index)) +
  geom_line(
    data = \(d) filter(d, !top3),
    colour = "grey70",
    linewidth = 0.15,
    alpha = 0.5
  ) +
  geom_line(
    data = \(d) filter(d, top3),
    aes(colour = "3 AF closest to BV510"),
    linewidth = 0.6
  ) +
  geom_line(
    data = fig1_bv,
    aes(colour = "BV510", group = 1),
    linewidth = 0.9
  ) +
  scale_colour_manual(
    values = c(`3 AF closest to BV510` = "#E69F00", BV510 = "#C0392B"),
    name = NULL
  ) +
  # every 8th detector named; 64 names do not fit
  scale_x_discrete(
    breaks = colnames(SPECTRA$sp188)[seq(1, ncol(SPECTRA$sp188), 8)]
  ) +
  coord_cartesian(ylim = c(NEG_CLIP, 1.05)) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "Detector", y = "Signal (row max = 1)") +
  theme(
    axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1),
    legend.position = "top"
  )
ggsave(
  plot = p1,
  filename = file.path(run_dir, "fig1_af_spectra.png"),
  width = 180,
  height = 170,
  units = "mm",
  dpi = 150
)

# Figure 1b: cosine to BV510 per AF spectrum; heavily used indices labelled ---------------
key_idx <- af_pooled |>
  filter(pop == FIG_POP, share >= KEY_SHARE) |>
  summarise(
    arms = str_c(sprintf("%s %.0f%%", arm, 100 * share), collapse = ", "),
    .by = c(af_set, af_index)
  ) |>
  inner_join(
    af_cos |> select(af_set, af_index, cos_bv510),
    join_by(af_set, af_index),
    relationship = "one-to-one",
    unmatched = c("error", "drop")
  ) |>
  mutate(tag = sprintf("AF%d", af_index)) |>
  inner_join(
    set_label |> select(af_set, label),
    join_by(af_set),
    relationship = "many-to-one",
    unmatched = c("error", "drop")
  )
fig1b_dat <- af_cos |>
  left_join(
    set_label |> select(af_set, label),
    join_by(af_set),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(label = factor(label, levels = rev(set_label$label)))
key_legend <- sprintf("AF used by >= %.0f%% of an arm's CD4T", 100 * KEY_SHARE)
p1b <- ggplot(fig1b_dat, aes(cos_bv510, label)) +
  geom_point(
    colour = "grey60",
    size = 0.6,
    alpha = 0.5,
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED)
  ) +
  geom_point(
    data = key_idx |>
      mutate(label = factor(label, levels = rev(set_label$label))),
    aes(colour = key_legend),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = key_idx |>
      mutate(label = factor(label, levels = rev(set_label$label))),
    aes(label = tag),
    size = 3,
    min.segment.length = 0,
    seed = JITTER_SEED,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = set_names("#C0392B", key_legend), name = NULL) +
  labs(x = "Cosine of AF spectrum to BV510", y = NULL) +
  scale_y_discrete(labels = scales::label_wrap(22))
ggsave(
  plot = p1b,
  filename = file.path(run_dir, "fig1b_cosine_by_set.png"),
  width = 150,
  height = 100,
  units = "mm",
  dpi = 150
)
write_csv(key_idx, file.path(run_dir, "key_af_indices.csv"))

# Figure 2: BV510 distribution by AF Index, one panel per arm -----------------------------
fig2_dat <- cells_fig |>
  mutate(class = factor(class, levels = rev(class_levels)))
p2 <- ggplot(fig2_dat, aes(bv_asinh, class, fill = source)) +
  geom_density_ridges(scale = 1.1, rel_min_height = 0.005, linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed") +
  scale_fill_manual(values = PALETTE, name = "Spectra/AF set") +
  coord_cartesian(xlim = XLIM) +
  facet_wrap(~arm, ncol = 3) +
  labs(
    x = sprintf("BV510-A in %s, arcsinh", FIG_POP),
    y = "AF Index rank in arm"
  )
ggsave(
  plot = p2,
  filename = file.path(run_dir, "fig2_bv510_by_af_index.png"),
  width = 200,
  height = 220,
  units = "mm",
  dpi = 150
)
write_csv(top_idx, file.path(run_dir, "fig2_source_top_indices.csv"))

# Figure 3: paired view, reference BV510 (x) vs arm BV510 (y), split by AF Index class ----
ref_bv <- cells_fig |>
  filter(arm == REF_ARM) |>
  select(sample, event, bv_ref = bv_asinh)
paired <- cells_fig |>
  filter(arm %in% PAIRED_ARMS) |>
  left_join(
    ref_bv,
    join_by(sample, event),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(arm = factor(as.character(arm), levels = PAIRED_ARMS))
fig3_ann <- top_idx |>
  filter(arm %in% PAIRED_ARMS) |>
  mutate(
    arm = factor(as.character(arm), levels = PAIRED_ARMS),
    class = factor(class, levels = class_levels),
    short = sprintf("AF%d", af_index)
  )
# column strips name the AF rank in the row arm (rank 1 = most cells)
rank_labels <- set_names(
  c(str_c("AF rank ", seq_len(TOP_INDEX)), "other AF"),
  class_levels
)
p3 <- ggplot(paired, aes(bv_ref, bv_asinh)) +
  geom_hex(bins = HEX_BINS) +
  geom_abline(linewidth = 0.2, colour = "grey30") +
  geom_text(
    data = fig3_ann,
    aes(x = XLIM[1], y = XLIM[2], label = short),
    hjust = 0,
    vjust = 1,
    size = 3.2,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(trans = "log10", name = "Cells") +
  coord_equal(xlim = XLIM, ylim = XLIM) +
  facet_grid(arm ~ class, labeller = labeller(class = rank_labels)) +
  labs(
    x = sprintf("BV510-A in %s, arcsinh", REF_ARM),
    y = "BV510-A in arm, arcsinh"
  )
ggsave(
  plot = p3,
  filename = file.path(run_dir, "fig3_paired_bv510.png"),
  width = 230,
  height = 160,
  units = "mm",
  dpi = 150
)
fig3_src <- paired |>
  summarise(
    n = n(),
    frac_near_zero = mean(near_zero),
    median_ref_bv510 = median(bv_ref),
    median_arm_bv510 = median(bv_asinh),
    .by = c(arm, class)
  )
write_csv(fig3_src, file.path(run_dir, "fig3_source_class_summary.csv"))

# Figure 4: AF-A abundance vs BV510 --------------------------------------------------------
fig4_dat <- cells_fig |>
  filter(arm %in% c(REF_ARM, PAIRED_ARMS)) |>
  mutate(arm = factor(as.character(arm), levels = c(REF_ARM, PAIRED_ARMS)))
p4 <- ggplot(fig4_dat, aes(afa_asinh, bv_asinh)) +
  geom_hex(bins = HEX_BINS) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", name = "Cells") +
  coord_cartesian(xlim = c(-0.5, 3.5), ylim = XLIM) +
  facet_wrap(~arm, ncol = 2) +
  labs(
    x = sprintf("AF-A abundance, arcsinh(x / %g)", COFACTOR),
    y = "BV510-A, arcsinh"
  )
ggsave(
  plot = p4,
  filename = file.path(run_dir, "fig4_afa_vs_bv510.png"),
  width = 150,
  height = 130,
  units = "mm",
  dpi = 150
)

# Figure 5: near-zero fraction by AF-A decile and by cosine of the assigned AF spectrum -----
by_decile <- cells_fig |>
  mutate(afa_decile = ntile(afa, N_DECILES), .by = arm) |>
  summarise(
    n = n(),
    frac_near_zero = mean(near_zero),
    median_afa = median(afa),
    .by = c(arm, source, sample, afa_decile)
  )
write_csv(by_decile, file.path(run_dir, "near_zero_by_afa_decile.csv"))
by_cos <- cells_fig |>
  mutate(
    cos_bin = floor(cos_bv510 / COS_BIN_WIDTH) *
      COS_BIN_WIDTH +
      COS_BIN_WIDTH / 2
  ) |>
  summarise(
    n = n(),
    frac_near_zero = mean(near_zero),
    .by = c(arm, source, sample, cos_bin)
  ) |>
  filter(n >= MIN_CELLS_BIN)
write_csv(by_cos, file.path(run_dir, "near_zero_by_cos_bin.csv"))
panel5 <- \(dat, xvar, xlab, brks) {
  ggplot(dat, aes({{ xvar }}, frac_near_zero, colour = source)) +
    geom_point(size = 1, alpha = 0.6) +
    stat_summary(
      aes(group = source),
      fun = mean,
      geom = "line",
      linewidth = 0.4
    ) +
    scale_colour_manual(values = PALETTE, name = "Spectra/AF set") +
    scale_x_continuous(breaks = brks) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    facet_wrap(~arm, nrow = 2) +
    labs(
      x = xlab,
      y = sprintf("%s near 0", FIG_POP)
    )
}
p5 <- cowplot::plot_grid(
  panel5(by_decile, afa_decile, "AF-A decile in arm", c(1, 5, 10)),
  panel5(by_cos, cos_bin, "Cosine of assigned AF to BV510", seq(0, 1, 0.2)),
  ncol = 1,
  labels = c("a", "b"),
  label_size = 14
)
ggsave(
  plot = p5,
  filename = file.path(run_dir, "fig5_near_zero_by_afa_cos.png"),
  width = 220,
  height = 190,
  units = "mm",
  dpi = 150
)

# Test 4 / Figure 6: BV510 content of the assigned AF spectrum x AF-A vs observed change ---
# bv510_content = BV510 coefficient when the AF spectrum (row max 1) is fitted on the 33 fluor signatures,
# i.e. the BV510 signal the fluor-only fit would read from that AF shape per unit AF-A.
absorb <- paired |>
  mutate(
    predicted = afa * bv510_content,
    observed = bv510 - sinh(bv_ref) * COFACTOR
  )
absorb_sum <- absorb |>
  summarise(
    n = n(),
    spearman = cor(predicted, observed, method = "spearman"),
    median_predicted = median(predicted),
    median_observed = median(observed),
    median_abs_predicted = median(abs(predicted)),
    .by = arm
  )
write_csv(absorb_sum, file.path(run_dir, "absorption_summary.csv"))
p6a <- ggplot(
  af_cos |>
    left_join(
      set_label |> select(af_set, label),
      join_by(af_set),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    mutate(label = factor(label, levels = rev(set_label$label))),
  aes(bv510_content, label)
) +
  geom_point(
    colour = "grey60",
    size = 0.6,
    alpha = 0.5,
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED)
  ) +
  geom_point(
    data = key_idx |>
      inner_join(
        af_cos |> select(af_set, af_index, bv510_content),
        join_by(af_set, af_index),
        relationship = "one-to-one",
        unmatched = c("error", "drop")
      ) |>
      mutate(label = factor(label, levels = rev(set_label$label))),
    aes(colour = key_legend),
    size = 2
  ) +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed") +
  scale_colour_manual(values = set_names("#C0392B", key_legend), name = NULL) +
  labs(x = "BV510 content of AF spectrum (per unit AF-A)", y = NULL) +
  coord_cartesian(xlim = c(-0.5, 1.2)) +
  scale_y_discrete(labels = scales::label_wrap(22)) +
  theme(plot.margin = margin(5, 5, 5, 25)) # room for the set names after cowplot
p6b <- ggplot(
  absorb,
  aes(asinh(predicted / COFACTOR), asinh(observed / COFACTOR))
) +
  geom_hex(bins = HEX_BINS) +
  geom_abline(linewidth = 0.2) +
  scale_fill_viridis_c(trans = "log10", name = "Cells") +
  coord_equal(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
  facet_wrap(~arm, nrow = 1) +
  labs(
    x = "Predicted BV510 from AF, arcsinh",
    y = sprintf("BV510 change vs %s, arcsinh", REF_ARM)
  )
p6 <- cowplot::plot_grid(
  p6a,
  p6b,
  ncol = 1,
  rel_heights = c(1, 1.1),
  labels = c("a", "b"),
  label_size = 14
)
ggsave(
  plot = p6,
  filename = file.path(run_dir, "fig6_absorption.png"),
  width = 180,
  height = 190,
  units = "mm",
  dpi = 150
)

# Session ---------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
