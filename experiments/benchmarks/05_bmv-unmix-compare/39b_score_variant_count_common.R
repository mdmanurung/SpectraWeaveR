# Purpose: score the arms of 39_variant_count.R against as188_joint, swap_var102 (37) and s102_c188.
#          Common-gates copy (plans/common_gates_plan.md, batch 3): the gates are the common SpectroFlo labels (clean3
#          rows) and the above-cut readout uses arcsinh 1.5 (G7a). Negatives for the rSD are below arcsinh 1.5 (G13).
#          Same method as 37b_score_swap.R. Also tabulates variant counts per dye and the angle spread of
#          the BV510 variants around the 1.8.8 BV510 signature.
#          On CD4T and CD8T cells defined by as188_joint's gates (rows paired across arms), per donor, for
#          BV510-A (TNF) and BV570-A (IFNg): Nguyen rSD of events below arcsinh 1.5, share of events with
#          |x| < 0.5 x as188_joint's rSD, and share of events above arcsinh 1.5 (G7a).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/39b_score_variant_count_common.R
#            from the repo root, after the four arms of 39_variant_count.R finished.
# Inputs:  run_in("common-gated")/gs/as188_joint        common-gate GatingSet [batch 0, 97d]
#          run_in("common-gated")/manifest.csv          clean3 FCS per arm [97c]
#          data_path("unmix_cells", <arm>, [<condition>,] <file>) for swap_var102 [37_swap_af_variants.R] and the four arms [39_variant_count.R]
#          data_path("autoSpectral_v1.8.8/AutoSpectral_unmixed", [<condition>,] <file>)   raw as188_joint
#          data_path("autoSpectral_v1.8.8/table_spectra", cfg("stamps")$as188_spectra), the variants .rds of 1.8.8 and 1.0.2
#          datasets/<SW_DATASET>.json (stamps, channels, targets, conditions, pos_rule, neg_rule), <ds>_samples.csv
#          The cuts come from cut_of(); one value, or per sample when they differ (3D: 1.5; helper98_common_scores.R labels both).
# Outputs: run_out("variant-count-common")/variant_counts.csv  variants per dye, 1.0.2 and 1.8.8
#          run_out("variant-count-common")/bv510_angles.csv    angle of each TNF-slot dye variant to the 1.8.8 signature
#          run_out("variant-count-common")/scores.csv           one row per arm, sample, population, channel
#          run_out("variant-count-common")/bv510_events.csv.gz  source data of the density figure
#          run_out("variant-count-common")/fig_bv510_density.{pdf,png,csv}
#          run_out("variant-count-common")/fig_scores.{pdf,png,csv}
#          run_out("variant-count-common")/sessionInfo.txt

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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ------------------------------------------------------------------------------
RUN_DIR <- run_out("variant-count-common")
GATED <- run_in("common-gated") # batch 0: gs/<arm> on SpectroFlo labels, manifest.csv
DATA <- data_path()
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# negatives (rSD, squeeze cut) are events below cut_of(arm, sample, channel, "neg"), positives above the "pos" cut
# (3D: arcsinh 1.5 for both, G13 and G7a, user 2026-11-11; was 1)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (user, 2026-11-10)
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
CHANNELS <- set_names(
  c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot),
  c(TNF_LABEL, IFNG_LABEL)
)
SAMPLES <- samples_tbl()
# 24H keeps arm and raw files in one subfolder per condition; 3D and EXV have none
cond_path <- \(dir, sample, file) {
  cond <- if (is.null(DS$conditions)) {
    NULL
  } else {
    SAMPLES$condition[SAMPLES$sample == sample]
  }
  do.call(file.path, as.list(c(dir, cond, file)))
}
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
VARIANT_DYE <- cfg("targets")$tnf_dye # dye whose variant angles are reported
stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_common_scores.R") # density_curves(), plot_*()

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
  file.path(d188, "table_spectra", cfg("stamps")$as188_spectra),
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
stopifnot(
  "manifest holds the samples of samples_tbl()" = setequal(
    ref$sample,
    SAMPLES$sample
  )
)
CHECK_SAMPLE <- sort(ref$sample)[[1]] # sample of the row-mapping check (3D: A3)

# Cuts per arm, sample and channel, read once (cut_of re-reads cuts.csv when the rule is not fixed)
cut_tbl <- tidyr::expand_grid(
  arm = names(ARM_COLORS),
  sample = ref$sample,
  channel = unname(CHANNELS)
) |>
  mutate(
    neg_cut = pmap_dbl(list(arm, sample, channel), \(a, s, ch) {
      cut_of(a, s, ch, "neg")
    }),
    pos_cut = pmap_dbl(list(arm, sample, channel), \(a, s, ch) {
      cut_of(a, s, ch, "pos")
    })
  )
cut_here <- \(a, s, ch, col) {
  v <- cut_tbl[[col]][
    cut_tbl$arm == a[1] & cut_tbl$sample == s[1] & cut_tbl$channel == ch[1]
  ]
  stopifnot("one cut per arm, sample and channel" = length(v) == 1)
  v
}
# helper98_common_scores.R takes one value or a per-sample cut; one value when all samples share it
PLOT_NEG_CUT <- cut_tbl |>
  filter(arm == "as188_joint", channel == CHANNELS[[TNF_LABEL]]) |>
  transmute(sample, cut = neg_cut)
if (n_distinct(PLOT_NEG_CUT$cut) == 1) {
  PLOT_NEG_CUT <- PLOT_NEG_CUT$cut[[1]]
}
PLOT_POS_CUT <- cut_tbl$pos_cut
keep <- map2(ref$file, ref$sample, \(clean_file, sample) {
  raw_key <- key_of(read_x(cond_path(
    raw_ref_dir,
    sample,
    basename(clean_file)
  )))
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
    cond_path(
      file.path(DATA, "unmix_cells", arm),
      sample,
      basename(clean_file)
    )
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
chk <- events |> filter(arm == "s102_c188", sample == CHECK_SAMPLE)
raw_s102 <- read_x(cond_path(
  file.path(DATA, "unmix_cells/s102_c188"),
  CHECK_SAMPLE,
  basename(chk$file)
))[keep[[CHECK_SAMPLE]], CHANNELS]
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
    rsd_neg = nguyen(x[
      asinh(x / COFACTOR) < cut_here(arm, sample, channel, "neg_cut")
    ]),
    pos_share = mean(
      asinh(x / COFACTOR) > cut_here(arm, sample, channel, "pos_cut")
    ),
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
  filter(channel == CHANNELS[[TNF_LABEL]]) |>
  mutate(arcsinh = round(asinh(x / COFACTOR), 4)) |>
  select(arm, sample, pop, arcsinh) |>
  write_csv(file.path(RUN_DIR, "bv510_events.csv.gz"))

# Figures: BV510 density and the three scores per arm (helper98_common_scores.R) --------------
curves <- long |>
  filter(channel == CHANNELS[[TNF_LABEL]]) |>
  transmute(arm, sample, pop, arcsinh = asinh(x / COFACTOR)) |>
  density_curves(from = -1.5, to = 4)
p_dens <- plot_density_curves(curves, ARM_COLORS, PLOT_NEG_CUT)
save_panel(p_dens, curves, RUN_DIR, "fig_bv510_density", 183, 90)
scores_plot <- plot_arm_scores(scores, ARM_COLORS, PLOT_POS_CUT)
save_panel(scores_plot$plot, scores_plot$frame, RUN_DIR, "fig_scores", 183, 160)
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
