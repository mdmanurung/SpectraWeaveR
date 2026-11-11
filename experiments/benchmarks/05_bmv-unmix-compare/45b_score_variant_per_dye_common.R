# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates and rows are the common SpectroFlo
#          labels (clean3), the squeeze cut is refitted as 0.5 x the negative rSD of as188_joint on those cells,
#          and the above-cut readout uses arcsinh 1.5 (G7a); negatives for the rSD stay below arcsinh 1.
#          Scores the per-dye variant arms of 45_variant_per_dye.R with the 39b metrics, against as188_joint
#          and var102_except_bv510 (39). On CD4T and CD8T cells defined by as188_joint's gates (rows paired
#          across arms), per donor, for BV510-A (TNF) and BV570-A (IFNg): Nguyen rSD of events below arcsinh 1
#          (ratio to as188_joint), share of events with |x| < 0.5 x as188_joint's rSD (squeezed), and share of
#          events above arcsinh 1.5. Draws one ranked dot plot of the CD8T TNF squeezed share per dye arm.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/45b_score_variant_per_dye_common.R
#            from the repo root, after the arms of 45_variant_per_dye.R finished.
# Inputs:  experiments/runs/2026-11-10_bmv-variant-per-dye/dyes.txt   dye per arm [written at submit; read only]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/as188_joint   common-gate GatingSet [batch 0, 97d]
#          experiments/runs/2026-11-11_bmv-common-gated/manifest.csv     clean3 FCS per arm [97c]
#          experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_<dye_id> [45_variant_per_dye.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/var102_except_bv510   [39_variant_count.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed   raw as188_joint
# Outputs: <RUN_DIR>/scores.csv                 one row per arm, donor, population, channel
#          <RUN_DIR>/fig_cd8t_tnf_squeezed.{pdf,png,csv}   ranked dot plot, CD8T TNF squeezed share per dye arm
#          <RUN_DIR>/sessionInfo.txt

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
RUN_DIR <- "experiments/runs/2026-11-11_bmv-variant-per-dye-common"
DYES_FILE <- "experiments/runs/2026-11-10_bmv-variant-per-dye/dyes.txt" # the 32 dyes of the old run, read only
GATED <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0: gs/<arm> on SpectroFlo labels, manifest.csv
DATA <- "experiments/bmv_test/compare_3D/data"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1.5 # negatives (rSD, squeeze cut) are events below arcsinh 1.5, as the positivity cut (G13, user 2026-11-11; was 1)
POS_CUT <- 1.5 # positives are events above arcsinh 1.5 (G7a, user 2026-11-11; was 1)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (user, 2026-11-10)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A")
POPS <- c("CD4T", "CD8T") # gates of as188_joint
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
REF_ARMS <- c(
  as188_joint = "1.8.8 (as188_joint)",
  var102_except_bv510 = "all but BV510 from 1.0.2"
)

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
dyes <- readLines(DYES_FILE)
dye_arms <- set_names(paste0("var102_only_", gsub("[ -]", "_", dyes)), dyes)

# Rows kept by the cleaning, from raw as188_joint vs its cleaned file ---------------------------
manifest <- read_csv(file.path(GATED, "manifest.csv"), show_col_types = FALSE)
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
gate_rows <- map(set_names(ref$sample), \(s) {
  map(set_names(POPS), \(p) gh_pop_get_indices(gs[[s]], p))
})
stopifnot(
  "gate rows match cleaned rows" = all(map2_lgl(
    keep,
    gate_rows,
    \(k, g) all(lengths(g) == sum(k))
  ))
)

# Events of one arm and donor on CD4T and CD8T, long ----------------------------------------
# as188_joint is read from its cleaned file; the unmixed arms from their raw files, masked to the cleaned rows
arm_events <- \(arm, sample, clean_file) {
  e <- if (arm == "as188_joint") {
    read_x(clean_file)
  } else {
    read_x(file.path(DATA, "unmix_cells", arm, basename(clean_file)))[
      keep[[sample]],
      ,
      drop = FALSE
    ]
  }
  stopifnot("row count equals cleaned rows" = nrow(e) == sum(keep[[sample]]))
  map(POPS, \(p) {
    as_tibble(e[gate_rows[[sample]][[p]], CHANNELS, drop = FALSE]) |>
      mutate(pop = p)
  }) |>
    list_rbind() |>
    pivot_longer(
      all_of(unname(CHANNELS)),
      names_to = "channel",
      values_to = "x"
    )
}
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84)))) # q84 - q50

# Reference rSD from as188_joint -----------------------------------------------------------
ref_events <- map2(ref$sample, ref$file, \(s, f) {
  arm_events("as188_joint", s, f) |> mutate(sample = s)
}) |>
  list_rbind()
ref_rsd <- ref_events |>
  summarise(
    ref_rsd = nguyen(x[asinh(x / COFACTOR) < NEG_CUT]),
    .by = c(sample, pop, channel)
  )

# Scores per arm, one arm at a time to keep memory low --------------------------------------
score_arm <- \(arm) {
  map2(ref$sample, ref$file, \(s, f) {
    ev <- if (arm == "as188_joint") {
      filter(ref_events, sample == s)
    } else {
      arm_events(arm, s, f) |> mutate(sample = s)
    }
    ev |>
      left_join(
        ref_rsd,
        by = c("sample", "pop", "channel"),
        relationship = "many-to-one"
      ) |>
      summarise(
        n = n(),
        rsd_neg = nguyen(x[asinh(x / COFACTOR) < NEG_CUT]),
        squeezed_share = mean(abs(x) < SQUEEZE_FRAC * ref_rsd[1]),
        pos_share = mean(asinh(x / COFACTOR) > POS_CUT),
        ref_rsd = ref_rsd[1],
        .by = c(sample, pop, channel)
      )
  }) |>
    list_rbind() |>
    mutate(arm = arm, rsd_ratio = rsd_neg / ref_rsd, .before = 1)
}
all_arms <- c(names(REF_ARMS), unname(dye_arms))
scores <- map(all_arms, score_arm) |>
  list_rbind() |>
  left_join(
    tibble(arm = unname(dye_arms), dye = names(dye_arms)),
    by = "arm",
    relationship = "many-to-one"
  )
write_csv(scores, file.path(RUN_DIR, "scores.csv"))

# Figure: CD8T TNF squeezed share per dye arm, ranked by donor median -----------------------
plot_frame <- scores |>
  filter(pop == "CD8T", channel == CHANNELS[["TNF"]], !is.na(dye)) |>
  mutate(median_share = median(squeezed_share), .by = dye) |>
  mutate(dye = reorder(dye, median_share))
ref_lines <- scores |>
  filter(
    pop == "CD8T",
    channel == CHANNELS[["TNF"]],
    arm %in% names(REF_ARMS)
  ) |>
  summarise(median_share = median(squeezed_share), .by = arm) |>
  mutate(label = REF_ARMS[arm])
p <- ggplot(plot_frame, aes(squeezed_share, dye)) +
  geom_vline(
    data = ref_lines,
    aes(xintercept = median_share, linetype = label),
    colour = "grey40",
    linewidth = 0.3
  ) +
  geom_point(aes(colour = sample), size = 1) +
  geom_point(aes(x = median_share), shape = 3, size = 1.5, stroke = 0.4) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_linetype_manual(
    values = c("dashed", "dotted"),
    name = "Reference (median)"
  ) +
  labs(x = "CD8T TNF share squeezed near 0", y = "Dye given its 1.0.2 variants")
save_panel(p, plot_frame, RUN_DIR, "fig_cd8t_tnf_squeezed", 183, 160)
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
