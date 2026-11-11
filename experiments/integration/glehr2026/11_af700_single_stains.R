# experiments/integration/glehr2026/11_af700_single_stains.R
# Purpose: find which part of the stain widens the AF700 negative on the Fortessa (10 showed the
#          empty detector accounts for only part). From the CD3 single stains of each donor and
#          instrument: (A) AF700 SD of CD3- cells in the CD3-AF700 stain (AF700 reagent, no other
#          dye), beside the unstained and panel-negative SD; (B) per source stain, AF700 SD in CD3+
#          over CD3- cells (that dye's spread into AF700).
# Run:     Rscript experiments/integration/glehr2026/11_af700_single_stains.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*_ungated_*_NN-CD3-<dye>.fcs  single stains
#          <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs       authors' gated panel
#          experiments/runs/2026-11-10_glehr-invert/authors_m_M.csv            Aurora maps, from 01
#          experiments/runs/2026-11-10_glehr-invert-fortessa/authors_m_M.csv   Fortessa maps, from 01
#          experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv  panel neg_sd, 08
#          experiments/runs/2026-11-10_glehr-unstained-noise-v2/unstained_sd.csv   unstained SD, 10
# Outputs: <run_dir>/single_stain_sd.csv     donor, instrument, stain, side, n, AF700 SD, kept share
#          <run_dir>/fig_af700_stain_sd.png, figure_source_data_af700_stain_sd.csv        (A)
#          <run_dir>/fig_af700_single_spread.png, figure_source_data_af700_single_spread.csv (B)
#          <run_dir>/session_info.txt
#
# Decisions G15, G16, G17 (plans/2026-11-10_diagnostics-next.md).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

TARGET <- "AF700-A" # G17
GATE_PROBS <- c(0.05, 0.95) # G15
SPLIT_COFACTOR <- 6000 # split as in 02 (G17)
SPLIT_SEED <- 1
SPLIT_NSTART <- 10
MIN_EVENTS <- 200 # a side with fewer events gets NA
LATER_SESSION <- paste0("d", 6:9) # Aurora donors 6-9 acquired later with other voltages
INSTRUMENTS <- c("aurora", "fortessa")

maps <- bind_rows(
  aurora = read.csv("experiments/runs/2026-11-10_glehr-invert/authors_m_M.csv"),
  fortessa = read.csv(
    "experiments/runs/2026-11-10_glehr-invert-fortessa/authors_m_M.csv"
  ),
  .id = "instrument"
) |>
  filter(channel == TARGET)

read_ev <- function(path) {
  exprs(read.FCS(
    path,
    truncate_max_range = FALSE,
    transformation = FALSE,
    emptyValue = FALSE
  ))
}
robust_sd <- function(v) {
  if (length(v) < MIN_EVENTS) NA_real_ else 1.4826 * mad(v, constant = 1)
}

# Gate, split and AF700 SD per donor, instrument and single stain ----------------------------
donor_dirs <- sort(basename(list.dirs(
  file.path(data_dir, "BAD", "aurora"),
  recursive = FALSE
)))
grid <- expand.grid(
  donor_i = seq_along(donor_dirs),
  instrument = INSTRUMENTS,
  stringsAsFactors = FALSE
)
single_sd <- lapply(seq_len(nrow(grid)), function(j) {
  dir <- file.path(
    data_dir,
    "BAD",
    grid$instrument[j],
    donor_dirs[grid$donor_i[j]]
  )
  donor <- paste0("d", grid$donor_i[j])
  fs <- list.files(dir, pattern = "\\.fcs$", full.names = TRUE)
  panel <- read_ev(fs[grepl("relativised.*12-panel", fs)])
  q <- function(x) quantile(x, GATE_PROBS)
  box <- list(
    fsc = q(panel[, "FSC-A"]),
    ssc = q(panel[, "SSC-A"]),
    ratio = q(panel[, "FSC-H"] / panel[, "FSC-A"])
  )
  ctrl_files <- fs[grepl("_ungated_.*_[0-9]{2}-CD3-", fs)]
  dyes <- sub(".*-CD3-(.*)\\.fcs$", "\\1", ctrl_files)
  n_total <- integer(length(ctrl_files))
  controls <- lapply(seq_along(ctrl_files), function(k) {
    e <- read_ev(ctrl_files[k])
    r <- e[, "FSC-H"] / e[, "FSC-A"]
    keep <- e[, "FSC-A"] >= box$fsc[[1]] &
      e[, "FSC-A"] <= box$fsc[[2]] &
      e[, "SSC-A"] >= box$ssc[[1]] &
      e[, "SSC-A"] <= box$ssc[[2]] &
      r >= box$ratio[[1]] &
      r <= box$ratio[[2]]
    n_total[k] <<- nrow(e)
    e[keep, , drop = FALSE]
  })
  names(controls) <- paste0(donor, "_", dyes)
  control_map <- data.frame(
    unit = donor,
    channel = paste0(dyes, "-A"),
    control = names(controls)
  )
  positive <- sw_relativise_split_kmeans(
    controls,
    control_map,
    split_transform = \(x) asinh(x / SPLIT_COFACTOR),
    seed = SPLIT_SEED,
    nstart = SPLIT_NSTART
  )
  m <- filter(maps, instrument == grid$instrument[j], .data$donor == .env$donor)
  stopifnot("one AF700 map per donor and instrument" = nrow(m) == 1)
  bind_rows(lapply(seq_along(controls), function(k) {
    target <- asinh((m$intercept + m$slope * controls[[k]][, TARGET]) / 1000)
    pos <- positive[[names(controls)[k]]]
    tibble(
      donor = donor,
      instrument = grid$instrument[j],
      stain = dyes[k],
      side = c("CD3-", "CD3+"),
      n = c(sum(!pos), sum(pos)),
      kept_share = nrow(controls[[k]]) / n_total[k],
      af700_sd = c(robust_sd(target[!pos]), robust_sd(target[pos]))
    )
  }))
}) |>
  bind_rows()
write.csv(
  single_sd,
  file.path(run_dir, "single_stain_sd.csv"),
  row.names = FALSE
)
print(as.data.frame(
  single_sd |>
    summarise(
      kept_min = round(min(kept_share), 3),
      kept_max = round(max(kept_share), 3),
      n_min = min(n),
      n_na = sum(is.na(af700_sd)),
      .by = instrument
    )
))

# (A) AF700 SD of negatives by sample type ------------------------------------------------------
# Unstained (10), CD3- cells of the CD3-AF700 stain, and the panel negative (08), same scale.
unstained <- read.csv(
  "experiments/runs/2026-11-10_glehr-unstained-noise-v2/unstained_sd.csv"
) |>
  filter(channel == TARGET) |>
  transmute(donor, instrument, sample_type = "unstained", sd = unstained_sd)
panel_neg <- read.csv(
  "experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv"
) |>
  filter(check == "neg_sd", channel == TARGET) |>
  separate_wider_delim(file, "_", names = c("donor", "instrument")) |>
  transmute(donor, instrument, sample_type = "panel negative", sd = value)
af700_only <- single_sd |>
  filter(stain == "AF700", side == "CD3-") |>
  transmute(
    donor,
    instrument,
    sample_type = "AF700 stain, CD3- cells",
    sd = af700_sd
  )
sd_frame <- bind_rows(unstained, af700_only, panel_neg) |>
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("unstained", "AF700 stain, CD3- cells", "panel negative")
    ),
    session = if_else(
      donor %in% LATER_SESSION,
      "Aurora later session (d6-d9)",
      "d1-d5"
    )
  )
stopifnot(
  "one SD per donor, instrument and sample type" = nrow(sd_frame) == 9 * 2 * 3
)
write.csv(
  sd_frame,
  file.path(run_dir, "figure_source_data_af700_stain_sd.csv"),
  row.names = FALSE
)
print(as.data.frame(
  sd_frame |>
    summarise(
      median_sd = round(median(sd), 3),
      min = round(min(sd), 3),
      max = round(max(sd), 3),
      .by = c(instrument, sample_type)
    )
))
p_a <- ggplot(sd_frame, aes(sample_type, sd, group = donor)) +
  geom_line(alpha = 0.4) +
  geom_point(aes(shape = session), size = 1.2) +
  facet_wrap(vars(instrument)) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "AF700 robust SD, asinh(x / 1000) of relativised values",
    caption = "Lines: donors; robust SD = 1.4826 x MAD"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(run_dir, "fig_af700_stain_sd.png"),
  p_a,
  width = 160,
  height = 90,
  units = "mm",
  dpi = 200
)

# (B) Spread of each source dye into AF700 ----------------------------------------------------
spread <- single_sd |>
  filter(stain != "AF700") |>
  select(donor, instrument, stain, side, af700_sd) |>
  pivot_wider(names_from = side, values_from = af700_sd) |>
  mutate(
    spread_ratio = `CD3+` / `CD3-`,
    session = if_else(
      donor %in% LATER_SESSION,
      "Aurora later session (d6-d9)",
      "d1-d5"
    )
  )
write.csv(
  spread,
  file.path(run_dir, "figure_source_data_af700_single_spread.csv"),
  row.names = FALSE
)
print(as.data.frame(
  spread |>
    summarise(
      median_ratio = round(median(spread_ratio, na.rm = TRUE), 2),
      .by = c(stain, instrument)
    ) |>
    pivot_wider(names_from = instrument, values_from = median_ratio)
))
p_b <- ggplot(spread, aes(stain, spread_ratio, colour = instrument)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(
    aes(shape = session),
    position = position_dodge(width = 0.6),
    size = 1.2
  ) +
  scale_y_log10() +
  labs(
    x = "Source single stain (CD3-<dye>)",
    y = "AF700 SD in CD3+ / CD3- cells",
    caption = "One point per donor; dashed line: no spread"
  )
ggsave(
  file.path(run_dir, "fig_af700_single_spread.png"),
  p_b,
  width = 160,
  height = 90,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
