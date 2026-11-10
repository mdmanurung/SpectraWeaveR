# experiments/integration/glehr2026/10_unstained_noise.R
# Purpose: test whether a wider detector noise floor on the Fortessa accounts for the wider AF700
#          and PB panel negatives (08). Per donor and channel, compare the Fortessa/Aurora robust
#          SD of unstained cells with the same ratio of the panel negatives, on the same scale.
# Run:     Rscript experiments/integration/glehr2026/10_unstained_noise.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*_11-none.fcs              unstained, ungated
#          <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs  authors' gated panel (gate)
#          experiments/runs/2026-11-10_glehr-invert/authors_m_M.csv           Aurora maps, from 01
#          experiments/runs/2026-11-10_glehr-invert-fortessa/authors_m_M.csv  Fortessa maps, from 01
#          experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv  panel neg_sd, 08
# Outputs: <run_dir>/unstained_sd.csv          donor, instrument, channel, n kept, kept share, SD
#          <run_dir>/ratios.csv                donor, channel, unstained and panel ratios, their ratio
#          <run_dir>/fig_unstained_noise.png, figure_source_data_unstained_noise.csv
#          <run_dir>/fig_unstained_gate.png, figure_source_data_unstained_gate.csv
#          <run_dir>/session_info.txt
#
# Decisions G15 and G16 (plans/2026-11-10_diagnostics-next.md).

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

CHANNELS <- c("AF700-A", "PB-A", "FITC-A", "APC-A", "PE-A", "PC7-A") # inherited from 08 (C1)
GATE_PROBS <- c(0.05, 0.95) # G15 (revised): 5th-95th; at 1st-99th the Fortessa box kept granulocytes
LATER_SESSION <- paste0("d", 6:9) # Aurora donors 6-9 acquired later with other voltages (header audit)
INSTRUMENTS <- c("aurora", "fortessa")

maps <- bind_rows(
  aurora = read.csv("experiments/runs/2026-11-10_glehr-invert/authors_m_M.csv"),
  fortessa = read.csv(
    "experiments/runs/2026-11-10_glehr-invert-fortessa/authors_m_M.csv"
  ),
  .id = "instrument"
)
panel_neg <- read.csv(
  "experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv"
) |>
  filter(check == "neg_sd") |>
  separate_wider_delim(file, "_", names = c("donor", "instrument")) |>
  select(donor, instrument, channel, panel_sd = value)

# Gate and robust SD per donor and instrument ------------------------------------------------
# The box comes from the authors' gated panel of the same donor and instrument (scatter is not
# rescaled by relativisation); unstained events are mapped with the authors' linear map, then
# asinh(x / 1000) as in 08.
donor_dirs <- sort(basename(list.dirs(
  file.path(data_dir, "BAD", "aurora"),
  recursive = FALSE
)))
read_ev <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  stopifnot("one file per pattern" = length(f) == 1)
  exprs(read.FCS(
    f,
    truncate_max_range = FALSE,
    transformation = FALSE,
    emptyValue = FALSE
  ))
}
grid <- expand.grid(
  donor_i = seq_along(donor_dirs),
  instrument = INSTRUMENTS,
  stringsAsFactors = FALSE
)
per_file <- lapply(seq_len(nrow(grid)), function(j) {
  dir <- file.path(
    data_dir,
    "BAD",
    grid$instrument[j],
    donor_dirs[grid$donor_i[j]]
  )
  donor <- paste0("d", grid$donor_i[j])
  panel <- read_ev(dir, "relativised.*12-panel\\.fcs$")
  unst <- read_ev(dir, "_11-none\\.fcs$")
  ratio_p <- panel[, "FSC-H"] / panel[, "FSC-A"]
  ratio_u <- unst[, "FSC-H"] / unst[, "FSC-A"]
  inside <- function(x, ref) {
    q <- quantile(ref, GATE_PROBS)
    x >= q[[1]] & x <= q[[2]]
  }
  keep <- inside(unst[, "FSC-A"], panel[, "FSC-A"]) &
    inside(unst[, "SSC-A"], panel[, "SSC-A"]) &
    inside(ratio_u, ratio_p)
  m <- maps |>
    filter(instrument == grid$instrument[j], .data$donor == .env$donor)
  sl <- setNames(m$slope, m$channel)[CHANNELS]
  ic <- setNames(m$intercept, m$channel)[CHANNELS]
  stopifnot("map for every channel" = !anyNA(sl) && !anyNA(ic))
  rel <- asinh(
    sweep(sweep(unst[keep, CHANNELS, drop = FALSE], 2, sl, `*`), 2, ic, `+`) /
      1000
  )
  list(
    sd = tibble(
      donor = donor,
      instrument = grid$instrument[j],
      channel = CHANNELS,
      n_total = nrow(unst),
      n_kept = sum(keep),
      kept_share = mean(keep),
      unstained_sd = apply(rel, 2, \(v) 1.4826 * mad(v, constant = 1))
    ),
    gate = tibble(
      donor = donor,
      instrument = grid$instrument[j],
      fsc_a = unst[, "FSC-A"],
      ssc_a = unst[, "SSC-A"],
      kept = keep
    )
  )
})
unstained_sd <- bind_rows(lapply(per_file, `[[`, "sd"))
write.csv(
  unstained_sd,
  file.path(run_dir, "unstained_sd.csv"),
  row.names = FALSE
)

# Fortessa over Aurora ratios ------------------------------------------------------------------
# One row per donor and channel; the estimand is the unstained ratio over the panel-negative ratio.
ratios <- unstained_sd |>
  select(donor, instrument, channel, unstained_sd) |>
  left_join(
    panel_neg,
    by = join_by(donor, instrument, channel),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  pivot_wider(
    names_from = instrument,
    values_from = c(unstained_sd, panel_sd)
  ) |>
  mutate(
    unstained_ratio = unstained_sd_fortessa / unstained_sd_aurora,
    panel_ratio = panel_sd_fortessa / panel_sd_aurora,
    unstained_over_panel = unstained_ratio / panel_ratio,
    session = if_else(
      donor %in% LATER_SESSION,
      "Aurora later session (d6-d9)",
      "d1-d5"
    )
  )
write.csv(ratios, file.path(run_dir, "ratios.csv"), row.names = FALSE)
print(as.data.frame(
  ratios |>
    summarise(
      panel_ratio = round(median(panel_ratio), 2),
      unstained_ratio = round(median(unstained_ratio), 2),
      over_median = round(median(unstained_over_panel), 2),
      over_min = round(min(unstained_over_panel), 2),
      over_max = round(max(unstained_over_panel), 2),
      .by = channel
    )
))
print(as.data.frame(
  unstained_sd |>
    summarise(
      kept_min = round(min(kept_share), 3),
      kept_max = round(max(kept_share), 3),
      .by = instrument
    )
))

# Figures ----------------------------------------------------------------------------------------
plot_frame <- ratios |>
  mutate(channel = factor(channel, levels = CHANNELS))
write.csv(
  plot_frame,
  file.path(run_dir, "figure_source_data_unstained_noise.csv"),
  row.names = FALSE
)
p <- ggplot(
  plot_frame,
  aes(panel_ratio, unstained_ratio, colour = channel, shape = session)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  theme(legend.position = "right") +
  labs(
    x = "Panel-negative SD, Fortessa / Aurora",
    y = "Unstained SD, Fortessa / Aurora",
    caption = "One point per donor and channel; dashed line: equal ratios\nasinh(x / 1000) of relativised values"
  )
ggsave(
  file.path(run_dir, "fig_unstained_noise.png"),
  p,
  width = 160,
  height = 110,
  units = "mm",
  dpi = 200
)

set.seed(1) # display subsample only
gate_frame <- bind_rows(lapply(per_file, `[[`, "gate")) |>
  filter(donor == "d1") |>
  slice_sample(n = 20000, by = instrument)
write.csv(
  gate_frame,
  file.path(run_dir, "figure_source_data_unstained_gate.csv"),
  row.names = FALSE
)
p_gate <- ggplot(gate_frame, aes(fsc_a, ssc_a, colour = kept)) +
  geom_point(size = 0.1, alpha = 0.3) +
  facet_wrap(vars(instrument), scales = "free") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75")) +
  labs(
    x = "FSC-A",
    y = "SSC-A",
    colour = "kept by box",
    caption = "Donor 1 unstained, 20,000 events per instrument"
  )
ggsave(
  file.path(run_dir, "fig_unstained_gate.png"),
  p_gate,
  width = 160,
  height = 80,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
