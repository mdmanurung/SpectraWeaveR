# experiments/integration/glehr2026/13_af700_compensation.R
# Purpose: find whether APC (CD4) or AA750 (CD3) to AF700 under-compensation on the Fortessa
#          makes the CD4+ T cells' AF700 wide and shifted (12). Under-compensation gives AF700 a
#          positive slope on APC or AA750 and a median shift in single stains; spread does not.
#          (E1) panel CD4+ T cells: Huber fit of AF700 on APC + AA750, and the AF700 SD before and
#          after removing the fit. (E2) APC and AA750 single stains: AF700 median, CD3+ minus CD3-.
# Run:     Rscript experiments/integration/glehr2026/13_af700_compensation.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs        authors' CD3+ panel
#          <data_dir>/BAD/<instrument>/<donor>/*_ungated_*_NN-CD3-{APC,AA750}.fcs single stains
#          experiments/runs/2026-11-10_glehr-invert/authors_m_M.csv              Aurora maps, from 01
#          experiments/runs/2026-11-10_glehr-invert-fortessa/authors_m_M.csv     Fortessa maps, from 01
# Outputs: <run_dir>/panel_fit.csv      donor, instrument, n, slopes, SD before and after
#          <run_dir>/single_shift.csv   donor, instrument, stain, n per side, AF700 median shift
#          <run_dir>/fig_af700_comp_slope.png, figure_source_data_af700_comp_slope.csv
#          <run_dir>/fig_af700_comp_shift.png, figure_source_data_af700_comp_shift.csv
#          <run_dir>/session_info.txt
#
# Decisions G15, G16, G17, G18, G19 (plans/2026-11-10_diagnostics-next.md).

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

TARGET <- "AF700-A" # CD8
CD4 <- "APC-A"
CD3 <- "AA750-A"
SOURCE_STAINS <- c("APC", "AA750") # E2
SPLIT_SEED <- 1 # G17, G18
SPLIT_NSTART <- 10
SPLIT_COFACTOR <- 6000 # single-stain split as in 02 (G17)
GATE_PROBS <- c(0.05, 0.95) # G15
MIN_EVENTS <- 200
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
robust_sd <- function(v) 1.4826 * mad(v, constant = 1)

donor_dirs <- sort(basename(list.dirs(
  file.path(data_dir, "BAD", "aurora"),
  recursive = FALSE
)))
grid <- expand.grid(
  donor_i = seq_along(donor_dirs),
  instrument = INSTRUMENTS,
  stringsAsFactors = FALSE
)

# E1: Huber fit in panel CD4+ T cells ----------------------------------------------------------
# Panel values are already relativised (linear). CD4+ T cells by k-means on asinh(CD4 / 1000) (G18).
panel_fit <- lapply(seq_len(nrow(grid)), function(j) {
  dir <- file.path(
    data_dir,
    "BAD",
    grid$instrument[j],
    donor_dirs[grid$donor_i[j]]
  )
  e <- read_ev(list.files(
    dir,
    pattern = "relativised.*12-panel\\.fcs$",
    full.names = TRUE
  ))
  set.seed(SPLIT_SEED)
  km <- kmeans(asinh(e[, CD4] / 1000), centers = 2, nstart = SPLIT_NSTART)
  cd4_pos <- km$cluster == which.max(km$centers[, 1])
  d <- data.frame(
    af700 = e[cd4_pos, TARGET],
    apc = e[cd4_pos, CD4],
    aa750 = e[cd4_pos, CD3]
  )
  fit <- MASS::rlm(
    af700 ~ apc + aa750,
    data = d,
    psi = MASS::psi.huber,
    maxit = 50
  )
  stopifnot("rlm converged" = fit$converged)
  # Removing the fitted slopes keeps the intercept, so the centre stays put.
  adjusted <- d$af700 -
    coef(fit)[["apc"]] * d$apc -
    coef(fit)[["aa750"]] * d$aa750
  tibble(
    donor = paste0("d", grid$donor_i[j]),
    instrument = grid$instrument[j],
    n = nrow(d),
    slope_apc = coef(fit)[["apc"]],
    slope_aa750 = coef(fit)[["aa750"]],
    sd_before = robust_sd(asinh(d$af700 / 1000)),
    sd_after = robust_sd(asinh(adjusted / 1000)),
    median_before = median(asinh(d$af700 / 1000))
  )
}) |>
  bind_rows()
write.csv(panel_fit, file.path(run_dir, "panel_fit.csv"), row.names = FALSE)
print(as.data.frame(
  panel_fit |>
    summarise(
      across(
        c(n, slope_apc, slope_aa750, sd_before, sd_after, median_before),
        \(v) signif(median(v), 3)
      ),
      .by = instrument
    )
))

# E2: AF700 median shift in the APC and AA750 single stains ----------------------------------
single_shift <- lapply(seq_len(nrow(grid)), function(j) {
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
  m <- filter(maps, instrument == grid$instrument[j], .data$donor == .env$donor)
  stopifnot("one AF700 map per donor and instrument" = nrow(m) == 1)
  bind_rows(lapply(SOURCE_STAINS, function(s) {
    f <- fs[grepl(paste0("_ungated_.*_[0-9]{2}-CD3-", s, "\\.fcs$"), fs)]
    stopifnot("one single stain per dye" = length(f) == 1)
    e <- read_ev(f)
    r <- e[, "FSC-H"] / e[, "FSC-A"]
    keep <- e[, "FSC-A"] >= box$fsc[[1]] &
      e[, "FSC-A"] <= box$fsc[[2]] &
      e[, "SSC-A"] >= box$ssc[[1]] &
      e[, "SSC-A"] <= box$ssc[[2]] &
      r >= box$ratio[[1]] &
      r <= box$ratio[[2]]
    e <- e[keep, , drop = FALSE]
    ch <- paste0(s, "-A")
    set.seed(SPLIT_SEED)
    km <- kmeans(
      asinh(e[, ch] / SPLIT_COFACTOR),
      centers = 2,
      nstart = SPLIT_NSTART
    )
    pos <- km$cluster == which.max(km$centers[, 1])
    af700 <- m$intercept + m$slope * e[, TARGET] # linear relativised AF700
    tibble(
      donor = donor,
      instrument = grid$instrument[j],
      stain = s,
      n_neg = sum(!pos),
      n_pos = sum(pos),
      shift = if (min(sum(pos), sum(!pos)) < MIN_EVENTS) {
        NA_real_
      } else {
        median(af700[pos]) - median(af700[!pos])
      }
    )
  }))
}) |>
  bind_rows()
write.csv(
  single_shift,
  file.path(run_dir, "single_shift.csv"),
  row.names = FALSE
)
print(as.data.frame(
  single_shift |>
    summarise(
      n_min = min(pmin(n_neg, n_pos)),
      shift_median = round(median(shift, na.rm = TRUE), 1),
      shift_min = round(min(shift, na.rm = TRUE), 1),
      shift_max = round(max(shift, na.rm = TRUE), 1),
      .by = c(stain, instrument)
    )
))

# Figures ----------------------------------------------------------------------------------------
session_of <- function(d) {
  if_else(d %in% LATER_SESSION, "donors d6-d9 (later Aurora session)", "d1-d5")
}
slope_frame <- panel_fit |>
  pivot_longer(
    c(slope_apc, slope_aa750, sd_before, sd_after),
    names_to = "quantity",
    values_to = "value"
  ) |>
  mutate(
    quantity = factor(
      quantity,
      levels = c("slope_apc", "slope_aa750", "sd_before", "sd_after"),
      labels = c(
        "slope on APC (CD4)",
        "slope on AA750 (CD3)",
        "AF700 SD before",
        "AF700 SD after removing fit"
      )
    ),
    session = session_of(donor)
  )
write.csv(
  slope_frame,
  file.path(run_dir, "figure_source_data_af700_comp_slope.csv"),
  row.names = FALSE
)
p_slope <- ggplot(slope_frame, aes(instrument, value)) +
  geom_hline(
    data = data.frame(quantity = levels(slope_frame$quantity)[1:2], y = 0) |>
      mutate(
        quantity = factor(quantity, levels = levels(slope_frame$quantity))
      ),
    aes(yintercept = y),
    linetype = 2,
    colour = "grey50"
  ) +
  geom_point(
    aes(shape = session),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    size = 1.2
  ) +
  facet_wrap(vars(quantity), scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = NULL,
    caption = "Panel CD4+ T cells, one point per donor. Slopes: Huber fit on linear relativised values (AF700 units per APC or AA750 unit).\nSDs: 1.4826 x MAD on asinh(x / 1000)."
  )
ggsave(
  file.path(run_dir, "fig_af700_comp_slope.png"),
  p_slope,
  width = 180,
  height = 80,
  units = "mm",
  dpi = 200
)

shift_frame <- single_shift |>
  mutate(session = session_of(donor))
write.csv(
  shift_frame,
  file.path(run_dir, "figure_source_data_af700_comp_shift.csv"),
  row.names = FALSE
)
p_shift <- ggplot(shift_frame, aes(instrument, shift)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(
    aes(shape = session),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    size = 1.2
  ) +
  facet_wrap(
    vars(stain),
    labeller = labeller(stain = \(s) paste0("CD3-", s, " single stain"))
  ) +
  labs(
    x = NULL,
    y = "AF700 median, CD3+ minus CD3-\n(relativised units)",
    caption = "One point per donor; CD3+ in the AF700 control = 1e4 units"
  )
ggsave(
  file.path(run_dir, "fig_af700_comp_shift.png"),
  p_shift,
  width = 120,
  height = 80,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
