# Purpose: redraw the spectra-angle figure (f6 in report.qmd) without the long title, with legible
#          fonts, from the saved table; the data are unchanged. The plot does not depend on the channel
#          set, so the "unpinned" copy is the same figure under its second name.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/44_replot_spectra_angle.R
#          from the repo root.
# Inputs:  experiments/runs/2026-11-09_bmv-unmix-compare-final/spectra_compare.csv  angle per dye and
#          1.8.8 run [04_spectra.R]
# Outputs: experiments/runs/2026-11-10_bmv-figures-v2/2026-11-09_bmv-unmix-compare-final/
#            f6_spectra_angle.png           angle between 1.0.2 and 1.8.8 spectra per dye
#            f6_spectra_angle_unpinned.png  the same figure, second name used by the report
#            sessionInfo_44.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

out_dir <- "experiments/runs/2026-11-10_bmv-figures-v2/2026-11-09_bmv-unmix-compare-final"
stopifnot(
  "44 outputs must not exist" = !file.exists(file.path(
    out_dir,
    "f6_spectra_angle.png"
  ))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Angle per dye, both versions present (as in report.qmd f6) ---------------------------------
spectra_angle <- read_csv(
  "experiments/runs/2026-11-09_bmv-unmix-compare-final/spectra_compare.csv",
  show_col_types = FALSE
) |>
  filter(status == "both") |>
  mutate(fluorophore = reorder(fluorophore, angle_deg)) # stats::reorder; forcats is not in the pixi R

# Plot -----------------------------------------------------------------------------------------
p_angle <- ggplot(
  spectra_angle,
  aes(angle_deg, fluorophore, colour = stamp_188)
) +
  geom_point(size = 2) +
  guides(colour = guide_legend(nrow = 2)) + # three run stamps do not fit on one line
  labs(
    title = "Spectra angle, 1.0.2 vs 1.8.8",
    x = "Angle (degrees)",
    y = NULL,
    colour = "1.8.8 run"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(size = 9), legend.position = "bottom")

ggsave(
  file.path(out_dir, "f6_spectra_angle.png"),
  p_angle,
  width = 160,
  height = 190,
  units = "mm",
  dpi = 150
)
ggsave(
  file.path(out_dir, "f6_spectra_angle_unpinned.png"),
  p_angle,
  width = 160,
  height = 190,
  units = "mm",
  dpi = 150
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_44.txt")
)
