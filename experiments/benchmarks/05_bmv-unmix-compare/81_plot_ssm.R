# Purpose: draw the AutoSpectral 1.8.8 spillover spread matrix (SSM) shipped in the zip, because the zip's own
#          heatmap is blank (its normalized SSM is all NA; the raw SSM has values) (user, 2026-11-10). Rows and
#          columns are as stored in the CSV; the direction (spreading dye vs receiving dye) is not documented in
#          the package code found, so the figure says "row, column" only. Display only, no new computation.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/81_plot_ssm.R   (repo root)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_similarity_heatmap/
#            autospectral_spillover_spread_matrix.csv  [user zip]
# Outputs: experiments/runs/2026-11-10_bmv-ssm-188/fig_ssm_188.png   heatmap with values, diagonal blank
#          experiments/runs/2026-11-10_bmv-ssm-188/ssm_long.csv       the plotted values, long format
#          experiments/runs/2026-11-10_bmv-ssm-188/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

# Parameters ---------------------------------------------------------------------------------------
SSM_CSV <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_similarity_heatmap/autospectral_spillover_spread_matrix.csv"
OUT_DIR <- "experiments/runs/2026-11-10_bmv-ssm-188"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Long table of the SSM, dye order as in the CSV ------------------------------------------------------
ssm <- as.matrix(read.csv(SSM_CSV, row.names = 1, check.names = FALSE))
dyes <- rownames(ssm)
stopifnot("square SSM with matching names" = identical(dyes, colnames(ssm)))
ssm_long <- as_tibble(ssm, rownames = "row_dye") |>
  pivot_longer(-row_dye, names_to = "col_dye", values_to = "ssm") |>
  mutate(
    row_dye = factor(row_dye, levels = rev(dyes)),
    col_dye = factor(col_dye, levels = dyes),
    label = if_else(is.na(ssm) | ssm < 0.5, "", sprintf("%.0f", ssm))
  )
write_csv(ssm_long, file.path(OUT_DIR, "ssm_long.csv"))

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(ssm_long, aes(col_dye, row_dye, fill = ssm)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = label), size = 1.8) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    trans = "sqrt",
    na.value = "grey85",
    name = "SSM"
  ) +
  coord_fixed() +
  labs(
    x = "column dye",
    y = "row dye",
    title = "AutoSpectral 1.8.8 spillover spread matrix (zip), values >= 0.5 printed; diagonal blank"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave(
  file.path(OUT_DIR, "fig_ssm_188.png"),
  p,
  width = 8,
  height = 7.5,
  dpi = 170,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
