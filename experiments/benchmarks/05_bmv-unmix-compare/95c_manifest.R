# Purpose: merge the per-script manifest parts of the 95c restyles into MANIFEST_restyle_common.csv and
#          list every certified figure (audit_freshness.csv plus the old restyle folders) that has no new figure,
#          with the reason: skipped as spike (batch 5), superseded, or no new source.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_manifest.R   (repo root)
# Inputs:  experiments/runs/2026-11-11_bmv-restyle-common/manifest_parts/95c_*.csv   one per script
#          experiments/runs/2026-11-11_bmv-restyle/audit_freshness.csv   the certified figure list
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/MANIFEST_restyle_common.csv
#          experiments/runs/2026-11-11_bmv-restyle-common/MANIFEST_coverage.csv   every certified figure and its status

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
})
OUT_ROOT <- "experiments/runs/2026-11-11_bmv-restyle-common"
OLD_ROOT <- "experiments/runs/2026-11-11_bmv-restyle"

# Merge -----------------------------------------------------------------------------------------------
parts <- list.files(
  file.path(OUT_ROOT, "manifest_parts"),
  pattern = "^95c_.*[.]csv$",
  full.names = TRUE
)
manifest <- parts |>
  map(\(f) read_csv(f, col_types = cols(.default = "c"))) |>
  list_rbind() |>
  mutate(across(c(old_figure, new_figure), \(x) na_if(x, "NA"))) |>
  distinct(old_figure, new_figure, .keep_all = TRUE)
# rows with no new figure (no source) carry NA; every other row must point at a file
drawn <- filter(manifest, !is.na(new_figure))
stopifnot("every new figure exists" = all(file.exists(drawn$new_figure)))
write_csv(manifest, file.path(OUT_ROOT, "MANIFEST_restyle_common.csv"))

# Coverage of the certified list -------------------------------------------------------------------------
audit <- read_csv(
  file.path(OLD_ROOT, "audit_freshness.csv"),
  show_col_types = FALSE
) |>
  filter(!str_detect(.data$used_in, "^dashboard v19 img"))
# old figure key as the manifest writes it: report figures by their s* name, restyle figures by path under runs/
key <- \(x) {
  str_remove(x, "^experiments/runs/") |>
    str_replace("^2026-11-1[01]_bmv-report/figures-v[0-9]+/", "report/")
}
manifest_key <- manifest |>
  mutate(
    k = key(.data$old_figure) |>
      str_replace("^2026-11-1[01]_bmv-report/figures-v[0-9]+/", "report/")
  )
coverage <- audit |>
  transmute(
    used_in,
    old_figure = .data$figure,
    k = key(.data$figure),
    superseded = .data$current %in% FALSE & str_detect(.data$note, "superseded")
  ) |>
  left_join(
    select(manifest_key, k, new_figure, status),
    by = "k",
    relationship = "many-to-many"
  ) |>
  mutate(
    status = case_when(
      !is.na(.data$status) ~ .data$status,
      .data$superseded ~ "superseded by per_dye_v8 / af_swap_v3",
      str_detect(
        .data$old_figure,
        "spike|18_|20_|22_|23_|21_raw_companions/spike_raw"
      ) ~ "skipped as spike (batch 5)",
      .default = "no new source"
    )
  ) |>
  select(-k, -superseded)
write_csv(coverage, file.path(OUT_ROOT, "MANIFEST_coverage.csv"))
print(count(coverage, status))
