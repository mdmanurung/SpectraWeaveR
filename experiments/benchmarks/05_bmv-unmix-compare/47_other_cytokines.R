#!/usr/bin/env Rscript
# Purpose: check whether the unmixing arms also differ for IL-2, IL-10 and Th2, as they do for TNF, or
#          agree, as they do for IFNg. Per arm, donor and parent (CD4T, CD8T): % positive at a fixed cut,
#          density shape near 0, co-expression with TNF and IFNg, Spearman rho for every cytokine pair
#          and for three expected links (IFNg-Tbet, IL-10-FoxP3, IL-2-CD25), and negative spread against
#          as188_joint on cells negative in every arm.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/47_other_cytokines.R \
#            experiments/runs/2026-11-10_bmv-gated-pooled experiments/runs/2026-11-10_bmv-other-cytokines
#          from the repo root.
# Inputs:  <gated_dir>/manifest.csv [01_manifest-gate.R]
#          <gated_dir>/gs/<arm>/ saved GatingSets with CD4T and CD8T [19a_gate_arm.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  ARM_LABELS, VERSION_ROWS, ZIP_ARMS
# Outputs: <run_dir>/pct_pos.csv                % above the cut per arm, donor, parent, cytokine
#          <run_dir>/il2_coexpression.csv       % IL-2+ among TNF+ and among IFNg+
#          <run_dir>/rho_pairs.csv              Spearman rho per arm, donor, parent, pair
#          <run_dir>/neg_spread.csv             Nguyen rSD on shared negatives and ratio to as188_joint
#          <run_dir>/fig_a_pct_pos.png          % positive, donors joined (+ fig_a_source.csv)
#          <run_dir>/fig_b_density_<cytokine>.png  density near 0 per donor and parent (+ fig_b_source.csv)
#          <run_dir>/fig_c_il2_coexpression.png % IL-2+ among TNF+ and IFNg+ (+ fig_c_source.csv)
#          <run_dir>/fig_d_rho_cytokine_pairs.png  rho for the 10 cytokine pairs (+ fig_d_source.csv)
#          <run_dir>/fig_e_rho_expected_links.png  rho for the 3 expected links (+ fig_e_source.csv)
#          <run_dir>/fig_f_neg_spread.png       rSD ratio to as188_joint (+ fig_f_source.csv)
#          <run_dir>/sessionInfo.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # same ggplot2 as 38/40/46

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS, ZIP_ARMS, ZIP_FILL

# Parameters ------------------------------------------------------------------------------
CUT <- 1 # fixed cut, arcsinh, for every cytokine, arm and donor (user, 2026-11-10: "Both")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# TNF and IFNg are stored as arcsinh in the GatingSets (19a); every other channel is stored linear.
CYTOKINES <- c(
  TNF = "BV510-A",
  IFNg = "BV570-A",
  `IL-2` = "BV605-A",
  `IL-10` = "BB700-A",
  Th2 = "eFluor 450-A"
)
STORED_ARCSINH <- c("BV510-A", "BV570-A")
NEW_CYTOKINES <- c("IL-2", "IL-10", "Th2") # the three this script adds
# expected biological links, used as reference: rank correlation only, no cut
LINKS <- tribble(
  ~link           , ~from   , ~to     , ~parents          ,
  "IFNg - Tbet"   , "IFNg"  , "Tbet"  , c("CD4T", "CD8T") ,
  "IL-10 - FoxP3" , "IL-10" , "FoxP3" , "CD4T"            ,
  "IL-2 - CD25"   , "IL-2"  , "CD25"  , c("CD4T", "CD8T")
)
MARKERS <- c(Tbet = "PE-Cy5-A", FoxP3 = "RB780-A", CD25 = "SBUV605-A")
POPULATIONS <- c("CD4T", "CD8T")
REF_ARM <- "as188_joint" # reference arm for the spread ratio (plan 2026-11-09)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)
GROUP_BREAKS <- c(1.5, 5.5, 9.5) # dashed lines between the version groups
VERSION_GROUP <- VERSION_ROWS |>
  imap(\(arms, group) set_names(rep(group, length(arms)), arms)) |>
  unname() |>
  unlist()
DENSITY_RANGE <- c(-1.5, 4) # arcsinh display range for densities; display only

cli_args <- commandArgs(trailingOnly = TRUE)
gated_dir <- cli_args[[1]]
run_dir <- cli_args[[2]]
stopifnot("47 run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir)
manifest <- read_csv(
  file.path(gated_dir, "manifest.csv"),
  show_col_types = FALSE
)
donors <- sort(unique(manifest$sample))
gating_sets <- ARM_ORDER |>
  set_names() |>
  map(\(a) load_gs(file.path(gated_dir, "gs", a)))

# Per donor: every arm's events on one arcsinh scale -------------------------------------------
# All arms come from the same raw events, so row i is the same cell in every arm. Linear channels get
# arcsinh(x / 6000) in memory; TNF and IFNg are already arcsinh.
channels <- c(CYTOKINES, MARKERS)
events_of <- function(gs, s) {
  x <- exprs(gh_pop_get_data(gs[[s]], "root"))[, channels]
  colnames(x) <- names(channels)
  linear <- !channels %in% STORED_ARCSINH
  x[, linear] <- asinh(x[, linear] / COFACTOR)
  x
}
by_donor <- donors |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(gs) events_of(gs, s))
    stopifnot("same events in every arm" = n_distinct(map_int(ev, nrow)) == 1)
    in_pop <- map(gating_sets, \(gs) {
      map(set_names(POPULATIONS), \(pop) gh_pop_get_indices(gs[[s]], pop))
    })
    list(ev = ev, in_pop = in_pop)
  })

# Per arm, donor, parent: % positive, IL-2 co-expression, rank correlations ---------------------
pair_names <- combn(names(CYTOKINES), 2, simplify = FALSE)
per_cell_metrics <- by_donor |>
  imap(\(d, s) {
    ARM_ORDER |>
      set_names() |>
      map(\(a) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], ]
            pos <- x[, names(CYTOKINES)] > CUT
            pairs <- pair_names |>
              map(\(p) {
                tibble(
                  pair = paste(p, collapse = " - "),
                  rho = cor(x[, p[[1]]], x[, p[[2]]], method = "spearman")
                )
              }) |>
              list_rbind()
            links <- LINKS |>
              dplyr::filter(map_lgl(parents, \(ps) pop %in% ps)) |>
              mutate(
                rho = map2_dbl(from, to, \(cx, cy) {
                  cor(x[, cx], x[, cy], method = "spearman")
                })
              ) |>
              select(pair = link, rho)
            list(
              pct = enframe(
                100 * colMeans(pos),
                name = "cytokine",
                value = "pct_pos"
              ),
              il2 = tibble(
                among = c("TNF+", "IFNg+"),
                n_among = c(sum(pos[, "TNF"]), sum(pos[, "IFNg"])),
                pct_il2 = 100 *
                  c(
                    mean(pos[pos[, "TNF"], "IL-2"]),
                    mean(pos[pos[, "IFNg"], "IL-2"])
                  )
              ),
              rho = bind_rows(
                mutate(pairs, kind = "cytokine pair"),
                mutate(links, kind = "expected link")
              )
            )
          })
      })
  })

# flatten donor -> arm -> parent into long tables
flatten_metric <- function(metric) {
  per_cell_metrics |>
    map(\(by_a) {
      map(by_a, \(by_p) list_rbind(map(by_p, metric), names_to = "parent")) |>
        list_rbind(names_to = "arm")
    }) |>
    list_rbind(names_to = "donor") |>
    mutate(arm = factor(arm, levels = ARM_ORDER))
}
pct_pos <- flatten_metric("pct")
il2_coexpression <- flatten_metric("il2")
rho_pairs <- flatten_metric("rho")
write_csv(pct_pos, file.path(run_dir, "pct_pos.csv"))
write_csv(il2_coexpression, file.path(run_dir, "il2_coexpression.csv"))
write_csv(rho_pairs, file.path(run_dir, "rho_pairs.csv"))

# Negative spread on cells negative in every arm ------------------------------------------------
# Shared negatives for a cytokine: in the parent in every arm and below the cut in every arm. Nguyen
# robust SD (84th - 50th percentile) on the linear scale, as in 34_rank_symmetric.R.
neg_spread <- by_donor |>
  imap(\(d, s) {
    POPULATIONS |>
      set_names() |>
      map(\(pop) {
        in_all <- reduce(map(d$in_pop, pop), `&`)
        names(CYTOKINES) |>
          set_names() |>
          map(\(cy) {
            neg_all <- reduce(map(d$ev, \(x) x[, cy] < CUT), `&`)
            keep <- in_all & neg_all
            map(d$ev, \(x) {
              lin <- COFACTOR * sinh(x[keep, cy])
              tibble(
                n_shared = sum(keep),
                rsd = unname(diff(quantile(lin, c(0.5, 0.84))))
              )
            }) |>
              list_rbind(names_to = "arm")
          }) |>
          list_rbind(names_to = "cytokine")
      }) |>
      list_rbind(names_to = "parent")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(
    rsd_ratio = rsd / rsd[arm == REF_ARM],
    .by = c(donor, parent, cytokine)
  ) |>
  mutate(arm = factor(arm, levels = ARM_ORDER))
write_csv(neg_spread, file.path(run_dir, "neg_spread.csv"))

# Shared figure pieces -------------------------------------------------------------------------
zip_shade <- annotate(
  "rect",
  xmin = match(ZIP_ARMS, ARM_ORDER) - 0.5,
  xmax = match(ZIP_ARMS, ARM_ORDER) + 0.5,
  ymin = -Inf,
  ymax = Inf,
  fill = ZIP_FILL
)
group_lines <- geom_vline(
  xintercept = GROUP_BREAKS,
  linetype = "dashed",
  colour = "grey60"
)
donor_layers <- list(
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40"),
  geom_line(aes(group = donor, colour = donor), alpha = 0.6, linewidth = 0.4),
  geom_point(aes(colour = donor), size = 1.6),
  scale_colour_brewer(palette = "Dark2", name = "Donor"),
  scale_x_discrete(labels = ARM_LABELS),
  theme_bw(base_size = 11),
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    plot.margin = margin(5, 5, 5, 45) # room for the first slanted arm label
  )
)
rho_theme <- list(
  geom_vline(xintercept = 0, colour = "grey60"),
  scale_y_discrete(labels = ARM_LABELS),
  theme_bw(base_size = 11),
  theme(
    axis.text.y = element_text(
      size = 9,
      face = if_else(rev(ARM_ORDER) %in% ZIP_ARMS, "bold", "plain")
    ),
    axis.text.x = element_text(size = 8),
    legend.position = "top"
  )
)

# Figure a: % positive at the cut, new cytokines ----------------------------------------------
fig_a <- pct_pos |>
  dplyr::filter(cytokine %in% NEW_CYTOKINES) |>
  mutate(cytokine = factor(cytokine, levels = NEW_CYTOKINES))
write_csv(fig_a, file.path(run_dir, "fig_a_source.csv"))
p_a <- ggplot(fig_a, aes(arm, pct_pos)) +
  zip_shade +
  group_lines +
  donor_layers +
  facet_grid(parent ~ cytokine, scales = "free_y") +
  labs(x = NULL, y = "% positive of parent (arcsinh > 1)")
ggsave(
  file.path(run_dir, "fig_a_pct_pos.png"),
  p_a,
  width = 13,
  height = 7.5,
  dpi = 130,
  bg = "white"
)

# Figure b: density near 0, one file per new cytokine ----------------------------------------
# lines per arm, coloured by version group; same events in every arm, so shape differences come from
# unmixing
fig_b <- by_donor |>
  imap(\(d, s) {
    ARM_ORDER |>
      set_names() |>
      map(\(a) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], NEW_CYTOKINES]
            NEW_CYTOKINES |>
              set_names() |>
              map(\(cy) {
                dens <- density(
                  x[, cy],
                  from = DENSITY_RANGE[[1]],
                  to = DENSITY_RANGE[[2]],
                  n = 256
                )
                tibble(value = dens$x, density = dens$y)
              }) |>
              list_rbind(names_to = "cytokine")
          }) |>
          list_rbind(names_to = "parent")
      }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(
    version = factor(VERSION_GROUP[arm], levels = names(VERSION_ROWS)),
    arm = factor(arm, levels = ARM_ORDER)
  )
write_csv(fig_b, file.path(run_dir, "fig_b_source.csv"))
version_colours <- c(
  SpectroFlo = "#666666",
  `AutoSpectral 1.0.2` = "#D55E00",
  `AutoSpectral 1.8.8` = "#0072B2",
  `Crossed: spectra x code` = "#009E73"
)
p_b_of <- function(cy) {
  ggplot(
    dplyr::filter(fig_b, cytokine == cy),
    aes(value, density, group = arm, colour = version)
  ) +
    geom_vline(xintercept = CUT, linetype = "dotted") +
    geom_line(linewidth = 0.45, alpha = 0.8) +
    scale_colour_manual(values = version_colours, name = NULL) +
    facet_grid(parent ~ donor, scales = "free_y") +
    labs(x = paste0(cy, " (arcsinh)"), y = "Density") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}
ggsave(
  file.path(run_dir, "fig_b_density_IL-2.png"),
  p_b_of("IL-2"),
  width = 13,
  height = 6,
  dpi = 130,
  bg = "white"
)
ggsave(
  file.path(run_dir, "fig_b_density_IL-10.png"),
  p_b_of("IL-10"),
  width = 13,
  height = 6,
  dpi = 130,
  bg = "white"
)
ggsave(
  file.path(run_dir, "fig_b_density_Th2.png"),
  p_b_of("Th2"),
  width = 13,
  height = 6,
  dpi = 130,
  bg = "white"
)

# Figure c: % IL-2+ among TNF+ and among IFNg+ ------------------------------------------------
fig_c <- il2_coexpression |>
  mutate(among = factor(among, levels = c("TNF+", "IFNg+")))
write_csv(fig_c, file.path(run_dir, "fig_c_source.csv"))
p_c <- ggplot(fig_c, aes(arm, pct_il2)) +
  zip_shade +
  group_lines +
  donor_layers +
  facet_grid(
    parent ~ among,
    labeller = labeller(among = \(v) paste("among", v))
  ) +
  labs(x = NULL, y = "% IL-2+ (arcsinh > 1)")
ggsave(
  file.path(run_dir, "fig_c_il2_coexpression.png"),
  p_c,
  width = 11,
  height = 7.5,
  dpi = 130,
  bg = "white"
)

# Figures d and e: Spearman rho, donor median over donor points -------------------------------
rho_frame <- rho_pairs |>
  mutate(arm = factor(arm, levels = rev(ARM_ORDER)))
rho_median <- rho_frame |>
  summarise(rho = median(rho), .by = c(arm, parent, pair, kind))
fig_d <- dplyr::filter(rho_frame, kind == "cytokine pair")
fig_e <- dplyr::filter(rho_frame, kind == "expected link")
write_csv(fig_d, file.path(run_dir, "fig_d_source.csv"))
write_csv(fig_e, file.path(run_dir, "fig_e_source.csv"))
p_d <- ggplot(fig_d, aes(rho, arm, colour = parent)) +
  geom_point(size = 0.8, alpha = 0.4, position = position_dodge(width = 0.6)) +
  geom_point(
    data = dplyr::filter(rho_median, kind == "cytokine pair"),
    shape = 18,
    size = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_colour_manual(
    values = c(CD4T = "#CC79A7", CD8T = "#0072B2"),
    name = NULL
  ) +
  facet_wrap(~pair, nrow = 2) +
  labs(x = "Spearman rho (diamond = donor median)", y = NULL) +
  rho_theme
ggsave(
  file.path(run_dir, "fig_d_rho_cytokine_pairs.png"),
  p_d,
  width = 15,
  height = 9,
  dpi = 130,
  bg = "white"
)
p_e <- ggplot(fig_e, aes(rho, arm, colour = parent)) +
  geom_point(size = 0.8, alpha = 0.4, position = position_dodge(width = 0.6)) +
  geom_point(
    data = dplyr::filter(rho_median, kind == "expected link"),
    shape = 18,
    size = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_colour_manual(
    values = c(CD4T = "#CC79A7", CD8T = "#0072B2"),
    name = NULL
  ) +
  facet_wrap(~pair, nrow = 1) +
  labs(x = "Spearman rho (diamond = donor median)", y = NULL) +
  rho_theme
ggsave(
  file.path(run_dir, "fig_e_rho_expected_links.png"),
  p_e,
  width = 11,
  height = 5.5,
  dpi = 130,
  bg = "white"
)

# Figure f: negative spread ratio to as188_joint --------------------------------------------
fig_f <- neg_spread |>
  mutate(cytokine = factor(cytokine, levels = names(CYTOKINES)))
write_csv(fig_f, file.path(run_dir, "fig_f_source.csv"))
p_f <- ggplot(fig_f, aes(arm, rsd_ratio)) +
  zip_shade +
  group_lines +
  geom_hline(yintercept = 1, colour = "grey50") +
  donor_layers +
  facet_grid(parent ~ cytokine) +
  labs(x = NULL, y = "Negative rSD / as188_joint")
ggsave(
  file.path(run_dir, "fig_f_neg_spread.png"),
  p_f,
  width = 16,
  height = 7.5,
  dpi = 130,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
