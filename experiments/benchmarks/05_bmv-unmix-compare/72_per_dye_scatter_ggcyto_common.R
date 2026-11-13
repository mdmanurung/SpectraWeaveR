# Purpose: copy of 72_per_dye_scatter_ggcyto.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          G2, G3, G7a, G11). Cytokine scatter per arm set with ggcyto: hexbin density (ggcyto's default
#          pseudocolor, reversed Brewer Spectral, sqrt counts), the positive cut (cut_of(): 3D arcsinh 1.5, G7a) drawn with
#          geom_gate on every panel and the % above it as text, column strips coloured by arm (user,
#          2026-11-10). CD8T cells carry the common SpectroFlo CD8T label in every arm, so each column holds
#          the same cells; the only line drawn is the cut. The script stops unless FSC-A and SSC-A of every
#          row are identical in every arm (G8). Rows are donors and columns arms, so one donor reads across
#          methods. Blocks side by side, the y marker against each x marker. Arm sets:
#            per_dye  y TNF; 1.0.2 variants swapped in one dye at a time, then SBUV510 + AF405 (45); x: IFNg, GLUT1, SSC
#            af_swap  y TNF; AF spectra and variants swapped between versions (37, 48); x: IFNg, AF, SSC
#            tnf_new  y TNF; arms from tests 3, Q4/Q9 and the downstream comparison; x: IFNg, GLUT1, SSC
#            th2      y Th2 (eFluor 450); arms that squeeze Th2 and the AF405 knock-out; x: GLUT1, IFNg, SSC
#          The set th2 needs the Th2 dye (config targets th2_dye); a set that needs a target the dataset lacks stops.
#          24H: rows are condition and sample; pct_above_cut.csv and the figure source csv gain condition.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/72_per_dye_scatter_ggcyto_common.R <set>
#            from the repo root. <set> is per_dye, af_swap, tnf_new or th2.
# Inputs:  registry key common-gated (datasets/<ds>_runs.csv): gs/<arm>/ common-gate GatingSets [batch 0], README.md
#          registry key common-cuts, only when the dataset's pos_rule is not fixed: cuts.csv [97e]
#          datasets/<ds>.json (channels, targets, lib_dirs, pos_rule) and datasets/<ds>_samples.csv [00_samples.R]
#          BV510/BV570 stored arcsinh; Alexa Fluor 405, AF-A and eFluor 450 stored linear
# Outputs: <out_dir>/fig_<y>_vs_<x>.{pdf,png,csv}      one figure per x marker, 183 mm wide, 300 dpi PNG
#          <out_dir>/pct_above_cut.csv                  % of CD8T above the cut, per arm and donor (sample)
#          <out_dir>/pairing.csv                        per arm and donor: FSC-A and SSC-A equal to SpectroFlo
#          <out_dir>/sessionInfo.txt
#          <out_dir> is run_out("common-72-<set>"), a new folder dated SW_RUN_DATE.

# dataset.R loads dplyr and readr only, so the library path can still change before ggplot2 loads
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of(), samples_tbl()
LIB_GGPLOT2 <- grep("lib_ggplot2_35$", cfg("lib_dirs"), value = TRUE)
stopifnot("one lib_ggplot2_35 in config lib_dirs" = length(LIB_GGPLOT2) == 1)
.libPaths(c(LIB_GGPLOT2, .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  # dataset.R attached dplyr below flowCore. Attach it again on top, with filter(), because the dataset.R
  # helpers call filter() and flowCore's filter would win (this script itself never calls filter()).
  detach("package:dplyr", unload = FALSE)
  library(dplyr)
  library(purrr)
  library(readr)
})
source(file.path(
  Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib"),
  "theme.R"
)) # save_panel, FIG_FONT

# Parameters ---------------------------------------------------------------------------------------
GLUT1_DYE_ID <- gsub("[ -]", "_", cfg("targets")$glut1_dye) # as in the arm folder names
GLUT1_ARM <- paste0("var102_only_", GLUT1_DYE_ID) # the GLUT1 dye with its 1.0.2 variants only
SET <- commandArgs(trailingOnly = TRUE)[1]
COMMON <- run_in("common-gated") # batch 0; gs/<arm> for all arms
# arms per set (label = arm folder), strip colour per arm, y marker, x blocks, output folder
SETS <- list(
  per_dye = list(
    arms = c(
      "1.8.8\nvariants" = "as188_joint",
      "opt188\n(optimized)" = "opt188_pregate_joint_p1af1",
      "+SBUV510\n1.0.2" = "var102_only_SBUV510",
      "+ AF405\n1.0.2" = GLUT1_ARM,
      "+SBUV510\n+AF405 1.0.2" = paste0("var102_only_SBUV510__", GLUT1_DYE_ID),
      "all 1.0.2\nvariants" = "swap_var102"
    ),
    # grey baseline, light to dark red as more of the TNF squeeze is swapped in
    fill = c("#d9d9d9", "#9ecae1", "#fdbb84", "#e34a33", "#b30000", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    out_key = "common-72-per_dye"
  ),
  af_swap = list(
    arms = c(
      "1.8.8\nall" = "as188_joint",
      "+AF\n1.0.2" = "swap_af102",
      "+variants\n1.0.2" = "swap_var102",
      "+AF+var\n1.0.2" = "swap_af102_var102",
      "1.0.2, AF\n1.8.8" = "rev_af188",
      "1.0.2\nall" = "s102_c188"
    ),
    # 1.8.8 spectra: grey to dark red as 1.0.2 parts go in; 1.0.2 spectra: blues
    fill = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000", "#9ecae1", "#08519c"),
    y = "TNF",
    blocks = c("IFNg", "AF", "SSC"),
    out_key = "common-72-af_swap"
  ),
  tnf_new = list(
    arms = c(
      "1.8.8\nzip" = "as188_joint",
      "opt188" = "opt188_pregate_joint_p1af1",
      "refined\nAF" = "af405_afrefine",
      "AF405 1.0.2\ncode, clean" = "var102code_clean_af405",
      "1.0.2 code\n1.8.8 spectra" = "s188_c102",
      "1.0.2\nzip" = "as102_var"
    ),
    # 1.8.8 arms grey/blue, 1.0.2-code or 1.0.2-variant arms red
    fill = c("#d9d9d9", "#9ecae1", "#08519c", "#fdbb84", "#e34a33", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    out_key = "common-72-tnf_new"
  ),
  th2 = list(
    arms = c(
      "1.8.8\nzip" = "as188_joint",
      "opt188" = "opt188_pregate_joint_p1af1",
      "+ AF405\n1.0.2" = GLUT1_ARM,
      "AF405 knock-out\n(no BV510 part)" = "af405proof_ko_bv510",
      "1.0.2 code\n1.8.8 spectra" = "s188_c102",
      "1.0.2\nzip" = "as102_var"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#e34a33", "#fdbb84", "#08519c", "#7f0000"),
    y = "Th2",
    blocks = c("GLUT1", "IFNg", "SSC"),
    out_key = "common-72-th2"
  )
)
# SpectroFlo 3.3.0 export, first column of every set as the vendor reference (user, 2026-11-10)
REFERENCE <- c("SpectroFlo\n(reference)" = "spectroflo")
REFERENCE_FILL <- "#252525"
DARK_FILLS <- c("#e34a33", "#b30000", "#7f0000", "#08519c", REFERENCE_FILL) # strips that take white text
POP <- "CD8T"
# positive cut of the y marker: cut_of(arm, sample, channel, "pos"), arcsinh, cofactor 6000 (G7a)
COFACTOR <- 6000 # arcsinh cofactor for linear channels, as for the cytokines (CLAUDE.md)
AXIS <- c(-1.5, 4) # display window only; the % above the cut uses every event
BINS <- 60 # hexagon bins per panel; display only, fewer for the ~27 mm panels
PANEL_MM <- 23 # panel column width; seven columns fill a 183 mm double column
LABEL_PT <- 7 # in-panel % and strip text; save_panel allows 5-7 pt
FIG_WIDTH_MM <- 183 # double column, the only width save_panel takes besides 89
PAIR_COLS <- c("FSC-A", "SSC-A") # scatter does not depend on unmixing; equal per row means the same cells
# Th2 only where the panel has the dye (config targets th2_dye); a set that needs a missing target stops below
CHANNELS <- c(
  TNF = cfg("channels")$TNF_slot,
  if (!is.null(cfg("targets")$th2_dye)) {
    c(Th2 = paste0(cfg("targets")$th2_dye, "-A"))
  },
  IFNg = cfg("channels")$IFNg_slot,
  GLUT1 = paste0(cfg("targets")$glut1_dye, "-A"),
  AF = "AF-A",
  SSC = "SSC-A"
)
LINEAR_TO_ARCSINH <- unname(CHANNELS[intersect(
  c("GLUT1", "AF", "Th2"),
  names(CHANNELS)
)])
# display names of the keys above: the dataset's marker labels for the TNF and IFNg dye slots
LABELS <- c(
  TNF = cfg("targets")$tnf_label,
  Th2 = "Th2",
  IFNg = cfg("targets")$ifng_label,
  GLUT1 = "GLUT1",
  AF = "AF",
  SSC = "SSC"
)
SCALE_LABEL <- c(
  TNF = "arcsinh",
  Th2 = "arcsinh, cofactor 6000",
  IFNg = "arcsinh, cofactor 6000",
  GLUT1 = "arcsinh, cofactor 6000",
  AF = "arcsinh, cofactor 6000",
  SSC = "linear, x 10^6"
)
# x display window per block; SSC 0-2.5e6 holds the CD8T cloud (~1e6) without the few outliers to 4e6
X_LIMITS <- list(IFNg = AXIS, GLUT1 = AXIS, AF = AXIS, SSC = c(0, 2.5e6))

stopifnot("set is one of the arm sets" = SET %in% names(SETS))
ARMS <- c(REFERENCE, SETS[[SET]]$arms)
ARM_FILL <- set_names(c(REFERENCE_FILL, SETS[[SET]]$fill), names(ARMS))
Y_NAME <- SETS[[SET]]$y
BLOCKS <- SETS[[SET]]$blocks
absent_targets <- setdiff(c(Y_NAME, BLOCKS), names(CHANNELS))
if (length(absent_targets) > 0) {
  stop(
    sprintf(
      "72_per_dye_scatter_ggcyto_common.R: set %s needs %s, missing in dataset %s",
      SET,
      paste(absent_targets, collapse = ", "),
      DS$dataset
    ),
    call. = FALSE
  )
}
Y_CH <- CHANNELS[[Y_NAME]]
Y_LABEL <- LABELS[[Y_NAME]]
OUT_DIR <- run_out(SETS[[SET]]$out_key)
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: rows by condition
# a fixed cut is the same for every sample; the table records it only when it differs by sample
CUT_COL <- if (identical(cfg("pos_rule")$type, "fixed")) {
  character()
} else {
  "pos_cut"
}

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "every arm has a strip colour" = setequal(names(ARMS), names(ARM_FILL))
)
dir.create(OUT_DIR, recursive = TRUE)

# Row pairing: FSC-A and SSC-A of every root row equal to SpectroFlo's, per arm and donor (G8) ---------
gs_list <- map(ARMS, \(a) load_gs(file.path(COMMON, "gs", a)))
root_scatter <- \(gs, s) exprs(gh_pop_get_data(gs[[s]], "root"))[, PAIR_COLS]
pairing <- gs_list |>
  imap(\(gs, label) {
    tibble(donor = sampleNames(gs)) |>
      mutate(
        arm = label,
        n_rows = map_int(donor, \(s) nrow(root_scatter(gs, s))),
        paired = map_lgl(donor, \(s) {
          identical(
            root_scatter(gs, s),
            root_scatter(gs_list[[names(REFERENCE)]], s)
          )
        })
      )
  }) |>
  list_rbind()
write_csv(
  mutate(pairing, arm = gsub("\n", " ", arm)),
  file.path(OUT_DIR, "pairing.csv")
)
stopifnot(
  "every arm pairs with SpectroFlo on FSC-A and SSC-A" = all(pairing$paired)
)

# One flowSet of CD8T events, arm x donor, linear channels on arcsinh --------------------------------
frames <- gs_list |>
  imap(\(gs, label) {
    sampleNames(gs) |>
      set_names(\(s) paste(label, s, sep = "|")) |>
      map(\(s) {
        cytoframe_to_flowFrame(gh_pop_get_data(gs[[s]], POP))[, unname(
          CHANNELS
        )]
      })
  }) |>
  list_flatten(name_spec = "{inner}")
fs <- flowSet(frames)
fs <- transform(
  fs,
  transformList(
    LINEAR_TO_ARCSINH,
    arcsinhTransform(a = 0, b = 1 / COFACTOR, c = 0)
  )
)
pd <- tibble(name = sampleNames(fs)) |>
  tidyr::separate_wider_delim(
    name,
    "|",
    names = c("arm", "donor"),
    cols_remove = FALSE
  ) |>
  mutate(arm = factor(arm, levels = names(ARMS)))
# positive cut of the y marker per arm and sample; 24H adds the condition of each sample
pd <- pd |>
  mutate(
    pos_cut = map2_dbl(ARMS[as.character(arm)], donor, \(a, s) {
      cut_of(a, s, Y_CH, "pos")
    })
  )
if (!is.null(DS$conditions)) {
  pd <- pd |>
    left_join(
      select(samples_tbl(), donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = "error"
    )
}
pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)

# % of CD8T above the cut per arm and donor, printed in each panel ----------------------------------
pct_pos <- pd |>
  mutate(
    pct_above_cut = map2_dbl(name, pos_cut, \(n, ct) {
      100 * mean(exprs(fs[[n]])[, Y_CH] >= ct)
    })
  ) |>
  select(name, pct_above_cut, arm, donor, all_of(COND_BY), all_of(CUT_COL))
write_csv(
  pct_pos |>
    mutate(
      arm_id = ARMS[as.character(arm)],
      arm = gsub("\n", " ", arm),
      marker = Y_LABEL
    ) |>
    select(
      arm_id,
      arm,
      donor,
      all_of(COND_BY),
      marker,
      pct_above_cut,
      all_of(CUT_COL)
    ),
  file.path(OUT_DIR, "pct_above_cut.csv")
)

# One ggcyto block per x marker: hexbin, fixed cut, %, donor x arm facets ----------------------------
ROW_FACET <- paste(c(COND_BY, "donor"), collapse = " + ") # 24H: condition, then sample
block <- function(x_name) {
  x_ch <- CHANNELS[[x_name]]
  # ggcyto on a flowSet takes one gate per sample, named by sample; the cut is that sample's positive cut
  gate <- filterList(
    set_names(
      map(pd$pos_cut, \(ct) {
        rectangleGate(
          filterId = "above cut",
          .gate = setNames(list(c(-Inf, Inf), c(ct, Inf)), c(x_ch, Y_CH))
        )
      }),
      pd$name
    )
  )
  # % from pct_pos, every event (geom_stats fails on the SSC block, so all blocks use this label)
  p <- ggcyto(fs, aes(x = !!rlang::sym(x_ch), y = !!rlang::sym(Y_CH))) +
    geom_hex(bins = BINS) +
    geom_gate(gate, colour = "red") +
    geom_label(
      data = pct_pos,
      aes(
        x = -Inf,
        y = Inf,
        label = sprintf("%s+ %.1f%%", Y_LABEL, pct_above_cut)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = LABEL_PT,
      size.unit = "pt",
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = unit(0.12, "lines")
    ) +
    facet_grid(as.formula(paste(ROW_FACET, "~ arm")), switch = "y") +
    ggcyto_par_set(limits = list(x = X_LIMITS[[x_name]], y = AXIS)) +
    # standard flow pseudocolor, as ggcyto's own default: blue (few events) to red (many)
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(
      x = sprintf("%s, %s (%s)", LABELS[[x_name]], x_ch, SCALE_LABEL[[x_name]]),
      y = sprintf("%s, %s (%s)", Y_LABEL, Y_CH, SCALE_LABEL[[Y_NAME]]),
      title = sprintf("CD8 T cells: %s vs %s", Y_LABEL, LABELS[[x_name]])
    ) +
    ggpubr::theme_pubr(base_size = 7, base_family = FIG_FONT, border = TRUE) +
    # donor strips left, outside, no box (house style); arm strips keep their fill
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = LABEL_PT),
      strip.text.y.left = element_text(
        size = LABEL_PT,
        face = "bold",
        angle = 0
      ),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      axis.text = element_text(size = LABEL_PT),
      panel.spacing = unit(1, "mm")
    )
  # SSC ticks in millions, so the panels' labels do not overlap
  if (x_name == "SSC") {
    p <- p +
      scale_x_continuous(breaks = c(0, 1e6, 2e6), labels = c("0", "1", "2"))
  }
  g <- ggplotGrob(as.ggplot(p))
  # colour each column strip by its arm; strip-t grobs run left to right in facet order
  strips <- grep("^strip-t", g$layout$name)
  strips <- strips[order(g$layout$l[strips])]
  for (i in seq_along(strips)) {
    fill <- ARM_FILL[[levels(pd$arm)[i]]]
    rect <- g$grobs[[strips[i]]]$grobs[[1]]$children[[1]]
    rect$gp$fill <- fill
    g$grobs[[strips[i]]]$grobs[[1]]$children[[1]] <- rect
    txt <- g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]
    txt$gp$col <- if (fill %in% DARK_FILLS) "white" else "black"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]] <- txt
  }
  g
}

# Figure ---------------------------------------------------------------------------------------------
# one figure per x marker; 183 mm double column (save_panel takes 89 or 183 mm)
fig_height_mm <- PANEL_MM * n_distinct(pd$donor) + 35
walk(BLOCKS, \(x_name) {
  save_panel(
    block(x_name),
    mutate(pct_pos, arm = gsub("\n", " ", arm), x_marker = LABELS[[x_name]]),
    OUT_DIR,
    sprintf("fig_%s_vs_%s", Y_LABEL, LABELS[[x_name]]),
    width_mm = FIG_WIDTH_MM,
    height_mm = fig_height_mm
  )
})

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
