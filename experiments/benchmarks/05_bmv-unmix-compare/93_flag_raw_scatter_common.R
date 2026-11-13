# Purpose: copy of 93_flag_raw_scatter.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G3, G7a).
#          Raw per-cell view behind the all-marker flag figures (grill Q4, README row "Raw companions", user
#          2026-11-11). For every marker flagged (revised Q2 rule) in as102_var, s188_c102 or opt188: a ggcyto
#          hexbin of that marker against GLUT1 (Alexa Fluor 405-A) in the population where it is flagged (CD8T,
#          else CD4T, else live), on the common SpectroFlo labels, as the flag scores were. Columns are arms
#          (SpectroFlo and 1.8.8 zip as references), rows donors; the positive cut (cut_of, kind "pos"; fixed
#          arcsinh 1.5 in 3D, G7a) is drawn with geom_gate and the % above it printed per panel. The script stops unless FSC-A and SSC-A of
#          every root row are identical in every arm (G8).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/93_flag_raw_scatter_common.R <flags.csv>
#            (repo root); <flags.csv> is the new 86b flags.csv from batch 3. The out dir is run_out("common-flag-raw-scatter").
# Inputs:  the dataset config and run registry (datasets/<SW_DATASET>.json, <SW_DATASET>_runs.csv); run folders come
#          from the registry keys below, the sample table from samples_tbl()
#          <flags.csv>                                                        as102_var, s188_c102 flags [86b, batch 3]
#          run_in("common-downstream-opt188")/marker_flags.csv                opt188 flags [89 common]
#          run_in("common-gated")/gs/{as188_joint,opt188_pregate_joint_p1af1,as102_var,
#            spectroflo,s188_c102}   common-gate GatingSets [batch 0], README.md
#          cut_of(arm, sample, "<channel>", "pos")   positive cut per arm, sample and channel (fixed 1.5 in 3D, G7a)
# Outputs: <out_dir>/fig_flag_<marker>_<pop>.{png,pdf,csv}   one figure per flagged marker, 183 mm wide, 300 dpi PNG;
#                                                         with conditions: one per condition, suffix _<condition>
#          <out_dir>/flagged_markers.csv                 marker, population shown, arms that flag it
#          <out_dir>/pct_above_cut.csv                   % above the cut per marker, arm and donor
#          <out_dir>/pairing.csv                         arm, donor, root rows, FSC-A and SSC-A equal to 1.8.8 zip
#                                                        (donor = sample id; with conditions, a condition column too)
#          <out_dir>/sessionInfo.txt

# dataset.R loads only dplyr and readr, so it can run before .libPaths puts ggplot2 3.5.2 (for ggcyto) first
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB_GGPLOT2 <- grep("lib_ggplot2", cfg("lib_dirs"), value = TRUE, fixed = TRUE)
stopifnot("one ggplot2 lib dir in lib_dirs" = length(LIB_GGPLOT2) == 1)
.libPaths(c(LIB_GGPLOT2, .libPaths()))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})

# dataset.R calls a bare filter(); here flowCore's filter generic comes first on the search path, so give
# its two helpers a scope where filter is dplyr's (dataset.R should say dplyr::filter itself)
.dataset_scope <- new.env(parent = globalenv())
.dataset_scope$filter <- dplyr::filter
environment(run_in) <- .dataset_scope
environment(cut_of) <- .dataset_scope

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- run_out("common-flag-raw-scatter")
FLAGS_86B <- commandArgs(trailingOnly = TRUE)[1] # new 86b flags.csv, batch 3
FLAGS_89 <- file.path(run_in("common-downstream-opt188"), "marker_flags.csv")
COMMON <- run_in("common-gated") # batch 0; gs/<arm> for all arms
HEIGHT_PER_ROW_MM <- 35 # figure height per sample row; 5 rows give the 175 mm of the 3D figure
FIG_WIDTH_MM <- 183 # double column, the width save_panel takes besides 89
PAIR_COLS <- c("FSC-A", "SSC-A") # equal per row means every arm holds the same cells (G8)
# column label = GatingSet; references first, then the flagged arms
ARMS <- file.path(
  COMMON,
  "gs",
  c(
    SpectroFlo = "spectroflo",
    "1.8.8 zip" = "as188_joint",
    opt188 = "opt188_pregate_joint_p1af1",
    "1.0.2 code,\n1.8.8 spectra" = "s188_c102",
    "1.0.2 zip" = "as102_var"
  )
) |>
  set_names(c(
    "SpectroFlo",
    "1.8.8 zip",
    "opt188",
    "1.0.2 code,\n1.8.8 spectra",
    "1.0.2 zip"
  ))
FLAG_ARM <- c(
  as102_var = "1.0.2 zip",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  opt188 = "opt188"
)
ARM_FILL <- c("#252525", "#d9d9d9", "#9ecae1", "#08519c", "#7f0000") # as the 72 scatters
DARK_FILLS <- c("#252525", "#08519c", "#7f0000") # strips that take white text
POP_PRIORITY <- c("CD8T", "CD4T", "live") # population shown when a marker is flagged in several
COFACTOR <- 6000 # README "Transform"
# 19a_gate_arm.R:113 writes these channels as arcsinh already; all other fluorescence channels are linear
STORED_ARCSINH <- unlist(
  cfg("channels")[c("CD3", "CD4", "CD8", "TNF_slot", "IFNg_slot")],
  use.names = FALSE
)
GLUT1 <- paste0(cfg("targets")$glut1_dye, "-A")
AXIS <- c(-1.5, 4) # display window only; the % above the cut uses every cell
BINS <- 60 # hexagons per panel; display only

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "pass the new 86b flags.csv as the first argument" = !is.na(FLAGS_86B) &&
    file.exists(FLAGS_86B),
  "89 common run finished" = file.exists(FLAGS_89)
)
dir.create(OUT_DIR)
source(file.path(
  Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib"),
  "theme.R"
)) # save_panel, FIG_FONT

# Flagged markers and the population each is shown in ----------------------------------------------------
flags <- bind_rows(
  read_csv(FLAGS_86B, show_col_types = FALSE) |>
    dplyr::filter(source == "screen", arm %in% names(FLAG_ARM)),
  read_csv(FLAGS_89, show_col_types = FALSE) |>
    mutate(arm = "opt188")
) |>
  dplyr::filter(flag)
flagged <- flags |>
  mutate(pop_rank = match(pop, POP_PRIORITY)) |>
  summarise(
    pop = POP_PRIORITY[min(pop_rank)],
    flagged_by = paste(sort(unique(arm)), collapse = "; "),
    .by = c(marker, channel)
  )
write_csv(flagged, file.path(OUT_DIR, "flagged_markers.csv"))
# a non-fixed positive rule has cuts only for its own channels; a flagged marker outside them has no cut,
# and which cut to use for it is the user's choice, so the run stops here and names the channels
POS_RULE <- cfg("pos_rule")
no_cut <- if (identical(POS_RULE$type, "fixed")) {
  character()
} else {
  setdiff(flagged$channel, POS_RULE$channels)
}
if (length(no_cut)) {
  stop(
    "93_flag_raw_scatter_common.R: no positive cut in dataset ",
    DS$dataset,
    " for flagged channel(s) ",
    paste(no_cut, collapse = ", "),
    "; pos_rule has cuts for ",
    paste(POS_RULE$channels, collapse = ", "),
    "; the user must decide how these markers are cut",
    call. = FALSE
  )
}

# Live cells of each arm and donor, rows paired on the common labels ------------------------------------
gs_list <- map(ARMS, load_gs)
gs_ref <- gs_list[["1.8.8 zip"]]
smp <- samples_tbl()
donors <- sampleNames(gs_ref)
stopifnot(
  "gating set samples are the sample table" = setequal(donors, smp$sample)
)
# arm key of each column, as cuts.csv names it
ARM_KEY <- setNames(basename(ARMS), names(ARMS))
# positive cut of every arm, sample and flagged channel, read once (cut_of re-reads cuts.csv per call)
pos_cuts <- tidyr::expand_grid(
  arm = names(ARMS),
  donor = donors,
  channel = unique(flagged$channel)
) |>
  mutate(
    cut = pmap_dbl(list(arm, donor, channel), \(a, d, ch) {
      cut_of(ARM_KEY[[a]], d, ch, "pos")
    })
  )
POS_CUT <- set_names(
  pos_cuts$cut,
  paste(pos_cuts$arm, pos_cuts$donor, pos_cuts$channel, sep = "|")
)
pos_cut <- \(arm, donor, channel) {
  POS_CUT[[paste(arm, donor, channel, sep = "|")]]
}
# with conditions, a condition column joins the sample-level tables; 3D and EXV keep their columns
add_condition <- \(d) {
  if (is.null(DS$conditions)) {
    return(d)
  }
  left_join(
    d,
    select(smp, donor = sample, condition),
    by = join_by(donor),
    relationship = "many-to-one",
    unmatched = "error"
  )
}
channels <- unique(c(GLUT1, flagged$channel))
live <- map(set_names(donors), \(d) {
  which(gh_pop_get_indices(gs_ref[[d]], "live"))
})
# population membership within the live rows, for each population shown
in_pop <- map(set_names(donors), \(d) {
  map(set_names(POP_PRIORITY), \(p) {
    gh_pop_get_indices(gs_ref[[d]], p)[live[[d]]]
  })
})
cells <- imap(gs_list, \(gs, arm) {
  map(set_names(donors), \(d) {
    m <- exprs(gh_pop_get_data(gs[[d]], "root"))[
      live[[d]],
      c("FSC-A", channels),
      drop = FALSE
    ]
    linear <- setdiff(channels, STORED_ARCSINH)
    m[, linear] <- asinh(m[, linear] / COFACTOR)
    m
  })
})
# scatter does not depend on unmixing, so equal FSC-A and SSC-A on every root row show every arm holds the same cells
root_scatter <- \(gs, d) exprs(gh_pop_get_data(gs[[d]], "root"))[, PAIR_COLS]
pairing <- gs_list |>
  imap(\(gs, arm) {
    tibble(
      donor = donors,
      n_rows = map_int(donors, \(d) nrow(root_scatter(gs, d))),
      paired = map_lgl(donors, \(d) {
        identical(root_scatter(gs, d), root_scatter(gs_ref, d))
      })
    )
  }) |>
  list_rbind(names_to = "arm")
write_csv(
  add_condition(mutate(pairing, arm = gsub("\n", " ", arm))),
  file.path(OUT_DIR, "pairing.csv")
)
stopifnot("rows pair across arms" = all(pairing$paired))

# % above the cut per flagged marker, arm and donor --------------------------------------------------------
pct_above <- flagged |>
  select(which_marker = marker, which_channel = channel, which_pop = pop) |>
  pmap(\(which_marker, which_channel, which_pop) {
    tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
      mutate(
        marker = which_marker,
        pop = which_pop,
        pct_above_cut = map2_dbl(arm, donor, \(a, d) {
          100 *
            mean(
              cells[[a]][[d]][in_pop[[d]][[which_pop]], which_channel] >=
                pos_cut(a, d, which_channel)
            )
        })
      )
  }) |>
  list_rbind()
write_csv(
  add_condition(mutate(pct_above, arm = gsub("\n", " ", arm))),
  file.path(OUT_DIR, "pct_above_cut.csv")
)

# One ggcyto figure per flagged marker ------------------------------------------------------------------
# with conditions, one figure per condition: its samples are the rows
flag_figure <- \(which_marker, channel, pop, flagged_by, set_donors, suffix) {
  frames <- tidyr::expand_grid(arm = names(ARMS), donor = set_donors) |>
    mutate(name = paste(arm, donor, sep = "|")) |>
    with(set_names(
      map2(arm, donor, \(a, d) {
        flowFrame(cells[[a]][[d]][in_pop[[d]][[pop]], c(GLUT1, channel)])
      }),
      name
    ))
  fs <- flowSet(frames)
  pd <- tibble(name = sampleNames(fs)) |>
    tidyr::separate_wider_delim(
      name,
      "|",
      names = c("arm", "donor"),
      cols_remove = FALSE
    ) |>
    mutate(arm = factor(arm, levels = names(ARMS)))
  pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
  # one gate per panel: the cut belongs to the arm and sample
  gate <- pd |>
    with(set_names(
      map2(as.character(arm), donor, \(a, d) {
        rectangleGate(
          filterId = "above cut",
          .gate = setNames(
            list(c(-Inf, Inf), c(pos_cut(a, d, channel), Inf)),
            c(GLUT1, channel)
          )
        )
      }),
      name
    )) |>
    filterList()
  labels <- dplyr::filter(
    pct_above,
    marker == which_marker,
    donor %in% set_donors
  ) |>
    mutate(arm = factor(arm, levels = names(ARMS)))
  flagged_labels <- FLAG_ARM[strsplit(flagged_by, "; ")[[1]]]
  p <- ggcyto(fs, aes(x = !!rlang::sym(GLUT1), y = !!rlang::sym(channel))) +
    geom_hex(bins = BINS) +
    geom_gate(gate, colour = "red") +
    geom_label(
      data = labels,
      aes(
        x = -Inf,
        y = Inf,
        label = sprintf("%s+ %.1f%%", which_marker, pct_above_cut)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 7,
      size.unit = "pt",
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0
    ) +
    facet_grid(donor ~ arm, switch = "y") +
    ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(
      x = sprintf("GLUT1, %s (arcsinh, cofactor 6000)", GLUT1),
      y = sprintf("%s, %s (arcsinh, cofactor 6000)", which_marker, channel),
      title = sprintf(
        "%s: %s vs GLUT1 (flagged in %s)",
        pop,
        which_marker,
        paste(gsub("\n", " ", flagged_labels), collapse = ", ")
      )
    ) +
    ggpubr::theme_pubr(base_size = 7, base_family = FIG_FONT, border = TRUE) +
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 7, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  g <- ggplotGrob(as.ggplot(p))
  # colour each column strip by its arm; strip-t grobs run left to right in facet order
  strips <- grep("^strip-t", g$layout$name)
  strips <- strips[order(g$layout$l[strips])]
  for (i in seq_along(strips)) {
    rect <- g$grobs[[strips[i]]]$grobs[[1]]$children[[1]]
    rect$gp$fill <- ARM_FILL[i]
    g$grobs[[strips[i]]]$grobs[[1]]$children[[1]] <- rect
    txt <- g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]
    txt$gp$col <- if (ARM_FILL[i] %in% DARK_FILLS) "white" else "black"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]] <- txt
  }
  save_panel(
    g,
    mutate(labels, arm = gsub("\n", " ", arm)),
    OUT_DIR,
    sprintf(
      "fig_flag_%s_%s%s",
      gsub("[^A-Za-z0-9]", "", which_marker),
      pop,
      suffix
    ),
    width_mm = FIG_WIDTH_MM,
    height_mm = HEIGHT_PER_ROW_MM * length(set_donors)
  )
}
# rows per figure: all samples without conditions, else the samples of each condition
sample_sets <- if (is.null(DS$conditions)) {
  list(all = donors)
} else {
  split(
    smp$sample[smp$sample %in% donors],
    smp$condition[smp$sample %in% donors]
  )
}
tidyr::expand_grid(
  flagged |> select(which_marker = marker, channel, pop, flagged_by),
  set = names(sample_sets)
) |>
  mutate(
    set_donors = sample_sets[set],
    suffix = if (is.null(DS$conditions)) "" else paste0("_", set)
  ) |>
  select(-set) |>
  pwalk(flag_figure)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
