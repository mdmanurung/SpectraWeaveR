# Synthetic references: the sw_synref class and the draw functions -----------------
#
# A synthetic reference is built from real study events. It is a diagnostic
# comparator, a null arm or a benchmark base. It is never a measured technical
# reference, and it never feeds the paired_shift or HEB bridge path (DR3, DR4).
# Cells are not replicates: donors and aliquots are.

.sw_synref_schema <- "sw_synref/1"

# Class -------------------------------------------------------------------------

# Build the object and stamp its digest.
.sw_new_synref <- function(events, provenance, recipe, realised, scope, floor) {
  obj <- structure(
    list(
      schema_version = .sw_synref_schema,
      events = events,
      provenance = provenance,
      recipe = recipe,
      realised = realised,
      scope = scope,
      floor = floor,
      certainty = "diagnostic_only"
    ),
    class = "sw_synref"
  )
  obj$digest <- .sw_model_digest(obj)
  obj
}

#' @export
print.sw_synref <- function(x, ...) {
  cat(sprintf(
    "<sw_synref> %d reference(s), %d events, scope %s, certainty %s\n",
    length(x[["events"]]),
    sum(vapply(x[["events"]], nrow, integer(1))),
    x[["scope"]],
    x[["certainty"]]
  ))
  if (!is.na(x[["floor"]])) {
    cat("floor: ", x[["floor"]], "\n", sep = "")
  }
  cat("Cells are not replicates; donors and aliquots are.\n")
  invisible(x)
}

# Shared argument checks --------------------------------------------------------

.sw_synref_check_data <- function(data) {
  ok <- is.list(data) &&
    !is.data.frame(data) &&
    length(data) > 0 &&
    !is.null(names(data)) &&
    !anyNA(names(data)) &&
    !anyDuplicated(names(data)) &&
    all(vapply(data, \(m) is.matrix(m) && is.numeric(m), logical(1)))
  if (!ok) {
    stop(
      "'data' must be a uniquely named list of numeric matrices, one per file.",
      call. = FALSE
    )
  }
  cn <- lapply(data, colnames)
  if (any(vapply(cn, is.null, logical(1))) || length(unique(cn)) != 1) {
    stop("Every matrix in 'data' needs the same column names.", call. = FALSE)
  }
  invisible(NULL)
}

.sw_synref_check_manifest <- function(manifest, ids) {
  if (
    !is.data.frame(manifest) ||
      !all(c("sample_id", "batch") %in% names(manifest))
  ) {
    stop("'manifest' needs columns 'sample_id' and 'batch'.", call. = FALSE)
  }
  if (anyDuplicated(manifest[["sample_id"]])) {
    stop("'manifest$sample_id' must be unique.", call. = FALSE)
  }
  absent <- setdiff(ids, as.character(manifest[["sample_id"]]))
  if (length(absent) > 0) {
    stop(
      "Sample(s) not in 'manifest': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.sw_synref_check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != round(x)) {
    stop(
      sprintf("'%s' must be one positive whole number.", name),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.sw_synref_check_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("'seed' must be one finite number.", call. = FALSE)
  }
  invisible(NULL)
}

# Per-event table of one group: sample, row, batch, donor, population.
.sw_synref_events_table <- function(data, manifest, ids, labels) {
  man <- as.data.frame(manifest)
  rownames(man) <- as.character(man[["sample_id"]])
  pieces <- lapply(ids, function(id) {
    n <- nrow(data[[id]])
    donor <- if ("donor" %in% names(man)) {
      as.character(man[id, "donor"])
    } else {
      NA_character_
    }
    pop <- if (is.null(labels)) {
      rep(NA_character_, n)
    } else {
      as.character(labels[[id]])
    }
    data.frame(
      source_sample = id,
      source_row = seq_len(n),
      batch = as.character(man[id, "batch"]),
      donor = donor,
      population = pop,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

# Rebuild the event matrix from provenance rows, in provenance order.
.sw_synref_matrix <- function(data, prov) {
  markers <- colnames(data[[1]])
  out <- matrix(
    NA_real_,
    nrow(prov),
    length(markers),
    dimnames = list(NULL, markers)
  )
  for (s in unique(prov[["source_sample"]])) {
    sel <- which(prov[["source_sample"]] == s)
    out[sel, ] <- data[[s]][prov[["source_row"]][sel], , drop = FALSE]
  }
  out
}

# Draw from real events -----------------------------------------------------------

#' Draw a synthetic reference from real study events
#'
#' Draws events without replacement from the files of each group, with an
#' equal number of events per stratum, and returns one synthetic reference per
#' group and draw. A group must sit inside one batch, because a reference drawn
#' across batches would carry the batch effect it is meant to compare against
#' (within-batch only). The function stops when a stratum has fewer events than
#' requested. It never filters silently. The result is a diagnostic comparator
#' with certainty `"diagnostic_only"`. Cells are not replicates; donors and
#' aliquots are.
#'
#' @param data Named list of numeric event matrices, one per file, with the
#'   same column names. The caller transforms them first; nothing is
#'   transformed here.
#' @param manifest Data frame with `sample_id` and `batch`, and `donor` when
#'   `strata = "donor"`.
#' @param groups Named list. Each element is a character vector of
#'   `sample_id`s that feed one reference. The name becomes the reference id.
#'   All files of a group must share one batch.
#' @param strata One of `"none"`, `"file"`, `"donor"` or `"population"`
#'   (uncalibrated). `"none"` is one pooled stratum per group: `n_events` are
#'   drawn from all its files together, so the share per file is random and
#'   unequal. Only `"file"` gives equal events per file. `"donor"` and
#'   `"population"` give one stratum per donor or population label, with
#'   `n_events` drawn from each.
#' @param labels Named list of population labels, one vector per file, same
#'   length as its rows. Required when `strata = "population"`. `NULL` is an
#'   explicit "no labels".
#' @param n_events Events drawn per stratum (uncalibrated). The function stops
#'   when a stratum is short.
#' @param n_draws Number of independent draws per group (uncalibrated).
#' @param seed Random seed (uncalibrated). The caller's RNG state is restored
#'   afterwards.
#'
#' @return An `sw_synref` object: `events` (named list of matrices),
#'   `provenance` (synref_id, draw, source_sample, source_row, batch, donor,
#'   population, stratum), `recipe` (every argument), `realised` (requested
#'   and realised events per stratum), `scope`, `certainty` and `digest`.
#'
#' @examples
#' set.seed(1)
#' d <- list(
#'   a = cbind(M1 = rnorm(60), M2 = rnorm(60)),
#'   b = cbind(M1 = rnorm(60), M2 = rnorm(60))
#' )
#' man <- data.frame(sample_id = c("a", "b"), batch = "b1")
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(ref1 = c("a", "b")), strata = "file",
#'   labels = NULL, n_events = 20, n_draws = 2, seed = 1
#' )
#' ref$realised
#'
#' @export
sw_synref_downsample <- function(
  data,
  manifest,
  groups,
  strata,
  labels,
  n_events,
  n_draws,
  seed
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    groups = missing(groups),
    strata = missing(strata),
    labels = missing(labels),
    n_events = missing(n_events),
    n_draws = missing(n_draws),
    seed = missing(seed)
  ))
  recipe <- list(
    groups = groups,
    strata = strata,
    n_events = n_events,
    n_draws = n_draws,
    seed = seed,
    labels_given = !is.null(labels),
    markers = if (is.list(data) && length(data) > 0) {
      colnames(data[[1]])
    } else {
      NULL
    }
  )
  .sw_synref_check_data(data)
  .sw_synref_check_count(n_events, "n_events")
  .sw_synref_check_count(n_draws, "n_draws")
  .sw_synref_check_seed(seed)
  if (
    !is.character(strata) ||
      length(strata) != 1 ||
      !strata %in% c("none", "file", "donor", "population")
  ) {
    stop(
      "'strata' must be one of none, file, donor, population.",
      call. = FALSE
    )
  }
  if (
    !is.list(groups) ||
      length(groups) == 0 ||
      is.null(names(groups)) ||
      anyNA(names(groups)) ||
      anyDuplicated(names(groups)) ||
      !all(vapply(groups, \(g) is.character(g) && length(g) > 0, logical(1)))
  ) {
    stop(
      "'groups' must be a uniquely named list of sample_id vectors.",
      call. = FALSE
    )
  }
  repeated <- names(groups)[vapply(groups, anyDuplicated, integer(1)) > 0]
  if (length(repeated) > 0) {
    stop(
      "Group(s) list a sample_id more than once: ",
      paste(repeated, collapse = ", "),
      ". A repeated file would let one event enter a reference twice.",
      call. = FALSE
    )
  }
  all_ids <- unique(unlist(groups, use.names = FALSE))
  absent <- setdiff(all_ids, names(data))
  if (length(absent) > 0) {
    stop(
      "Group sample(s) not in 'data': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  .sw_synref_check_manifest(manifest, all_ids)
  if (strata == "population") {
    if (is.null(labels)) {
      stop("strata = 'population' needs 'labels'.", call. = FALSE)
    }
    ok <- all(vapply(
      all_ids,
      \(id) length(labels[[id]]) == nrow(data[[id]]),
      logical(1)
    ))
    if (!ok) {
      stop(
        "'labels' needs one vector per file with one value per event.",
        call. = FALSE
      )
    }
  }
  if (strata == "donor" && !"donor" %in% names(manifest)) {
    stop(
      "strata = 'donor' needs a 'donor' column in 'manifest'.",
      call. = FALSE
    )
  }
  man_batch <- stats::setNames(
    as.character(manifest[["batch"]]),
    as.character(manifest[["sample_id"]])
  )
  for (g in names(groups)) {
    b <- unique(man_batch[groups[[g]]])
    if (length(b) > 1) {
      stop(
        sprintf(
          "Group '%s' spans batches %s. A synthetic reference is within-batch only.",
          g,
          paste(b, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  drawn <- .sw_with_seed(seed, {
    out <- list()
    for (g in names(groups)) {
      ev <- .sw_synref_events_table(
        data,
        manifest,
        groups[[g]],
        if (strata == "population") labels else NULL
      )
      if (strata == "donor" && anyNA(ev[["donor"]])) {
        stop(
          sprintf("Group '%s' has files with a missing donor.", g),
          call. = FALSE
        )
      }
      if (strata == "population" && anyNA(ev[["population"]])) {
        stop(
          sprintf("Group '%s' has events with a missing population label.", g),
          call. = FALSE
        )
      }
      ev[["stratum"]] <- switch(
        strata,
        none = "all",
        file = ev[["source_sample"]],
        donor = ev[["donor"]],
        population = ev[["population"]]
      )
      for (k in seq_len(n_draws)) {
        key <- sprintf("%s_d%02d", g, k)
        picks <- lapply(sort(unique(ev[["stratum"]])), function(s) {
          idx <- which(ev[["stratum"]] == s)
          if (length(idx) < n_events) {
            stop(
              sprintf(
                "Stratum '%s' of group '%s' has %d events; %d requested.",
                s,
                g,
                length(idx),
                n_events
              ),
              call. = FALSE
            )
          }
          idx[sort(sample.int(length(idx), n_events))]
        })
        prov <- ev[unlist(picks), , drop = FALSE]
        prov <- cbind(synref_id = key, draw = k, prov)
        rownames(prov) <- NULL
        out[[key]] <- prov
      }
    }
    out
  })

  prov <- tibble::as_tibble(do.call(rbind, unname(drawn)))
  prov <- prov[c(
    "synref_id",
    "draw",
    "source_sample",
    "source_row",
    "batch",
    "donor",
    "population",
    "stratum"
  )]
  events <- lapply(drawn, \(p) .sw_synref_matrix(data, p))
  realised <- prov |>
    dplyr::summarise(
      realised = dplyr::n(),
      .by = c("synref_id", "draw", "stratum")
    ) |>
    dplyr::mutate(requested = as.integer(n_events), .after = "stratum")
  .sw_new_synref(
    events,
    prov,
    recipe,
    realised,
    scope = "within_batch",
    floor = NA_character_
  )
}

# Disjoint pseudo-aliquots ---------------------------------------------------------

#' Cut one real sample into disjoint pseudo-aliquots
#'
#' Splits the events of one file into disjoint parts of equal size. The parts
#' have identical biology, so differences between them are sampling noise
#' only. The result is labelled `floor = "sampling only"`: it is a floor for
#' event-level noise, not for aliquot-to-aliquot or run-to-run noise. Cells
#' are not replicates; donors and aliquots are.
#'
#' @param data Named list of numeric event matrices (already transformed by the
#'   caller).
#' @param manifest Data frame with `sample_id` and `batch`, and optionally
#'   `donor`.
#' @param source The `sample_id` to split.
#' @param n_parts Number of disjoint parts (uncalibrated).
#' @param events_per_part Events per part (uncalibrated). The function stops
#'   when `n_parts * events_per_part` exceeds the events of the source.
#' @param seed Random seed (uncalibrated). The caller's RNG state is restored.
#'
#' @return An `sw_synref` with one reference per part, named
#'   `<source>_part<k>`.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(200), M2 = rnorm(200)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' sw_synref_split(d, man, source = "a", n_parts = 4, events_per_part = 40, seed = 1)
#'
#' @export
sw_synref_split <- function(
  data,
  manifest,
  source,
  n_parts,
  events_per_part,
  seed
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    source = missing(source),
    n_parts = missing(n_parts),
    events_per_part = missing(events_per_part),
    seed = missing(seed)
  ))
  recipe <- list(
    source = source,
    n_parts = n_parts,
    events_per_part = events_per_part,
    seed = seed,
    markers = if (is.list(data) && length(data) > 0) {
      colnames(data[[1]])
    } else {
      NULL
    }
  )
  .sw_synref_check_data(data)
  .sw_synref_check_count(n_parts, "n_parts")
  .sw_synref_check_count(events_per_part, "events_per_part")
  .sw_synref_check_seed(seed)
  if (
    !is.character(source) || length(source) != 1 || !source %in% names(data)
  ) {
    stop("'source' must be one name from 'data'.", call. = FALSE)
  }
  .sw_synref_check_manifest(manifest, source)
  n_have <- nrow(data[[source]])
  if (n_parts * events_per_part > n_have) {
    stop(
      sprintf(
        "Source '%s' has %d events; %d parts of %d need %d.",
        source,
        n_have,
        n_parts,
        events_per_part,
        n_parts * events_per_part
      ),
      call. = FALSE
    )
  }
  ev <- .sw_synref_events_table(data, manifest, source, NULL)
  perm <- .sw_with_seed(seed, sample.int(n_have, n_parts * events_per_part))
  parts <- split(perm, rep(seq_len(n_parts), each = events_per_part))
  drawn <- lapply(seq_len(n_parts), function(k) {
    p <- ev[sort(parts[[k]]), , drop = FALSE]
    p <- cbind(
      synref_id = sprintf("%s_part%02d", source, k),
      draw = k,
      p,
      stratum = source
    )
    rownames(p) <- NULL
    p
  })
  prov <- tibble::as_tibble(do.call(rbind, drawn))
  prov <- prov[c(
    "synref_id",
    "draw",
    "source_sample",
    "source_row",
    "batch",
    "donor",
    "population",
    "stratum"
  )]
  events <- stats::setNames(
    lapply(drawn, \(p) .sw_synref_matrix(data, p)),
    vapply(drawn, \(p) p[["synref_id"]][1], "")
  )
  realised <- prov |>
    dplyr::summarise(
      realised = dplyr::n(),
      .by = c("synref_id", "draw", "stratum")
    ) |>
    dplyr::mutate(requested = as.integer(events_per_part), .after = "stratum")
  .sw_new_synref(
    events,
    prov,
    recipe,
    realised,
    scope = "within_batch",
    floor = "sampling only"
  )
}

# Mix real populations at set proportions ------------------------------------------

# Whole-number counts per population from proportions, largest remainder.
.sw_synref_counts <- function(proportions, n_events) {
  raw <- proportions * n_events
  counts <- floor(raw)
  short <- n_events - sum(counts)
  if (short > 0) {
    bump <- order(raw - counts, decreasing = TRUE)[seq_len(short)]
    counts[bump] <- counts[bump] + 1
  }
  counts
}

#' Mix real populations at set proportions
#'
#' Draws events without replacement from the labelled populations of one
#' within-batch group and returns a synthetic reference whose population
#' shares follow `proportions`. The requested and the realised counts and
#' proportions are reported for every population and draw. When a population
#' has fewer events than requested, `short_stratum` decides: `"stop"` stops
#' the call, `"cap"` takes every event of that population and reports the
#' shortfall in `realised`. The result is a diagnostic comparator with
#' certainty `"diagnostic_only"`. Cells are not replicates; donors and
#' aliquots are.
#'
#' @param data Named list of numeric event matrices, one per file, with the
#'   same column names. The caller transforms them first.
#' @param manifest Data frame with `sample_id` and `batch`.
#' @param groups Named list. Each element is a character vector of
#'   `sample_id`s that feed one reference. All files of a group must share one
#'   batch. Populations are pooled over the files of a group.
#' @param labels Named list of population labels, one vector per file, same
#'   length as its rows. It cannot be `NULL`: a mix needs labels.
#' @param proportions Named numeric vector of target shares (uncalibrated),
#'   names are population labels, values are positive and sum to 1.
#' @param n_events Events in each draw before any cap (uncalibrated). Counts
#'   per population follow the largest-remainder rule.
#' @param n_draws Number of independent draws per group (uncalibrated).
#' @param seed Random seed (uncalibrated). The caller's RNG state is restored.
#' @param short_stratum `"stop"` or `"cap"` (uncalibrated). What to do when a
#'   population has fewer events than its requested count.
#'
#' @return An `sw_synref` object. `realised` has `synref_id`, `draw`,
#'   `stratum` (the population), `requested`, `realised`, `requested_prop`
#'   and `realised_prop`.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(200), M2 = rnorm(200)))
#' lab <- list(a = rep(c("T", "B"), c(150, 50)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' ref <- sw_synref_mix(
#'   d, man, groups = list(r = "a"), labels = lab,
#'   proportions = c(T = 0.6, B = 0.4), n_events = 50, n_draws = 1,
#'   seed = 1, short_stratum = "stop"
#' )
#' ref$realised
#'
#' @export
sw_synref_mix <- function(
  data,
  manifest,
  groups,
  labels,
  proportions,
  n_events,
  n_draws,
  seed,
  short_stratum
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    groups = missing(groups),
    labels = missing(labels),
    proportions = missing(proportions),
    n_events = missing(n_events),
    n_draws = missing(n_draws),
    seed = missing(seed),
    short_stratum = missing(short_stratum)
  ))
  recipe <- list(
    kind = "mix",
    groups = groups,
    proportions = proportions,
    n_events = n_events,
    n_draws = n_draws,
    seed = seed,
    short_stratum = short_stratum,
    labels_given = !is.null(labels),
    markers = if (is.list(data) && length(data) > 0) {
      colnames(data[[1]])
    } else {
      NULL
    }
  )
  .sw_synref_check_data(data)
  .sw_synref_check_count(n_events, "n_events")
  .sw_synref_check_count(n_draws, "n_draws")
  .sw_synref_check_seed(seed)
  if (
    !is.character(short_stratum) ||
      length(short_stratum) != 1 ||
      !short_stratum %in% c("stop", "cap")
  ) {
    stop("'short_stratum' must be one of stop, cap.", call. = FALSE)
  }
  if (is.null(labels)) {
    stop("sw_synref_mix() needs 'labels'; NULL is not allowed.", call. = FALSE)
  }
  if (
    !is.numeric(proportions) ||
      is.null(names(proportions)) ||
      anyNA(names(proportions)) ||
      anyDuplicated(names(proportions)) ||
      anyNA(proportions) ||
      any(proportions <= 0) ||
      abs(sum(proportions) - 1) > 1e-8
  ) {
    stop(
      "'proportions' must be a uniquely named vector of positive shares that sum to 1.",
      call. = FALSE
    )
  }
  if (
    !is.list(groups) ||
      length(groups) == 0 ||
      is.null(names(groups)) ||
      anyNA(names(groups)) ||
      anyDuplicated(names(groups)) ||
      !all(vapply(groups, \(g) is.character(g) && length(g) > 0, logical(1)))
  ) {
    stop(
      "'groups' must be a uniquely named list of sample_id vectors.",
      call. = FALSE
    )
  }
  repeated <- names(groups)[vapply(groups, anyDuplicated, integer(1)) > 0]
  if (length(repeated) > 0) {
    stop(
      "Group(s) list a sample_id more than once: ",
      paste(repeated, collapse = ", "),
      ". A repeated file would let one event enter a reference twice.",
      call. = FALSE
    )
  }
  all_ids <- unique(unlist(groups, use.names = FALSE))
  absent <- setdiff(all_ids, names(data))
  if (length(absent) > 0) {
    stop(
      "Group sample(s) not in 'data': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  .sw_synref_check_manifest(manifest, all_ids)
  ok <- all(vapply(
    all_ids,
    \(id) length(labels[[id]]) == nrow(data[[id]]),
    logical(1)
  ))
  if (!ok) {
    stop(
      "'labels' needs one vector per file with one value per event.",
      call. = FALSE
    )
  }
  man_batch <- stats::setNames(
    as.character(manifest[["batch"]]),
    as.character(manifest[["sample_id"]])
  )
  for (g in names(groups)) {
    b <- unique(man_batch[groups[[g]]])
    if (length(b) > 1) {
      stop(
        sprintf(
          "Group '%s' spans batches %s. A synthetic reference is within-batch only.",
          g,
          paste(b, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  counts <- .sw_synref_counts(proportions, n_events)

  drawn <- .sw_with_seed(seed, {
    out <- list()
    for (g in names(groups)) {
      ev <- .sw_synref_events_table(data, manifest, groups[[g]], labels)
      if (anyNA(ev[["population"]])) {
        stop(
          sprintf("Group '%s' has events with a missing population label.", g),
          call. = FALSE
        )
      }
      ev[["stratum"]] <- ev[["population"]]
      unknown <- setdiff(names(counts), ev[["stratum"]])
      if (length(unknown) > 0 && short_stratum == "stop") {
        stop(
          sprintf(
            "Population(s) %s of 'proportions' not in the labels of group '%s'.",
            paste(unknown, collapse = ", "),
            g
          ),
          call. = FALSE
        )
      }
      if (length(unknown) > 0) {
        warning(
          sprintf(
            "Population(s) %s of 'proportions' have no events in group '%s'; 'cap' draws none.",
            paste(unknown, collapse = ", "),
            g
          ),
          call. = FALSE
        )
      }
      if (length(unknown) == length(counts)) {
        stop(
          sprintf("Group '%s' has none of the requested populations.", g),
          call. = FALSE
        )
      }
      for (k in seq_len(n_draws)) {
        key <- sprintf("%s_d%02d", g, k)
        picks <- lapply(names(counts), function(s) {
          idx <- which(ev[["stratum"]] == s)
          want <- counts[[s]]
          if (length(idx) < want) {
            if (short_stratum == "stop") {
              stop(
                sprintf(
                  "Population '%s' of group '%s' has %d events; %d requested.",
                  s,
                  g,
                  length(idx),
                  want
                ),
                call. = FALSE
              )
            }
            want <- length(idx)
          }
          if (want == 0) {
            return(integer(0))
          }
          idx[sort(sample.int(length(idx), want))]
        })
        prov <- ev[unlist(picks), , drop = FALSE]
        prov <- cbind(synref_id = key, draw = k, prov)
        rownames(prov) <- NULL
        out[[key]] <- prov
      }
    }
    out
  })

  prov <- tibble::as_tibble(do.call(rbind, unname(drawn)))
  prov <- prov[c(
    "synref_id",
    "draw",
    "source_sample",
    "source_row",
    "batch",
    "donor",
    "population",
    "stratum"
  )]
  events <- lapply(drawn, \(p) .sw_synref_matrix(data, p))
  # One row per draw and requested population, also when nothing was drawn.
  realised <- dplyr::cross_join(
    tibble::tibble(
      synref_id = names(drawn),
      draw = as.integer(sub("^.*_d", "", names(drawn)))
    ),
    tibble::tibble(
      stratum = names(counts),
      requested = as.integer(counts),
      requested_prop = as.numeric(proportions[names(counts)])
    )
  ) |>
    dplyr::left_join(
      prov |>
        dplyr::summarise(
          n_drawn = dplyr::n(),
          .by = c("synref_id", "stratum")
        ),
      by = c("synref_id", "stratum"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      realised = dplyr::coalesce(.data[["n_drawn"]], 0L),
      realised_prop = .data[["realised"]] / sum(.data[["realised"]]),
      .by = "synref_id"
    ) |>
    dplyr::select(
      "synref_id",
      "draw",
      "stratum",
      "requested",
      "realised",
      "requested_prop",
      "realised_prop"
    )
  .sw_new_synref(
    events,
    prov,
    recipe,
    realised,
    scope = "within_batch",
    floor = NA_character_
  )
}

# Quantile goal across batches -----------------------------------------------------

#' Quantile goal table across batches
#'
#' Draws `n_events` events from the pooled files of each batch, takes the
#' quantiles `probs` of every marker per batch, and sets the goal to the mean
#' or the median of the batch quantiles across batches. The result is a goal
#' for quantile alignment. It is never a reference: its role is `"goal"`, it
#' blends batches, and it cannot feed the paired shift or HEB bridge path.
#' The batch is the unit, not the cell. The table reports the spread of each
#' quantile across batches without a roll-up over markers.
#'
#' @param data Named list of numeric event matrices, one per file, with the
#'   same column names. The caller transforms them first.
#' @param manifest Data frame with `sample_id` and `batch`. At least two
#'   batches must appear.
#' @param probs Quantile probabilities in (0, 1) (uncalibrated).
#' @param goal_statistic `"mean"` or `"median"` (uncalibrated). How the batch
#'   quantiles combine into the goal. Sources disagree, so the choice is yours.
#' @param n_events Events drawn per batch from its pooled files
#'   (uncalibrated). The call stops when a batch has fewer.
#' @param seed Random seed (uncalibrated). The caller's RNG state is restored.
#'
#' @return An `sw_synref_goal` list: `goal` (marker, prob, goal, spread_sd,
#'   spread_min, spread_max, n_batches), `batch_quantiles` (batch, marker,
#'   prob, value, delta from goal, n_events), `recipe`, `role = "goal"`,
#'   `certainty = "diagnostic_only"` and `digest`.
#'
#' @examples
#' set.seed(1)
#' d <- list(
#'   a = cbind(M1 = rnorm(100)), b = cbind(M1 = rnorm(100, 1)),
#'   c = cbind(M1 = rnorm(100, 2))
#' )
#' man <- data.frame(sample_id = c("a", "b", "c"), batch = c("b1", "b2", "b3"))
#' sw_synref_goal(d, man, probs = c(0.1, 0.5, 0.9), goal_statistic = "median",
#'   n_events = 80, seed = 1)$goal
#'
#' @export
sw_synref_goal <- function(
  data,
  manifest,
  probs,
  goal_statistic,
  n_events,
  seed
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    probs = missing(probs),
    goal_statistic = missing(goal_statistic),
    n_events = missing(n_events),
    seed = missing(seed)
  ))
  .sw_synref_check_data(data)
  .sw_synref_check_count(n_events, "n_events")
  .sw_synref_check_seed(seed)
  if (
    !is.character(goal_statistic) ||
      length(goal_statistic) != 1 ||
      !goal_statistic %in% c("mean", "median")
  ) {
    stop("'goal_statistic' must be one of mean, median.", call. = FALSE)
  }
  if (
    !is.numeric(probs) ||
      length(probs) == 0 ||
      anyNA(probs) ||
      any(probs <= 0 | probs >= 1) ||
      anyDuplicated(probs)
  ) {
    stop(
      "'probs' must be unique numbers strictly between 0 and 1.",
      call. = FALSE
    )
  }
  .sw_synref_check_manifest(manifest, names(data))
  batch_of <- stats::setNames(
    as.character(manifest[["batch"]]),
    as.character(manifest[["sample_id"]])
  )[names(data)]
  if (anyNA(batch_of)) {
    stop("Every file in 'data' needs a batch.", call. = FALSE)
  }
  batches <- sort(unique(batch_of))
  if (length(batches) < 2) {
    stop("A goal needs at least two batches.", call. = FALSE)
  }
  markers <- colnames(data[[1]])
  fn <- if (goal_statistic == "mean") base::mean else stats::median

  per_batch <- .sw_with_seed(seed, {
    lapply(batches, function(b) {
      pooled <- do.call(rbind, unname(data[batch_of == b]))
      if (nrow(pooled) < n_events) {
        stop(
          sprintf(
            "Batch '%s' has %d events; %d requested.",
            b,
            nrow(pooled),
            n_events
          ),
          call. = FALSE
        )
      }
      pooled[sort(sample.int(nrow(pooled), n_events)), , drop = FALSE]
    })
  })
  names(per_batch) <- batches
  bq <- dplyr::bind_rows(lapply(batches, function(b) {
    q <- apply(per_batch[[b]], 2, stats::quantile, probs = probs, names = FALSE)
    q <- matrix(q, nrow = length(probs), dimnames = list(NULL, markers))
    tibble::tibble(
      batch = b,
      marker = rep(markers, each = length(probs)),
      prob = rep(probs, times = length(markers)),
      value = as.vector(q),
      n_events = as.integer(n_events)
    )
  }))
  goal <- bq |>
    dplyr::summarise(
      goal = fn(.data[["value"]]),
      spread_sd = stats::sd(.data[["value"]]),
      spread_min = min(.data[["value"]]),
      spread_max = max(.data[["value"]]),
      n_batches = dplyr::n(),
      .by = c("marker", "prob")
    )
  bq <- bq |>
    dplyr::left_join(
      goal[c("marker", "prob", "goal")],
      by = c("marker", "prob"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(delta = .data[["value"]] - .data[["goal"]]) |>
    dplyr::select("batch", "marker", "prob", "value", "delta", "n_events")
  obj <- structure(
    list(
      goal = goal,
      batch_quantiles = bq,
      recipe = list(
        probs = probs,
        goal_statistic = goal_statistic,
        n_events = n_events,
        seed = seed,
        batches = batches,
        markers = markers
      ),
      role = "goal",
      certainty = "diagnostic_only"
    ),
    class = "sw_synref_goal"
  )
  obj$digest <- .sw_model_digest(obj)
  obj
}

#' @export
print.sw_synref_goal <- function(x, ...) {
  cat(sprintf(
    "<sw_synref_goal> %d batches, %d markers, %d quantiles, goal = %s of batch quantiles\n",
    length(x[["recipe"]][["batches"]]),
    length(x[["recipe"]][["markers"]]),
    length(x[["recipe"]][["probs"]]),
    x[["recipe"]][["goal_statistic"]]
  ))
  cat("role: goal. Not a reference; it blends batches.\n")
  invisible(x)
}
