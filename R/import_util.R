# R/import_util.R
# Shared helper for the import-step functions.

# flowCore::read.FCSheader() can leave a connection open when the file is not a
# valid FCS file, and R then warns "closing unused connection" at garbage
# collection. Close any connection the read opened, whether it succeeded or not.
.sw_read_fcs_header <- function(path, ...) {
  conns_before <- getAllConnections()
  on.exit(
    for (id in setdiff(getAllConnections(), conns_before)) {
      try(close(getConnection(id)), silent = TRUE)
    },
    add = TRUE
  )
  flowCore::read.FCSheader(path, ...)
}

# Shared provenance block for the import step. Every import result carries it
# as attribute `provenance`, with the same leading keys in the same order:
# package_version, flowcore_version, r_version, platform, datetime, seed,
# session, args. Function-specific keys follow in `extra`. `seed` is stored as
# given (NULL or NA when nothing is random). Versions come from .sw_pkg_info()
# (R/unmix_report.R).
.sw_import_provenance <- function(args, seed = NULL, extra = list()) {
  si <- utils::sessionInfo()
  c(
    list(
      package_version = .sw_pkg_info("SpectraWeaveR")$version,
      flowcore_version = .sw_pkg_info("flowCore")$version,
      r_version = R.version.string,
      platform = R.version$platform,
      datetime = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      seed = seed,
      session = list(
        running = if (is.null(si$running)) NA_character_ else si$running,
        locale = si$locale,
        attached = unname(vapply(
          si$otherPkgs,
          \(p) paste0(p$Package, " ", p$Version),
          character(1)
        ))
      ),
      args = args
    ),
    extra
  )
}

# Package version table kept as an extra key by results that listed it before.
.sw_import_packages <- function() {
  dplyr::bind_rows(lapply(c("SpectraWeaveR", "flowCore"), .sw_pkg_info))
}
