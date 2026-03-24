# =============================================================================
# zzz.R — Package load / attach hooks
# =============================================================================

.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage(
      "INLAocc requires the INLA package. Install with:\n",
      '  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)'
    )
  }
}

#' @noRd
check_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "INLA is required but not installed. Install with:\n",
      '  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)',
      call. = FALSE
    )
  }
}

#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
