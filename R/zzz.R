# =============================================================================
# zzz.R — Package load / attach hooks
# =============================================================================

#' @importFrom stats as.formula binomial cor model.matrix pnorm prcomp
#'   quantile rbinom rnorm runif sd terms var

.onLoad <- function(libname, pkgname) {
  # Soft-check for INLA at load time (no error, just a message)
  if (!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage(
      "INLAocc requires the INLA package. Install with:\n",
      '  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)'
    )
  }
}

#' Check that INLA is available (internal)
#' @keywords internal
check_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "INLA is required but not installed. Install with:\n",
      '  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)',
      call. = FALSE
    )
  }
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
