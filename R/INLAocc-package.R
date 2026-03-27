#' @keywords internal
"_PACKAGE"

#' @importFrom stats as.formula binom.test binomial coef confint cor dist
#'   fitted glm ks.test lm logLik lowess median model.matrix nobs pnorm
#'   prcomp predict qnorm quantile rbeta rbinom residuals rnorm runif sd
#'   setNames simulate terms update update.formula var vcov
#' @importFrom Matrix Matrix
#' @importFrom graphics abline arrows axis barplot hist image legend lines
#'   mtext par plot points polygon rect segments text
#' @importFrom grDevices adjustcolor colorRampPalette hcl.colors rgb
#' @importFrom utils tail
NULL

#' Tidy model output into a data.frame
#'
#' Generic for producing tidy coefficient tables. See
#' \code{\link{tidy.occu_inla}} for the occupancy model method.
#'
#' @param x a model object
#' @param ... additional arguments passed to methods
#' @return A data.frame.
#' @export
tidy <- function(x, ...) UseMethod("tidy")

#' Extract random effects
#'
#' Generic for extracting random effect estimates. See
#' \code{\link{ranef.occu_inla}} for the occupancy model method.
#'
#' @param object a model object
#' @param ... additional arguments passed to methods
#' @return A list of data.frames.
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' Glance at model-level statistics
#'
#' Generic for producing single-row model summaries. See
#' \code{\link{glance.occu_inla}} for the occupancy model method.
#'
#' @param x a model object
#' @param ... additional arguments passed to methods
#' @return A one-row data.frame.
#' @export
glance <- function(x, ...) UseMethod("glance")
