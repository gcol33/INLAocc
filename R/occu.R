# =============================================================================
# occu.R — Unified entry point for all occupancy models
# =============================================================================

#' Fit occupancy models using INLA
#'
#' Single entry point for all occupancy model types. The model structure is
#' controlled by flags: \code{spatial}, \code{temporal}, \code{multispecies},
#' \code{integrated}, \code{n.factors}, and \code{svc}. Covers all model types
#' from \pkg{spOccupancy} (PGOcc, spPGOcc, tPGOcc, msPGOcc, etc.) via
#' flag combinations.
#'
#' @param occ.formula RHS formula for occupancy (psi). Supports mixed-model
#'   random effects syntax: \code{~ elev + forest + (1 | region)}.
#' @param det.formula RHS formula for detection (p). NULL for JSDM models
#'   (set \code{multispecies = "jsdm"}).
#' @param data An \code{occu_data} or \code{occu_data_ms} object, or a raw
#'   list with components \code{y}, \code{occ.covs}, \code{det.covs}
#'   (spOccupancy-compatible).
#' @param spatial N x 2 coordinate matrix or \code{occu_spatial()} object.
#'   Enables spatial SPDE models.
#' @param temporal \code{"ar1"} or \code{"iid"} for multi-season models.
#'   NULL (default) for single-season.
#' @param multispecies \code{FALSE} (default), \code{TRUE} for community
#'   models, or \code{"jsdm"} for joint species distribution models
#'   (no detection process).
#' @param integrated \code{TRUE} for multi-data-source models. Requires
#'   \code{data$y} to be a list of detection matrices.
#' @param n.factors Integer number of latent factors. Enables latent factor
#'   variants of multi-species models.
#' @param svc Integer vector of occupancy design matrix columns that get
#'   spatially-varying coefficients (1 = intercept). Requires \code{spatial}.
#' @param spde.args List of arguments passed to \code{occu_spatial()} when
#'   \code{spatial} is a coordinate matrix.
#' @param priors \code{occu_priors()} object or spOccupancy-compatible named
#'   list.
#' @param occ.re Explicit list of \code{occu_re()} specs for occupancy.
#' @param det.re Explicit list of \code{occu_re()} specs for detection.
#' @param max.iter Maximum EM iterations (default 50).
#' @param tol Convergence tolerance (default 1e-4).
#' @param damping EM damping factor 0--1 (default 0.3).
#' @param k.fold Number of cross-validation folds (default 0, no CV).
#' @param ensemble Logical; if `TRUE`, use ensemble averaging across
#'   multiple imputation chains (default `FALSE`).
#' @param num.threads Thread specification for INLA, as `"A:B"` where A is
#'   outer and B is inner OpenMP threads (default `"1:1"`).
#' @param verbose 0 = silent, 1 = iteration summaries, 2 = full INLA output.
#'
#' @return An S3 object whose class depends on the model type (e.g.,
#'   \code{"occu_inla"}, \code{"occu_inla_spatial"}, \code{"occu_inla_ms"}).
#'
#' @examples
#' \donttest{
#' # Single-species (cf. PGOcc)
#' # fit <- occu(~ elev, ~ effort, data)
#'
#' # Spatial (cf. spPGOcc)
#' # fit <- occu(~ elev, ~ effort, data, spatial = coords)
#'
#' # Multi-species (cf. msPGOcc)
#' # fit <- occu(~ elev, ~ effort, ms_data, multispecies = TRUE)
#'
#' # Temporal (cf. tPGOcc)
#' # fit <- occu(~ elev, ~ effort, data, temporal = "ar1")
#' }
#'
#' @export
occu <- function(occ.formula, det.formula = NULL, data,
                 spatial      = NULL,
                 temporal     = NULL,
                 multispecies = FALSE,
                 integrated   = FALSE,
                 n.factors    = NULL,
                 svc          = NULL,
                 spde.args    = list(),
                 priors       = NULL,
                 occ.re       = NULL,
                 det.re       = NULL,
                 max.iter     = 50L,
                 tol          = 1e-4,
                 damping      = 0.3,
                 correction   = "auto",
                 ensemble     = FALSE,
                 num.threads  = "1:1",
                 k.fold       = 0L,
                 verbose      = 1L) {

  check_inla()

  # --- Default detection formula ---
  if (is.null(det.formula) && !identical(multispecies, "jsdm")) {
    warning("No detection formula supplied; defaulting to ~ 1 (intercept-only detection)")
    det.formula <- ~ 1
  }

  # --- 1. Resolve model type ---
  model_key <- resolve_model_type(
    spatial, temporal, multispecies, integrated, n.factors, svc, det.formula
  )
  spec <- MODEL_REGISTRY[[model_key]]
  if (is.null(spec))
    stop(sprintf("Unsupported model combination (key: %s). Please file an issue.", model_key))

  # --- 2. Shared prep: spatial ---
  spatial_obj <- prepare_spatial(spatial, spde.args)

  # --- 3. Shared prep: data coercion ---
  data <- prepare_data(data, multispecies, integrated, temporal)

  # --- 4. Shared prep: formulas ---
  is_jsdm <- identical(multispecies, "jsdm")
  occ_parsed <- parse_re_formula(occ.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)

  det_parsed <- NULL
  det_re_all <- NULL
  if (!is_jsdm && !is.null(det.formula)) {
    det_parsed <- parse_re_formula(det.formula)
    det_re_all <- merge_re(det_parsed$re_list, det.re)
  }

  # --- 5. Shared prep: priors ---
  if (!is.null(priors) && !inherits(priors, "occu_priors")) {
    priors <- do.call(occu_priors, priors)
  }

  # --- 6. Validate correction ---
  valid_corrections <- c("auto", "mi", "gibbs", "none")
  if (!correction %in% valid_corrections) {
    stop(sprintf("correction must be one of [%s], got '%s'",
                 paste(valid_corrections, collapse = ", "), correction))
  }

  # --- 7. Build common args ---
  args <- list(
    data        = data,
    occ_parsed  = occ_parsed,
    det_parsed  = det_parsed,
    occ_re      = occ_re_all,
    det_re      = det_re_all,
    spatial     = spatial_obj,
    temporal    = temporal,
    n.factors   = n.factors,
    svc         = svc,
    priors      = priors,
    max.iter    = max.iter,
    tol         = tol,
    damping     = damping,
    correction  = correction,
    ensemble    = ensemble,
    num.threads = num.threads,
    k.fold      = k.fold,
    verbose     = verbose,
    # Keep raw formulas for storage in result
    occ.formula = occ.formula,
    det.formula = det.formula
  )

  # --- 7. Dispatch to engine ---
  t0 <- proc.time()
  engine_fn <- get(spec$engine, envir = asNamespace("INLAocc"), mode = "function")
  result <- engine_fn(args)

  # --- 8. Store metadata and assign class ---
  result$call       <- match.call()
  result$model_type <- model_key
  result$run_time   <- proc.time() - t0
  result$occ_fixed_formula <- occ_parsed$fixed
  result$det_fixed_formula <- if (!is.null(det_parsed)) det_parsed$fixed else NULL
  class(result) <- spec$class
  result
}
