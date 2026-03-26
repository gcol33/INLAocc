# =============================================================================
# occu_fit.R — Internal helpers and deprecated legacy wrappers
# =============================================================================


# ---------------------------------------------------------------------------
# Internal helpers (used by engines)
# ---------------------------------------------------------------------------

#' @noRd
merge_re <- function(formula_re, explicit_re) {
  if (inherits(explicit_re, "occu_re")) explicit_re <- list(explicit_re)
  combined <- c(formula_re, explicit_re)
  if (length(combined) == 0) return(NULL)
  combined
}

#' @noRd
coerce_data <- function(data) {
  if (inherits(data, "occu_data")) return(data)
  occu_format(
    y        = data$y,
    occ.covs = data$occ.covs,
    det.covs = data$det.covs,
    coords   = data$coords
  )
}

#' @noRd
pool_community_effects <- function(beta_list) {
  beta_list <- Filter(function(x) !is.null(x), beta_list)
  if (length(beta_list) == 0) return(NULL)

  coef_names <- rownames(beta_list[[1]])
  community <- data.frame(
    parameter      = coef_names,
    community_mean = NA_real_,
    community_sd   = NA_real_,
    species_sd     = NA_real_,
    n_species      = length(beta_list)
  )

  for (k in seq_along(coef_names)) {
    sp_means <- vapply(beta_list, function(b) b$mean[k], numeric(1))
    sp_sds   <- vapply(beta_list, function(b) b$sd[k], numeric(1))
    community$community_mean[k] <- mean(sp_means)
    community$community_sd[k]   <- sqrt(mean(sp_sds^2) + var(sp_means))
    community$species_sd[k]     <- sd(sp_means)
  }
  community
}


#' Bootstrap-ensemble community pooling for robustness against outlier species
#'
#' Fits community estimates on B random subsets (each of size prop * n_species),
#' then pools across subsets. Species that distort the community mean in the
#' full-data estimate get diluted because they appear in only some subsets.
#'
#' @param beta_list Named list of summary.fixed data frames (one per species)
#' @param B Number of bootstrap subsets (default 50)
#' @param prop Proportion of species per subset (default 0.7)
#' @param seed Random seed for reproducibility
#' @return Community data frame (same format as pool_community_effects)
#' @noRd
pool_community_ensemble <- function(beta_list, B = 50L, prop = 0.7,
                                     seed = NULL) {
  beta_list <- Filter(function(x) !is.null(x), beta_list)
  n_sp <- length(beta_list)
  if (n_sp < 3) return(pool_community_effects(beta_list))

  subset_size <- max(3L, round(n_sp * prop))
  coef_names <- rownames(beta_list[[1]])
  n_coef <- length(coef_names)

  if (!is.null(seed)) set.seed(seed)

  # Each subset produces a community estimate
  ensemble_means <- matrix(NA, n_coef, B)
  ensemble_sds   <- matrix(NA, n_coef, B)

  for (b in seq_len(B)) {
    idx <- sample(n_sp, subset_size, replace = FALSE)
    sub_pool <- pool_community_effects(beta_list[idx])
    ensemble_means[, b] <- sub_pool$community_mean
    ensemble_sds[, b]   <- sub_pool$community_sd
  }

  # Pool across subsets: mean of means, inflated SD
  data.frame(
    parameter      = coef_names,
    community_mean = rowMeans(ensemble_means),
    community_sd   = sqrt(rowMeans(ensemble_sds^2) +
                            apply(ensemble_means, 1, var)),
    species_sd     = apply(ensemble_means, 1, sd),
    n_species      = n_sp,
    n_subsets       = B
  )
}


# =============================================================================
# Deprecated wrappers — use occu() instead
# =============================================================================

#' @rdname occu
#' @usage NULL
#' @export
occu_inla <- function(occ.formula, det.formula, data, priors = NULL,
                      occ.re = NULL, det.re = NULL, k.fold = 0,
                      max.iter = 50, tol = 1e-4, damping = 0.3,
                      verbose = 1) {
  .Deprecated("occu")
  occu(occ.formula, det.formula, data, priors = priors,
       occ.re = occ.re, det.re = det.re, k.fold = k.fold,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
spatial_occu_inla <- function(occ.formula, det.formula, data,
                              coords = NULL, priors = NULL,
                              occ.re = NULL, det.re = NULL,
                              spde.args = list(), k.fold = 0,
                              max.iter = 50, tol = 1e-4, damping = 0.3,
                              verbose = 1) {
  .Deprecated("occu")
  spatial <- coords %||% data$coords
  occu(occ.formula, det.formula, data, spatial = spatial,
       spde.args = spde.args, priors = priors,
       occ.re = occ.re, det.re = det.re, k.fold = k.fold,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
intOccu_inla <- function(occ.formula, det.formula, data,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         max.iter = 50, tol = 1e-4, damping = 0.3,
                         verbose = 1) {
  .Deprecated("occu")
  occu(occ.formula, det.formula, data, integrated = TRUE,
       priors = priors, occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
spIntOccu_inla <- function(occ.formula, det.formula, data,
                           coords = NULL, priors = NULL,
                           occ.re = NULL, det.re = NULL,
                           spde.args = list(),
                           max.iter = 50, tol = 1e-4, damping = 0.3,
                           verbose = 1) {
  .Deprecated("occu")
  spatial <- coords %||% data$coords
  occu(occ.formula, det.formula, data, spatial = spatial,
       integrated = TRUE, spde.args = spde.args, priors = priors,
       occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
ms_occu_inla <- function(occ.formula, det.formula, data,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         k.fold = 0,
                         max.iter = 50, tol = 1e-4, damping = 0.3,
                         verbose = 1) {
  .Deprecated("occu")
  occu(occ.formula, det.formula, data, multispecies = TRUE,
       priors = priors, occ.re = occ.re, det.re = det.re, k.fold = k.fold,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
lfMsOccu_inla <- function(occ.formula, det.formula, data,
                          n.factors = NULL,
                          priors = NULL, occ.re = NULL, det.re = NULL,
                          max.iter = 30, tol = 1e-3, damping = 0.3,
                          verbose = 1) {
  .Deprecated("occu")
  occu(occ.formula, det.formula, data, multispecies = TRUE,
       n.factors = n.factors, priors = priors,
       occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
sfMsOccu_inla <- function(occ.formula, det.formula, data,
                          coords = NULL, n.factors = NULL,
                          priors = NULL, occ.re = NULL, det.re = NULL,
                          spde.args = list(),
                          max.iter = 30, tol = 1e-3, damping = 0.3,
                          verbose = 1) {
  .Deprecated("occu")
  spatial <- coords %||% data$coords
  occu(occ.formula, det.formula, data, spatial = spatial,
       multispecies = TRUE, n.factors = n.factors,
       spde.args = spde.args, priors = priors,
       occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
svcOccu_inla <- function(occ.formula, det.formula, data,
                         coords = NULL, svc.cols = 1L,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         spde.args = list(),
                         max.iter = 50, tol = 1e-4, damping = 0.3,
                         verbose = 1) {
  .Deprecated("occu")
  spatial <- coords %||% data$coords
  occu(occ.formula, det.formula, data, spatial = spatial,
       svc = svc.cols, spde.args = spde.args, priors = priors,
       occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
temporal_occu_inla <- function(occ.formula, det.formula, data,
                               ar1 = TRUE, priors = NULL,
                               occ.re = NULL, det.re = NULL,
                               coords = NULL, spde.args = list(),
                               max.iter = 30, tol = 1e-3, damping = 0.3,
                               verbose = 1) {
  .Deprecated("occu")
  spatial <- coords %||% data$coords
  temporal <- if (ar1) "ar1" else "iid"
  occu(occ.formula, det.formula, data, spatial = spatial,
       temporal = temporal, spde.args = spde.args, priors = priors,
       occ.re = occ.re, det.re = det.re,
       max.iter = max.iter, tol = tol, damping = damping, verbose = verbose)
}

#' @rdname occu
#' @usage NULL
#' @export
lfJSDM_inla <- function(formula, data, n.factors = NULL,
                        priors = NULL, verbose = 1) {
  .Deprecated("occu")
  occu(formula, data = data, multispecies = "jsdm",
       n.factors = n.factors, priors = priors, verbose = verbose)
}
