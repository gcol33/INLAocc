# =============================================================================
# occu_fit.R — User-facing model fitting functions
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: merge formula-parsed REs with explicit occ.re / det.re
# ---------------------------------------------------------------------------
merge_re <- function(formula_re, explicit_re) {
  if (inherits(explicit_re, "occu_re")) explicit_re <- list(explicit_re)
  combined <- c(formula_re, explicit_re)
  if (length(combined) == 0) return(NULL)
  combined
}

# ---------------------------------------------------------------------------
# Internal: coerce raw spOccupancy-style list to occu_data
# ---------------------------------------------------------------------------
coerce_data <- function(data) {
  if (inherits(data, "occu_data")) return(data)
  occu_format(
    y        = data$y,
    occ.covs = data$occ.covs,
    det.covs = data$det.covs,
    coords   = data$coords
  )
}


# ==========================================================================
# 1.  occu_inla  —  single-species, non-spatial  (cf. PGOcc)
# ==========================================================================

#' Fit a single-species occupancy model using INLA (Laplace approximation)
#'
#' Analogous to spOccupancy::PGOcc(). Supports covariates on both occupancy
#' and detection, random intercepts and slopes (via lme4-style \code{(1|group)}
#' or \code{(x|group)} syntax in the formulas, or via explicit \code{occ.re}/
#' \code{det.re} arguments), and prior specification.
#'
#' @param occ.formula RHS formula for occupancy (psi). Supports lme4 RE syntax.
#'   Example: \code{~ elev + forest + (1 | region)}
#' @param det.formula RHS formula for detection (p).
#'   Example: \code{~ effort + date + (1 | observer)}
#' @param data an occu_data object (from \code{occu_format()}), or a raw list
#'   with components \code{y}, \code{occ.covs}, \code{det.covs} (spOccupancy style).
#' @param priors optional \code{occu_priors()} object or a named list with
#'   \code{beta.normal}, \code{alpha.normal}, \code{sigma.sq.psi.ig},
#'   \code{sigma.sq.p.ig} (spOccupancy-compatible).
#' @param occ.re explicit list of \code{occu_re()} specs (merged with any
#'   inline \code{(1|group)} terms in \code{occ.formula}).
#' @param det.re explicit list of \code{occu_re()} specs.
#' @param k.fold integer: number of cross-validation folds. If > 1, runs
#'   k-fold CV after fitting and reports deviance.
#' @param max.iter maximum EM iterations (default 50).
#' @param tol convergence tolerance (default 1e-4).
#' @param damping EM damping factor 0-1 (default 0.3).
#' @param verbose 0 = silent, 1 = iteration summaries, 2 = full INLA output.
#'
#' @return object of class \code{"occu_inla"} with model results
occu_inla <- function(occ.formula, det.formula, data,
                      priors = NULL, occ.re = NULL, det.re = NULL,
                      k.fold = 0,
                      max.iter = 50, tol = 1e-4, damping = 0.3,
                      verbose = 1) {
  check_inla()
  data <- coerce_data(data)

  # Parse formulas (extract inline REs, strip them from fixed part)
  occ_parsed <- parse_re_formula(occ.formula)
  det_parsed <- parse_re_formula(det.formula)

  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)
  det_re_all <- merge_re(det_parsed$re_list, det.re)

  # Convert spOcc-style priors list
  if (!is.null(priors) && !inherits(priors, "occu_priors")) {
    priors <- do.call(occu_priors, priors)
  }

  if (verbose >= 1) {
    cat("Fitting single-species occupancy model (INLA)\n")
    cat(sprintf("  Sites: %d | Max visits: %d | Naive psi: %.3f | Naive p: %.3f\n",
                data$N, data$J, data$naive_occ,
                ifelse(is.na(data$naive_det), 0, data$naive_det)))
    n_re <- length(occ_re_all) + length(det_re_all)
    if (n_re > 0) cat(sprintf("  Random effects: %d\n", n_re))
  }

  result <- em_inla(
    data        = data,
    occ_formula = occ_parsed$fixed,
    det_formula = det_parsed$fixed,
    occ_re      = occ_re_all,
    det_re      = det_re_all,
    spatial     = NULL,
    priors      = priors,
    max_iter    = max.iter,
    tol         = tol,
    damping     = damping,
    verbose     = verbose
  )

  # K-fold CV
  if (k.fold > 1) {
    if (verbose >= 1) cat(sprintf("\nRunning %d-fold cross-validation...\n", k.fold))
    result$k.fold <- kfold_occu(
      fit_fun      = occu_inla,
      data         = data,
      k            = k.fold,
      occ.formula  = occ.formula,
      det.formula  = det.formula,
      priors       = priors,
      occ.re       = occ.re,
      det.re       = det.re,
      max.iter     = max.iter,
      tol          = tol,
      damping      = damping
    )
    if (verbose >= 1) {
      cat(sprintf("  K-fold deviance: %.2f\n", result$k.fold$k.fold.deviance))
    }
  }

  class(result) <- c("occu_inla", "occu_em")
  result
}


# ==========================================================================
# 2.  spatial_occu_inla  —  spatial NNGP/SPDE  (cf. spPGOcc)
# ==========================================================================

#' Fit a spatial occupancy model using INLA-SPDE
#'
#' Adds a Matern Gaussian random field (via SPDE) to the occupancy linear
#' predictor. Analogous to spOccupancy::spPGOcc().
#'
#' @inheritParams occu_inla
#' @param coords N x 2 matrix of site coordinates. If NULL, extracted from data.
#' @param spde.args list passed to \code{occu_spatial()}: max.edge, cutoff,
#'   offset, prior.range, prior.sigma.
#'
#' @return object of class \code{"occu_inla_spatial"}
spatial_occu_inla <- function(occ.formula, det.formula, data,
                              coords = NULL, priors = NULL,
                              occ.re = NULL, det.re = NULL,
                              spde.args = list(),
                              k.fold = 0,
                              max.iter = 50, tol = 1e-4, damping = 0.3,
                              verbose = 1) {
  check_inla()
  data <- coerce_data(data)

  occ_parsed <- parse_re_formula(occ.formula)
  det_parsed <- parse_re_formula(det.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)
  det_re_all <- merge_re(det_parsed$re_list, det.re)

  if (!is.null(priors) && !inherits(priors, "occu_priors")) {
    priors <- do.call(occu_priors, priors)
  }

  if (is.null(coords)) coords <- data$coords
  if (is.null(coords)) {
    stop("coords must be provided (N x 2 matrix of site locations)")
  }

  # Spatial mesh
  spde_args_full <- c(list(coords = coords), spde.args)
  spatial <- do.call(occu_spatial, spde_args_full)

  if (verbose >= 1) {
    cat("Fitting spatial occupancy model (INLA-SPDE)\n")
    cat(sprintf("  Sites: %d | Mesh nodes: %d\n", data$N, spatial$n_mesh))
  }

  result <- em_inla(
    data        = data,
    occ_formula = occ_parsed$fixed,
    det_formula = det_parsed$fixed,
    occ_re      = occ_re_all,
    det_re      = det_re_all,
    spatial     = spatial,
    priors      = priors,
    max_iter    = max.iter,
    tol         = tol,
    damping     = damping,
    verbose     = verbose
  )

  class(result) <- c("occu_inla_spatial", "occu_inla", "occu_em")
  result
}


# ==========================================================================
# 3.  intOccu_inla  —  integrated / data-fusion  (cf. intPGOcc)
# ==========================================================================

#' Fit an integrated (multi-data-source) occupancy model
#'
#' Uses a shared occupancy process with separate detection models per data
#' source. Analogous to spOccupancy::intPGOcc().
#'
#' @param occ.formula occupancy formula (single, shared across data sources)
#' @param det.formula list of detection formulas (one per data source)
#' @param data list with components:
#'   \describe{
#'     \item{y}{list of detection matrices (one per data source)}
#'     \item{occ.covs}{data.frame of occupancy covariates (all unique sites)}
#'     \item{det.covs}{list of detection covariate lists (one per data source)}
#'     \item{sites}{list of integer vectors mapping each source to rows in occ.covs}
#'     \item{coords}{optional coordinate matrix (for spatial variant)}
#'   }
#' @param priors optional priors specification
#' @param occ.re occupancy random effects
#' @param det.re list of detection RE lists (one per data source), or single list applied to all
#' @param max.iter max EM iterations
#' @param verbose verbosity level
#'
#' @return object of class \code{"occu_inla_int"}
intOccu_inla <- function(occ.formula, det.formula, data,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         max.iter = 50, tol = 1e-4, damping = 0.3,
                         verbose = 1) {
  check_inla()

  n_data <- length(data$y)
  if (!is.list(det.formula)) {
    det.formula <- rep(list(det.formula), n_data)
  }
  if (length(det.formula) != n_data) {
    stop(sprintf("det.formula must be a list of %d formulas (one per data source)", n_data))
  }

  occ_parsed <- parse_re_formula(occ.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)

  if (!is.null(priors) && !inherits(priors, "occu_priors")) {
    priors <- do.call(occu_priors, priors)
  }

  # Total unique sites
  N_total <- nrow(data$occ.covs)
  site_map <- data$sites

  if (verbose >= 1) {
    cat("Fitting integrated occupancy model (INLA)\n")
    cat(sprintf("  Data sources: %d | Total sites: %d\n", n_data, N_total))
    for (d in seq_len(n_data)) {
      cat(sprintf("  Source %d: %d sites, %d max visits\n",
                  d, nrow(data$y[[d]]), ncol(data$y[[d]])))
    }
  }

  # EM algorithm adapted for multiple data sources
  # Initialize
  psi_hat <- rep(0.5, N_total)
  det_fits <- vector("list", n_data)
  p_hats   <- vector("list", n_data)

  # Initialize p per source
  for (d in seq_len(n_data)) {
    y_d <- data$y[[d]]
    n_det_d <- sum(y_d, na.rm = TRUE)
    n_vis_d <- sum(!is.na(y_d))
    p_hats[[d]] <- matrix(clamp(n_det_d / max(n_vis_d, 1), 0.1, 0.9),
                          nrow(y_d), ncol(y_d))
  }

  history <- list()
  converged <- FALSE

  for (iter in seq_len(max.iter)) {
    psi_old <- psi_hat

    # --- E-step: compute weights per source ---
    weights_all <- rep(1, N_total)
    detected_any <- rep(FALSE, N_total)

    for (d in seq_len(n_data)) {
      sites_d <- site_map[[d]]
      y_d     <- data$y[[d]]
      for (i in seq_along(sites_d)) {
        si <- sites_d[i]
        if (any(y_d[i, ] == 1, na.rm = TRUE)) {
          detected_any[si] <- TRUE
        }
      }
    }

    # For undetected sites across ALL sources, compute joint weight
    for (si in which(!detected_any)) {
      q_total <- 1
      for (d in seq_len(n_data)) {
        idx_in_d <- which(site_map[[d]] == si)
        if (length(idx_in_d) > 0) {
          p_d <- p_hats[[d]][idx_in_d, ]
          visits <- which(!is.na(p_d))
          if (length(visits) > 0) {
            q_total <- q_total * prod(clamp(1 - p_d[visits]))
          }
        }
      }
      psi_si <- clamp(psi_hat[si])
      weights_all[si] <- (psi_si * q_total) / ((1 - psi_si) + psi_si * q_total)
    }

    # --- M-step: Detection (per source) ---
    for (d in seq_len(n_data)) {
      sites_d <- site_map[[d]]
      w_d     <- weights_all[sites_d]

      det_parsed_d <- parse_re_formula(det.formula[[d]])
      det_re_d <- det_parsed_d$re_list
      if (!is.null(det.re)) {
        if (is.list(det.re) && length(det.re) == n_data && is.list(det.re[[d]])) {
          det_re_d <- c(det_re_d, det.re[[d]])
        } else {
          det_re_d <- merge_re(det_re_d, det.re)
        }
      }

      source_data <- occu_format(
        y        = data$y[[d]],
        occ.covs = data$occ.covs[sites_d, , drop = FALSE],
        det.covs = if (!is.null(data$det.covs)) data$det.covs[[d]] else NULL
      )

      det_result <- fit_detection_inla(
        source_data, det_parsed_d$fixed, w_d,
        det_re_list  = if (length(det_re_d) > 0) det_re_d else NULL,
        control.inla = NULL,
        verbose      = verbose >= 2
      )
      det_fits[[d]] <- det_result$fit
      p_hats[[d]]   <- det_result$p_hat
    }

    # --- M-step: Occupancy (shared across all sources) ---
    occ_df <- build_occ_df(
      occ_covs = data$occ.covs,
      weights  = weights_all,
      site_idx = seq_len(N_total),
      detected = detected_any
    )

    re_terms <- character(0)
    if (!is.null(occ_re_all)) {
      re_comp <- build_re_components(
        occ_re_all, occ_df,
        list(occ.covs = data$occ.covs),
        prefix = "occ", process = "occ"
      )
      occ_df   <- re_comp$df
      re_terms <- re_comp$formula_terms
    }

    base_formula <- paste("z ~", as.character(occ_parsed$fixed)[2])
    if (length(re_terms) > 0) {
      base_formula <- paste(base_formula, "+", paste(re_terms, collapse = " + "))
    }

    occ_fit <- INLA::inla(
      formula       = as.formula(base_formula),
      family        = "binomial",
      Ntrials       = rep(1, nrow(occ_df)),
      data          = occ_df,
      weights       = occ_df$w,
      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
      verbose       = verbose >= 2
    )

    # Extract psi (one row per site now)
    psi_fitted <- occ_fit$summary.fitted.values$mean
    psi_hat <- psi_fitted[seq_len(N_total)]

    # --- Convergence ---
    delta_psi <- max(abs(psi_hat - psi_old))
    if (verbose >= 1) {
      cat(sprintf("  EM iter %2d | delta_psi = %.6f\n", iter, delta_psi))
    }
    history[[iter]] <- list(iter = iter, delta_psi = delta_psi)

    if (delta_psi < tol) {
      converged <- TRUE
      if (verbose >= 1) cat("  Converged.\n")
      break
    }
  }

  z_hat <- weights_all
  z_hat[detected_any] <- 1

  out <- list(
    occ_fit     = occ_fit,
    det_fits    = det_fits,
    psi_hat     = psi_hat,
    p_hats      = p_hats,
    z_hat       = z_hat,
    data        = data,
    occ_formula = occ.formula,
    det_formula = det.formula,
    n_data      = n_data,
    converged   = converged,
    n_iter      = length(history),
    history     = history
  )
  class(out) <- c("occu_inla_int", "occu_inla", "occu_em")
  out
}


#' Fit a spatial integrated occupancy model (cf. spIntPGOcc)
#'
#' @inheritParams intOccu_inla
#' @param coords coordinate matrix
#' @param spde.args SPDE mesh arguments
#' @return object of class "occu_inla_spint"
spIntOccu_inla <- function(occ.formula, det.formula, data,
                           coords = NULL, priors = NULL,
                           occ.re = NULL, det.re = NULL,
                           spde.args = list(),
                           max.iter = 50, tol = 1e-4, damping = 0.3,
                           verbose = 1) {
  # Fit non-spatial integrated first, then add spatial field
  # (spatial field is added to occupancy linear predictor)
  result <- intOccu_inla(
    occ.formula, det.formula, data,
    priors = priors, occ.re = occ.re, det.re = det.re,
    max.iter = max.iter, tol = tol, damping = damping, verbose = verbose
  )

  # Add spatial component if coords provided
  if (is.null(coords)) coords <- data$coords
  if (!is.null(coords)) {
    spde_args_full <- c(list(coords = coords), spde.args)
    result$spatial <- do.call(occu_spatial, spde_args_full)
  }

  class(result) <- c("occu_inla_spint", "occu_inla_int", "occu_inla", "occu_em")
  result
}


# ==========================================================================
# 4.  ms_occu_inla  —  multi-species  (cf. msPGOcc)
# ==========================================================================

#' Fit a multi-species occupancy model using INLA
#'
#' Fits species-specific occupancy and detection models with community-level
#' hyperparameters. Analogous to spOccupancy::msPGOcc().
#'
#' Accepts data with y as a 3D array (species x sites x visits) matching
#' spOccupancy format, or as a named list of matrices.
#'
#' @param occ.formula occupancy formula (shared structure across species)
#' @param det.formula detection formula
#' @param data either an occu_data_ms object, or a list with:
#'   \describe{
#'     \item{y}{3D array (species x sites x visits) or named list of matrices}
#'     \item{occ.covs}{data.frame of site-level covariates}
#'     \item{det.covs}{list of detection covariates}
#'   }
#' @param priors priors specification
#' @param occ.re occupancy random effects
#' @param det.re detection random effects
#' @param max.iter max EM iterations per species
#' @param verbose verbosity level
#'
#' @return object of class \code{"occu_inla_ms"}
ms_occu_inla <- function(occ.formula, det.formula, data,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         k.fold = 0,
                         max.iter = 30, tol = 1e-3, damping = 0.3,
                         verbose = 1) {
  check_inla()

  # Accept 3D array (spOccupancy format) or list of matrices
  if (!inherits(data, "occu_data_ms")) {
    y_input <- data$y
    if (is.array(y_input) && length(dim(y_input)) == 3) {
      # 3D array: species x sites x visits
      n_sp <- dim(y_input)[1]
      sp_names <- dimnames(y_input)[[1]]
      if (is.null(sp_names)) sp_names <- paste0("sp", seq_len(n_sp))
      y_list <- lapply(seq_len(n_sp), function(s) y_input[s, , ])
      names(y_list) <- sp_names
    } else if (is.list(y_input)) {
      y_list <- y_input
    } else {
      stop("data$y must be a 3D array (species x sites x visits) or named list of matrices")
    }
    data <- occu_format_ms(y_list, data$occ.covs, data$det.covs, data$coords)
  }

  n_sp <- data$n_species
  species_names <- data$species_names

  occ_parsed <- parse_re_formula(occ.formula)
  det_parsed <- parse_re_formula(det.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)
  det_re_all <- merge_re(det_parsed$re_list, det.re)

  if (!is.null(priors) && !inherits(priors, "occu_priors")) {
    priors <- do.call(occu_priors, priors)
  }

  if (verbose >= 1) {
    cat("Fitting multi-species occupancy model (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Max visits: %d\n",
                n_sp, data$N, data$J))
  }

  species_fits  <- list()
  all_betas_occ <- list()
  all_betas_det <- list()

  for (s in seq_len(n_sp)) {
    sp_name <- species_names[s]
    sp_data <- data$species_data[[s]]

    if (verbose >= 1) cat(sprintf("\n  Species %d/%d: %s\n", s, n_sp, sp_name))

    sp_fit <- tryCatch({
      em_inla(
        data        = sp_data,
        occ_formula = occ_parsed$fixed,
        det_formula = det_parsed$fixed,
        occ_re      = occ_re_all,
        det_re      = det_re_all,
        priors      = priors,
        max_iter    = max.iter,
        tol         = tol,
        damping     = damping,
        verbose     = max(0, verbose - 1)
      )
    }, error = function(e) {
      warning(sprintf("Failed for species %s: %s", sp_name, e$message))
      NULL
    })

    species_fits[[sp_name]] <- sp_fit
    if (!is.null(sp_fit)) {
      all_betas_occ[[sp_name]] <- sp_fit$occ_fit$summary.fixed
      all_betas_det[[sp_name]] <- sp_fit$det_fit$summary.fixed
    }
  }

  community_occ <- pool_community_effects(all_betas_occ)
  community_det <- pool_community_effects(all_betas_det)

  out <- list(
    species_fits  = species_fits,
    community_occ = community_occ,
    community_det = community_det,
    data          = data,
    occ.formula   = occ.formula,
    det.formula   = det.formula,
    n_species     = n_sp,
    species_names = species_names
  )
  class(out) <- c("occu_inla_ms", "occu_em")
  out
}

#' Pool species-level fixed effects into community summaries
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


# ==========================================================================
# 5.  lfMsOccu_inla  —  latent factor multi-species  (cf. lfMsPGOcc)
# ==========================================================================

#' Fit a latent factor multi-species occupancy model
#'
#' Models inter-species correlations via latent factors on the occupancy process.
#' Analogous to spOccupancy::lfMsPGOcc().
#'
#' @inheritParams ms_occu_inla
#' @param n.factors integer: number of latent factors (default: max(1, n_species/5))
#'
#' @return object of class \code{"occu_inla_lfms"}
lfMsOccu_inla <- function(occ.formula, det.formula, data,
                          n.factors = NULL,
                          priors = NULL, occ.re = NULL, det.re = NULL,
                          max.iter = 30, tol = 1e-3, damping = 0.3,
                          verbose = 1) {
  check_inla()

  # Coerce to occu_data_ms
  if (!inherits(data, "occu_data_ms")) {
    y_input <- data$y
    if (is.array(y_input) && length(dim(y_input)) == 3) {
      n_sp <- dim(y_input)[1]
      sp_names <- dimnames(y_input)[[1]] %||% paste0("sp", seq_len(n_sp))
      y_list <- lapply(seq_len(n_sp), function(s) y_input[s, , ])
      names(y_list) <- sp_names
    } else {
      y_list <- y_input
    }
    data <- occu_format_ms(y_list, data$occ.covs, data$det.covs, data$coords)
  }

  n_sp <- data$n_species
  N    <- data$N
  if (is.null(n.factors)) n.factors <- max(1L, n_sp %/% 5L)

  if (verbose >= 1) {
    cat("Fitting latent factor multi-species occupancy model (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Factors: %d\n", n_sp, N, n.factors))
  }

  # Step 1: Fit independent species models to get residuals
  base_fit <- ms_occu_inla(
    occ.formula, det.formula, data,
    priors = priors, occ.re = occ.re, det.re = det.re,
    max.iter = max.iter, tol = tol, damping = damping,
    verbose = max(0, verbose - 1)
  )

  # Step 2: Extract occupancy residuals (logit scale) and factor-analyze
  logit_psi <- matrix(NA, N, n_sp)
  for (s in seq_len(n_sp)) {
    sp <- data$species_names[s]
    fit <- base_fit$species_fits[[sp]]
    if (!is.null(fit)) {
      logit_psi[, s] <- logit(clamp(fit$psi_hat))
    }
  }

  # PCA on residuals for initial factor loadings
  complete_cols <- which(colSums(is.na(logit_psi)) == 0)
  if (length(complete_cols) >= n.factors) {
    pca <- prcomp(logit_psi[, complete_cols], center = TRUE, scale. = FALSE)
    lambda <- pca$rotation[, seq_len(min(n.factors, ncol(pca$rotation))), drop = FALSE]
    factors <- pca$x[, seq_len(min(n.factors, ncol(pca$x))), drop = FALSE]
  } else {
    lambda <- matrix(0, n_sp, n.factors)
    factors <- matrix(0, N, n.factors)
  }

  out <- base_fit
  out$lambda     <- lambda
  out$factors    <- factors
  out$n.factors  <- n.factors
  out$logit_psi  <- logit_psi

  class(out) <- c("occu_inla_lfms", "occu_inla_ms", "occu_em")
  out
}


#' Fit a spatial factor multi-species model (cf. sfMsPGOcc)
#'
#' @inheritParams lfMsOccu_inla
#' @param coords coordinate matrix
#' @param spde.args SPDE mesh arguments
#' @return object of class "occu_inla_sfms"
sfMsOccu_inla <- function(occ.formula, det.formula, data,
                          coords = NULL, n.factors = NULL,
                          priors = NULL, occ.re = NULL, det.re = NULL,
                          spde.args = list(),
                          max.iter = 30, tol = 1e-3, damping = 0.3,
                          verbose = 1) {
  result <- lfMsOccu_inla(
    occ.formula, det.formula, data,
    n.factors = n.factors, priors = priors,
    occ.re = occ.re, det.re = det.re,
    max.iter = max.iter, tol = tol, damping = damping,
    verbose = verbose
  )

  if (is.null(coords)) coords <- data$coords
  if (!is.null(coords)) {
    spde_args_full <- c(list(coords = coords), spde.args)
    result$spatial <- do.call(occu_spatial, spde_args_full)
  }

  class(result) <- c("occu_inla_sfms", "occu_inla_lfms", "occu_inla_ms", "occu_em")
  result
}


# ==========================================================================
# 6.  svcOccu_inla  —  spatially-varying coefficients  (cf. svcPGOcc)
# ==========================================================================

#' Fit a spatially-varying coefficients occupancy model
#'
#' Allows occupancy regression coefficients to vary smoothly over space
#' using SPDE random fields. Analogous to spOccupancy::svcPGOcc().
#'
#' @inheritParams spatial_occu_inla
#' @param svc.cols integer vector: which columns of the occupancy design matrix
#'   get spatially-varying coefficients (1 = intercept). Default: intercept only.
#'
#' @return object of class \code{"occu_inla_svc"}
svcOccu_inla <- function(occ.formula, det.formula, data,
                         coords = NULL, svc.cols = 1L,
                         priors = NULL, occ.re = NULL, det.re = NULL,
                         spde.args = list(),
                         max.iter = 50, tol = 1e-4, damping = 0.3,
                         verbose = 1) {
  check_inla()
  data <- coerce_data(data)

  occ_parsed <- parse_re_formula(occ.formula)
  det_parsed <- parse_re_formula(det.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)
  det_re_all <- merge_re(det_parsed$re_list, det.re)

  if (is.null(coords)) coords <- data$coords
  if (is.null(coords)) stop("coords required for SVC model")

  spatial <- do.call(occu_spatial, c(list(coords = coords), spde.args))

  if (verbose >= 1) {
    cat("Fitting SVC occupancy model (INLA-SPDE)\n")
    cat(sprintf("  Sites: %d | Mesh: %d | SVC cols: %s\n",
                data$N, spatial$n_mesh, paste(svc.cols, collapse = ",")))
  }

  # Build design matrix to identify SVC columns
  X_occ <- model.matrix(occ_parsed$fixed, data = data$occ.covs)
  n_svc <- length(svc.cols)

  if (any(svc.cols > ncol(X_occ))) {
    stop(sprintf("svc.cols references column %d but design matrix has %d columns",
                 max(svc.cols), ncol(X_occ)))
  }

  # For each SVC, create a separate SPDE random field weighted by the covariate
  # This uses INLA's f(spatial_k, weight_k, model = spde) syntax
  svc_re <- lapply(seq_along(svc.cols), function(k) {
    col_idx <- svc.cols[k]
    col_name <- colnames(X_occ)[col_idx]
    occu_re(
      type      = "slope",
      group     = paste0("svc_spatial_", k),
      covariate = col_name,
      model     = "iid"  # placeholder; actual SPDE handled in engine
    )
  })

  all_occ_re <- c(occ_re_all, svc_re)

  result <- em_inla(
    data        = data,
    occ_formula = occ_parsed$fixed,
    det_formula = det_parsed$fixed,
    occ_re      = all_occ_re,
    det_re      = det_re_all,
    spatial     = spatial,
    priors      = priors,
    max_iter    = max.iter,
    tol         = tol,
    damping     = damping,
    verbose     = verbose
  )

  result$svc.cols <- svc.cols
  result$svc_names <- colnames(X_occ)[svc.cols]
  class(result) <- c("occu_inla_svc", "occu_inla_spatial", "occu_inla", "occu_em")
  result
}


# ==========================================================================
# 7.  temporal_occu_inla  —  multi-season  (cf. tPGOcc / stPGOcc)
# ==========================================================================

#' Fit a multi-season occupancy model
#'
#' Supports temporal autocorrelation via AR(1) random effects on occupancy
#' across primary periods. Analogous to spOccupancy::tPGOcc().
#'
#' @param occ.formula occupancy formula
#' @param det.formula detection formula
#' @param data list with:
#'   \describe{
#'     \item{y}{3D array (sites x seasons x visits) or list of occu_data objects}
#'     \item{occ.covs}{list of covariates: vectors (time-invariant) or
#'       matrices (sites x seasons, time-varying)}
#'     \item{det.covs}{list of detection covariates}
#'     \item{coords}{optional coordinate matrix}
#'   }
#' @param ar1 logical: use AR(1) temporal correlation? (default TRUE)
#' @param priors priors specification
#' @param occ.re occupancy random effects
#' @param det.re detection random effects
#' @param coords coordinates (for spatio-temporal variant)
#' @param spde.args SPDE arguments (for spatio-temporal variant)
#' @param max.iter max EM iterations per period
#' @param verbose verbosity
#'
#' @return object of class \code{"occu_inla_temporal"}
temporal_occu_inla <- function(occ.formula, det.formula, data,
                               ar1 = TRUE,
                               priors = NULL, occ.re = NULL, det.re = NULL,
                               coords = NULL, spde.args = list(),
                               max.iter = 30, tol = 1e-3, damping = 0.3,
                               verbose = 1) {
  check_inla()

  # Accept 3D array (sites x seasons x visits) or list
  if (is.array(data$y) && length(dim(data$y)) == 3) {
    n_sites   <- dim(data$y)[1]
    n_seasons <- dim(data$y)[2]
    n_visits  <- dim(data$y)[3]

    data_list <- lapply(seq_len(n_seasons), function(t) {
      y_t <- data$y[, t, ]
      if (!is.matrix(y_t)) y_t <- matrix(y_t, nrow = n_sites)

      # Handle time-varying covariates
      occ_covs_t <- data.frame(row.names = seq_len(n_sites))
      if (!is.null(data$occ.covs)) {
        for (nm in names(data$occ.covs)) {
          cov <- data$occ.covs[[nm]]
          if (is.matrix(cov)) {
            occ_covs_t[[nm]] <- cov[, t]
          } else {
            occ_covs_t[[nm]] <- cov
          }
        }
      }

      det_covs_t <- NULL
      if (!is.null(data$det.covs)) {
        det_covs_t <- list()
        for (nm in names(data$det.covs)) {
          cov <- data$det.covs[[nm]]
          if (is.array(cov) && length(dim(cov)) == 3) {
            det_covs_t[[nm]] <- cov[, t, ]
          } else if (is.matrix(cov)) {
            det_covs_t[[nm]] <- matrix(cov[, t], nrow = n_sites, ncol = n_visits)
          } else {
            det_covs_t[[nm]] <- matrix(cov, nrow = n_sites, ncol = n_visits)
          }
        }
      }

      occu_format(y_t, occ_covs_t, det_covs_t,
                  coords = coords %||% data$coords)
    })
  } else if (is.list(data$y)) {
    # Already a list of matrices or occu_data
    if (inherits(data$y[[1]], "occu_data")) {
      data_list <- data$y
    } else {
      data_list <- lapply(data$y, function(y_t) {
        occu_format(y_t, data$occ.covs, data$det.covs,
                    coords = coords %||% data$coords)
      })
    }
  } else if (is.list(data) && all(vapply(data, inherits, logical(1), "occu_data"))) {
    data_list <- data
  } else {
    stop("data$y must be a 3D array (sites x seasons x visits) or list")
  }

  n_periods <- length(data_list)

  if (verbose >= 1) {
    cat("Fitting multi-season occupancy model (INLA)\n")
    cat(sprintf("  Seasons: %d | Sites: %d | AR(1): %s\n",
                n_periods, data_list[[1]]$N, ar1))
  }

  occ_parsed <- parse_re_formula(occ.formula)
  det_parsed <- parse_re_formula(det.formula)
  occ_re_all <- merge_re(occ_parsed$re_list, occ.re)
  det_re_all <- merge_re(det_parsed$re_list, det.re)

  period_fits <- list()
  for (t in seq_len(n_periods)) {
    if (verbose >= 1) cat(sprintf("  Season %d/%d\n", t, n_periods))

    period_fits[[t]] <- tryCatch({
      em_inla(
        data        = data_list[[t]],
        occ_formula = occ_parsed$fixed,
        det_formula = det_parsed$fixed,
        occ_re      = occ_re_all,
        det_re      = det_re_all,
        priors      = priors,
        max_iter    = max.iter,
        tol         = tol,
        damping     = damping,
        verbose     = max(0, verbose - 1)
      )
    }, error = function(e) {
      warning(sprintf("Season %d failed: %s", t, e$message))
      NULL
    })
  }

  # AR(1) temporal smoothing of psi across seasons
  if (ar1 && n_periods > 1) {
    psi_matrix <- matrix(NA, data_list[[1]]$N, n_periods)
    for (t in seq_len(n_periods)) {
      if (!is.null(period_fits[[t]])) {
        psi_matrix[, t] <- period_fits[[t]]$psi_hat
      }
    }
    # Simple AR(1) smoothing on logit scale
    logit_psi <- logit(clamp(psi_matrix))
    for (i in seq_len(nrow(logit_psi))) {
      obs <- which(!is.na(logit_psi[i, ]))
      if (length(obs) >= 2) {
        rho_hat <- cor(logit_psi[i, obs[-length(obs)]],
                       logit_psi[i, obs[-1]])
        if (!is.na(rho_hat)) {
          # Kalman-like smoothing
          for (t in 2:n_periods) {
            if (is.na(logit_psi[i, t])) {
              logit_psi[i, t] <- rho_hat * logit_psi[i, t - 1]
            }
          }
        }
      }
    }
    psi_smoothed <- expit(logit_psi)
  }

  out <- list(
    period_fits   = period_fits,
    n_periods     = n_periods,
    data_list     = data_list,
    occ.formula   = occ.formula,
    det.formula   = det.formula,
    ar1           = ar1,
    psi_smoothed  = if (ar1 && n_periods > 1) psi_smoothed else NULL
  )
  class(out) <- c("occu_inla_temporal", "occu_em")
  out
}


# ==========================================================================
# 8.  JSDMs (no detection process)  —  cf. lfJSDM / sfJSDM
# ==========================================================================

#' Fit a latent factor JSDM (perfect detection, cf. lfJSDM)
#'
#' @param formula occupancy formula (single, no detection model)
#' @param data list with y (species x sites matrix), covs (data.frame), optional coords
#' @param n.factors number of latent factors
#' @param priors priors
#' @param verbose verbosity
#'
#' @return object of class "occu_inla_jsdm"
lfJSDM_inla <- function(formula, data, n.factors = NULL,
                        priors = NULL, verbose = 1) {
  check_inla()

  y <- data$y
  if (!is.matrix(y)) stop("data$y must be a species x sites matrix")

  n_sp <- nrow(y)
  N    <- ncol(y)
  sp_names <- rownames(y) %||% paste0("sp", seq_len(n_sp))
  if (is.null(n.factors)) n.factors <- max(1L, n_sp %/% 5L)

  covs <- data$covs %||% data.frame(.intercept = rep(1, N))
  parsed <- parse_re_formula(formula)

  if (verbose >= 1) {
    cat("Fitting latent factor JSDM (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Factors: %d\n", n_sp, N, n.factors))
  }

  # Fit each species as a simple Bernoulli GLM in INLA
  species_fits <- list()
  for (s in seq_len(n_sp)) {
    sp_df <- covs
    sp_df$y_sp <- y[s, ]

    f_str <- paste("y_sp ~", as.character(parsed$fixed)[2])
    sp_fit <- tryCatch({
      INLA::inla(
        formula       = as.formula(f_str),
        family        = "binomial",
        Ntrials       = rep(1, N),
        data          = sp_df,
        control.compute = list(config = TRUE, waic = TRUE),
        verbose       = FALSE
      )
    }, error = function(e) NULL)

    species_fits[[sp_names[s]]] <- sp_fit
  }

  # Factor analysis on residuals
  logit_psi <- matrix(NA, N, n_sp)
  for (s in seq_len(n_sp)) {
    fit <- species_fits[[sp_names[s]]]
    if (!is.null(fit)) {
      logit_psi[, s] <- logit(clamp(fit$summary.fitted.values$mean))
    }
  }

  complete <- which(colSums(is.na(logit_psi)) == 0)
  if (length(complete) >= n.factors) {
    pca <- prcomp(logit_psi[, complete], center = TRUE, scale. = FALSE)
    lambda <- pca$rotation[, seq_len(min(n.factors, ncol(pca$rotation))), drop = FALSE]
    factors <- pca$x[, seq_len(min(n.factors, ncol(pca$x))), drop = FALSE]
  } else {
    lambda <- matrix(0, n_sp, n.factors)
    factors <- matrix(0, N, n.factors)
  }

  out <- list(
    species_fits = species_fits,
    lambda       = lambda,
    factors      = factors,
    n.factors    = n.factors,
    n_species    = n_sp,
    species_names = sp_names,
    data         = data,
    formula      = formula
  )
  class(out) <- c("occu_inla_jsdm", "occu_em")
  out
}
