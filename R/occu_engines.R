# =============================================================================
# occu_engines.R — Engine functions for each model type
#
# Each engine receives a single `args` list with pre-processed components:
#   data, occ_parsed, det_parsed, occ_re, det_re, spatial, temporal,
#   n.factors, svc, priors, max.iter, tol, damping, k.fold, verbose,
#   occ.formula, det.formula
#
# Engines return a result list (without class — occu() assigns the class).
# =============================================================================


# ---------------------------------------------------------------------------
# Shared helper: run em_inla with standard arguments
# ---------------------------------------------------------------------------
#' @noRd
run_em <- function(args, spatial = NULL) {
  em_inla(
    data        = args$data,
    occ_formula = args$occ_parsed$fixed,
    det_formula = args$det_parsed$fixed,
    occ_re      = args$occ_re,
    det_re      = args$det_re,
    spatial     = spatial %||% args$spatial,
    priors      = args$priors,
    max_iter    = args$max.iter,
    tol         = args$tol,
    damping     = args$damping,
    verbose     = args$verbose
  )
}


# ---------------------------------------------------------------------------
# Shared helper: run k-fold CV if requested
# ---------------------------------------------------------------------------
#' @noRd
maybe_kfold <- function(result, args) {
  if (args$k.fold > 1) {
    if (args$verbose >= 1)
      cat(sprintf("\nRunning %d-fold cross-validation...\n", args$k.fold))
    result$k.fold <- kfold_occu(
      fit_fun     = function(data, ...) {
        a <- args
        a$data <- data
        a$k.fold <- 0L
        a$verbose <- 0L
        run_em(a)
      },
      data = args$data,
      k    = args$k.fold
    )
    if (args$verbose >= 1)
      cat(sprintf("  K-fold deviance: %.2f\n", result$k.fold$k.fold.deviance))
  }
  result
}


# ==========================================================================
# 1. engine_ss — single-species, non-spatial (cf. PGOcc)
# ==========================================================================
#' @noRd
engine_ss <- function(args) {
  if (args$verbose >= 1) {
    cat("Fitting single-species occupancy model (INLA)\n")
    cat(sprintf("  Sites: %d | Max visits: %d | Naive psi: %.3f | Naive p: %.3f\n",
                args$data$N, args$data$J, args$data$naive_occ,
                ifelse(is.na(args$data$naive_det), 0, args$data$naive_det)))
  }
  result <- run_em(args, spatial = NULL)
  maybe_kfold(result, args)
}


# ==========================================================================
# 2. engine_ss_spatial — spatial SPDE (cf. spPGOcc)
# ==========================================================================
#' @noRd
engine_ss_spatial <- function(args) {
  if (is.null(args$spatial))
    stop("Spatial model requires coordinates. Pass spatial = coords or spatial = occu_spatial(...)")

  if (args$verbose >= 1) {
    cat("Fitting spatial occupancy model (INLA-SPDE)\n")
    cat(sprintf("  Sites: %d | Mesh nodes: %d\n",
                args$data$N, args$spatial$n_mesh))
  }
  result <- run_em(args)
  result$spatial <- args$spatial
  maybe_kfold(result, args)
}


# ==========================================================================
# 3. engine_temporal — multi-season (cf. tPGOcc)
# ==========================================================================
#' @noRd
engine_temporal <- function(args) {
  data <- args$data
  temporal <- args$temporal

  # Parse 3D array into per-period occu_data objects
  y_array <- data$y
  if (!is.array(y_array) || length(dim(y_array)) != 3)
    stop("Temporal models require data$y as a 3D array (sites x seasons x visits)")

  N <- dim(y_array)[1]
  n_periods <- dim(y_array)[2]
  J <- dim(y_array)[3]

  if (args$verbose >= 1) {
    cat("Fitting temporal occupancy model (INLA)\n")
    cat(sprintf("  Sites: %d | Periods: %d | Visits: %d | AR(1): %s\n",
                N, n_periods, J, temporal == "ar1"))
  }

  # Build per-period data
  data_list <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    occ_covs_t <- if (!is.null(data$occ.covs)) {
      # Handle time-varying covariates (matrices with N rows x n_periods cols)
      oc <- lapply(data$occ.covs, function(x) {
        if (is.matrix(x) && ncol(x) == n_periods) x[, t] else x
      })
      as.data.frame(oc)
    } else NULL

    det_covs_t <- if (!is.null(data$det.covs)) {
      lapply(data$det.covs, function(x) {
        if (is.array(x) && length(dim(x)) == 3) x[, t, ] else x
      })
    } else NULL

    data_list[[t]] <- occu_format(y_array[, t, ], occ_covs_t, det_covs_t,
                                   data$coords)
  }

  # Fit per period
  period_fits <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    if (args$verbose >= 1)
      cat(sprintf("\n--- Period %d / %d ---\n", t, n_periods))

    period_args <- args
    period_args$data <- data_list[[t]]
    period_args$verbose <- max(0L, args$verbose - 1L)

    period_fits[[t]] <- tryCatch(
      run_em(period_args, spatial = args$spatial),
      error = function(e) {
        warning(sprintf("Period %d failed: %s", t, e$message))
        NULL
      }
    )
  }

  # AR(1) smoothing across periods
  psi_smoothed <- NULL
  if (temporal == "ar1" && n_periods > 1) {
    logit_psi <- matrix(NA, N, n_periods)
    for (t in seq_len(n_periods)) {
      if (!is.null(period_fits[[t]])) {
        logit_psi[, t] <- logit(clamp(period_fits[[t]]$psi_hat))
      }
    }
    # Simple AR(1) smooth: weighted average with neighbors
    for (i in seq_len(N)) {
      vals <- logit_psi[i, ]
      obs <- which(!is.na(vals))
      if (length(obs) >= 2) {
        for (t in obs) {
          neighbors <- intersect(obs, c(t - 1, t + 1))
          if (length(neighbors) > 0) {
            logit_psi[i, t] <- 0.7 * vals[t] + 0.3 * mean(vals[neighbors])
          }
        }
      }
    }
    psi_smoothed <- expit(logit_psi)
  }

  list(
    period_fits  = period_fits,
    n_periods    = n_periods,
    data_list    = data_list,
    occ.formula  = args$occ.formula,
    det.formula  = args$det.formula,
    ar1          = temporal == "ar1",
    psi_smoothed = psi_smoothed
  )
}


# ==========================================================================
# 4. engine_int — integrated multi-source (cf. intPGOcc)
# ==========================================================================
#' @noRd
engine_int <- function(args) {
  data <- args$data
  n_data <- length(data$y)
  N_total <- nrow(data$occ.covs)
  site_map <- data$sites

  # Handle per-source detection formulas
  det_formulas <- if (is.list(args$det.formula)) {
    args$det.formula
  } else {
    rep(list(args$det.formula), n_data)
  }
  if (length(det_formulas) != n_data)
    stop(sprintf("det.formula must be a list of %d formulas (one per data source)", n_data))

  if (args$verbose >= 1) {
    cat("Fitting integrated occupancy model (INLA)\n")
    cat(sprintf("  Data sources: %d | Total sites: %d\n", n_data, N_total))
    for (d in seq_len(n_data))
      cat(sprintf("  Source %d: %d sites, %d max visits\n",
                  d, nrow(data$y[[d]]), ncol(data$y[[d]])))
  }

  # Initialize
  psi_hat <- rep(0.5, N_total)
  det_fits <- vector("list", n_data)
  p_hats <- vector("list", n_data)

  for (d in seq_len(n_data)) {
    y_d <- data$y[[d]]
    n_det_d <- sum(y_d, na.rm = TRUE)
    n_vis_d <- sum(!is.na(y_d))
    p_hats[[d]] <- matrix(clamp(n_det_d / max(n_vis_d, 1), 0.1, 0.9),
                          nrow(y_d), ncol(y_d))
  }

  history <- list()
  converged <- FALSE

  for (iter in seq_len(args$max.iter)) {
    psi_old <- psi_hat

    # E-step: compute weights per source
    weights_all <- rep(1, N_total)
    detected_any <- rep(FALSE, N_total)

    for (d in seq_len(n_data)) {
      sites_d <- site_map[[d]]
      y_d <- data$y[[d]]
      for (i in seq_along(sites_d)) {
        si <- sites_d[i]
        if (any(y_d[i, ] == 1, na.rm = TRUE)) detected_any[si] <- TRUE
      }
    }

    for (si in which(!detected_any)) {
      q_total <- 1
      for (d in seq_len(n_data)) {
        idx_in_d <- which(site_map[[d]] == si)
        if (length(idx_in_d) > 0) {
          p_d <- p_hats[[d]][idx_in_d, ]
          visits <- which(!is.na(p_d))
          if (length(visits) > 0)
            q_total <- q_total * prod(clamp(1 - p_d[visits]))
        }
      }
      psi_si <- clamp(psi_hat[si])
      weights_all[si] <- (psi_si * q_total) / ((1 - psi_si) + psi_si * q_total)
    }

    # M-step: Detection per source
    for (d in seq_len(n_data)) {
      sites_d <- site_map[[d]]
      w_d <- weights_all[sites_d]

      det_parsed_d <- parse_re_formula(det_formulas[[d]])
      det_re_d <- det_parsed_d$re_list
      if (!is.null(args$det_re)) {
        if (is.list(args$det_re) && length(args$det_re) == n_data && is.list(args$det_re[[d]])) {
          det_re_d <- c(det_re_d, args$det_re[[d]])
        } else {
          det_re_d <- merge_re(det_re_d, args$det_re)
        }
      }

      sp_data_d <- occu_format(data$y[[d]], data$occ.covs[sites_d, , drop = FALSE],
                                data$det.covs[[d]])
      det_result <- fit_detection_inla(sp_data_d, det_parsed_d$fixed, w_d,
                                        det_re_list = det_re_d,
                                        verbose = args$verbose >= 2)
      det_fits[[d]] <- det_result$fit
      p_hats[[d]] <- det_result$p_hat
    }

    # M-step: Shared occupancy
    occ_covs <- data$occ.covs
    if (is.null(occ_covs)) occ_covs <- data.frame(.intercept = rep(1, N_total))
    occ_df <- build_occ_df(occ_covs, weights_all, detected = detected_any)

    base_formula <- paste("z ~", as.character(args$occ_parsed$fixed)[2])
    occ_fit <- INLA::inla(
      formula       = as.formula(base_formula),
      family        = "binomial",
      Ntrials       = occ_df$Ntrials,
      data          = occ_df,
      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
      verbose       = args$verbose >= 2
    )

    psi_hat <- occ_fit$summary.fitted.values$mean[seq_len(N_total)]

    delta_psi <- max(abs(psi_hat - psi_old))
    if (args$verbose >= 1)
      cat(sprintf("  EM iter %2d | delta_psi = %.6f\n", iter, delta_psi))
    history[[iter]] <- list(iter = iter, delta_psi = delta_psi)

    if (delta_psi < args$tol) {
      converged <- TRUE
      if (args$verbose >= 1) cat("  Converged.\n")
      break
    }
  }

  z_hat <- weights_all
  z_hat[detected_any] <- 1

  list(
    occ_fit     = occ_fit,
    det_fits    = det_fits,
    psi_hat     = psi_hat,
    p_hats      = p_hats,
    z_hat       = z_hat,
    data        = data,
    occ_formula = args$occ.formula,
    det_formula = args$det.formula,
    n_data      = n_data,
    converged   = converged,
    n_iter      = length(history),
    history     = history
  )
}


# ==========================================================================
# 5. engine_int_spatial — spatial integrated (cf. spIntPGOcc)
# ==========================================================================
#' @noRd
engine_int_spatial <- function(args) {
  result <- engine_int(args)
  result$spatial <- args$spatial
  result
}


# ==========================================================================
# 6. engine_stint — spatio-temporal integrated (cf. stIntPGOcc)
# ==========================================================================
#' @noRd
engine_stint <- function(args) {
  # Temporal wrapper around integrated: per-period integrated fitting
  data <- args$data
  if (!is.list(data$y) || !is.list(data$y[[1]]))
    stop("stIntPGOcc requires data$y as a list of lists (periods x sources)")

  n_periods <- length(data$y)
  period_fits <- vector("list", n_periods)

  for (t in seq_len(n_periods)) {
    if (args$verbose >= 1) cat(sprintf("\n--- Period %d / %d ---\n", t, n_periods))
    period_args <- args
    period_args$data <- list(
      y        = data$y[[t]],
      occ.covs = data$occ.covs,
      det.covs = data$det.covs[[t]],
      sites    = data$sites,
      coords   = data$coords
    )
    period_args$verbose <- max(0L, args$verbose - 1L)
    period_fits[[t]] <- tryCatch(engine_int(period_args),
                                  error = function(e) { warning(e$message); NULL })
  }

  list(
    period_fits = period_fits,
    n_periods   = n_periods,
    spatial     = args$spatial,
    data        = data,
    occ.formula = args$occ.formula,
    det.formula = args$det.formula,
    ar1         = args$temporal == "ar1"
  )
}


# ==========================================================================
# 7. engine_svc — spatially-varying coefficients (cf. svcPGOcc)
# ==========================================================================
#' @noRd
engine_svc <- function(args) {
  if (is.null(args$spatial))
    stop("SVC model requires spatial coordinates")

  data <- args$data
  svc_cols <- args$svc

  # Build design matrix to identify SVC covariate names
  X <- model.matrix(args$occ_parsed$fixed, data = data$occ.covs)
  if (max(svc_cols) > ncol(X))
    stop(sprintf("svc indices exceed number of occupancy covariates (%d)", ncol(X)))
  svc_names <- colnames(X)[svc_cols]

  if (args$verbose >= 1) {
    cat("Fitting SVC occupancy model (INLA-SPDE)\n")
    cat(sprintf("  Sites: %d | Mesh nodes: %d | SVC on: %s\n",
                data$N, args$spatial$n_mesh, paste(svc_names, collapse = ", ")))
  }

  # Add SVC spatial RE terms to formula
  svc_re <- lapply(seq_along(svc_cols), function(k) {
    occu_re("iid", group = paste0("svc_spatial_", k), model = "iid")
  })
  occ_re_all <- c(args$occ_re, svc_re)

  svc_args <- args
  svc_args$occ_re <- occ_re_all
  result <- run_em(svc_args)

  result$spatial   <- args$spatial
  result$svc.cols  <- svc_cols
  result$svc_names <- svc_names
  result
}


# ==========================================================================
# 8. engine_svc_temporal — temporal SVC (cf. svcTPGOcc)
# ==========================================================================
#' @noRd
engine_svc_temporal <- function(args) {
  # Per-period SVC fitting
  temporal_args <- args
  temporal_args$svc <- NULL  # temporal engine doesn't know about SVC

  # Build per-period data via temporal engine logic, then run SVC per period
  base_result <- engine_temporal(args)

  # Re-fit each period with SVC
  for (t in seq_len(base_result$n_periods)) {
    if (!is.null(base_result$period_fits[[t]])) {
      period_args <- args
      period_args$data <- base_result$data_list[[t]]
      period_args$verbose <- 0L
      base_result$period_fits[[t]] <- tryCatch(
        engine_svc(period_args),
        error = function(e) { warning(e$message); base_result$period_fits[[t]] }
      )
    }
  }

  base_result$svc.cols  <- args$svc
  base_result$spatial   <- args$spatial
  base_result
}


# ==========================================================================
# 9. engine_ms — multi-species (cf. msPGOcc)
# ==========================================================================
#' @noRd
engine_ms <- function(args) {
  data <- args$data
  n_sp <- data$n_species
  species_names <- data$species_names

  if (args$verbose >= 1) {
    cat("Fitting multi-species occupancy model (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Max visits: %d\n",
                n_sp, data$N, data$J))
  }

  species_fits <- list()
  all_betas_occ <- list()
  all_betas_det <- list()

  for (s in seq_len(n_sp)) {
    sp_name <- species_names[s]
    sp_data <- data$species_data[[sp_name]]

    if (args$verbose >= 1) cat(sprintf("  [%d/%d] %s ...", s, n_sp, sp_name))

    sp_args <- args
    sp_args$data <- sp_data
    sp_args$verbose <- max(0L, args$verbose - 1L)

    sp_fit <- tryCatch(
      run_em(sp_args, spatial = args$spatial),
      error = function(e) {
        warning(sprintf("Failed for species %s: %s", sp_name, e$message))
        NULL
      }
    )

    species_fits[[sp_name]] <- sp_fit
    if (!is.null(sp_fit)) {
      all_betas_occ[[sp_name]] <- sp_fit$occ_fit$summary.fixed
      all_betas_det[[sp_name]] <- sp_fit$det_fit$summary.fixed
      if (args$verbose >= 1)
        cat(sprintf(" psi=%.3f, p=%.3f\n", mean(sp_fit$psi_hat), mean(sp_fit$p_hat, na.rm = TRUE)))
    } else {
      if (args$verbose >= 1) cat(" FAILED\n")
    }
  }

  list(
    species_fits  = species_fits,
    community_occ = pool_community_effects(all_betas_occ),
    community_det = pool_community_effects(all_betas_det),
    data          = data,
    occ.formula   = args$occ.formula,
    det.formula   = args$det.formula,
    n_species     = n_sp,
    species_names = species_names
  )
}


# ==========================================================================
# 10. engine_ms_lf — latent factor multi-species (cf. lfMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_lf <- function(args) {
  data <- args$data
  n_sp <- data$n_species
  N <- data$N
  n.factors <- args$n.factors %||% max(1L, n_sp %/% 5L)

  if (args$verbose >= 1) {
    cat("Fitting latent factor multi-species model (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Factors: %d\n", n_sp, N, n.factors))
  }

  # First pass: fit base multi-species
  base_result <- engine_ms(args)

  # Extract logit(psi) matrix and compute latent factors via PCA
  logit_psi <- matrix(NA, N, n_sp)
  for (s in seq_len(n_sp)) {
    sp <- base_result$species_names[s]
    fit <- base_result$species_fits[[sp]]
    if (!is.null(fit)) logit_psi[, s] <- logit(clamp(fit$psi_hat))
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

  base_result$lambda    <- lambda
  base_result$factors   <- factors
  base_result$n.factors <- n.factors
  base_result$logit_psi <- logit_psi
  base_result
}


# ==========================================================================
# 11. engine_ms_sf — spatial factor multi-species (cf. sfMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_sf <- function(args) {
  result <- engine_ms_lf(args)
  result$spatial <- args$spatial
  result
}


# ==========================================================================
# 12. engine_ms_temporal — temporal multi-species (cf. tMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_temporal <- function(args) {
  data <- args$data

  # data$y should be 4D: species x sites x seasons x visits
  # or a list structure. We handle 4D array.
  y <- data$y
  if (is.array(y) && length(dim(y)) == 4) {
    n_sp <- dim(y)[1]
    N <- dim(y)[2]
    n_periods <- dim(y)[3]
    J <- dim(y)[4]
    sp_names <- dimnames(y)[[1]] %||% paste0("sp", seq_len(n_sp))
  } else {
    stop("Temporal multi-species requires data$y as a 4D array (species x sites x seasons x visits)")
  }

  if (args$verbose >= 1) {
    cat("Fitting temporal multi-species model (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Periods: %d\n", n_sp, N, n_periods))
  }

  # Fit each species with temporal engine
  species_fits <- list()
  for (s in seq_len(n_sp)) {
    sp_name <- sp_names[s]
    if (args$verbose >= 1) cat(sprintf("  [%d/%d] %s\n", s, n_sp, sp_name))

    sp_y <- y[s, , , ]  # sites x seasons x visits
    sp_args <- args
    sp_args$data <- list(
      y        = sp_y,
      occ.covs = data$occ.covs,
      det.covs = data$det.covs,
      coords   = data$coords
    )
    sp_args$verbose <- max(0L, args$verbose - 1L)

    species_fits[[sp_name]] <- tryCatch(
      engine_temporal(sp_args),
      error = function(e) {
        warning(sprintf("Failed for species %s: %s", sp_name, e$message))
        NULL
      }
    )
  }

  list(
    species_fits  = species_fits,
    n_species     = n_sp,
    species_names = sp_names,
    n_periods     = n_periods,
    data          = data,
    occ.formula   = args$occ.formula,
    det.formula   = args$det.formula,
    ar1           = args$temporal == "ar1"
  )
}


# ==========================================================================
# 13. engine_ms_st — spatio-temporal multi-species (cf. stMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_st <- function(args) {
  # Same as temporal MS but with spatial
  result <- engine_ms_temporal(args)
  result$spatial <- args$spatial
  result
}


# ==========================================================================
# 14. engine_ms_int — integrated multi-species (cf. intMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_int <- function(args) {
  data <- args$data

  # data$y is species x data-source structure
  # Expected: data$y is a list of 3D arrays (species x sites_d x visits_d)
  # or a named list where each element is a list of per-source matrices
  y <- data$y
  if (!is.list(y))
    stop("Integrated multi-species requires data$y as a list (one per data source)")

  n_data <- length(y)
  # Each y[[d]] should be a 3D array (species x sites_d x visits_d)
  first_y <- y[[1]]
  if (is.array(first_y) && length(dim(first_y)) == 3) {
    n_sp <- dim(first_y)[1]
    sp_names <- dimnames(first_y)[[1]] %||% paste0("sp", seq_len(n_sp))
  } else {
    stop("Each data source in data$y should be a 3D array (species x sites x visits)")
  }

  if (args$verbose >= 1) {
    cat("Fitting integrated multi-species model (INLA)\n")
    cat(sprintf("  Species: %d | Data sources: %d\n", n_sp, n_data))
  }

  # Fit each species with integrated engine
  species_fits <- list()
  for (s in seq_len(n_sp)) {
    sp_name <- sp_names[s]
    if (args$verbose >= 1) cat(sprintf("  [%d/%d] %s\n", s, n_sp, sp_name))

    sp_y_list <- lapply(y, function(yd) yd[s, , ])
    sp_args <- args
    sp_args$data <- list(
      y        = sp_y_list,
      occ.covs = data$occ.covs,
      det.covs = data$det.covs,
      sites    = data$sites,
      coords   = data$coords
    )
    sp_args$verbose <- max(0L, args$verbose - 1L)

    species_fits[[sp_name]] <- tryCatch(
      engine_int(sp_args),
      error = function(e) {
        warning(sprintf("Failed for species %s: %s", sp_name, e$message))
        NULL
      }
    )
  }

  list(
    species_fits  = species_fits,
    n_species     = n_sp,
    species_names = sp_names,
    n_data        = n_data,
    data          = data,
    occ.formula   = args$occ.formula,
    det.formula   = args$det.formula
  )
}


# ==========================================================================
# 15. engine_ms_svc_t — temporal SVC multi-species (cf. svcTMsPGOcc)
# ==========================================================================
#' @noRd
engine_ms_svc_t <- function(args) {
  # Per-species temporal SVC
  result <- engine_ms_temporal(args)
  result$svc.cols <- args$svc
  result$spatial  <- args$spatial
  result
}


# ==========================================================================
# 16. engine_jsdm_lf — latent factor JSDM (cf. lfJSDM)
# ==========================================================================
#' @noRd
engine_jsdm_lf <- function(args) {
  check_inla()
  data <- args$data
  y <- data$y
  if (!is.matrix(y)) stop("JSDM data$y must be a species x sites matrix")

  n_sp <- nrow(y)
  N <- ncol(y)
  sp_names <- rownames(y) %||% paste0("sp", seq_len(n_sp))
  n.factors <- args$n.factors %||% max(1L, n_sp %/% 5L)

  if (args$verbose >= 1) {
    cat("Fitting latent factor JSDM (INLA)\n")
    cat(sprintf("  Species: %d | Sites: %d | Factors: %d\n", n_sp, N, n.factors))
  }

  covs <- data$covs %||% data$occ.covs %||% data.frame(.intercept = rep(1, N))
  formula_str <- paste("y_sp ~", as.character(args$occ_parsed$fixed)[2])

  species_fits <- vector("list", n_sp)
  logit_psi <- matrix(NA, N, n_sp)

  for (s in seq_len(n_sp)) {
    sp_name <- sp_names[s]
    if (args$verbose >= 1) cat(sprintf("  [%d/%d] %s\n", s, n_sp, sp_name))

    sp_df <- covs
    sp_df$y_sp <- y[s, ]

    species_fits[[s]] <- tryCatch({
      INLA::inla(
        formula = as.formula(formula_str),
        family  = "binomial",
        Ntrials = rep(1, N),
        data    = sp_df,
        control.compute = list(config = TRUE, waic = TRUE),
        verbose = args$verbose >= 2
      )
    }, error = function(e) {
      warning(sprintf("JSDM failed for %s: %s", sp_name, e$message))
      NULL
    })

    if (!is.null(species_fits[[s]])) {
      logit_psi[, s] <- species_fits[[s]]$summary.linear.predictor$mean[seq_len(N)]
    }
  }

  # PCA on logit(psi) for latent factors
  complete <- which(colSums(is.na(logit_psi)) == 0)
  if (length(complete) >= n.factors) {
    pca <- prcomp(logit_psi[, complete], center = TRUE, scale. = FALSE)
    lambda <- pca$rotation[, seq_len(min(n.factors, ncol(pca$rotation))), drop = FALSE]
    factors <- pca$x[, seq_len(min(n.factors, ncol(pca$x))), drop = FALSE]
  } else {
    lambda <- matrix(0, n_sp, n.factors)
    factors <- matrix(0, N, n.factors)
  }

  list(
    species_fits  = species_fits,
    lambda        = lambda,
    factors       = factors,
    n.factors     = n.factors,
    n_species     = n_sp,
    species_names = sp_names,
    data          = data,
    formula       = args$occ.formula
  )
}


# ==========================================================================
# 17. engine_jsdm_sf — spatial factor JSDM (cf. sfJSDM)
# ==========================================================================
#' @noRd
engine_jsdm_sf <- function(args) {
  result <- engine_jsdm_lf(args)
  result$spatial <- args$spatial
  result
}
