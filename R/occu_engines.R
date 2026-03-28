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
    correction  = args$correction %||% "auto",
    num.threads = args$num.threads %||% "1:1",
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
#
# Joint EM across all periods:
#   - Detection: fit per-period (different data per decade)
#   - Occupancy: single stacked INLA call with f(period, model="ar1")
#     so temporal information is shared across decades
# ==========================================================================
#' @noRd
engine_temporal <- function(args) {
  check_inla()
  data <- args$data
  temporal <- args$temporal

  y_array <- data$y
  if (!is.array(y_array) || length(dim(y_array)) != 3)
    stop("Temporal models require data$y as a 3D array (sites x seasons x visits)")

  N <- dim(y_array)[1]
  n_periods <- dim(y_array)[2]
  J <- dim(y_array)[3]
  use_ar1 <- temporal == "ar1"

  if (args$verbose >= 1) {
    cat("Fitting temporal occupancy model (INLA, joint)\n")
    cat(sprintf("  Sites: %d | Periods: %d | Visits: %d | AR(1): %s\n",
                N, n_periods, J, use_ar1))
  }

  # --- Build per-period data for detection ---
  data_list <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    occ_covs_t <- if (!is.null(data$occ.covs)) {
      oc <- lapply(data$occ.covs, function(x) {
        if (is.matrix(x) && ncol(x) == n_periods) x[, t] else x
      })
      as.data.frame(oc)
    } else data.frame(.intercept = rep(1, N))

    det_covs_t <- if (!is.null(data$det.covs)) {
      lapply(data$det.covs, function(x) {
        if (is.array(x) && length(dim(x)) == 3) x[, t, ] else x
      })
    } else NULL

    data_list[[t]] <- occu_format(y_array[, t, ], occ_covs_t, det_covs_t,
                                   data$coords)
  }

  # --- Initialize ---
  psi_mat <- matrix(0.5, N, n_periods)  # psi[site, period]
  p_list  <- vector("list", n_periods)  # p_hat per period
  for (t in seq_len(n_periods)) {
    naive_p <- clamp(data_list[[t]]$naive_det, 0.1, 0.9)
    p_list[[t]] <- matrix(naive_p, N, J)
  }

  # Prep detection for each period
  det_preps <- lapply(seq_len(n_periods), function(t) {
    prep_detection(data_list[[t]], args$det_parsed$fixed, args$det_re)
  })

  # --- EM loop ---
  ctrl_fast <- list(config = TRUE, dic = FALSE, waic = FALSE)
  ctrl_full <- list(config = TRUE, dic = TRUE, waic = TRUE)
  converged <- FALSE
  history <- list()

  for (iter in seq_len(args$max.iter)) {
    psi_old <- psi_mat

    # --- M-step: Detection per period ---
    for (t in seq_len(n_periods)) {
      weights_t <- compute_weights(data_list[[t]]$y, psi_mat[, t], p_list[[t]])
      det_result <- fit_detection_inla(
        det_preps[[t]], data_list[[t]], weights_t,
        control.inla = NULL, control.compute = ctrl_fast,
        verbose = FALSE
      )
      p_new <- det_result$p_hat
      p_new[is.na(p_new)] <- p_list[[t]][is.na(p_new)]
      p_list[[t]] <- p_new
    }

    # --- M-step: Joint occupancy with AR1 ---
    # Stack all site-periods into one INLA call
    weights_all <- matrix(NA, N, n_periods)
    for (t in seq_len(n_periods)) {
      weights_all[, t] <- compute_weights(
        data_list[[t]]$y, psi_mat[, t], p_list[[t]]
      )
    }

    psi_mat <- fit_joint_occupancy(
      data_list, weights_all, args$occ_parsed$fixed,
      n_periods, use_ar1, ctrl_fast, args$verbose >= 2
    )

    # --- Convergence ---
    delta_psi <- max(abs(psi_mat - psi_old))

    history[[iter]] <- list(iter = iter, delta_psi = delta_psi)

    if (args$verbose >= 1) {
      cat(sprintf("  EM iter %2d | delta_psi = %.6f | mean_psi = %.3f\n",
                  iter, delta_psi, mean(psi_mat)))
    }

    if (delta_psi < args$tol) {
      converged <- TRUE
      if (args$verbose >= 1) cat("  Converged.\n")
      break
    }
  }

  if (!converged && args$verbose >= 1) {
    cat(sprintf("  Warning: EM did not converge after %d iterations.\n", args$max.iter))
  }

  # --- Final full-quality fit ---
  occ_fit <- fit_joint_occupancy(
    data_list, weights_all, args$occ_parsed$fixed,
    n_periods, use_ar1, ctrl_full, args$verbose >= 2,
    return_fit = TRUE
  )

  # Build per-period results for compatibility
  period_fits <- lapply(seq_len(n_periods), function(t) {
    list(psi_hat = psi_mat[, t], p_hat = p_list[[t]])
  })

  list(
    occ_fit      = occ_fit$fit,
    period_fits  = period_fits,
    psi_mat      = psi_mat,
    n_periods    = n_periods,
    data_list    = data_list,
    occ.formula  = args$occ.formula,
    det.formula  = args$det.formula,
    ar1          = use_ar1,
    converged    = converged,
    n_iter       = length(history),
    history      = history
  )
}


#' Fit joint occupancy model across all periods with AR1 temporal effect
#' @noRd
fit_joint_occupancy <- function(data_list, weights_all, occ_formula,
                                 n_periods, use_ar1, control.compute,
                                 verbose = FALSE, return_fit = FALSE) {
  check_inla()
  N <- data_list[[1]]$N

  # Stack all site-period rows using M=1000 binomial encoding.
  # MI debiasing is applied at the per-species level (engine_ms_temporal).
  rows <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    occ_df_t <- build_occ_df(
      occ_covs = data_list[[t]]$occ.covs,
      weights  = weights_all[, t],
      site_idx = seq_len(N),
      detected = data_list[[t]]$detected
    )
    occ_df_t$period <- t
    rows[[t]] <- occ_df_t
  }
  stacked <- do.call(rbind, rows)

  # Formula: covariates + AR1 temporal effect
  base_formula <- paste("z ~", as.character(occ_formula)[2])
  if (use_ar1) {
    base_formula <- paste(base_formula, "+ f(period, model = 'ar1')")
  } else {
    base_formula <- paste(base_formula, "+ f(period, model = 'iid')")
  }

  n_stack <- nrow(stacked)
  fit <- INLA::inla(
    formula         = as.formula(base_formula),
    family          = "binomial",
    Ntrials         = stacked$Ntrials,
    data            = stacked,
    num.threads     = "1:1",
    control.inla    = INLA::inla.set.control.inla.default(),
    control.compute = control.compute,
    verbose         = verbose
  )

  # Extract psi per site-period (direct indexing, no which())
  psi_fitted <- fit$summary.fitted.values$mean
  psi_mat <- matrix(NA, N, n_periods)
  psi_mat[cbind(stacked$site, stacked$period)] <- psi_fitted

  if (return_fit) {
    return(list(fit = fit, psi_mat = psi_mat))
  }
  psi_mat
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
      det_prep_d <- prep_detection(sp_data_d, det_parsed_d$fixed, det_re_d)
      det_result <- fit_detection_inla(
        det_prep_d, sp_data_d, w_d,
        control.compute = list(config = FALSE, dic = FALSE, waic = FALSE),
        verbose = args$verbose >= 2
      )
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

  # Add SVC spatial RE terms: one iid RE per SVC covariate, indexed by site.
  # The group column must exist in occ.covs before the RE machinery resolves it.
  svc_re <- vector("list", length(svc_cols))
  for (k in seq_along(svc_cols)) {
    col_name <- paste0("svc_spatial_", k)
    data$occ.covs[[col_name]] <- seq_len(data$N)
    svc_re[[k]] <- occu_re("iid", group = col_name, model = "iid")
  }
  occ_re_all <- c(args$occ_re, svc_re)

  svc_args <- args
  svc_args$data <- data
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

  # Per-species fitting function
  fit_one_species <- function(s) {
    sp_name <- species_names[s]
    sp_data <- data$species_data[[sp_name]]
    sp_args <- args
    sp_args$data <- sp_data
    sp_args$verbose <- 0L
    tryCatch(
      run_em(sp_args, spatial = args$spatial),
      error = function(e) {
        warning(sprintf("Failed for species %s: %s", sp_name, e$message))
        NULL
      }
    )
  }

  # Parallel if available and enough species to amortize cluster startup cost.
  # Socket clusters on Windows take ~3-5s to spawn + load INLA on each worker,

  # so parallel is only faster when there are enough species per core.
  n_cores <- getOption("INLAocc.cores", 1L)
  use_parallel <- n_cores > 1 &&
    n_sp >= 2L * n_cores &&
    requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    if (args$verbose >= 1) cat(sprintf("  Using %d cores\n", min(n_cores, n_sp)))
    cl <- parallel::makeCluster(min(n_cores, n_sp),
                                    port = sample(49152:65535, 1))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("args", "data", "species_names"),
                             envir = environment())
    parallel::clusterEvalQ(cl, suppressMessages(library(INLAocc)))
    species_fits_raw <- parallel::parLapply(cl, seq_len(n_sp), fit_one_species)
  } else {
    species_fits_raw <- lapply(seq_len(n_sp), function(s) {
      if (args$verbose >= 1) cat(sprintf("  [%d/%d] %s ...", s, n_sp, species_names[s]))
      result <- fit_one_species(s)
      if (args$verbose >= 1) {
        if (!is.null(result))
          cat(sprintf(" psi=%.3f, p=%.3f\n", mean(result$psi_hat), mean(result$p_hat, na.rm = TRUE)))
        else
          cat(" FAILED\n")
      }
      result
    })
  }

  species_fits <- list()
  all_betas_occ <- list()
  all_betas_det <- list()

  for (s in seq_len(n_sp)) {
    sp_name <- species_names[s]
    sp_fit <- species_fits_raw[[s]]

    species_fits[[sp_name]] <- sp_fit
    if (!is.null(sp_fit)) {
      all_betas_occ[[sp_name]] <- sp_fit$occ_fit$summary.fixed
      all_betas_det[[sp_name]] <- sp_fit$det_fit$summary.fixed
    }
  }

  # Community pooling: ensemble for robustness if requested
  pool_fn <- if (isTRUE(args$ensemble)) pool_community_ensemble
             else pool_community_effects

  list(
    species_fits  = species_fits,
    community_occ = pool_fn(all_betas_occ),
    community_det = pool_fn(all_betas_det),
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
# 12. engine_ms_temporal — parallel per-species + empirical Bayes shrinkage
#
# Architecture:
#   1. Fit each species independently via engine_temporal (parallelizable)
#   2. Pool betas across species → community mean & variance
#   3. Shrink each species' betas toward community mean (EB shrinkage)
#   4. Refit with shrunk priors (one more EM round per species)
#
# This is O(n_species / n_cores) instead of O(n_species^2) for stacked.
# ==========================================================================
#' @noRd
engine_ms_temporal <- function(args) {
  check_inla()
  data <- args$data

  y <- data$y
  if (!is.array(y) || length(dim(y)) != 4)
    stop("Temporal multi-species requires data$y as a 4D array (species x sites x seasons x visits)")

  n_sp <- dim(y)[1]
  N <- dim(y)[2]
  n_periods <- dim(y)[3]
  J <- dim(y)[4]
  sp_names <- dimnames(y)[[1]] %||% paste0("sp", seq_len(n_sp))
  use_ar1 <- args$temporal == "ar1"

  if (args$verbose >= 1) {
    cat("Fitting multi-species temporal model (parallel + EB shrinkage)\n")
    cat(sprintf("  Species: %d | Sites: %d | Periods: %d | Visits: %d | AR1: %s\n",
                n_sp, N, n_periods, J, use_ar1))
  }

  # --- Stage 1: Fit each species independently ---
  if (args$verbose >= 1) cat("\n  Stage 1: Per-species temporal fits\n")

  fit_one_species <- function(s) {
    sp_args <- args
    sp_args$data <- list(
      y        = y[s, , , ],
      occ.covs = data$occ.covs,
      det.covs = data$det.covs,
      coords   = data$coords
    )
    sp_args$verbose <- 0L

    tryCatch(
      engine_temporal(sp_args),
      error = function(e) {
        warning(sprintf("Species %s failed: %s", sp_names[s], e$message))
        NULL
      }
    )
  }

  # Parallel if enough species to amortize cluster startup cost
  n_cores <- getOption("INLAocc.cores", 1L)
  use_parallel <- n_cores > 1 &&
    n_sp >= 2L * n_cores &&
    requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    if (args$verbose >= 1)
      cat(sprintf("  Using %d cores\n", min(n_cores, n_sp)))
    cl <- parallel::makeCluster(min(n_cores, n_sp),
                                    port = sample(49152:65535, 1))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressMessages(library(INLAocc)))
    parallel::clusterExport(cl, c("args", "y", "data", "sp_names"),
                             envir = environment())
    species_fits <- parallel::parLapply(cl, seq_len(n_sp), fit_one_species)
  } else {
    species_fits <- lapply(seq_len(n_sp), function(s) {
      if (args$verbose >= 1)
        cat(sprintf("    [%d/%d] %s ... ", s, n_sp, sp_names[s]))
      t0 <- proc.time()
      result <- fit_one_species(s)
      if (args$verbose >= 1)
        cat(sprintf("%.1fs\n", (proc.time() - t0)["elapsed"]))
      result
    })
  }
  names(species_fits) <- sp_names

  # --- Stage 2: Empirical Bayes community shrinkage ---
  if (args$verbose >= 1) cat("\n  Stage 2: Community shrinkage\n")

  # Collect betas and their SEs from per-species fits
  occ_betas <- list()
  occ_ses   <- list()
  for (s in seq_len(n_sp)) {
    fit <- species_fits[[s]]
    if (!is.null(fit) && !is.null(fit$occ_fit)) {
      summ <- fit$occ_fit$summary.fixed
      occ_betas[[sp_names[s]]] <- summ$mean
      occ_ses[[sp_names[s]]]   <- summ$sd
    }
  }

  community <- NULL
  if (length(occ_betas) >= 2) {
    beta_mat <- do.call(rbind, occ_betas)  # species x covariates
    se_mat   <- do.call(rbind, occ_ses)
    n_coef   <- ncol(beta_mat)
    coef_names <- colnames(beta_mat) %||%
      rownames(species_fits[[1]]$occ_fit$summary.fixed)

    # Community parameters
    mu_comm    <- colMeans(beta_mat)
    tau2_comm  <- apply(beta_mat, 2, var)  # between-species variance

    # EB shrinkage: beta_s* = mu + (tau2 / (tau2 + sigma_s^2)) * (beta_s - mu)
    # Vectorized: broadcast tau2_comm across species rows
    sigma2_mat <- se_mat^2
    tau2_row <- matrix(tau2_comm, nrow = nrow(beta_mat), ncol = n_coef, byrow = TRUE)
    shrink_mat <- tau2_row / (tau2_row + sigma2_mat)
    shrink_mat[is.na(shrink_mat) | tau2_row < 1e-10] <- 0
    mu_row <- matrix(mu_comm, nrow = nrow(beta_mat), ncol = n_coef, byrow = TRUE)
    beta_shrunk <- mu_row + shrink_mat * (beta_mat - mu_row)

    community <- data.frame(
      parameter      = coef_names,
      community_mean = mu_comm,
      community_sd   = sqrt(tau2_comm),
      n_species      = nrow(beta_mat)
    )

    if (args$verbose >= 1) {
      cat("    Community means: [",
          paste(round(mu_comm, 3), collapse = ", "), "]\n")
      cat("    Between-species SD: [",
          paste(round(sqrt(tau2_comm), 3), collapse = ", "), "]\n")

      # Report shrinkage magnitude
      shrink_pct <- mean(abs(beta_shrunk - beta_mat) /
                           (abs(beta_mat) + 0.01)) * 100
      cat(sprintf("    Mean shrinkage: %.1f%%\n", shrink_pct))
    }

    # --- Stage 3: Apply shrinkage to psi estimates ---
    # Build design matrix once (invariant across species)
    occ_covs <- if (!is.null(data$occ.covs)) {
      oc <- lapply(data$occ.covs, function(x) {
        if (is.matrix(x) && ncol(x) == n_periods) x[, 1] else x
      })
      as.data.frame(oc)
    } else data.frame(.intercept = rep(1, N))
    X <- model.matrix(args$occ_parsed$fixed, data = occ_covs)

    for (s in seq_len(n_sp)) {
      sp <- sp_names[s]
      fit <- species_fits[[sp]]
      if (is.null(fit) || is.null(fit$psi_mat)) next

      if (sp %in% rownames(beta_shrunk)) {
        beta_diff <- beta_shrunk[sp, ] - beta_mat[sp, ]
        logit_shift <- as.vector(X %*% beta_diff)

        # Vectorized across periods: apply shift to psi_mat columns
        psi_mat_old <- clamp(fit$psi_mat)
        fit$psi_mat <- expit(logit(psi_mat_old) + logit_shift)
        for (t in seq_len(n_periods)) {
          fit$period_fits[[t]]$psi_hat <- fit$psi_mat[, t]
        }
      }
      species_fits[[sp]] <- fit
    }
  }

  # --- Build output ---
  # Reshape period_fits to match expected format: list[[species]][[period]]
  period_fits <- lapply(sp_names, function(sp) {
    fit <- species_fits[[sp]]
    if (is.null(fit)) return(NULL)
    lapply(seq_len(n_periods), function(t) fit$period_fits[[t]])
  })
  names(period_fits) <- sp_names

  # Build psi array
  psi_arr <- array(NA, dim = c(n_sp, N, n_periods))
  for (s in seq_len(n_sp)) {
    fit <- species_fits[[sp_names[s]]]
    if (!is.null(fit) && !is.null(fit$psi_mat)) {
      psi_arr[s, , ] <- fit$psi_mat
    }
  }

  list(
    species_fits  = species_fits,
    period_fits   = period_fits,
    psi_arr       = psi_arr,
    community     = community,
    n_species     = n_sp,
    species_names = sp_names,
    n_periods     = n_periods,
    data          = data,
    occ.formula   = args$occ.formula,
    det.formula   = args$det.formula,
    ar1           = use_ar1
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
