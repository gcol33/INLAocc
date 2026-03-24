# =============================================================================
# occu_engine.R — EM-INLA engine for occupancy models
# =============================================================================
# Core algorithm: Expectation-Maximization with INLA at each M-step.
#
# E-step:  Compute P(z_i = 1 | y, theta) for undetected sites
# M-step:  Fit detection and occupancy sub-models via INLA
#
# This approach:
#   - Handles arbitrary covariates on both processes
#   - Supports random intercepts and slopes via INLA's f()
#   - Can incorporate spatial effects via SPDE
#   - Is computationally efficient (INLA is fast per iteration)
# =============================================================================

#' @noRd
fit_detection_inla <- function(data, det_formula, weights,
                               det_re_list = NULL,
                               control.inla = NULL,
                               verbose = FALSE) {
  check_inla()

  # Build long-format detection data
  det_df <- build_det_df(
    y        = data$y,
    det_covs = data$det.covs,
    site_idx = seq_len(data$N),
    weights  = weights
  )

  # Add random effect columns
  re_terms <- character(0)
  if (!is.null(det_re_list)) {
    re_comp <- build_re_components(
      det_re_list, det_df, data,
      prefix = "det", process = "det"
    )
    det_df   <- re_comp$df
    re_terms <- re_comp$formula_terms
  }

  # Build formula
  base_formula <- paste("y_det ~", as.character(det_formula)[2])
  if (length(re_terms) > 0) {
    base_formula <- paste(base_formula, "+", paste(re_terms, collapse = " + "))
  }
  full_formula <- as.formula(base_formula)

  # Default INLA control
  if (is.null(control.inla)) {
    control.inla <- INLA::inla.set.control.inla.default()
  }

  # Fit
  fit <- INLA::inla(
    formula       = full_formula,
    family        = "binomial",
    Ntrials       = rep(1, nrow(det_df)),
    data          = det_df,
    weights       = det_df$w,
    control.inla  = control.inla,
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
    verbose       = verbose
  )

  # Extract fitted detection probabilities as N x J matrix
  p_fitted <- fit$summary.fitted.values$mean
  p_hat <- matrix(NA, data$N, data$J)

  # Map back from long format to matrix
  for (r in seq_len(nrow(det_df))) {
    i <- det_df$site[r]
    j <- det_df$visit[r]
    p_hat[i, j] <- p_fitted[r]
  }

  list(fit = fit, p_hat = p_hat, det_df = det_df)
}


#' @noRd
fit_occupancy_inla <- function(data, occ_formula, weights,
                               occ_re_list = NULL,
                               spatial = NULL,
                               control.inla = NULL,
                               verbose = FALSE) {
  check_inla()

  # Build augmented occupancy data
  occ_df <- build_occ_df(
    occ_covs = data$occ.covs,
    weights  = weights,
    site_idx = seq_len(data$N),
    detected = data$detected
  )

  # Add random effect columns
  re_terms <- character(0)
  if (!is.null(occ_re_list)) {
    re_comp <- build_re_components(
      occ_re_list, occ_df, data,
      prefix = "occ", process = "occ"
    )
    occ_df   <- re_comp$df
    re_terms <- re_comp$formula_terms
  }

  # Build formula
  base_formula <- paste("z ~", as.character(occ_formula)[2])
  if (length(re_terms) > 0) {
    base_formula <- paste(base_formula, "+", paste(re_terms, collapse = " + "))
  }

  # Spatial component
  use_spatial <- !is.null(spatial)
  if (use_spatial) {
    base_formula <- paste(base_formula,
                          '+ f(spatial, model = spde)')
    # Add spatial index to data
    occ_df$spatial <- occ_df$site
  }

  full_formula <- as.formula(base_formula)

  if (is.null(control.inla)) {
    control.inla <- INLA::inla.set.control.inla.default()
  }

  # Fit: with or without spatial stack
  if (use_spatial) {
    # Use INLA stack for spatial model
    A_obs <- spatial$A[occ_df$site, , drop = FALSE]

    stack <- INLA::inla.stack(
      data    = list(z = occ_df$z),
      A       = list(A_obs, 1),
      effects = list(
        list(spatial = seq_len(spatial$n_mesh)),
        occ_df[, setdiff(names(occ_df), c("z", "spatial")), drop = FALSE]
      ),
      tag     = "occ"
    )

    base_formula_sp <- paste("z ~", as.character(occ_formula)[2])
    if (length(re_terms) > 0) {
      base_formula_sp <- paste(base_formula_sp, "+",
                               paste(re_terms, collapse = " + "))
    }
    base_formula_sp <- paste(base_formula_sp,
                             "+ f(spatial, model = spde)")

    fit <- INLA::inla(
      formula         = as.formula(base_formula_sp),
      family          = "binomial",
      Ntrials         = occ_df$Ntrials,
      data            = INLA::inla.stack.data(stack),
      control.predictor = list(
        A       = INLA::inla.stack.A(stack),
        compute = TRUE
      ),
      control.inla    = control.inla,
      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
      verbose         = verbose
    )
  } else {
    fit <- INLA::inla(
      formula       = full_formula,
      family        = "binomial",
      Ntrials       = occ_df$Ntrials,
      data          = occ_df,
      control.inla  = control.inla,
      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
      verbose       = verbose
    )
  }

  # Extract fitted occupancy probabilities per site (one row per site now)
  psi_fitted <- fit$summary.fitted.values$mean
  psi_hat <- psi_fitted[seq_len(data$N)]

  list(fit = fit, psi_hat = psi_hat, occ_df = occ_df)
}


#' @noRd
em_inla <- function(data, occ_formula, det_formula,
                    occ_re = NULL, det_re = NULL,
                    spatial = NULL, priors = NULL,
                    max_iter = 50, tol = 1e-4, damping = 0.3,
                    control.inla = NULL,
                    verbose = 1) {

  check_inla()
  if (!inherits(data, "occu_data")) {
    stop("data must be an occu_data object (from occu_format())")
  }

  N <- data$N
  J <- data$J

  # --- Initialize ---
  psi_hat <- rep(clamp(data$naive_occ, 0.1, 0.9), N)
  p_hat   <- matrix(clamp(data$naive_det, 0.1, 0.9), N, J)

  # Better initialization: use detected sites for initial p
  if (sum(data$detected) > 5) {
    init_p <- sum(data$y[data$detected, ], na.rm = TRUE) /
      sum(!is.na(data$y[data$detected, ]))
    p_hat[] <- clamp(init_p, 0.05, 0.95)
  }

  weights <- compute_weights(data$y, psi_hat, p_hat)

  history <- list()
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    psi_old <- psi_hat
    p_old   <- p_hat

    # --- M-step: Detection ---
    det_result <- fit_detection_inla(
      data, det_formula, weights,
      det_re_list = det_re,
      control.inla = control.inla,
      verbose = verbose >= 2
    )
    p_hat <- det_result$p_hat
    # Fill NAs with previous values
    p_hat[is.na(p_hat)] <- p_old[is.na(p_hat)]

    # --- M-step: Occupancy ---
    occ_result <- fit_occupancy_inla(
      data, occ_formula, weights,
      occ_re_list = occ_re,
      spatial = spatial,
      control.inla = control.inla,
      verbose = verbose >= 2
    )
    psi_hat <- occ_result$psi_hat

    # --- E-step: Update weights ---
    weights_new <- compute_weights(data$y, psi_hat, p_hat)

    # Damping
    if (damping > 0) {
      weights <- damping * weights + (1 - damping) * weights_new
    } else {
      weights <- weights_new
    }

    # --- Convergence check ---
    delta_psi <- max(abs(psi_hat - psi_old))
    delta_p   <- max(abs(p_hat - p_old), na.rm = TRUE)
    ll        <- occu_loglik(data$y, psi_hat, p_hat)

    history[[iter]] <- list(
      iter      = iter,
      delta_psi = delta_psi,
      delta_p   = delta_p,
      loglik    = ll,
      mean_psi  = mean(psi_hat),
      mean_p    = mean(p_hat, na.rm = TRUE)
    )

    if (verbose >= 1) {
      cat(sprintf(
        "  EM iter %2d | delta_psi = %.6f | delta_p = %.6f | loglik = %.2f\n",
        iter, delta_psi, delta_p, ll
      ))
    }

    if (delta_psi < tol && delta_p < tol) {
      converged <- TRUE
      if (verbose >= 1) cat("  Converged.\n")
      break
    }
  }

  if (!converged && verbose >= 1) {
    cat(sprintf("  Warning: EM did not converge after %d iterations.\n", max_iter))
  }

  # Final E-step weights
  z_hat <- compute_weights(data$y, psi_hat, p_hat)

  out <- list(
    occ_fit     = occ_result$fit,
    det_fit     = det_result$fit,
    psi_hat     = psi_hat,
    p_hat       = p_hat,
    z_hat       = z_hat,
    weights     = weights,
    data        = data,
    occ_formula = occ_formula,
    det_formula = det_formula,
    occ_re      = occ_re,
    det_re      = det_re,
    spatial     = spatial,
    converged   = converged,
    n_iter      = length(history),
    history     = history,
    occ_df      = occ_result$occ_df,
    det_df      = det_result$det_df
  )
  class(out) <- "occu_em"
  out
}
