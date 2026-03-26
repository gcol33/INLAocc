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


# ---------------------------------------------------------------------------
# Prep functions: build once, reuse across EM iterations
# ---------------------------------------------------------------------------

#' @noRd
prep_detection <- function(data, det_formula, det_re_list = NULL) {
  det_df <- build_det_df(
    y        = data$y,
    det_covs = data$det.covs,
    site_idx = seq_len(data$N),
    weights  = rep(1, data$N)  # placeholder, updated per iteration
  )

  re_terms <- character(0)
  if (!is.null(det_re_list)) {
    re_comp <- build_re_components(
      det_re_list, det_df, data,
      prefix = "det", process = "det"
    )
    det_df   <- re_comp$df
    re_terms <- re_comp$formula_terms
  }

  base_formula <- paste("y_det ~", as.character(det_formula)[2])
  if (length(re_terms) > 0) {
    base_formula <- paste(base_formula, "+", paste(re_terms, collapse = " + "))
  }

  list(
    det_df   = det_df,
    formula  = as.formula(base_formula),
    Ntrials  = rep(1L, nrow(det_df))
  )
}


#' @noRd
prep_occupancy <- function(data, occ_formula, occ_re_list = NULL,
                           spatial = NULL) {
  use_spatial <- !is.null(spatial)
  spatial_type <- if (inherits(spatial, "occu_areal")) "areal"
                 else if (use_spatial) "spde"
                 else "none"

  # Build base occ_df.
  # Spatial models use weighted Bernoulli (M=1) to avoid overwhelming spatial priors.
  # Non-spatial models use M=1000 binomial encoding for better EM convergence.
  # The M=1000 attenuation bias is corrected by MI post-processing.
  if (use_spatial) {
    occ_df <- build_occ_df_weighted(
      occ_covs = data$occ.covs,
      weights  = rep(0.5, data$N),
      site_idx = seq_len(data$N),
      detected = data$detected
    )
  } else {
    occ_df <- build_occ_df(
      occ_covs = data$occ.covs,
      weights  = rep(0.5, data$N),
      site_idx = seq_len(data$N),
      detected = data$detected
    )
  }

  re_terms <- character(0)
  if (!is.null(occ_re_list)) {
    re_comp <- build_re_components(
      occ_re_list, occ_df, data,
      prefix = "occ", process = "occ"
    )
    occ_df   <- re_comp$df
    re_terms <- re_comp$formula_terms
  }

  base_formula <- paste("z ~", as.character(occ_formula)[2])
  if (length(re_terms) > 0) {
    base_formula <- paste(base_formula, "+", paste(re_terms, collapse = " + "))
  }

  # Spatial components (built once, data updated per iteration)
  stack_info <- NULL
  if (spatial_type == "spde") {
    base_formula <- paste(base_formula, "+ f(spatial, model = spde)")

    drop_cols <- c("z", "spatial", "Ntrials", "w", "site")
    cov_cols <- setdiff(names(occ_df), drop_cols)

    stack_info <- list(
      cov_cols  = cov_cols,
      n_mesh    = spatial$n_mesh,
      spde      = spatial$spde,
      A         = spatial$A
    )
  } else if (spatial_type == "areal") {
    # CAR/BYM2: region index directly, no mesh/stack needed
    occ_df$region <- occ_df$site  # region = site index for areal models
    areal_term <- sprintf(
      'f(region, model = "%s", graph = areal_graph, scale.model = %s, hyper = list(prec = list(prior = "%s", param = c(%s))))',
      spatial$model, as.character(spatial$scale.model),
      spatial$hyper$prec$prior,
      paste(spatial$hyper$prec$param, collapse = ", ")
    )
    base_formula <- paste(base_formula, "+", areal_term)
  }

  list(
    occ_df       = occ_df,
    formula      = as.formula(base_formula),
    spatial      = use_spatial,
    spatial_type = spatial_type,
    stack_info   = stack_info
  )
}


# ---------------------------------------------------------------------------
# Fit functions: fast per-iteration INLA calls
# ---------------------------------------------------------------------------

#' @noRd
fit_detection_inla <- function(det_prep, data, weights,
                               control.inla = NULL,
                               control.compute = NULL,
                               num.threads = "1:1",
                               verbose = FALSE) {
  det_df <- det_prep$det_df
  # Update weights in-place (site -> long format)
  det_df$w <- weights[det_df$site]

  if (is.null(control.inla)) {
    control.inla <- INLA::inla.set.control.inla.default()
  }
  if (is.null(control.compute)) {
    control.compute <- list(config = FALSE, dic = FALSE, waic = FALSE)
  }

  det_df$.Ntrials <- det_prep$Ntrials

  fit <- INLA::inla(
    formula         = det_prep$formula,
    family          = "binomial",
    Ntrials         = det_df$.Ntrials,
    data            = det_df,
    weights         = det_df$w,
    num.threads     = num.threads,
    control.inla    = control.inla,
    control.compute = control.compute,
    verbose         = verbose
  )

  p_fitted <- fit$summary.fitted.values$mean
  p_hat <- matrix(NA, data$N, data$J)
  p_hat[cbind(det_df$site, det_df$visit)] <- p_fitted

  list(fit = fit, p_hat = p_hat, det_df = det_df)
}


#' @noRd
fit_occupancy_inla <- function(occ_prep, data, weights, spatial_obj = NULL,
                               control.inla = NULL,
                               control.compute = NULL,
                               num.threads = "1:1",
                               verbose = FALSE) {
  occ_df <- occ_prep$occ_df
  use_spatial <- occ_prep$spatial

  # Update weights in the data frame
  if (use_spatial) {
    # Weighted Bernoulli for spatial models
    N <- data$N
    n_undet <- sum(!data$detected)
    w_clamped <- clamp(weights, 1e-6, 1 - 1e-6)
    occ_df$z[seq_len(N)]  <- 1L
    occ_df$w[seq_len(N)]  <- ifelse(data$detected, 1.0, w_clamped)
    if (n_undet > 0) {
      extra_idx <- (N + 1L):(N + n_undet)
      occ_df$z[extra_idx] <- 0L
      occ_df$w[extra_idx] <- 1 - w_clamped[!data$detected]
    }
  } else {
    # M=1000 binomial for non-spatial (better EM convergence, MI debiases after)
    M <- occ_df$Ntrials[1]
    occ_df$z <- ifelse(data$detected, M, round(weights * M))
  }

  if (is.null(control.inla)) {
    control.inla <- INLA::inla.set.control.inla.default()
  }
  if (is.null(control.compute)) {
    control.compute <- list(config = FALSE, dic = FALSE, waic = FALSE)
  }

  spatial_type <- occ_prep$spatial_type %||% (if (use_spatial) "spde" else "none")

  if (spatial_type == "spde") {
    si <- occ_prep$stack_info
    A_obs <- si$A[occ_df$site, , drop = FALSE]
    cov_df <- occ_df[, si$cov_cols, drop = FALSE]

    stack <- INLA::inla.stack(
      data    = list(z = occ_df$z, Ntrials = occ_df$Ntrials,
                     w_obs = occ_df$w),
      A       = list(A_obs, 1),
      effects = list(
        list(spatial = seq_len(si$n_mesh)),
        cov_df
      ),
      tag     = "occ"
    )

    stack_data <- INLA::inla.stack.data(stack)
    stack_data$spde <- si$spde
    n_stack <- nrow(stack_data)

    fit <- suppressWarnings(INLA::inla(
      formula         = occ_prep$formula,
      family          = "binomial",
      Ntrials         = rep(1L, n_stack),
      data            = stack_data,
      weights         = stack_data$w_obs,
      num.threads     = num.threads,
      control.predictor = list(
        A       = INLA::inla.stack.A(stack),
        compute = TRUE
      ),
      control.inla    = control.inla,
      control.compute = control.compute,
      verbose         = verbose
    ))

    idx <- INLA::inla.stack.index(stack, "occ")$data
    psi_all <- fit$summary.fitted.values$mean[idx]
    first_per_site <- !duplicated(occ_df$site)
    psi_hat <- psi_all[first_per_site]
  } else if (spatial_type == "areal") {
    # CAR/BYM2: no stack needed, region index is direct
    occ_df$areal_graph <- spatial_obj$graph
    inla_data <- occ_df
    inla_data$areal_graph <- spatial_obj$graph

    fit <- suppressWarnings(INLA::inla(
      formula         = occ_prep$formula,
      family          = "binomial",
      Ntrials         = rep(1L, nrow(occ_df)),
      data            = inla_data,
      weights         = occ_df$w,
      num.threads     = num.threads,
      control.inla    = control.inla,
      control.compute = control.compute,
      verbose         = verbose
    ))

    psi_all <- fit$summary.fitted.values$mean
    first_per_site <- !duplicated(occ_df$site)
    psi_hat <- psi_all[first_per_site]
  } else {
    fit <- INLA::inla(
      formula         = occ_prep$formula,
      family          = "binomial",
      Ntrials         = occ_df$Ntrials,
      data            = occ_df,
      num.threads     = num.threads,
      control.inla    = control.inla,
      control.compute = control.compute,
      verbose         = verbose
    )
    psi_hat <- fit$summary.fitted.values$mean[seq_len(data$N)]
  }

  list(fit = fit, psi_hat = psi_hat, occ_df = occ_df)
}


# ---------------------------------------------------------------------------
# Main EM loop
# ---------------------------------------------------------------------------

#' @noRd
em_inla <- function(data, occ_formula, det_formula,
                    occ_re = NULL, det_re = NULL,
                    spatial = NULL, priors = NULL,
                    max_iter = 50, tol = 1e-4, damping = 0.3,
                    num.threads = "1:1",
                    control.inla = NULL,
                    verbose = 1) {

  check_inla()
  if (!inherits(data, "occu_data")) {
    stop("data must be an occu_data object (from occu_format())")
  }

  N <- data$N
  J <- data$J

  # --- Prep: build data structures once ---
  det_prep <- prep_detection(data, det_formula, det_re)
  occ_prep <- prep_occupancy(data, occ_formula, occ_re, spatial)

  # INLA control: skip DIC/WAIC for intermediate iters, full for final
  ctrl_fast <- list(config = TRUE, dic = FALSE, waic = FALSE)
  ctrl_full <- list(config = TRUE,  dic = TRUE,  waic = TRUE)

  # --- Initialize with GLM warm start ---
  # A fast GLM on collapsed data (ever_detected ~ covariates) gives much
  # better starting psi than naive rates, reducing EM iterations by 30-50%.
  ever_detected <- rowSums(data$y == 1, na.rm = TRUE) > 0
  psi_hat <- tryCatch({
    glm_df <- data$occ.covs
    glm_df$.y <- as.integer(ever_detected)
    glm_fit <- glm(.y ~ ., data = glm_df, family = binomial)
    clamp(predict(glm_fit, type = "response"), 0.01, 0.99)
  }, error = function(e) {
    rep(clamp(data$naive_occ, 0.1, 0.9), N)
  })

  p_hat <- matrix(clamp(data$naive_det, 0.1, 0.9), N, J)
  if (sum(data$detected) > 5) {
    init_p <- sum(data$y[data$detected, ], na.rm = TRUE) /
      sum(!is.na(data$y[data$detected, ]))
    p_hat[] <- clamp(init_p, 0.05, 0.95)
  }

  weights <- compute_weights(data$y, psi_hat, p_hat)

  history <- list()
  converged <- FALSE
  det_result <- NULL
  occ_result <- NULL
  det_frozen <- FALSE
  delta_p <- Inf

  # --- Adaptive damping ---
  # damping = how much of OLD weights to keep (high = conservative, low = aggressive)
  # Start aggressive (low damping), back off when oscillating.
  # User-supplied damping sets the floor; we start at max(damping, 0.1).
  damp_floor <- damping
  damp_current <- max(damping, 0.1)  # start aggressive
  prev_delta_psi <- Inf

  for (iter in seq_len(max_iter)) {
    psi_old <- psi_hat
    p_old   <- p_hat

    # --- M-step: Detection ---
    # Skip when detection has stabilized (delta_p < tol for 2 consecutive iters)
    if (!det_frozen) {
      det_result <- fit_detection_inla(
        det_prep, data, weights,
        control.inla = control.inla,
        control.compute = ctrl_fast,
        num.threads = num.threads,
        verbose = verbose >= 2
      )
      p_hat <- det_result$p_hat
      p_hat[is.na(p_hat)] <- p_old[is.na(p_hat)]
    }

    # --- M-step: Occupancy ---
    occ_result <- fit_occupancy_inla(
      occ_prep, data, weights, spatial_obj = spatial,
      control.inla = control.inla,
      control.compute = ctrl_fast,
      num.threads = num.threads,
      verbose = verbose >= 2
    )
    psi_hat <- occ_result$psi_hat

    # --- E-step: Update weights ---
    weights_new <- compute_weights(data$y, psi_hat, p_hat)

    if (damp_current > 0) {
      weights <- damp_current * weights + (1 - damp_current) * weights_new
    } else {
      weights <- weights_new
    }

    # --- Convergence check ---
    delta_psi <- max(abs(psi_hat - psi_old))
    delta_p_new <- max(abs(p_hat - p_old), na.rm = TRUE)

    # --- Adapt damping ---
    # Oscillation = delta_psi increased from previous iteration
    if (iter >= 2) {
      if (delta_psi > prev_delta_psi * 1.05) {
        # Oscillating: increase damping (more conservative)
        damp_current <- min(0.9, damp_current + 0.15)
      } else if (delta_psi < prev_delta_psi * 0.5) {
        # Converging fast: decrease damping (more aggressive)
        damp_current <- max(damp_floor, damp_current - 0.1)
      }
    }
    prev_delta_psi <- delta_psi

    # Freeze detection after 2 consecutive sub-tol iterations
    if (delta_p < tol && delta_p_new < tol && iter >= 2) {
      det_frozen <- TRUE
    }
    delta_p <- delta_p_new

    history[[iter]] <- list(
      iter      = iter,
      delta_psi = delta_psi,
      delta_p   = delta_p,
      damping   = damp_current,
      mean_psi  = mean(psi_hat),
      mean_p    = mean(p_hat, na.rm = TRUE)
    )

    if (verbose >= 1) {
      ll <- occu_loglik(data$y, psi_hat, p_hat)
      history[[iter]]$loglik <- ll
      cat(sprintf(
        "  EM iter %2d | delta_psi = %.6f | delta_p = %.6f | damp = %.2f | loglik = %.2f%s\n",
        iter, delta_psi, delta_p, damp_current, ll,
        if (det_frozen) " [det frozen]" else ""
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

  # --- Final full-quality INLA fits (with DIC/WAIC/config) ---
  det_result <- fit_detection_inla(
    det_prep, data, weights,
    control.inla = control.inla,
    control.compute = ctrl_full,
    num.threads = num.threads,
    verbose = verbose >= 2
  )
  occ_result <- fit_occupancy_inla(
    occ_prep, data, weights, spatial_obj = spatial,
    control.inla = control.inla,
    control.compute = ctrl_full,
    num.threads = num.threads,
    verbose = verbose >= 2
  )
  psi_hat <- occ_result$psi_hat
  p_hat   <- det_result$p_hat
  p_hat[is.na(p_hat)] <- p_old[is.na(p_hat)]

  z_hat <- compute_weights(data$y, psi_hat, p_hat)

  # --- Multiple imputation to debias occupancy betas ---
  # The EM's soft weights attenuate betas toward zero. MI samples hard 0/1
  # from the converged weights and refits INLA on proper binary data.
  mi_result <- mi_occupancy(
    occ_prep, data, z_hat, spatial,
    control.inla = control.inla, K = 10L, verbose = verbose
  )
  if (!is.null(mi_result)) {
    occ_result$fit$summary.fixed <- mi_result$summary.fixed
    psi_hat <- mi_result$psi_hat
  }

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
    det_df      = det_result$det_df,
    mi          = mi_result
  )
  class(out) <- "occu_em"
  out
}


# ---------------------------------------------------------------------------
# Multiple imputation for occupancy beta debiasing
# ---------------------------------------------------------------------------
#' @noRd
mi_occupancy <- function(occ_prep, data, z_hat, spatial,
                          control.inla = NULL, K = 20L,
                          K_min = 5L, mi_tol = 0.05,
                          verbose = 1) {
  check_inla()
  N <- data$N
  detected <- data$detected
  undetected_idx <- which(!detected)

  if (length(undetected_idx) == 0) return(NULL)
  if (verbose >= 1) cat(sprintf("  MI debiasing (up to %d imputations)...\n", K))

  use_spatial <- occ_prep$spatial
  ctrl <- list(config = FALSE, dic = FALSE, waic = FALSE)

  # Build a clean formula for MI (no weighted Bernoulli extra rows)
  # Use the base formula without spatial terms — add spatial via stack
  cov_names <- setdiff(names(data$occ.covs), c("z", "site", "Ntrials", "w"))
  if (length(cov_names) > 0) {
    mi_formula_str <- paste("z ~", paste(cov_names, collapse = " + "))
  } else {
    mi_formula_str <- "z ~ 1"
  }
  if (use_spatial) {
    mi_formula_str <- paste(mi_formula_str, "+ f(spatial, model = spde)")
  }
  mi_formula <- as.formula(mi_formula_str)

  # --- Single imputation function ---
  run_one_mi <- function(seed) {
    set.seed(seed)
    z_imp <- as.integer(detected)
    z_imp[undetected_idx] <- rbinom(
      length(undetected_idx), 1, clamp(z_hat[undetected_idx])
    )

    occ_df <- data$occ.covs
    occ_df$site <- seq_len(N)
    occ_df$z <- z_imp
    nt <- rep(1L, N)

    if (use_spatial) {
      si <- occ_prep$stack_info
      A_obs <- si$A[seq_len(N), , drop = FALSE]
      drop_cols <- c("z", "site")
      cov_df <- occ_df[, setdiff(names(occ_df), drop_cols), drop = FALSE]

      stack <- INLA::inla.stack(
        data = list(z = z_imp),
        A = list(A_obs, 1),
        effects = list(
          list(spatial = seq_len(si$n_mesh)),
          cov_df
        ),
        tag = "occ"
      )
      stack_data <- INLA::inla.stack.data(stack)
      stack_data$spde <- si$spde
      n_stack <- nrow(stack_data)

      fit <- suppressWarnings(INLA::inla(
        formula = mi_formula,
        family = "binomial",
        Ntrials = rep(1L, n_stack),
        data = stack_data,
        control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
        num.threads = "1:1",
        control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
        control.compute = ctrl,
        verbose = FALSE
      ))
      idx <- INLA::inla.stack.index(stack, "occ")$data
      psi_k <- fit$summary.fitted.values$mean[idx]
    } else {
      fit <- suppressWarnings(INLA::inla(
        formula = mi_formula,
        family = "binomial",
        Ntrials = nt,
        data = occ_df,
        num.threads = "1:1",
        control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
        control.compute = ctrl,
        verbose = FALSE
      ))
      psi_k <- fit$summary.fitted.values$mean[seq_len(N)]
    }

    list(beta = fit$summary.fixed, psi = psi_k)
  }

  # --- Run K imputations (adaptive stopping) ---
  seeds <- sample.int(1e7, K)
  beta_samples <- list()
  psi_samples  <- matrix(NA, N, K)
  prev_V_between <- NULL
  K_used <- K

  # Running accumulators for incremental variance (Welford's algorithm)
  beta_sum <- NULL
  beta_sum_sq <- NULL

  for (k in seq_len(K)) {
    res <- run_one_mi(seeds[k])
    beta_samples[[k]] <- res$beta
    psi_samples[, k]  <- res$psi

    # Adaptive stopping after K_min imputations
    if (k >= K_min) {
      # Incremental mean/variance via running sums (O(K) total, not O(K²))
      bm <- res$beta$mean
      if (is.null(beta_sum)) {
        # First time reaching K_min: initialize from all samples so far
        all_means <- sapply(beta_samples[seq_len(k)], function(b) b$mean)
        if (!is.matrix(all_means)) all_means <- matrix(all_means, nrow = 1)
        beta_sum <- rowSums(all_means)
        beta_sum_sq <- rowSums(all_means^2)
      } else {
        beta_sum <- beta_sum + bm
        beta_sum_sq <- beta_sum_sq + bm^2
      }
      V_between <- (beta_sum_sq / k) - (beta_sum / k)^2

      if (!is.null(prev_V_between)) {
        denom <- pmax(prev_V_between, 1e-10)
        rel_change <- max(abs(V_between - prev_V_between) / denom)
        if (rel_change < mi_tol) {
          K_used <- k
          if (verbose >= 1)
            cat(sprintf("    MI converged at K=%d (rel_change=%.4f < %.2f)\n",
                        k, rel_change, mi_tol))
          break
        }
      }
      prev_V_between <- V_between
    }
  }

  # --- Pool with Rubin's rules ---
  K_actual <- K_used
  beta_samples <- beta_samples[seq_len(K_actual)]
  psi_samples  <- psi_samples[, seq_len(K_actual), drop = FALSE]

  n_coef <- nrow(beta_samples[[1]])
  coef_names <- rownames(beta_samples[[1]])

  beta_means <- sapply(beta_samples, function(b) b$mean)  # coef x K
  beta_vars  <- sapply(beta_samples, function(b) b$sd^2)
  if (!is.matrix(beta_means)) {
    beta_means <- matrix(beta_means, nrow = 1)
    beta_vars  <- matrix(beta_vars, nrow = 1)
  }

  pooled_mean <- rowMeans(beta_means)
  V_within    <- rowMeans(beta_vars)
  V_between   <- apply(beta_means, 1, var)
  V_total     <- V_within + (1 + 1/K_actual) * V_between
  pooled_sd   <- sqrt(V_total)

  summary_fixed <- data.frame(
    mean    = pooled_mean,
    sd      = pooled_sd,
    `0.025quant` = pooled_mean - 1.96 * pooled_sd,
    `0.5quant`   = pooled_mean,
    `0.975quant` = pooled_mean + 1.96 * pooled_sd,
    row.names = coef_names,
    check.names = FALSE
  )

  psi_hat <- rowMeans(psi_samples)

  if (verbose >= 1) {
    cat(sprintf("    MI betas (%d imputations): [%s]\n", K_actual,
                paste(round(pooled_mean, 3), collapse = ", ")))
    cat(sprintf("    MI psi range: [%.3f, %.3f] mean=%.3f\n",
                min(psi_hat), max(psi_hat), mean(psi_hat)))
  }

  list(
    summary.fixed = summary_fixed,
    psi_hat       = psi_hat,
    beta_samples  = beta_means,
    K             = K_actual
  )
}
