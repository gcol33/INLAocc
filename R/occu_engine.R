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
# Fit functions: GLM fast path for EM iterations (no random effects / spatial)
# ---------------------------------------------------------------------------

#' @noRd
fit_detection_glm <- function(det_prep, data, weights) {
  det_df <- det_prep$det_df
  det_df$w <- weights[det_df$site]

  fit <- glm(det_prep$formula, family = binomial, data = det_df,
             weights = det_df$w, control = glm.control(maxit = 50))

  p_fitted <- predict(fit, type = "response")
  p_hat <- matrix(NA, data$N, data$J)
  p_hat[cbind(det_df$site, det_df$visit)] <- p_fitted

  list(fit = fit, p_hat = p_hat, det_df = det_df)
}

#' @noRd
fit_occupancy_glm <- function(occ_prep, data, weights) {
  occ_df <- occ_prep$occ_df
  M <- occ_df$Ntrials[1]
  occ_df$z <- ifelse(data$detected, M, round(weights * M))
  occ_df$.failures <- M - occ_df$z

  # glm needs cbind(success, failure) for binomial with Ntrials
  rhs <- as.character(occ_prep$formula)[3]
  glm_formula <- as.formula(paste("cbind(z, .failures) ~", rhs))

  fit <- suppressWarnings(glm(glm_formula, family = binomial, data = occ_df,
             control = glm.control(maxit = 50)))

  psi_hat <- predict(fit, type = "response")[seq_len(data$N)]
  list(fit = fit, psi_hat = psi_hat, occ_df = occ_df)
}


# ---------------------------------------------------------------------------
# Fit functions: INLA calls (for final fit, Gibbs, and RE/spatial models)
# ---------------------------------------------------------------------------

#' @noRd
fit_detection_inla <- function(det_prep, data, weights,
                               control.inla = NULL,
                               control.compute = NULL,
                               control.mode = NULL,
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
    control.mode    = control.mode,
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
                               control.mode = NULL,
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
      control.mode    = control.mode,
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
      control.mode    = control.mode,
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
      control.mode    = control.mode,
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
                    correction = "auto",
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

  # GLM fast path: use glm() for EM iterations when the model is simple
  # (no random effects, no spatial). INLA is only needed for the final fit
  # and the Gibbs loop. This avoids INLA's ~0.3s/call process spawn overhead.
  use_glm_em <- is.null(occ_re) && is.null(det_re) && is.null(spatial)

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

  # Warm-start: pass previous mode to INLA to skip Newton iterations.
  # Structure is constant across EM iterations (same covariates, same N),
  # so both theta and x transfer for occupancy; detection has constant
  # structure too (same long-format data frame, only weights change).
  prev_occ_mode <- NULL
  prev_det_mode <- NULL

  extract_mode <- function(fit) {
    if (is.null(fit)) return(NULL)
    list(theta = fit$mode$theta, x = fit$mode$x)
  }

  mode_to_ctrl <- function(m) {
    if (is.null(m)) return(NULL)
    has_theta <- length(m$theta) > 0
    cm <- list(restart = has_theta)
    if (has_theta) cm$theta <- m$theta
    if (length(m$x) > 0) cm$x <- m$x
    cm
  }

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
      if (use_glm_em) {
        det_result <- fit_detection_glm(det_prep, data, weights)
      } else {
        det_result <- fit_detection_inla(
          det_prep, data, weights,
          control.inla = control.inla,
          control.compute = ctrl_fast,
          control.mode = mode_to_ctrl(prev_det_mode),
          num.threads = num.threads,
          verbose = verbose >= 2
        )
        prev_det_mode <- extract_mode(det_result$fit)
      }
      p_hat <- det_result$p_hat
      p_hat[is.na(p_hat)] <- p_old[is.na(p_hat)]
    }

    # --- M-step: Occupancy ---
    if (use_glm_em) {
      occ_result <- fit_occupancy_glm(occ_prep, data, weights)
    } else {
      occ_result <- fit_occupancy_inla(
        occ_prep, data, weights, spatial_obj = spatial,
        control.inla = control.inla,
        control.compute = ctrl_fast,
        control.mode = mode_to_ctrl(prev_occ_mode),
        num.threads = num.threads,
        verbose = verbose >= 2
      )
      prev_occ_mode <- extract_mode(occ_result$fit)
    }
    psi_hat <- occ_result$psi_hat

    # --- E-step: Update weights ---
    weights_new <- compute_weights(data$y, psi_hat, p_hat)

    weights <- damp_current * weights + (1 - damp_current) * weights_new

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

  # --- Post-EM debiasing ---
  # The EM's soft weights attenuate coefficients. Two correction methods:
  #   "mi"    — independent draws of z, refit both submodels, pool (Rubin 1987)
  #   "gibbs" — sequential Gibbs chain: z|params → occ|z → det|z (Tanner & Wong 1987)
  #   "auto"  — Gibbs for now; will benchmark MI vs Gibbs to decide
  #   "none"  — skip correction entirely
  # "auto": MI-glm for simple models (fast, no INLA overhead in correction),
  # Gibbs for spatial/RE models (needs INLA at each step for the spatial field).
  has_re_or_spatial <- !is.null(spatial) || !is.null(occ_re) || !is.null(det_re)
  correction_method <- if (correction == "auto") {
    if (has_re_or_spatial) "gibbs" else "mi"
  } else {
    correction
  }

  mi_result <- if (correction_method == "none") {
    NULL
  } else if (correction_method == "mi") {
    mi_joint(occ_prep, det_prep, data, z_hat, spatial,
             control.inla = control.inla,
             n_cores = getOption("INLAocc.mi.cores", 1L),
             verbose = verbose)
  } else {
    mi_data_augmentation(occ_prep, det_prep, data, z_hat, spatial,
                         control.inla = control.inla, verbose = verbose)
  }
  mi_det_result <- NULL
  if (!is.null(mi_result)) {
    occ_result$fit$summary.fixed <- mi_result$occ_summary
    det_result$fit$summary.fixed <- mi_result$det_summary
    psi_hat <- mi_result$psi_hat
    p_hat   <- mi_result$p_hat
    z_hat   <- compute_weights(data$y, psi_hat, p_hat)
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
# Joint MI: independent draws for both occ and det debiasing
# ---------------------------------------------------------------------------
#' Each imputation independently draws z from converged weights, then refits
#' both submodels conditional on it. Draws are independent (not a chain),
#' Within each round, draws are independent (Rubin's rules exact).
#' Across rounds, z_hat is updated from pooled estimates, fixing the bias
#' that comes from drawing z from the EM's attenuated weights.
#' @noRd
mi_joint <- function(occ_prep, det_prep, data, z_hat, spatial,
                     control.inla = NULL, K = 20L,
                     K_min = 5L, mi_tol = 0.05,
                     n_rounds = 3L,
                     n_cores = 1L,
                     verbose = 1) {
  check_inla()
  N <- data$N
  J <- data$J
  detected <- data$detected
  undetected_idx <- which(!detected)

  if (length(undetected_idx) == 0) return(NULL)
  if (verbose >= 1)
    cat(sprintf("  MI joint debiasing (%d rounds x %d imputations)...\n", n_rounds, K))

  use_spatial <- occ_prep$spatial
  ctrl <- list(config = FALSE, dic = FALSE, waic = FALSE)

  # Build formulas once
  cov_names <- setdiff(names(data$occ.covs), c("z", "site", "Ntrials", "w"))
  occ_fstr <- if (length(cov_names) > 0) paste("z ~", paste(cov_names, collapse = " + ")) else "z ~ 1"
  if (use_spatial) occ_fstr <- paste(occ_fstr, "+ f(spatial, model = spde)")
  occ_formula <- as.formula(occ_fstr)
  det_fstr <- paste("y_det ~", as.character(det_prep$formula)[3])
  det_formula <- as.formula(det_fstr, env = globalenv())

  # Warm-start from EM mode (same starting point for all imputations —
  # valid because draws are independent; we just reuse the optimization hint)
  em_occ_mode <- NULL
  em_det_mode <- NULL

  build_ctrl_mode <- function(prev_mode, transfer_x = TRUE) {
    if (is.null(prev_mode)) return(NULL)
    has_theta <- length(prev_mode$theta) > 0
    cm <- list(restart = has_theta)
    if (has_theta) cm$theta <- prev_mode$theta
    if (transfer_x && length(prev_mode$x) > 0) cm$x <- prev_mode$x
    cm
  }

  # --- Helper: pack glm summary into INLA-compatible format ---
  glm_to_summary <- function(fit) {
    cc <- summary(fit)$coefficients
    data.frame(
      mean = cc[, "Estimate"],
      sd   = cc[, "Std. Error"],
      row.names = rownames(cc),
      check.names = FALSE
    )
  }

  # --- Single imputation: draw z, fit occ|z, fit det|z=1 ---
  # Uses glm() for non-spatial models (milliseconds per call);
  # falls back to INLA only when spatial structure (SPDE) is needed.
  run_one <- function(seed) {
    set.seed(seed)
    z_imp <- as.integer(detected)
    z_imp[undetected_idx] <- rbinom(
      length(undetected_idx), 1, clamp(z_hat[undetected_idx])
    )

    # -- Occupancy | z --
    occ_df <- data$occ.covs
    occ_df$site <- seq_len(N)
    occ_df$z <- z_imp

    if (use_spatial) {
      si <- occ_prep$stack_info
      A_obs <- si$A[seq_len(N), , drop = FALSE]
      cov_df <- occ_df[, setdiff(names(occ_df), c("z", "site")), drop = FALSE]
      stack <- INLA::inla.stack(
        data = list(z = z_imp), A = list(A_obs, 1),
        effects = list(list(spatial = seq_len(si$n_mesh)), cov_df), tag = "occ"
      )
      stack_data <- INLA::inla.stack.data(stack)
      stack_data$spde <- si$spde
      occ_fit <- suppressWarnings(INLA::inla(
        formula = occ_formula, family = "binomial",
        Ntrials = rep(1L, nrow(stack_data)), data = stack_data,
        control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
        num.threads = "1:1",
        control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
        control.mode = build_ctrl_mode(em_occ_mode, transfer_x = TRUE),
        control.compute = ctrl, verbose = FALSE
      ))
      idx <- INLA::inla.stack.index(stack, "occ")$data
      psi_k <- occ_fit$summary.fitted.values$mean[idx]
      occ_beta <- occ_fit$summary.fixed
      # Cache INLA mode for warm-start
      if (is.null(em_occ_mode)) {
        em_occ_mode <<- list(theta = occ_fit$mode$theta, x = occ_fit$mode$x)
      }
    } else {
      # GLM: plain binomial, no INLA overhead
      occ_glm <- glm(occ_formula, family = binomial, data = occ_df,
                      control = glm.control(maxit = 50))
      psi_k <- predict(occ_glm, type = "response")
      occ_beta <- glm_to_summary(occ_glm)
    }

    # -- Detection | z=1 subset --
    occ_sites <- which(z_imp == 1L)
    y_sub <- data$y[occ_sites, , drop = FALSE]
    det_covs_sub <- if (!is.null(data$det.covs)) {
      lapply(data$det.covs, function(m) m[occ_sites, , drop = FALSE])
    }
    det_df <- build_det_df(y = y_sub, det_covs = det_covs_sub, site_idx = occ_sites)

    # GLM for detection (always non-spatial in MI correction)
    det_glm <- glm(det_formula, family = binomial, data = det_df,
                    control = glm.control(maxit = 50))
    p_fitted <- predict(det_glm, type = "response")
    p_k <- matrix(NA, N, J)
    p_k[cbind(det_df$site, det_df$visit)] <- p_fitted
    if (any(is.na(p_k))) {
      p_fill <- 1 / (1 + exp(-coef(det_glm)[1]))
      p_k[is.na(p_k)] <- p_fill
    }
    det_beta <- glm_to_summary(det_glm)

    list(occ_beta = occ_beta, det_beta = det_beta, psi = psi_k, p = p_k)
  }

  # --- Rubin's pooling ---
  pool_rubins <- function(draws) {
    Kn <- length(draws)
    coef_names <- rownames(draws[[1]])
    beta_means <- sapply(draws, function(b) b$mean)
    beta_vars  <- sapply(draws, function(b) b$sd^2)
    if (!is.matrix(beta_means)) {
      beta_means <- matrix(beta_means, nrow = 1)
      beta_vars  <- matrix(beta_vars, nrow = 1)
    }
    pooled_mean <- rowMeans(beta_means)
    V_within    <- rowMeans(beta_vars)
    V_between   <- apply(beta_means, 1, var)
    V_total     <- V_within + (1 + 1/Kn) * V_between
    pooled_sd   <- sqrt(V_total)
    data.frame(
      mean         = pooled_mean,
      sd           = pooled_sd,
      `0.025quant` = pooled_mean - 1.96 * pooled_sd,
      `0.5quant`   = pooled_mean,
      `0.975quant` = pooled_mean + 1.96 * pooled_sd,
      row.names    = coef_names,
      check.names  = FALSE
    )
  }

  # --- Run n_rounds of MI, updating z_hat between rounds ---
  # Within each round: K independent draws → pool with Rubin's rules.
  # Between rounds: recompute z_hat from pooled psi/p so next round
  # draws from a better distribution.
  # Early rounds just need to update z_hat (5 draws suffice).
  # Final round uses full K for accurate Rubin's pooling.
  K_early <- min(5L, K)
  current_z_hat <- z_hat
  use_parallel <- n_cores > 1L
  total_K <- 0L

  for (round in seq_len(n_rounds)) {
    is_final <- round == n_rounds
    K_round <- if (is_final) K else K_early
    if (verbose >= 1) cat(sprintf("    Round %d/%d (K=%d)\n", round, n_rounds, K_round))
    seeds <- sample.int(1e7, K_round)

    # run_one closure captures current_z_hat from this environment,
    # but it reads z_hat — update the binding
    z_hat <- current_z_hat

    if (use_parallel) {
      if (verbose >= 1) cat(sprintf("      %d imputations on %d cores...\n", K_round, n_cores))
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c(
        "run_one", "occ_formula", "det_formula", "occ_prep", "det_prep",
        "data", "z_hat", "detected", "undetected_idx", "use_spatial",
        "ctrl", "control.inla", "N", "J", "em_occ_mode", "em_det_mode",
        "build_ctrl_mode", "clamp", "build_det_df"
      ), envir = environment())
      parallel::clusterEvalQ(cl, { library(INLA); NULL })
      results <- parallel::parLapply(cl, seeds, run_one)
      parallel::stopCluster(cl)
      on.exit(NULL)  # clear the on.exit since we stopped manually
      K_used <- K_round
    } else {
      # Sequential with adaptive stopping
      results <- vector("list", K_round)
      occ_sum <- NULL; occ_sum_sq <- NULL
      prev_V <- NULL
      K_used <- K_round

      for (k in seq_len(K_round)) {
        results[[k]] <- run_one(seeds[k])

        if (verbose >= 1 && round == n_rounds) {
          cat(sprintf("      MI %2d: occ=[%s] det=[%s]\n", k,
                      paste(round(results[[k]]$occ_beta$mean, 3), collapse = ", "),
                      paste(round(results[[k]]$det_beta$mean, 3), collapse = ", ")))
        }

        # Adaptive stopping after K_min (only on final round for speed)
        if (round == n_rounds && k >= K_min) {
          combined <- c(results[[k]]$occ_beta$mean, results[[k]]$det_beta$mean)
          if (is.null(occ_sum)) {
            all_combined <- sapply(results[seq_len(k)], function(r)
              c(r$occ_beta$mean, r$det_beta$mean))
            if (!is.matrix(all_combined)) all_combined <- matrix(all_combined, nrow = 1)
            occ_sum <- rowSums(all_combined)
            occ_sum_sq <- rowSums(all_combined^2)
          } else {
            occ_sum <- occ_sum + combined
            occ_sum_sq <- occ_sum_sq + combined^2
          }
          V_between <- (occ_sum_sq / k) - (occ_sum / k)^2

          if (!is.null(prev_V)) {
            denom <- pmax(prev_V, 1e-10)
            rel_change <- max(abs(V_between - prev_V) / denom)
            if (rel_change < mi_tol) {
              K_used <- k
              if (verbose >= 1)
                cat(sprintf("      MI converged at K=%d (rel_change=%.4f)\n", k, rel_change))
              break
            }
          }
          prev_V <- V_between
        }
      }
    }

    # Unpack results
    K_actual <- K_used
    results <- results[seq_len(K_actual)]
    total_K <- total_K + K_actual
    occ_draws <- lapply(results, `[[`, "occ_beta")
    det_draws <- lapply(results, `[[`, "det_beta")
    psi_draws <- sapply(results, `[[`, "psi")
    if (!is.matrix(psi_draws)) psi_draws <- matrix(psi_draws, nrow = N)
    p_draws <- array(NA, dim = c(N, J, K_actual))
    for (k in seq_len(K_actual)) p_draws[, , k] <- results[[k]]$p

    # Pool this round
    psi_hat <- rowMeans(psi_draws)
    p_hat   <- apply(p_draws, c(1, 2), mean, na.rm = TRUE)
    p_hat[is.nan(p_hat)] <- NA

    if (verbose >= 1) {
      occ_pooled <- pool_rubins(occ_draws)$mean
      det_pooled <- pool_rubins(det_draws)$mean
      cat(sprintf("      Pooled: occ=[%s] det=[%s]\n",
                  paste(round(occ_pooled, 3), collapse = ", "),
                  paste(round(det_pooled, 3), collapse = ", ")))
    }

    # Update z_hat for next round
    current_z_hat <- compute_weights(data$y, psi_hat, p_hat)
  }

  list(
    occ_summary = pool_rubins(occ_draws),
    det_summary = pool_rubins(det_draws),
    psi_hat     = psi_hat,
    p_hat       = p_hat,
    n_iter      = total_K,
    n_burn      = 0L
  )
}


# ---------------------------------------------------------------------------
# Gibbs-style data augmentation for joint occ + det debiasing
# ---------------------------------------------------------------------------
#' Alternate between sampling z and refitting both submodels.
#' Each iteration: z|params → occ|z → det|z. Collect post-burn-in draws.
#' Follows tulpa's Gibbs pattern of component-wise conditional updates.
#' @noRd
mi_data_augmentation <- function(occ_prep, det_prep, data, z_hat, spatial,
                                 control.inla = NULL,
                                 n_iter = NULL, n_burn = NULL,
                                 verbose = 1) {
  check_inla()
  N <- data$N
  J <- data$J
  detected <- data$detected
  undetected_idx <- which(!detected)

  if (length(undetected_idx) == 0) return(NULL)


  # Adaptive iteration count: fewer Gibbs draws at large N where the EM
  # already converges well.  Floor of 10 (3 burn + 7 keep) ensures enough
  # post-burn-in draws even for species in hard parameter regimes.
  if (is.null(n_iter)) {
    n_iter <- max(10L, as.integer(16 - 2 * log10(N)))
  }
  if (is.null(n_burn)) {
    n_burn <- max(3L, n_iter %/% 3L)
  }
  if (verbose >= 1) cat(sprintf("  Gibbs MI (%d iter, %d burn-in)...\n", n_iter, n_burn))

  use_spatial <- occ_prep$spatial
  ctrl <- list(config = FALSE, dic = FALSE, waic = FALSE)

  # Build formulas once (global env to avoid INLA caching)
  cov_names <- setdiff(names(data$occ.covs), c("z", "site", "Ntrials", "w"))
  occ_fstr <- if (length(cov_names) > 0) paste("z ~", paste(cov_names, collapse = " + ")) else "z ~ 1"
  if (use_spatial) occ_fstr <- paste(occ_fstr, "+ f(spatial, model = spde)")
  occ_formula <- as.formula(occ_fstr)
  det_fstr <- paste("y_det ~", as.character(det_prep$formula)[3])
  det_formula <- as.formula(det_fstr, env = globalenv())

  # Current state: start from EM's z_hat
  current_z_hat <- z_hat

  # Warm-start: reuse previous iteration's mode to skip Newton iterations.
  # Occupancy model has fixed structure so both theta and x carry over.
  # Detection model changes subset size, so only theta transfers.
  prev_occ_mode_g <- NULL
  prev_det_mode_g <- NULL

  gibbs_extract_mode <- function(fit) {
    if (is.null(fit)) return(NULL)
    list(theta = fit$mode$theta, x = fit$mode$x)
  }

  # Storage for post-burn-in draws
  n_keep <- n_iter - n_burn
  occ_draws <- vector("list", n_keep)
  det_draws <- vector("list", n_keep)
  psi_draws <- matrix(NA, N, n_keep)
  p_draws   <- array(NA, dim = c(N, J, n_keep))

  gibbs_warm_start <- function(prev_mode, transfer_x = TRUE) {
    if (is.null(prev_mode)) return(list())
    has_theta <- length(prev_mode$theta) > 0
    cm <- list(restart = has_theta)
    if (has_theta) cm$theta <- prev_mode$theta
    if (transfer_x && length(prev_mode$x) > 0) cm$x <- prev_mode$x
    list(control.mode = cm)
  }

  for (iter in seq_len(n_iter)) {
    # --- Step 1: Sample z | (psi, p, y) ---
    z_imp <- as.integer(detected)
    z_imp[undetected_idx] <- rbinom(
      length(undetected_idx), 1, clamp(current_z_hat[undetected_idx])
    )

    # --- Step 2: Fit occupancy | z ---
    occ_df <- data$occ.covs
    occ_df$site <- seq_len(N)
    occ_df$z <- z_imp

    if (use_spatial) {
      si <- occ_prep$stack_info
      A_obs <- si$A[seq_len(N), , drop = FALSE]
      cov_df <- occ_df[, setdiff(names(occ_df), c("z", "site")), drop = FALSE]
      stack <- INLA::inla.stack(
        data = list(z = z_imp), A = list(A_obs, 1),
        effects = list(list(spatial = seq_len(si$n_mesh)), cov_df), tag = "occ"
      )
      stack_data <- INLA::inla.stack.data(stack)
      stack_data$spde <- si$spde
      occ_ws <- gibbs_warm_start(prev_occ_mode_g, transfer_x = TRUE)
      occ_fit <- suppressWarnings(INLA::inla(
        formula = occ_formula, family = "binomial",
        Ntrials = rep(1L, nrow(stack_data)), data = stack_data,
        control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
        num.threads = "1:1",
        control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
        control.mode = occ_ws$control.mode,
        control.compute = ctrl, verbose = FALSE
      ))
      idx <- INLA::inla.stack.index(stack, "occ")$data
      psi_iter <- occ_fit$summary.fitted.values$mean[idx]
    } else {
      occ_ws <- gibbs_warm_start(prev_occ_mode_g, transfer_x = TRUE)
      occ_fit <- suppressWarnings(INLA::inla(
        formula = occ_formula, family = "binomial",
        Ntrials = rep(1L, N), data = occ_df, num.threads = "1:1",
        control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
        control.mode = occ_ws$control.mode,
        control.compute = ctrl, verbose = FALSE
      ))
      psi_iter <- occ_fit$summary.fitted.values$mean[seq_len(N)]
    }

    # --- Step 3: Fit detection | z=1 subset ---
    occ_sites <- which(z_imp == 1L)
    y_sub <- data$y[occ_sites, , drop = FALSE]
    det_covs_sub <- if (!is.null(data$det.covs)) {
      lapply(data$det.covs, function(m) m[occ_sites, , drop = FALSE])
    }
    det_df <- build_det_df(y = y_sub, det_covs = det_covs_sub, site_idx = occ_sites)
    det_df$.Ntrials <- rep(1L, nrow(det_df))

    det_ws <- gibbs_warm_start(prev_det_mode_g, transfer_x = FALSE)
    det_fit <- suppressWarnings(INLA::inla(
      formula = det_formula, family = "binomial",
      Ntrials = det_df$.Ntrials, data = det_df, num.threads = "1:1",
      control.inla = control.inla %||% INLA::inla.set.control.inla.default(),
      control.mode = det_ws$control.mode,
      control.compute = ctrl, verbose = FALSE
    ))

    p_fitted <- det_fit$summary.fitted.values$mean
    p_iter <- matrix(NA, N, J)
    p_iter[cbind(det_df$site, det_df$visit)] <- p_fitted
    # Fill p for non-occupied sites using intercept-only prediction
    if (any(is.na(p_iter))) {
      p_fill <- 1 / (1 + exp(-det_fit$summary.fixed$mean[1]))
      p_iter[is.na(p_iter)] <- p_fill
    }

    # --- Cache fits for warm-starting next iteration ---
    prev_occ_mode_g <- gibbs_extract_mode(occ_fit)
    prev_det_mode_g <- gibbs_extract_mode(det_fit)

    # --- Update z_hat for next iteration ---
    current_z_hat <- compute_weights(data$y, psi_iter, p_iter)

    # Store post-burn-in draws
    if (iter > n_burn) {
      k <- iter - n_burn
      occ_draws[[k]] <- occ_fit$summary.fixed
      det_draws[[k]] <- det_fit$summary.fixed
      psi_draws[, k] <- psi_iter
      p_draws[, , k] <- p_iter

      if (verbose >= 1) {
        cat(sprintf("    iter %2d: occ=[%s] det=[%s]\n", iter,
                    paste(round(occ_fit$summary.fixed$mean, 3), collapse = ", "),
                    paste(round(det_fit$summary.fixed$mean, 3), collapse = ", ")))
      }

      # Early stopping: if we have >= 3 post-burn-in draws and the running
      # mean has stabilized (max change < 0.01), stop the chain.
      if (k >= 3) {
        safe_mean <- function(draws, idx) {
          m <- sapply(draws[idx], function(b) b$mean)
          if (!is.matrix(m)) m <- matrix(m, nrow = 1)
          rowMeans(m)
        }
        prev_occ <- safe_mean(occ_draws, seq_len(k - 1))
        curr_occ <- safe_mean(occ_draws, seq_len(k))
        prev_det <- safe_mean(det_draws, seq_len(k - 1))
        curr_det <- safe_mean(det_draws, seq_len(k))
        max_shift <- max(abs(c(curr_occ - prev_occ, curr_det - prev_det)))
        if (max_shift < 0.01) {
          if (verbose >= 1) cat(sprintf("    Converged at iter %d (shift=%.4f)\n", iter, max_shift))
          n_keep <- k
          break
        }
      }
    }
  }

  # Truncate draws to actual n_keep (early stopping may reduce it)
  occ_draws <- occ_draws[seq_len(n_keep)]
  det_draws <- det_draws[seq_len(n_keep)]
  psi_draws <- psi_draws[, seq_len(n_keep), drop = FALSE]
  p_draws   <- p_draws[, , seq_len(n_keep), drop = FALSE]

  # --- Pool with Rubin's rules ---
  pool_rubins <- function(draws) {
    K <- length(draws)
    coef_names <- rownames(draws[[1]])
    beta_means <- sapply(draws, function(b) b$mean)
    beta_vars  <- sapply(draws, function(b) b$sd^2)
    if (!is.matrix(beta_means)) {
      beta_means <- matrix(beta_means, nrow = 1)
      beta_vars  <- matrix(beta_vars, nrow = 1)
    }
    pooled_mean <- rowMeans(beta_means)
    V_within    <- rowMeans(beta_vars)
    V_between   <- apply(beta_means, 1, var)
    V_total     <- V_within + (1 + 1/K) * V_between
    pooled_sd   <- sqrt(V_total)
    data.frame(
      mean         = pooled_mean,
      sd           = pooled_sd,
      `0.025quant` = pooled_mean - 1.96 * pooled_sd,
      `0.5quant`   = pooled_mean,
      `0.975quant` = pooled_mean + 1.96 * pooled_sd,
      row.names    = coef_names,
      check.names  = FALSE
    )
  }

  psi_hat <- rowMeans(psi_draws)
  p_hat   <- apply(p_draws, c(1, 2), mean, na.rm = TRUE)
  p_hat[is.nan(p_hat)] <- NA

  if (verbose >= 1) {
    occ_pooled <- pool_rubins(occ_draws)$mean
    det_pooled <- pool_rubins(det_draws)$mean
    cat(sprintf("    Pooled: occ=[%s] det=[%s]\n",
                paste(round(occ_pooled, 3), collapse = ", "),
                paste(round(det_pooled, 3), collapse = ", ")))
  }

  list(
    occ_summary = pool_rubins(occ_draws),
    det_summary = pool_rubins(det_draws),
    psi_hat     = psi_hat,
    p_hat       = p_hat,
    n_iter      = n_iter,
    n_burn      = n_burn
  )
}
