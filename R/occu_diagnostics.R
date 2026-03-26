# =============================================================================
# occu_diagnostics.R — fitted, residuals, GOF, WAIC for INLA occupancy models
# =============================================================================

# ---------------------------------------------------------------------------
#  fitted.occu_inla  —  cf. fitted.PGOcc
# ---------------------------------------------------------------------------

#' Extract fitted values from an occupancy model
#'
#' Returns fitted detection-level values (y.rep) and detection probabilities.
#' Analogous to spOccupancy's fitted() method.
#'
#' @param object fitted occu_inla object
#' @param ... ignored
#'
#' @return list with:
#'   \describe{
#'     \item{y.rep}{N x J matrix of expected detection values (z_hat * p_hat)}
#'     \item{p}{N x J matrix of estimated detection probabilities}
#'     \item{psi}{length-N vector of estimated occupancy probabilities}
#'     \item{z}{length-N vector of posterior P(z=1)}
#'   }
#' @export
fitted.occu_inla <- function(object, ...) {
  y_rep <- object$z_hat * object$p_hat
  # For unvisited cells, set to NA
  y_rep[is.na(object$data$y)] <- NA

  list(
    y.rep = y_rep,
    p     = object$p_hat,
    psi   = object$psi_hat,
    z     = object$z_hat
  )
}


# ---------------------------------------------------------------------------
#  residuals.occu_inla  —  cf. residuals.PGOcc
# ---------------------------------------------------------------------------

#' Compute residuals for an occupancy model
#'
#' Occupancy residuals (site-level) and detection residuals (visit-level).
#'
#' @param object fitted occu_inla object
#' @param type "deviance" (default), "pearson", or "response"
#' @param ... ignored
#'
#' @return list with:
#'   \describe{
#'     \item{occ.resids}{length-N occupancy residuals}
#'     \item{det.resids}{N x J detection residuals}
#'   }
#' @export
residuals.occu_inla <- function(object, type = c("deviance", "pearson", "response"),
                                ...) {
  type <- match.arg(type)

  y     <- object$data$y
  psi   <- clamp(object$psi_hat)
  p_hat <- clamp(object$p_hat)
  z_hat <- object$z_hat
  N     <- object$data$N
  J     <- object$data$J

  # Detection residuals
  y_exp <- z_hat * p_hat
  if (type == "response") {
    det_resid <- y - y_exp
    occ_resid <- as.integer(object$data$detected) - psi
  } else if (type == "pearson") {
    det_var <- y_exp * (1 - p_hat)
    det_resid <- (y - y_exp) / sqrt(clamp(det_var, lo = 1e-6))
    occ_var <- psi * (1 - psi)
    occ_resid <- (as.integer(object$data$detected) - psi) / sqrt(clamp(occ_var, lo = 1e-6))
  } else {
    # Deviance residuals (vectorized)
    mu <- clamp(y_exp)
    det_resid <- ifelse(y == 1, sqrt(-2 * log(mu)), -sqrt(-2 * log(1 - mu)))
    det_resid[is.na(y)] <- NA

    d <- as.integer(object$data$detected)
    mu_occ <- clamp(psi)
    occ_resid <- ifelse(d == 1, sqrt(-2 * log(mu_occ)), -sqrt(-2 * log(1 - mu_occ)))
  }

  list(
    occ.resids = occ_resid,
    det.resids = det_resid
  )
}


# ---------------------------------------------------------------------------
#  fitted / residuals for integrated models
# ---------------------------------------------------------------------------

#' Extract fitted values from an integrated occupancy model
#'
#' @param object fitted occu_inla_int object
#' @param ... ignored
#' @return list with per-source y.rep, p, and shared psi, z
#' @export
fitted.occu_inla_int <- function(object, ...) {
  n_data <- object$n_data
  y_reps <- vector("list", n_data)
  for (d in seq_len(n_data)) {
    y_d <- object$data$y[[d]]
    sites_d <- object$data$sites[[d]]
    z_d <- object$z_hat[sites_d]
    p_d <- object$p_hats[[d]]
    y_rep <- z_d * p_d
    y_rep[is.na(y_d)] <- NA
    y_reps[[d]] <- y_rep
  }

  list(
    y.rep = y_reps,
    p     = object$p_hats,
    psi   = object$psi_hat,
    z     = object$z_hat
  )
}


#' Compute residuals for an integrated occupancy model
#'
#' @param object fitted occu_inla_int object
#' @param type "deviance" (default), "pearson", or "response"
#' @param ... ignored
#' @return list with per-source detection residuals and shared occupancy residuals
#' @export
residuals.occu_inla_int <- function(object,
                                    type = c("deviance", "pearson", "response"),
                                    ...) {
  type <- match.arg(type)
  N_total <- length(object$psi_hat)
  psi <- clamp(object$psi_hat)
  z_hat <- object$z_hat
  detected_any <- z_hat == 1

  # Occupancy residuals (shared)
  if (type == "response") {
    occ_resid <- as.integer(detected_any) - psi
  } else if (type == "pearson") {
    occ_var <- psi * (1 - psi)
    occ_resid <- (as.integer(detected_any) - psi) / sqrt(clamp(occ_var, lo = 1e-6))
  } else {
    d_vec <- as.integer(detected_any)
    mu_occ <- clamp(psi)
    occ_resid <- ifelse(d_vec == 1, sqrt(-2 * log(mu_occ)), -sqrt(-2 * log(1 - mu_occ)))
  }

  # Per-source detection residuals
  det_resids <- vector("list", object$n_data)
  for (d in seq_len(object$n_data)) {
    y_d <- object$data$y[[d]]
    sites_d <- object$data$sites[[d]]
    z_d <- z_hat[sites_d]
    p_d <- clamp(object$p_hats[[d]])
    y_exp <- z_d * p_d

    if (type == "response") {
      det_resids[[d]] <- y_d - y_exp
    } else if (type == "pearson") {
      det_var <- y_exp * (1 - p_d)
      det_resids[[d]] <- (y_d - y_exp) / sqrt(clamp(det_var, lo = 1e-6))
    } else {
      mu_d <- clamp(y_exp)
      dr <- ifelse(y_d == 1, sqrt(-2 * log(mu_d)), -sqrt(-2 * log(1 - mu_d)))
      dr[is.na(y_d)] <- NA
      det_resids[[d]] <- dr
    }
  }

  list(occ.resids = occ_resid, det.resids = det_resids)
}


# ---------------------------------------------------------------------------
#  fitted / residuals for multi-species models
# ---------------------------------------------------------------------------

#' Extract fitted values from a multi-species occupancy model
#'
#' @param object fitted occu_inla_ms object
#' @param ... ignored
#' @return named list of per-species fitted value lists
#' @export
fitted.occu_inla_ms <- function(object, ...) {
  results <- list()
  for (sp in object$species_names) {
    fit <- object$species_fits[[sp]]
    if (!is.null(fit)) {
      results[[sp]] <- fitted.occu_inla(fit, ...)
    }
  }
  results
}


#' Compute residuals for a multi-species occupancy model
#'
#' @param object fitted occu_inla_ms object
#' @param type "deviance" (default), "pearson", or "response"
#' @param ... ignored
#' @return named list of per-species residual lists
#' @export
residuals.occu_inla_ms <- function(object,
                                   type = c("deviance", "pearson", "response"),
                                   ...) {
  type <- match.arg(type)
  results <- list()
  for (sp in object$species_names) {
    fit <- object$species_fits[[sp]]
    if (!is.null(fit)) {
      results[[sp]] <- residuals.occu_inla(fit, type = type, ...)
    }
  }
  results
}


# ---------------------------------------------------------------------------
#  fitted / residuals for temporal models
# ---------------------------------------------------------------------------

#' Extract fitted values from a temporal occupancy model
#'
#' @param object fitted occu_inla_temporal object
#' @param ... ignored
#' @return list of per-period fitted value lists
#' @export
fitted.occu_inla_temporal <- function(object, ...) {
  results <- vector("list", object$n_periods)
  for (t in seq_len(object$n_periods)) {
    fit <- object$period_fits[[t]]
    if (!is.null(fit)) {
      results[[t]] <- fitted.occu_inla(fit, ...)
    }
  }
  results
}


#' Compute residuals for a temporal occupancy model
#'
#' @param object fitted occu_inla_temporal object
#' @param type "deviance" (default), "pearson", or "response"
#' @param ... ignored
#' @return list of per-period residual lists
#' @export
residuals.occu_inla_temporal <- function(object,
                                         type = c("deviance", "pearson", "response"),
                                         ...) {
  type <- match.arg(type)
  results <- vector("list", object$n_periods)
  for (t in seq_len(object$n_periods)) {
    fit <- object$period_fits[[t]]
    if (!is.null(fit)) {
      results[[t]] <- residuals.occu_inla(fit, type = type, ...)
    }
  }
  results
}


# ---------------------------------------------------------------------------
#  fitted for JSDM (no detection process)
# ---------------------------------------------------------------------------

#' Extract fitted values from a JSDM
#'
#' @param object fitted occu_inla_jsdm object
#' @param ... ignored
#' @return named list of per-species fitted probability vectors
#' @export
fitted.occu_inla_jsdm <- function(object, ...) {
  results <- list()
  for (s in seq_along(object$species_fits)) {
    fit <- object$species_fits[[s]]
    sp <- object$species_names[s]
    if (!is.null(fit)) {
      results[[sp]] <- expit(fit$summary.fitted.values$mean)
    }
  }
  results
}


# ---------------------------------------------------------------------------
#  ppcOccu  —  posterior predictive checks  (cf. ppcOcc)
# ---------------------------------------------------------------------------

#' Posterior predictive checks for occupancy models
#'
#' Computes a goodness-of-fit statistic for observed and replicated data.
#' Analogous to spOccupancy::ppcOcc().
#'
#' @param object fitted occu_inla object
#' @param fit.stat "freeman-tukey" (default) or "chi-squared"
#' @param group 1 = aggregate by site, 2 = aggregate by visit
#' @param n.samples number of replicated datasets to generate (default 500)
#'
#' @return list with:
#'   \describe{
#'     \item{fit.y}{vector of fit statistic for observed data across samples}
#'     \item{fit.y.rep}{vector of fit statistic for replicated data}
#'     \item{bayesian.p}{Bayesian p-value (proportion fit.y.rep > fit.y)}
#'   }
#' @export
ppcOccu <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                    group = 1, n.samples = 500) {
  fit.stat <- match.arg(fit.stat)

  y     <- object$data$y
  p_hat <- clamp(object$p_hat)
  z_hat <- object$z_hat
  N     <- object$data$N
  J     <- object$data$J

  # Expected values under the model
  E_y <- z_hat * p_hat
  E_y[is.na(y)] <- NA

  fit_y     <- numeric(n.samples)
  fit_y_rep <- numeric(n.samples)

  # Pre-compute NA mask and observed aggregates (invariant across samples)
  not_na <- !is.na(y)
  n_obs <- sum(not_na)
  psi_clamped <- clamp(object$psi_hat)
  if (group == 1) {
    obs_stat <- rowSums(y, na.rm = TRUE)
    exp_stat <- rowSums(E_y, na.rm = TRUE)
  } else {
    obs_stat <- colSums(y, na.rm = TRUE)
    exp_stat <- colSums(E_y, na.rm = TRUE)
  }

  for (s in seq_len(n.samples)) {
    # Generate replicated data (vectorized)
    z_rep <- rbinom(N, 1, psi_clamped)
    prob_rep <- z_rep * p_hat
    y_rep <- matrix(NA_integer_, N, J)
    y_rep[not_na] <- rbinom(n_obs, 1, prob_rep[not_na])

    E_y_rep <- z_rep * p_hat
    E_y_rep[!not_na] <- NA

    if (group == 1) {
      rep_stat   <- rowSums(y_rep, na.rm = TRUE)
      exp_r_stat <- rowSums(E_y_rep, na.rm = TRUE)
    } else {
      rep_stat   <- colSums(y_rep, na.rm = TRUE)
      exp_r_stat <- colSums(E_y_rep, na.rm = TRUE)
    }

    if (fit.stat == "freeman-tukey") {
      fit_y[s]     <- sum((sqrt(obs_stat) - sqrt(clamp(exp_stat, lo = 0)))^2)
      fit_y_rep[s] <- sum((sqrt(rep_stat) - sqrt(clamp(exp_r_stat, lo = 0)))^2)
    } else {
      fit_y[s]     <- sum((obs_stat - exp_stat)^2 / clamp(exp_stat, lo = 1e-6))
      fit_y_rep[s] <- sum((rep_stat - exp_r_stat)^2 / clamp(exp_r_stat, lo = 1e-6))
    }
  }

  # Group-level quantiles (spOccupancy compatibility)
  fit.y.group.quants     <- apply(matrix(fit_y, ncol = 1), 2,
                                   quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  fit.y.rep.group.quants <- apply(matrix(fit_y_rep, ncol = 1), 2,
                                   quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

  list(
    fit.y                  = fit_y,
    fit.y.rep              = fit_y_rep,
    fit.y.group.quants     = fit.y.group.quants,
    fit.y.rep.group.quants = fit.y.rep.group.quants,
    bayesian.p             = mean(fit_y_rep > fit_y),
    fit.stat               = fit.stat,
    group                  = group,
    n.samples              = n.samples
  )
}


# ---------------------------------------------------------------------------
#  waicOccu  —  WAIC computation  (cf. waicOcc)
# ---------------------------------------------------------------------------

#' Compute WAIC for occupancy models
#'
#' Returns WAIC for the occupancy and detection components separately
#' and combined. Analogous to spOccupancy::waicOcc().
#'
#' @param object fitted occu_inla object
#' @param by.sp logical: if TRUE and object is multi-species, return per-species WAIC
#'
#' @return data.frame with elpd, pD, and WAIC columns
#' @export
waicOccu <- function(object, by.sp = FALSE) {

  if (inherits(object, "occu_inla_ms") && by.sp) {
    # Per-species WAIC
    results <- lapply(object$species_names, function(sp) {
      fit <- object$species_fits[[sp]]
      if (is.null(fit)) return(data.frame(species = sp, elpd = NA, pD = NA, WAIC = NA))

      occ_waic <- fit$occ_fit$waic
      det_waic <- fit$det_fit$waic

      data.frame(
        species = sp,
        elpd    = sum(occ_waic$local.waic / -2, na.rm = TRUE) +
                  sum(det_waic$local.waic / -2, na.rm = TRUE),
        pD      = (occ_waic$p.eff %||% NA) + (det_waic$p.eff %||% NA),
        WAIC    = (occ_waic$waic %||% NA) + (det_waic$waic %||% NA)
      )
    })
    return(do.call(rbind, results))
  }

  if (inherits(object, "occu_inla_int")) {
    # Integrated model: per-source WAIC for detection
    occ_waic <- object$occ_fit$waic
    parts <- list(data.frame(
      component = "occupancy",
      elpd      = sum(occ_waic$local.waic / -2, na.rm = TRUE),
      pD        = occ_waic$p.eff %||% NA,
      WAIC      = occ_waic$waic %||% NA
    ))

    for (d in seq_along(object$det_fits)) {
      det_waic <- object$det_fits[[d]]$waic
      if (!is.null(det_waic)) {
        parts[[length(parts) + 1L]] <- data.frame(
          component = paste0("detection_", d),
          elpd      = sum(det_waic$local.waic / -2, na.rm = TRUE),
          pD        = det_waic$p.eff %||% NA,
          WAIC      = det_waic$waic %||% NA
        )
      }
    }

    results <- do.call(rbind, parts)
    # Total
    results <- rbind(results, data.frame(
      component = "total",
      elpd      = sum(results$elpd, na.rm = TRUE),
      pD        = sum(results$pD, na.rm = TRUE),
      WAIC      = sum(results$WAIC, na.rm = TRUE)
    ))
    return(results)
  }

  # Standard single-species model
  occ_waic <- object$occ_fit$waic
  det_waic <- object$det_fit$waic

  data.frame(
    component = c("occupancy", "detection", "total"),
    elpd = c(
      sum(occ_waic$local.waic / -2, na.rm = TRUE),
      sum(det_waic$local.waic / -2, na.rm = TRUE),
      sum(occ_waic$local.waic / -2, na.rm = TRUE) +
        sum(det_waic$local.waic / -2, na.rm = TRUE)
    ),
    pD = c(
      occ_waic$p.eff %||% NA,
      det_waic$p.eff %||% NA,
      (occ_waic$p.eff %||% 0) + (det_waic$p.eff %||% 0)
    ),
    WAIC = c(
      occ_waic$waic %||% NA,
      det_waic$waic %||% NA,
      (occ_waic$waic %||% 0) + (det_waic$waic %||% 0)
    )
  )
}


# ---------------------------------------------------------------------------
#  getSVCSamples  —  extract SVC posterior  (cf. getSVCSamples)
# ---------------------------------------------------------------------------

#' Extract spatially-varying coefficient summaries
#'
#' @param object fitted occu_inla_svc object
#' @return list of data.frames (one per SVC) with site-level mean, sd, quantiles
#' @export
getSVCSamples <- function(object) {
  if (!inherits(object, "occu_inla_svc")) {
    stop("getSVCSamples requires an occu_inla_svc object")
  }

  svc_names <- object$svc_names
  svc_results <- list()

  for (k in seq_along(svc_names)) {
    re_name <- paste0("occ_re_svc_spatial_", k)
    if (re_name %in% names(object$occ_fit$summary.random)) {
      re_summary <- object$occ_fit$summary.random[[re_name]]
      fixed_beta <- object$occ_fit$summary.fixed$mean[
        which(rownames(object$occ_fit$summary.fixed) == svc_names[k])
      ]
      if (length(fixed_beta) == 0) fixed_beta <- 0

      svc_results[[svc_names[k]]] <- data.frame(
        site  = seq_len(nrow(re_summary)),
        mean  = fixed_beta + re_summary$mean,
        sd    = re_summary$sd,
        q025  = fixed_beta + re_summary$`0.025quant`,
        q975  = fixed_beta + re_summary$`0.975quant`
      )
    }
  }

  svc_results
}


# ---------------------------------------------------------------------------
#  postHocLM  —  post-hoc linear model on parameters  (cf. postHocLM)
# ---------------------------------------------------------------------------

#' Fit a post-hoc linear model relating covariates to estimated parameters
#'
#' Useful for exploring drivers of occupancy/detection variation.
#' Supports both Bayesian (via INLA) and frequentist (\code{\link{lm}})
#' fitting.
#'
#' @param formula model formula (e.g., psi_hat ~ trait1 + trait2)
#' @param data data.frame with response and predictors
#' @param weights optional weights (e.g., inverse of psi standard errors)
#' @param method \code{"bayes"} (default, uses INLA) or \code{"freq"} (uses lm).
#'   Falls back to \code{"freq"} if INLA is not available.
#' @param n.samples number of posterior samples (default 1000, Bayesian only)
#'
#' @return A list of class \code{"postHocLM"} with:
#'   \describe{
#'     \item{beta.samples}{matrix of posterior coefficient samples (Bayesian only)}
#'     \item{tau.sq.samples}{posterior residual variance samples (Bayesian only)}
#'     \item{y.hat.samples}{posterior fitted value samples (Bayesian only)}
#'     \item{bayes.R2}{posterior samples of Bayesian R-squared (Bayesian only)}
#'     \item{summary}{data.frame of coefficient summaries}
#'     \item{lm.fit}{frequentist lm object (always available)}
#'     \item{method}{character: which method was used}
#'   }
#' @export
postHocLM <- function(formula, data, weights = NULL,
                       method = c("bayes", "freq"),
                       n.samples = 1000L) {
  method <- match.arg(method)

  # Always fit frequentist model
  lm_fit <- if (is.null(weights)) {
    lm(formula, data = data)
  } else {
    lm(formula, data = data, weights = weights)
  }

  # Try Bayesian version via INLA
  if (method == "bayes" && requireNamespace("INLA", quietly = TRUE)) {
    inla_fit <- tryCatch({
      INLA::inla(
        formula = formula,
        family  = "gaussian",
        data    = data,
        control.compute = list(config = TRUE),
        verbose = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(inla_fit)) {
      # Draw posterior samples
      samples <- tryCatch(
        INLA::inla.posterior.sample(n.samples, inla_fit),
        error = function(e) NULL
      )

      if (!is.null(samples)) {
        fix_names <- rownames(inla_fit$summary.fixed)
        n_coef <- length(fix_names)

        # Build name→index map once from first sample
        latent_names <- rownames(samples[[1]]$latent)
        fix_idx <- vapply(fix_names, function(nm) {
          idx <- grep(paste0("^", nm, ":"), latent_names)
          if (length(idx) > 0) idx[1] else NA_integer_
        }, integer(1))

        beta_samples <- matrix(NA_real_, n.samples, n_coef)
        colnames(beta_samples) <- fix_names
        tau_sq_samples <- numeric(n.samples)

        for (s in seq_len(n.samples)) {
          latent <- samples[[s]]$latent
          valid <- !is.na(fix_idx)
          beta_samples[s, valid] <- latent[fix_idx[valid]]
          hyp <- samples[[s]]$hyperpar
          if (length(hyp) > 0) {
            tau_sq_samples[s] <- 1 / max(hyp[1], 1e-8)
          }
        }

        # Vectorized fitted values: (n.samples x n_coef) %*% t(X) → n.samples x n
        X <- model.matrix(formula, data = data)
        y_hat_samples <- beta_samples %*% t(X)

        # Bayesian R-squared (vectorized)
        var_fit <- apply(y_hat_samples, 1, var)
        bayes_r2 <- var_fit / (var_fit + tau_sq_samples)

        out <- list(
          beta.samples   = beta_samples,
          tau.sq.samples = tau_sq_samples,
          y.hat.samples  = y_hat_samples,
          bayes.R2       = bayes_r2,
          summary        = inla_fit$summary.fixed,
          lm.fit         = lm_fit,
          method         = "bayes"
        )
        class(out) <- "postHocLM"
        return(out)
      }
    }
  }

  # Frequentist (explicit choice or INLA fallback)
  out <- list(
    beta.samples   = NULL,
    tau.sq.samples = NULL,
    y.hat.samples  = NULL,
    bayes.R2       = NULL,
    summary        = summary(lm_fit)$coefficients,
    lm.fit         = lm_fit,
    method         = "freq"
  )
  class(out) <- "postHocLM"
  out
}


# =============================================================================
# simulate / autocorrelation / variogram / DHARMa
# =============================================================================

# ---------------------------------------------------------------------------
#  simulate.occu_inla  —  posterior predictive simulation
# ---------------------------------------------------------------------------

#' Simulate replicate datasets from a fitted occupancy model
#'
#' Generates posterior predictive datasets by sampling the latent occupancy
#' state z from Bernoulli(psi) and observations y from Bernoulli(z * p)
#' at each visit.  Returns a site-level summary (number of detections per
#' site) suitable for DHARMa's \code{createDHARMa()}.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of replicate datasets (default 250)
#' @param seed random seed for reproducibility (default 123)
#' @param level \code{"site"} (default) returns an N x nsim matrix of
#'   per-site detection counts; \code{"obs"} returns an (N*J) x nsim matrix
#'   of individual visit outcomes.
#' @param ... ignored
#'
#' @return A matrix (see \code{level}).
#' @export
simulate.occu_inla <- function(object, nsim = 250L, seed = 123L,
                                level = c("site", "obs"), ...) {
  level <- match.arg(level)
  set.seed(seed)

  psi <- clamp(object$psi_hat)
  p   <- clamp(object$p_hat)
  N   <- object$data$N
  J   <- object$data$J
  not_na <- !is.na(object$data$y)

  if (level == "site") {
    out <- matrix(NA_integer_, N, nsim)
    for (s in seq_len(nsim)) {
      z <- rbinom(N, 1L, psi)
      y_sim <- matrix(0L, N, J)
      y_sim[not_na] <- rbinom(sum(not_na), 1L, (z * p)[not_na])
      out[, s] <- rowSums(y_sim)
    }
  } else {
    # obs-level: flatten N x J → (N*J) vector, only non-NA cells
    obs_idx <- which(not_na)
    n_obs <- length(obs_idx)
    out <- matrix(NA_integer_, n_obs, nsim)
    for (s in seq_len(nsim)) {
      z <- rbinom(N, 1L, psi)
      prob <- (z * p)[obs_idx]
      out[, s] <- rbinom(n_obs, 1L, prob)
    }
  }
  out
}


# ---------------------------------------------------------------------------
#  moranI  —  Moran's I test for spatial autocorrelation in residuals
# ---------------------------------------------------------------------------

#' Moran's I test for spatial autocorrelation
#'
#' Computes Moran's I on occupancy residuals using an inverse-distance
#' weight matrix (k nearest neighbours).  No external dependencies —
#' uses the normal approximation under the randomisation assumption.
#'
#' @param object fitted \code{occu_inla} object, or a numeric vector of
#'   residuals (in which case \code{coords} must be supplied)
#' @param coords optional N x 2 coordinate matrix.
#'   If \code{object} is an \code{occu_inla} fit, defaults to
#'   \code{object$data$coords}.
#' @param weights weight scheme: \code{"inverse"} (default) uses
#'   inverse-distance weights on all pairs (matches DHARMa),
#'   \code{"knn"} uses k nearest-neighbour binary weights
#' @param k number of nearest neighbours (only used when
#'   \code{weights = "knn"}, default 10)
#' @param resid.type residual type passed to \code{residuals()} if
#'   \code{object} is a fit (default \code{"deviance"})
#' @param alternative \code{"two.sided"} (default), \code{"greater"}, or
#'   \code{"less"}
#'
#' @return A list of class \code{"htest"} with components \code{statistic},
#'   \code{p.value}, \code{parameter} (expected I), and \code{method}.
#' @export
moranI <- function(object, coords = NULL,
                   weights = c("inverse", "knn"), k = 10L,
                   resid.type = "deviance",
                   alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  weights <- match.arg(weights)

  # Extract residuals and coords
  if (inherits(object, "occu_em")) {
    x <- residuals(object, type = resid.type)$occ.resids
    if (is.null(coords)) coords <- object$data$coords
  } else if (is.numeric(object)) {
    x <- object
  } else {
    stop("object must be an occu_inla fit or a numeric vector of residuals")
  }
  if (is.null(coords)) stop("coords required for Moran's I")
  coords <- as.matrix(coords)

  N <- length(x)
  if (nrow(coords) != N) {
    stop(sprintf("coords has %d rows but residuals has length %d",
                 nrow(coords), N))
  }

  # Build weight matrix
  D <- as.matrix(dist(coords))
  if (weights == "inverse") {
    W <- 1 / D
    diag(W) <- 0
    # Zero out infinite values (identical points)
    W[!is.finite(W)] <- 0
    method_str <- "Moran's I (inverse-distance weights)"
  } else {
    k <- min(k, N - 1L)
    W <- matrix(0, N, N)
    for (i in seq_len(N)) {
      nn <- order(D[i, ])[2:(k + 1L)]
      W[i, nn] <- 1
    }
    method_str <- sprintf("Moran's I (k=%d nearest neighbours)", k)
  }

  # Moran's I statistic
  xbar <- mean(x)
  dx <- x - xbar
  ss <- sum(dx^2)
  S0 <- sum(W)
  I <- (N / S0) * (as.numeric(dx %*% W %*% dx) / ss)

  # Expected value and variance under randomisation
  EI <- -1 / (N - 1)
  S1 <- 0.5 * sum((W + t(W))^2)
  S2 <- sum((rowSums(W) + colSums(W))^2)
  n2 <- N^2
  k2 <- (sum(dx^4) / N) / (ss / N)^2  # kurtosis
  VI <- (N * ((n2 - 3 * N + 3) * S1 - N * S2 + 3 * S0^2) -
         k2 * (N * (N - 1) * S1 - 2 * N * S2 + 6 * S0^2)) /
        ((N - 1) * (N - 2) * (N - 3) * S0^2) - EI^2

  Z <- (I - EI) / sqrt(VI)
  p_val <- switch(alternative,
    two.sided = 2 * pnorm(abs(Z), lower.tail = FALSE),
    greater   = pnorm(Z, lower.tail = FALSE),
    less      = pnorm(Z, lower.tail = TRUE)
  )

  structure(
    list(
      statistic = c("Moran's I" = I),
      parameter = c("Expected I" = EI),
      p.value   = p_val,
      alternative = alternative,
      method    = method_str,
      data.name = deparse(substitute(object))
    ),
    class = "htest"
  )
}


# ---------------------------------------------------------------------------
#  durbinWatson  —  Durbin-Watson test for temporal autocorrelation
# ---------------------------------------------------------------------------

#' Durbin-Watson test for temporal autocorrelation in residuals
#'
#' Tests whether occupancy residuals across time periods exhibit first-order
#' autocorrelation.  Designed for temporal occupancy models where each period
#' yields a site-averaged residual.
#'
#' @param object fitted \code{occu_inla_temporal} object, or a numeric vector
#'   of temporally-ordered residuals
#' @param resid.type residual type (default \code{"deviance"})
#' @param alternative \code{"two.sided"} (default), \code{"greater"} (positive
#'   autocorrelation), or \code{"less"} (negative autocorrelation)
#'
#' @return A list of class \code{"htest"} with the DW statistic, approximate
#'   p-value (normal approximation), and the lag-1 autocorrelation.
#' @export
durbinWatson <- function(object, resid.type = "deviance",
                          alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)

  if (inherits(object, "occu_inla_temporal")) {
    # Average occupancy residual per period
    r_list <- residuals(object, type = resid.type)
    x <- vapply(r_list, function(r) mean(r$occ.resids, na.rm = TRUE),
                numeric(1))
  } else if (is.numeric(object)) {
    x <- object
  } else {
    stop("object must be an occu_inla_temporal fit or a numeric vector")
  }

  n <- length(x)
  if (n < 3L) stop("need at least 3 time points for Durbin-Watson test")

  # DW statistic
  e <- x - mean(x)
  dw <- sum(diff(e)^2) / sum(e^2)

  # Under H0: E[DW] ≈ 2, r1 ≈ 0
  # Approximate: DW ≈ 2(1 - r1), so r1 ≈ 1 - DW/2
  r1 <- 1 - dw / 2

  # Normal approximation (Pan, 1968): E[DW] ≈ 2, Var[DW] ≈ 4/n
  Z <- (dw - 2) / sqrt(4 / n)
  p_val <- switch(alternative,
    two.sided = 2 * pnorm(abs(Z), lower.tail = FALSE),
    greater   = pnorm(Z, lower.tail = TRUE),   # DW < 2 → positive autocorr
    less      = pnorm(Z, lower.tail = FALSE)    # DW > 2 → negative autocorr
  )

  structure(
    list(
      statistic   = c("DW" = dw),
      parameter   = c("lag-1 r" = r1),
      p.value     = p_val,
      alternative = alternative,
      method      = "Durbin-Watson test for temporal autocorrelation",
      data.name   = deparse(substitute(object))
    ),
    class = "htest"
  )
}


# ---------------------------------------------------------------------------
#  variogram  —  empirical semivariogram of residuals
# ---------------------------------------------------------------------------

#' Empirical semivariogram of occupancy residuals
#'
#' Computes the empirical semivariogram in distance bins.  Useful for
#' visually assessing whether spatial structure remains in the residuals
#' after model fitting.
#'
#' @param object fitted \code{occu_inla} object, or a numeric vector of
#'   residuals (in which case \code{coords} must be supplied)
#' @param coords optional N x 2 coordinate matrix (defaults to
#'   \code{object$data$coords})
#' @param n.bins number of distance bins (default 15)
#' @param max.dist maximum distance to consider (default: half the max
#'   pairwise distance)
#' @param resid.type residual type (default \code{"deviance"})
#'
#' @return A data.frame of class \code{"occu_variogram"} with columns
#'   \code{dist} (bin midpoint), \code{gamma} (semivariance), and
#'   \code{n.pairs} (number of pairs in bin).  Has a \code{plot()} method.
#' @export
variogram <- function(object, coords = NULL, n.bins = 15L,
                       max.dist = NULL, resid.type = "deviance") {

  if (inherits(object, "occu_em")) {
    x <- residuals(object, type = resid.type)$occ.resids
    if (is.null(coords)) coords <- object$data$coords
  } else if (is.numeric(object)) {
    x <- object
  } else {
    stop("object must be an occu_inla fit or a numeric vector of residuals")
  }
  if (is.null(coords)) stop("coords required for variogram")
  coords <- as.matrix(coords)
  N <- length(x)

  # Pairwise distances (lower triangle only)
  D <- dist(coords)
  d_vec <- as.numeric(D)
  if (is.null(max.dist)) max.dist <- max(d_vec) / 2

  # Pairwise squared differences
  idx <- which(lower.tri(matrix(0, N, N)), arr.ind = TRUE)
  sq_diff <- (x[idx[, 1]] - x[idx[, 2]])^2

  # Bin edges
  breaks <- seq(0, max.dist, length.out = n.bins + 1L)
  mids   <- (breaks[-1] + breaks[-(n.bins + 1L)]) / 2

  gamma   <- numeric(n.bins)
  n_pairs <- integer(n.bins)
  for (b in seq_len(n.bins)) {
    in_bin <- d_vec >= breaks[b] & d_vec < breaks[b + 1L]
    n_pairs[b] <- sum(in_bin)
    if (n_pairs[b] > 0) {
      gamma[b] <- mean(sq_diff[in_bin]) / 2
    }
  }

  out <- data.frame(dist = mids, gamma = gamma, n.pairs = n_pairs)
  out <- out[out$n.pairs > 0, ]
  class(out) <- c("occu_variogram", "data.frame")
  out
}


#' @export
plot.occu_variogram <- function(x, ...) {
  plot(x$dist, x$gamma,
       xlab = "Distance", ylab = "Semivariance",
       main = "Empirical Semivariogram of Residuals",
       pch = 19, cex = sqrt(x$n.pairs / max(x$n.pairs)) * 2,
       ...)
  lines(x$dist, x$gamma, lty = 2, col = "grey50")
  invisible(x)
}


# ---------------------------------------------------------------------------
#  dharma  —  DHARMa convenience wrapper (optional dependency)
# ---------------------------------------------------------------------------

#' Create a DHARMa residuals object from a fitted occupancy model
#'
#' Thin convenience wrapper that calls \code{simulate()} on the fitted model
#' and passes the result to \code{DHARMa::createDHARMa()}.  All DHARMa
#' tests (\code{testSpatialAutocorrelation}, \code{testDispersion}, etc.)
#' then work on the returned object.
#'
#' Requires DHARMa to be installed.  For a dependency-free alternative,
#' use \code{\link{moranI}}, \code{\link{durbinWatson}}, or
#' \code{\link{variogram}} directly.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations (default 250)
#' @param seed random seed (default 123)
#' @param ... passed to \code{DHARMa::createDHARMa()}
#'
#' @return A \code{DHARMa} object (see \code{\link[DHARMa]{createDHARMa}})
#' @export
dharma <- function(object, nsim = 250L, seed = 123L, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("DHARMa package required. Install with: install.packages('DHARMa')")
  }

  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)
  fv   <- rowSums(fitted(object)$y.rep, na.rm = TRUE)

  DHARMa::createDHARMa(
    simulatedResponse       = sims,
    observedResponse        = obs,
    fittedPredictedResponse = fv,
    integerResponse         = TRUE,
    ...
  )
}


# =============================================================================
# Simulation-based diagnostic tests (native DHARMa equivalents)
# =============================================================================

# ---------------------------------------------------------------------------
#  pitResiduals  —  probability integral transform residuals
# ---------------------------------------------------------------------------

#' Compute PIT (scaled) residuals for an occupancy model
#'
#' For each site, computes the quantile of the observed detection count
#' within the posterior predictive distribution from \code{\link{simulate}}.
#' If the model is correct, these residuals are Uniform(0, 1).
#'
#' For integer-valued responses, a randomisation step avoids discrete
#' artefacts: the residual is drawn uniformly between
#' P(sim < obs) and P(sim <= obs).
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations (default 250)
#' @param seed random seed (default 123)
#'
#' @return A numeric vector of length N with values between 0 and 1.
#' @export
pitResiduals <- function(object, nsim = 250L, seed = 123L) {
  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)
  N    <- length(obs)

  # P(sim < obs) and P(sim <= obs) per site
  lower <- rowMeans(sims < obs)
  upper <- rowMeans(sims <= obs)

  # Randomise between lower and upper for integer data
  set.seed(seed + 1L)
  lower + runif(N) * (upper - lower)
}


# ---------------------------------------------------------------------------
#  testUniformity  —  KS test on PIT residuals
# ---------------------------------------------------------------------------

#' Test uniformity of PIT residuals (KS test)
#'
#' If the model is correct, PIT residuals should follow Uniform(0, 1).
#' This test applies a Kolmogorov-Smirnov test against that null.
#' Native equivalent of \code{DHARMa::testUniformity()}.
#'
#' @param object fitted \code{occu_inla} object, or a numeric vector of
#'   PIT residuals (from \code{\link{pitResiduals}})
#' @param nsim number of simulations if \code{object} is a fit (default 250)
#' @param seed random seed (default 123)
#' @param plot if \code{TRUE}, draws a QQ plot of residuals vs Uniform
#'
#' @return A list of class \code{"htest"} (KS test result).
#' @export
testUniformity <- function(object, nsim = 250L, seed = 123L, plot = FALSE) {
  if (inherits(object, "occu_em")) {
    r <- pitResiduals(object, nsim = nsim, seed = seed)
  } else if (is.numeric(object)) {
    r <- object
  } else {
    stop("object must be an occu_inla fit or a numeric vector of PIT residuals")
  }

  # Suppress ties warning — expected for integer data even after PIT randomisation
  result <- suppressWarnings(ks.test(r, "punif"))

  if (plot) {
    n <- length(r)
    expected <- (seq_len(n) - 0.5) / n
    plot(sort(r), expected,
         xlab = "Observed PIT residuals", ylab = "Expected Uniform",
         main = sprintf("QQ Uniform (KS p = %.3f)", result$p.value),
         pch = 19, cex = 0.6)
    abline(0, 1, col = "red", lty = 2)
  }

  result
}


# ---------------------------------------------------------------------------
#  testDispersion  —  over/underdispersion
# ---------------------------------------------------------------------------

#' Test for over- or underdispersion
#'
#' Compares the variance of observed site-level detection counts to the
#' variance expected under the fitted model (via simulation).  A ratio > 1
#' indicates overdispersion; < 1 indicates underdispersion.
#' Native equivalent of \code{DHARMa::testDispersion()}.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations (default 250)
#' @param seed random seed (default 123)
#' @param alternative \code{"two.sided"} (default), \code{"greater"}
#'   (overdispersion), or \code{"less"} (underdispersion)
#'
#' @return A list of class \code{"htest"} with the dispersion ratio and
#'   simulation-based p-value.
#' @export
testDispersion <- function(object, nsim = 250L, seed = 123L,
                            alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)

  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)

  var_obs <- var(obs)
  var_sim <- apply(sims, 2, var)
  ratio   <- var_obs / median(var_sim)

  # Simulation-based p-value
  p_greater <- mean(var_sim >= var_obs)
  p_less    <- mean(var_sim <= var_obs)
  p_val <- switch(alternative,
    two.sided = 2 * min(p_greater, p_less),
    greater   = p_greater,
    less      = p_less
  )
  p_val <- min(p_val, 1)

  structure(
    list(
      statistic   = c("dispersion ratio" = ratio),
      p.value     = p_val,
      alternative = alternative,
      method      = "Simulation-based dispersion test",
      data.name   = deparse(substitute(object))
    ),
    class = "htest"
  )
}


# ---------------------------------------------------------------------------
#  testOutliers  —  observations outside simulation envelope
# ---------------------------------------------------------------------------

#' Test for outliers (simulation envelope)
#'
#' Counts how many sites have observed detection counts outside the
#' min-to-max range of all simulations.  Under a correct model, the
#' expected number of such outliers is approximately
#' \code{2 * N * (1 / (nsim + 1))}.  A binomial test assesses whether
#' more outliers than expected are present.
#' Native equivalent of \code{DHARMa::testOutliers()}.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations (default 250)
#' @param seed random seed (default 123)
#'
#' @return A list of class \code{"htest"} (binomial test result).
#' @export
testOutliers <- function(object, nsim = 250L, seed = 123L) {
  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)
  N    <- length(obs)

  sim_min <- apply(sims, 1, min)
  sim_max <- apply(sims, 1, max)
  n_outliers <- sum(obs < sim_min | obs > sim_max)

  # Expected outlier probability under correct model
  p_expected <- 2 / (nsim + 1)

  result <- binom.test(n_outliers, N, p = p_expected,
                       alternative = "greater")
  result$method <- sprintf(
    "Simulation outlier test (%d/%d sites outside envelope)", n_outliers, N
  )
  result
}


# ---------------------------------------------------------------------------
#  testZeroInflation  —  excess zeros
# ---------------------------------------------------------------------------

#' Test for zero-inflation
#'
#' Compares the number of all-zero sites (no detections) in the observed
#' data to the distribution expected under the fitted model.
#' Native equivalent of \code{DHARMa::testZeroInflation()}.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations (default 250)
#' @param seed random seed (default 123)
#'
#' @return A list of class \code{"htest"} with the zero-inflation ratio
#'   and simulation-based p-value.
#' @export
testZeroInflation <- function(object, nsim = 250L, seed = 123L) {
  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)

  n_zero_obs <- sum(obs == 0)
  n_zero_sim <- colSums(sims == 0)
  ratio <- n_zero_obs / max(median(n_zero_sim), 1)

  # Two-sided simulation p-value
  p_greater <- mean(n_zero_sim >= n_zero_obs)
  p_less    <- mean(n_zero_sim <= n_zero_obs)
  p_val     <- min(2 * min(p_greater, p_less), 1)

  structure(
    list(
      statistic   = c("zero-inflation ratio" = ratio),
      parameter   = c("observed zeros" = n_zero_obs,
                       "expected zeros" = median(n_zero_sim)),
      p.value     = p_val,
      alternative = "two.sided",
      method      = "Simulation-based zero-inflation test",
      data.name   = deparse(substitute(object))
    ),
    class = "htest"
  )
}


# =============================================================================
# checkModel  —  diagnostic panel plot
# =============================================================================

#' Visual diagnostic panel for a fitted occupancy model
#'
#' Produces a 2x2 (or 2x3) panel of diagnostic plots:
#' \enumerate{
#'   \item QQ plot of PIT residuals vs Uniform(0, 1)
#'   \item Residuals vs fitted values
#'   \item Dispersion histogram (observed variance vs simulated)
#'   \item Moran's I correlogram (if coordinates are available)
#' }
#'
#' Analogous to \code{performance::check_model()} or
#' \code{DHARMa::plot.DHARMa()}.
#'
#' @param object fitted \code{occu_inla} object
#' @param nsim number of simulations for PIT and dispersion (default 250)
#' @param seed random seed (default 123)
#'
#' @return Invisible list of test results.
#' @export
checkModel <- function(object, nsim = 250L, seed = 123L) {
  sims <- simulate(object, nsim = nsim, seed = seed, level = "site")
  obs  <- rowSums(object$data$y, na.rm = TRUE)
  pit  <- pitResiduals(object, nsim = nsim, seed = seed)

  has_coords <- !is.null(object$data$coords)
  n_panels <- if (has_coords) 4L else 3L
  layout_mat <- if (has_coords) matrix(1:4, 2, 2) else matrix(1:3, 1, 3)

  old_par <- par(mfrow = if (has_coords) c(2, 2) else c(1, 3),
                 mar = c(4, 4, 2.5, 1))
  on.exit(par(old_par))

  # --- Panel 1: QQ Uniform ---
  n <- length(pit)
  expected <- (seq_len(n) - 0.5) / n
  ks_p <- suppressWarnings(ks.test(pit, "punif"))$p.value
  plot(sort(pit), expected,
       xlab = "PIT residuals", ylab = "Expected Uniform",
       main = sprintf("QQ Uniform (KS p = %.3f)", ks_p),
       pch = 19, cex = 0.5,
       col = adjustcolor("black", 0.6))
  abline(0, 1, col = "red", lty = 2, lwd = 1.5)

  # --- Panel 2: Residuals vs fitted ---
  fv <- rowSums(fitted(object)$y.rep, na.rm = TRUE)
  r  <- residuals(object, type = "deviance")$occ.resids
  plot(fv, r,
       xlab = "Fitted (site detection count)", ylab = "Deviance residuals",
       main = "Residuals vs Fitted",
       pch = 19, cex = 0.5,
       col = adjustcolor("black", 0.6))
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
  lo <- tryCatch(lowess(fv, r), error = function(e) NULL)
  if (!is.null(lo)) lines(lo, col = "blue", lwd = 1.5)

  # --- Panel 3: Dispersion ---
  var_sim <- apply(sims, 2, var)
  var_obs <- var(obs)
  ratio <- var_obs / median(var_sim)
  hist(var_sim, breaks = 20,
       main = sprintf("Dispersion (ratio = %.2f)", ratio),
       xlab = "Simulated variance", col = "grey80", border = "grey50")
  abline(v = var_obs, col = "red", lwd = 2)

  # --- Panel 4: Moran's I correlogram (if spatial) ---
  moran_result <- NULL
  if (has_coords) {
    coords <- object$data$coords
    ks <- c(3, 5, 8, 12, 20)
    ks <- ks[ks < nrow(coords)]
    I_vals <- numeric(length(ks))
    p_vals <- numeric(length(ks))
    for (j in seq_along(ks)) {
      mi <- moranI(r, coords = coords, weights = "knn", k = ks[j])
      I_vals[j] <- mi$statistic
      p_vals[j] <- mi$p.value
    }
    sig <- p_vals < 0.05
    plot(ks, I_vals, type = "b",
         xlab = "k neighbours", ylab = "Moran's I",
         main = "Spatial Correlogram",
         pch = ifelse(sig, 19, 1),
         col = ifelse(sig, "red", "black"))
    abline(h = -1 / (nrow(coords) - 1), col = "grey50", lty = 2)
    moran_result <- data.frame(k = ks, I = I_vals, p = p_vals)
  }

  invisible(list(
    ks_p     = ks_p,
    disp_ratio = ratio,
    moran    = moran_result
  ))
}
