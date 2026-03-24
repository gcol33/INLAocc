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
    # Deviance residuals
    det_resid <- matrix(NA, N, J)
    for (i in seq_len(N)) {
      for (j in seq_len(J)) {
        if (is.na(y[i, j])) next
        mu <- clamp(y_exp[i, j])
        if (y[i, j] == 1) {
          det_resid[i, j] <- sqrt(-2 * log(mu))
        } else {
          det_resid[i, j] <- -sqrt(-2 * log(1 - mu))
        }
      }
    }

    occ_resid <- numeric(N)
    for (i in seq_len(N)) {
      d <- as.integer(object$data$detected[i])
      mu <- clamp(psi[i])
      if (d == 1) {
        occ_resid[i] <- sqrt(-2 * log(mu))
      } else {
        occ_resid[i] <- -sqrt(-2 * log(1 - mu))
      }
    }
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
    occ_resid <- numeric(N_total)
    for (i in seq_len(N_total)) {
      d_i <- as.integer(detected_any[i])
      mu <- clamp(psi[i])
      occ_resid[i] <- if (d_i == 1) sqrt(-2 * log(mu)) else -sqrt(-2 * log(1 - mu))
    }
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
      N_d <- nrow(y_d)
      J_d <- ncol(y_d)
      dr <- matrix(NA, N_d, J_d)
      for (i in seq_len(N_d)) {
        for (j in seq_len(J_d)) {
          if (is.na(y_d[i, j])) next
          mu <- clamp(y_exp[i, j])
          dr[i, j] <- if (y_d[i, j] == 1) sqrt(-2 * log(mu)) else -sqrt(-2 * log(1 - mu))
        }
      }
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

  for (s in seq_len(n.samples)) {
    # Generate replicated data
    z_rep <- rbinom(N, 1, clamp(object$psi_hat))
    y_rep <- matrix(NA, N, J)
    for (i in seq_len(N)) {
      for (j in seq_len(J)) {
        if (!is.na(y[i, j])) {
          y_rep[i, j] <- rbinom(1, 1, z_rep[i] * p_hat[i, j])
        }
      }
    }

    E_y_rep <- z_rep * p_hat
    E_y_rep[is.na(y)] <- NA

    # Compute fit statistic
    if (group == 1) {
      # Aggregate by site
      obs_stat   <- rowSums(y, na.rm = TRUE)
      exp_stat   <- rowSums(E_y, na.rm = TRUE)
      rep_stat   <- rowSums(y_rep, na.rm = TRUE)
      exp_r_stat <- rowSums(E_y_rep, na.rm = TRUE)
    } else {
      # Aggregate by visit
      obs_stat   <- colSums(y, na.rm = TRUE)
      exp_stat   <- colSums(E_y, na.rm = TRUE)
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
    results <- data.frame(
      component = "occupancy",
      elpd      = sum(occ_waic$local.waic / -2, na.rm = TRUE),
      pD        = occ_waic$p.eff %||% NA,
      WAIC      = occ_waic$waic %||% NA
    )

    for (d in seq_along(object$det_fits)) {
      det_waic <- object$det_fits[[d]]$waic
      if (!is.null(det_waic)) {
        results <- rbind(results, data.frame(
          component = paste0("detection_", d),
          elpd      = sum(det_waic$local.waic / -2, na.rm = TRUE),
          pD        = det_waic$p.eff %||% NA,
          WAIC      = det_waic$waic %||% NA
        ))
      }
    }

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

        beta_samples <- matrix(NA_real_, n.samples, n_coef)
        colnames(beta_samples) <- fix_names
        tau_sq_samples <- numeric(n.samples)

        X <- model.matrix(formula, data = data)
        y_hat_samples <- matrix(NA_real_, n.samples, nrow(data))

        for (s in seq_len(n.samples)) {
          latent <- samples[[s]]$latent
          latent_names <- rownames(latent)
          for (k in seq_along(fix_names)) {
            idx <- grep(paste0("^", fix_names[k], ":"), latent_names)
            if (length(idx) > 0) beta_samples[s, k] <- latent[idx[1]]
          }
          # Residual precision → variance
          hyp <- samples[[s]]$hyperpar
          if (length(hyp) > 0) {
            tau_sq_samples[s] <- 1 / max(hyp[1], 1e-8)
          }
          y_hat_samples[s, ] <- as.vector(X %*% beta_samples[s, ])
        }

        # Bayesian R-squared: var(y_hat) / (var(y_hat) + sigma^2)
        resp_name <- as.character(formula)[2]
        y_obs <- data[[resp_name]]
        bayes_r2 <- numeric(n.samples)
        for (s in seq_len(n.samples)) {
          var_fit <- var(y_hat_samples[s, ])
          bayes_r2[s] <- var_fit / (var_fit + tau_sq_samples[s])
        }

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
