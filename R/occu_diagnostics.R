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

  list(
    fit.y       = fit_y,
    fit.y.rep   = fit_y_rep,
    bayesian.p  = mean(fit_y_rep > fit_y),
    fit.stat    = fit.stat,
    group       = group,
    n.samples   = n.samples
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
#' Accounts for uncertainty in the response (occupancy estimates) by
#' using weighted regression.
#'
#' @param formula model formula (e.g., psi_hat ~ trait1 + trait2)
#' @param data data.frame with response and predictors
#' @param weights optional weights (e.g., inverse of psi standard errors)
#'
#' @return lm or INLA fit object
postHocLM <- function(formula, data, weights = NULL) {
  if (is.null(weights)) {
    lm(formula, data = data)
  } else {
    lm(formula, data = data, weights = weights)
  }
}
