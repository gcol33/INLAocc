# =============================================================================
# occu_predict.R — Prediction and derived quantities
# =============================================================================

#' Predict from a fitted occupancy model
#'
#' Supports both in-sample and out-of-sample prediction. Out-of-sample
#' prediction uses a design matrix (X.0) matching spOccupancy's interface.
#'
#' @param object fitted occu_inla object
#' @param X.0 optional: design matrix for occupancy prediction at new locations.
#'   Must include intercept column if model has intercept.
#' @param X.p.0 optional: design matrix for detection prediction.
#' @param coords.0 optional: coordinates for spatial prediction (n_pred x 2).
#' @param ignore.RE logical: set random effects to zero? (default FALSE)
#' @param type "occupancy", "detection", or "both"
#' @param quantiles quantile levels for credible intervals
#' @param n.samples number of posterior samples for uncertainty
#' @param ... additional arguments
#'
#' @return list with psi.0 (occupancy predictions), p.0 (detection predictions),
#'   z.0 (latent state predictions)
predict.occu_inla <- function(object, X.0 = NULL, X.p.0 = NULL,
                              coords.0 = NULL,
                              ignore.RE = FALSE,
                              type = c("both", "occupancy", "detection"),
                              quantiles = c(0.025, 0.5, 0.975),
                              n.samples = 500, ...) {
  type <- match.arg(type)
  check_inla()

  out <- list()

  # =========================================================================
  # Occupancy predictions
  # =========================================================================
  if (type %in% c("both", "occupancy")) {
    if (is.null(X.0)) {
      # --- In-sample ---
      psi_mean <- object$psi_hat
      psi_samples <- NULL

      if (n.samples > 0) {
        samples <- tryCatch({
          INLA::inla.posterior.sample(n.samples, object$occ_fit)
        }, error = function(e) NULL)

        if (!is.null(samples)) {
          N <- object$data$N
          psi_samples <- matrix(NA, N, n.samples)
          for (s in seq_len(n.samples)) {
            psi_samples[, s] <- expit(
              samples[[s]]$latent[seq_len(N)]
            )
          }
        }
      }

      out$psi.0 <- list(
        mean      = psi_mean,
        sd        = if (!is.null(psi_samples)) apply(psi_samples, 1, sd) else NULL,
        quantiles = if (!is.null(psi_samples)) {
          t(apply(psi_samples, 1, quantile, probs = quantiles))
        } else NULL,
        samples   = psi_samples
      )
      out$z.0 <- list(mean = object$z_hat)

    } else {
      # --- Out-of-sample via design matrix X.0 ---
      if (!is.matrix(X.0)) X.0 <- as.matrix(X.0)

      beta <- object$occ_fit$summary.fixed$mean
      beta_sd <- object$occ_fit$summary.fixed$sd

      if (ncol(X.0) != length(beta)) {
        stop(sprintf(
          "X.0 has %d columns but model has %d occupancy coefficients. Include intercept column.",
          ncol(X.0), length(beta)
        ))
      }

      eta <- as.vector(X.0 %*% beta)

      # Spatial component at new locations
      if (!is.null(coords.0) && inherits(object, "occu_inla_spatial")) {
        A_pred <- predict_spatial_A(object$spatial, coords.0)
        sf <- extract_spatial_field(object$occ_fit, object$spatial)
        spatial_pred <- as.vector(A_pred %*% sf$mean)
        eta <- eta + spatial_pred
      }

      psi_pred <- expit(eta)

      # Uncertainty via posterior samples
      psi_samples <- NULL
      if (n.samples > 0) {
        samples <- tryCatch({
          INLA::inla.posterior.sample(n.samples, object$occ_fit)
        }, error = function(e) NULL)

        if (!is.null(samples)) {
          n_pred <- nrow(X.0)
          psi_samples <- matrix(NA, n_pred, n.samples)

          for (s in seq_len(n.samples)) {
            # Extract fixed effects from sample
            fix_names <- rownames(object$occ_fit$summary.fixed)
            beta_s <- vapply(fix_names, function(nm) {
              idx <- which(names(samples[[s]]$latent) == nm)
              if (length(idx) > 0) samples[[s]]$latent[idx[1]] else 0
            }, numeric(1))

            eta_s <- as.vector(X.0 %*% beta_s)

            # Add spatial if available
            if (!is.null(coords.0) && inherits(object, "occu_inla_spatial")) {
              sp_idx <- grep("spatial", names(samples[[s]]$latent))
              if (length(sp_idx) > 0) {
                w_s <- samples[[s]]$latent[sp_idx]
                eta_s <- eta_s + as.vector(A_pred %*% w_s)
              }
            }

            if (ignore.RE) {
              psi_samples[, s] <- expit(eta_s)
            } else {
              psi_samples[, s] <- expit(eta_s)
            }
          }
        }
      }

      out$psi.0 <- list(
        mean      = psi_pred,
        sd        = if (!is.null(psi_samples)) apply(psi_samples, 1, sd) else NULL,
        quantiles = if (!is.null(psi_samples)) {
          t(apply(psi_samples, 1, quantile, probs = quantiles))
        } else NULL,
        samples   = psi_samples
      )

      # Latent occupancy state prediction
      out$z.0 <- list(
        mean = ifelse(psi_pred > 0.5, 1, 0)
      )
    }
  }

  # =========================================================================
  # Detection predictions
  # =========================================================================
  if (type %in% c("both", "detection")) {
    if (is.null(X.p.0)) {
      # In-sample
      out$p.0 <- list(
        mean   = object$p_hat,
        fitted = object$det_fit$summary.fitted.values
      )
    } else {
      if (!is.matrix(X.p.0)) X.p.0 <- as.matrix(X.p.0)

      alpha <- object$det_fit$summary.fixed$mean
      if (ncol(X.p.0) != length(alpha)) {
        stop(sprintf(
          "X.p.0 has %d columns but model has %d detection coefficients",
          ncol(X.p.0), length(alpha)
        ))
      }

      eta_det <- as.vector(X.p.0 %*% alpha)
      out$p.0 <- list(mean = expit(eta_det))
    }
  }

  out
}


#' Predict occupancy at new spatial locations
#'
#' Convenience wrapper for spatial models. Constructs X.0 from a covariate
#' data.frame and handles coordinate projection.
#'
#' @param object fitted occu_inla_spatial object
#' @param newcoords n_pred x 2 matrix of prediction coordinates
#' @param newocc.covs data.frame of occupancy covariates at prediction locations
#' @param ignore.RE ignore random effects?
#' @param n.samples number of posterior samples
#'
#' @return list with predicted psi at new locations
predict_spatial <- function(object, newcoords, newocc.covs = NULL,
                            ignore.RE = FALSE, n.samples = 500) {
  check_inla()

  if (!inherits(object, "occu_inla_spatial")) {
    stop("predict_spatial requires an occu_inla_spatial object")
  }

  # Build design matrix from covariates
  X.0 <- NULL
  if (!is.null(newocc.covs)) {
    X.0 <- model.matrix(object$occ_formula, data = newocc.covs)
  } else {
    # Intercept only
    X.0 <- matrix(1, nrow = nrow(newcoords), ncol = 1)
  }

  predict.occu_inla(
    object, X.0 = X.0, coords.0 = newcoords,
    ignore.RE = ignore.RE,
    type = "occupancy", n.samples = n.samples
  )
}


#' Compute species richness from multi-species model
#'
#' @param object fitted occu_inla_ms object
#' @return data.frame with site-level richness estimates and uncertainty
richness <- function(object) {
  if (!inherits(object, "occu_inla_ms")) {
    stop("richness() requires an occu_inla_ms object")
  }

  N    <- object$data$N
  n_sp <- object$n_species

  z_matrix <- matrix(NA, N, n_sp)
  for (s in seq_len(n_sp)) {
    sp_name <- object$species_names[s]
    if (!is.null(object$species_fits[[sp_name]])) {
      z_matrix[, s] <- object$species_fits[[sp_name]]$z_hat
    }
  }

  rich_mean <- rowSums(z_matrix, na.rm = TRUE)
  rich_sd   <- sqrt(rowSums(z_matrix * (1 - z_matrix), na.rm = TRUE))

  data.frame(
    site       = seq_len(N),
    richness   = rich_mean,
    sd         = rich_sd,
    n_detected = rowSums(z_matrix > 0.5, na.rm = TRUE)
  )
}


#' Compute marginal effects for a covariate
#'
#' @param object fitted occu_inla object
#' @param covariate character: name of covariate
#' @param process "occupancy" or "detection"
#' @param values optional: vector of covariate values
#' @param n_points number of prediction points (if values is NULL)
#' @param other.means named list of values for other covariates (default: means)
#'
#' @return data.frame with covariate values and predicted probability + CI
marginal_effect <- function(object, covariate,
                            process = c("occupancy", "detection"),
                            values = NULL, n_points = 100,
                            other.means = NULL) {
  process <- match.arg(process)

  if (process == "occupancy") {
    fit  <- object$occ_fit
    covs <- object$data$occ.covs
  } else {
    fit  <- object$det_fit
    covs <- object$det_df
  }

  fixed      <- fit$summary.fixed
  coef_names <- rownames(fixed)

  # Find the covariate in the model
  idx <- which(coef_names == covariate)
  if (length(idx) == 0) {
    stop(sprintf("Coefficient '%s' not found. Available: %s",
                 covariate, paste(coef_names, collapse = ", ")))
  }

  # Covariate range
  cov_vals <- covs[[covariate]]
  if (is.null(cov_vals)) {
    stop(sprintf("Covariate '%s' not found in %s data", covariate, process))
  }
  if (is.null(values)) {
    values <- seq(min(cov_vals, na.rm = TRUE),
                  max(cov_vals, na.rm = TRUE),
                  length.out = n_points)
  }

  # Build prediction design matrix with other covariates at their means
  n_coef <- length(coef_names)
  X_pred <- matrix(0, length(values), n_coef)
  colnames(X_pred) <- coef_names

  # Intercept
  X_pred[, 1] <- 1

  # Target covariate
  X_pred[, idx] <- values

  # Other covariates at mean
  for (k in seq_along(coef_names)) {
    if (k == 1 || k == idx) next
    nm <- coef_names[k]
    if (!is.null(other.means) && nm %in% names(other.means)) {
      X_pred[, k] <- other.means[[nm]]
    } else if (!is.null(covs[[nm]])) {
      X_pred[, k] <- mean(covs[[nm]], na.rm = TRUE)
    }
  }

  eta     <- as.vector(X_pred %*% fixed$mean)
  eta_sd  <- sqrt(diag(X_pred %*% diag(fixed$sd^2) %*% t(X_pred)))
  eta_lo  <- eta - 1.96 * eta_sd
  eta_hi  <- eta + 1.96 * eta_sd

  data.frame(
    covariate = covariate,
    value     = values,
    estimate  = expit(eta),
    lower     = expit(eta_lo),
    upper     = expit(eta_hi),
    eta       = eta,
    eta_sd    = eta_sd
  )
}
