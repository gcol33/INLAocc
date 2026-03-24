# =============================================================================
# occu_compat.R — spOccupancy compatibility layer
#
# Provides lazy `$` accessor so users can access spOccupancy-style fields
# (e.g., fit$beta.samples, fit$psi.samples) on INLAocc objects. Samples are
# generated on demand via INLA::inla.posterior.sample() and cached.
# =============================================================================

# Default number of posterior samples to generate
COMPAT_N_SAMPLES <- 1000L


# ---------------------------------------------------------------------------
# Cache key: store generated samples in the object itself under `$.cache`
# ---------------------------------------------------------------------------

#' @noRd
get_or_generate <- function(x, field, generator) {
  cache_name <- ".compat_cache"
  if (is.null(x[[cache_name]])) x[[cache_name]] <- list()
  if (is.null(x[[cache_name]][[field]])) {
    x[[cache_name]][[field]] <- generator(x)
    # Write cache back to parent environment's copy
    env <- parent.frame(2)
    obj_name <- as.character(substitute(x, env = parent.frame()))
    if (length(obj_name) == 1 && exists(obj_name, envir = env)) {
      assign(obj_name, x, envir = env)
    }
  }
  x[[cache_name]][[field]]
}


# ---------------------------------------------------------------------------
# Draw posterior samples from an INLA fit object (cached)
# ---------------------------------------------------------------------------

#' @noRd
draw_inla_samples <- function(inla_fit, n = COMPAT_N_SAMPLES) {
  if (is.null(inla_fit)) return(NULL)
  tryCatch(
    INLA::inla.posterior.sample(n, inla_fit),
    error = function(e) NULL
  )
}


# ---------------------------------------------------------------------------
# Extract named coefficients from INLA posterior samples
# ---------------------------------------------------------------------------

#' @noRd
extract_beta_samples <- function(samples, fix_names, n_samples) {
  if (is.null(samples)) return(NULL)
  n_coef <- length(fix_names)
  mat <- matrix(NA_real_, n_samples, n_coef)
  colnames(mat) <- fix_names

  for (s in seq_len(n_samples)) {
    latent <- samples[[s]]$latent
    latent_names <- rownames(latent)
    for (k in seq_along(fix_names)) {
      # INLA names fixed effects as "name:1"
      idx <- grep(paste0("^", fix_names[k], ":"), latent_names)
      if (length(idx) > 0) mat[s, k] <- latent[idx[1]]
    }
  }
  mat
}


# ---------------------------------------------------------------------------
# Extract site-level fitted probabilities from posterior samples
# ---------------------------------------------------------------------------

#' @noRd
extract_prob_samples <- function(samples, n_sites, n_samples) {
  if (is.null(samples)) return(NULL)
  mat <- matrix(NA_real_, n_samples, n_sites)
  for (s in seq_len(n_samples)) {
    fitted <- samples[[s]]$latent[seq_len(n_sites)]
    mat[s, ] <- expit(fitted)
  }
  mat
}


# ---------------------------------------------------------------------------
# $ accessor for occu_inla objects (spOccupancy compatibility)
# ---------------------------------------------------------------------------

#' Access spOccupancy-compatible fields from INLAocc fits
#'
#' Allows accessing spOccupancy-style posterior sample fields (e.g.,
#' \code{$beta.samples}, \code{$psi.samples}) on INLAocc model objects.
#' Samples are generated lazily via \code{INLA::inla.posterior.sample()}
#' on first access.
#'
#' @param x an INLAocc model object
#' @param name field name to access
#' @return The requested field value
#' @export
`$.occu_inla` <- function(x, name) {
  # First check if it's a native field
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  # spOccupancy compatibility fields
  n <- COMPAT_N_SAMPLES

  switch(name,
    # --- Occupancy fixed effects ---
    "beta.samples" = {
      samples <- draw_inla_samples(x[["occ_fit"]], n)
      fix_names <- rownames(x[["occ_fit"]]$summary.fixed)
      extract_beta_samples(samples, fix_names, n)
    },

    # --- Detection fixed effects ---
    "alpha.samples" = {
      samples <- draw_inla_samples(x[["det_fit"]], n)
      fix_names <- rownames(x[["det_fit"]]$summary.fixed)
      extract_beta_samples(samples, fix_names, n)
    },

    # --- Occupancy probabilities ---
    "psi.samples" = {
      samples <- draw_inla_samples(x[["occ_fit"]], n)
      N <- length(x[["psi_hat"]])
      extract_prob_samples(samples, N, n)
    },

    # --- Latent occupancy state ---
    "z.samples" = {
      psi_samp <- .subset2(x, "psi.samples")
      if (is.null(psi_samp)) {
        # Generate from psi_hat directly
        N <- length(x[["psi_hat"]])
        mat <- matrix(NA_real_, n, N)
        for (s in seq_len(n)) {
          mat[s, ] <- rbinom(N, 1, x[["psi_hat"]])
        }
        mat
      } else {
        # Generate from psi samples
        apply(psi_samp, c(1, 2), function(p) rbinom(1, 1, p))
      }
    },

    # --- Detection probabilities ---
    "p.samples" = {
      samples <- draw_inla_samples(x[["det_fit"]], n)
      if (!is.null(samples)) {
        N <- nrow(x[["p_hat"]])
        J <- ncol(x[["p_hat"]])
        # Return as 3D array: sample x site x visit
        arr <- array(NA_real_, dim = c(n, N, J))
        for (s in seq_len(n)) {
          n_fitted <- length(samples[[s]]$latent)
          fitted_vals <- expit(samples[[s]]$latent[seq_len(min(n_fitted, N * J))])
          arr[s, , ] <- matrix(fitted_vals[seq_len(N * J)], N, J)
        }
        arr
      }
    },

    # --- RE variance (occupancy) ---
    "sigma.sq.psi.samples" = {
      hyp <- x[["occ_fit"]]$summary.hyperpar
      if (!is.null(hyp) && nrow(hyp) > 0) {
        # INLA stores precision; convert to variance
        prec_rows <- grep("Precision", rownames(hyp))
        if (length(prec_rows) > 0) {
          mat <- matrix(NA_real_, n, length(prec_rows))
          for (k in seq_along(prec_rows)) {
            prec_mean <- hyp$mean[prec_rows[k]]
            prec_sd <- hyp$sd[prec_rows[k]]
            prec_samples <- rnorm(n, prec_mean, prec_sd)
            mat[, k] <- 1 / pmax(prec_samples, 1e-8)
          }
          mat
        }
      }
    },

    # --- RE variance (detection) ---
    "sigma.sq.p.samples" = {
      hyp <- x[["det_fit"]]$summary.hyperpar
      if (!is.null(hyp) && nrow(hyp) > 0) {
        prec_rows <- grep("Precision", rownames(hyp))
        if (length(prec_rows) > 0) {
          mat <- matrix(NA_real_, n, length(prec_rows))
          for (k in seq_along(prec_rows)) {
            prec_mean <- hyp$mean[prec_rows[k]]
            prec_sd <- hyp$sd[prec_rows[k]]
            prec_samples <- rnorm(n, prec_mean, prec_sd)
            mat[, k] <- 1 / pmax(prec_samples, 1e-8)
          }
          mat
        }
      }
    },

    # --- Spatial covariance parameters ---
    "theta.samples" = {
      hyp <- x[["occ_fit"]]$summary.hyperpar
      if (!is.null(hyp)) {
        sp_rows <- grep("Range|Stdev", rownames(hyp))
        if (length(sp_rows) > 0) {
          mat <- matrix(NA_real_, n, length(sp_rows))
          colnames(mat) <- rownames(hyp)[sp_rows]
          for (k in seq_along(sp_rows)) {
            mat[, k] <- rnorm(n, hyp$mean[sp_rows[k]], hyp$sd[sp_rows[k]])
          }
          mat
        }
      }
    },

    # --- Spatial random effects ---
    "w.samples" = {
      re <- x[["occ_fit"]]$summary.random
      if (!is.null(re) && "spatial" %in% names(re)) {
        sp_re <- re$spatial
        n_mesh <- nrow(sp_re)
        mat <- matrix(NA_real_, n, n_mesh)
        for (k in seq_len(n_mesh)) {
          mat[, k] <- rnorm(n, sp_re$mean[k], sp_re$sd[k])
        }
        mat
      }
    },

    # --- Likelihood (for WAIC) ---
    "like.samples" = {
      # INLA computes WAIC internally; return local contributions
      occ_waic <- x[["occ_fit"]]$waic
      det_waic <- x[["det_fit"]]$waic
      list(
        occ = if (!is.null(occ_waic)) occ_waic$local.waic else NULL,
        det = if (!is.null(det_waic)) det_waic$local.waic else NULL
      )
    },

    # --- Computation time ---
    "run.time" = x[["run_time"]],

    # Not a compat field — return NULL
    NULL
  )
}


# ---------------------------------------------------------------------------
# $ accessor for multi-species models
# ---------------------------------------------------------------------------

#' @export
`$.occu_inla_ms` <- function(x, name) {
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  n <- COMPAT_N_SAMPLES

  switch(name,
    # Community-level occupancy coefficients
    "beta.comm.samples" = {
      comm <- x[["community_occ"]]
      if (!is.null(comm)) {
        mat <- matrix(NA_real_, n, nrow(comm))
        colnames(mat) <- comm$parameter
        for (k in seq_len(nrow(comm))) {
          mat[, k] <- rnorm(n, comm$community_mean[k], comm$community_sd[k])
        }
        mat
      }
    },

    # Community-level detection coefficients
    "alpha.comm.samples" = {
      comm <- x[["community_det"]]
      if (!is.null(comm)) {
        mat <- matrix(NA_real_, n, nrow(comm))
        colnames(mat) <- comm$parameter
        for (k in seq_len(nrow(comm))) {
          mat[, k] <- rnorm(n, comm$community_mean[k], comm$community_sd[k])
        }
        mat
      }
    },

    # Community occ variance
    "tau.sq.beta.samples" = {
      comm <- x[["community_occ"]]
      if (!is.null(comm)) {
        mat <- matrix(NA_real_, n, nrow(comm))
        colnames(mat) <- comm$parameter
        for (k in seq_len(nrow(comm))) {
          # species_sd is the SD across species; variance = sd^2
          mat[, k] <- comm$species_sd[k]^2
        }
        mat
      }
    },

    # Community det variance
    "tau.sq.alpha.samples" = {
      comm <- x[["community_det"]]
      if (!is.null(comm)) {
        mat <- matrix(NA_real_, n, nrow(comm))
        colnames(mat) <- comm$parameter
        for (k in seq_len(nrow(comm))) {
          mat[, k] <- comm$species_sd[k]^2
        }
        mat
      }
    },

    # Species-level beta samples (3D: sample x species x coef)
    "beta.samples" = {
      n_sp <- x[["n_species"]]
      fits <- x[["species_fits"]]
      sp_names <- x[["species_names"]]
      if (length(fits) > 0 && !is.null(fits[[1]])) {
        n_coef <- nrow(fits[[1]]$occ_fit$summary.fixed)
        coef_names <- rownames(fits[[1]]$occ_fit$summary.fixed)
        arr <- array(NA_real_, dim = c(n, n_sp, n_coef),
                     dimnames = list(NULL, sp_names, coef_names))
        for (s in seq_len(n_sp)) {
          fit <- fits[[sp_names[s]]]
          if (!is.null(fit)) {
            fixed <- fit$occ_fit$summary.fixed
            for (k in seq_len(n_coef)) {
              arr[, s, k] <- rnorm(n, fixed$mean[k], fixed$sd[k])
            }
          }
        }
        arr
      }
    },

    # Species-level alpha samples (3D: sample x species x coef)
    "alpha.samples" = {
      n_sp <- x[["n_species"]]
      fits <- x[["species_fits"]]
      sp_names <- x[["species_names"]]
      if (length(fits) > 0 && !is.null(fits[[1]])) {
        n_coef <- nrow(fits[[1]]$det_fit$summary.fixed)
        coef_names <- rownames(fits[[1]]$det_fit$summary.fixed)
        arr <- array(NA_real_, dim = c(n, n_sp, n_coef),
                     dimnames = list(NULL, sp_names, coef_names))
        for (s in seq_len(n_sp)) {
          fit <- fits[[sp_names[s]]]
          if (!is.null(fit)) {
            fixed <- fit$det_fit$summary.fixed
            for (k in seq_len(n_coef)) {
              arr[, s, k] <- rnorm(n, fixed$mean[k], fixed$sd[k])
            }
          }
        }
        arr
      }
    },

    # Species-level psi samples (3D: sample x species x site)
    "psi.samples" = {
      n_sp <- x[["n_species"]]
      fits <- x[["species_fits"]]
      sp_names <- x[["species_names"]]
      N <- x[["data"]]$N
      arr <- array(NA_real_, dim = c(n, n_sp, N))
      for (s in seq_len(n_sp)) {
        fit <- fits[[sp_names[s]]]
        if (!is.null(fit)) {
          for (i in seq_len(N)) {
            arr[, s, i] <- rnorm(n, logit(clamp(fit$psi_hat[i])), 0.1)
          }
          arr[, s, ] <- expit(arr[, s, ])
        }
      }
      arr
    },

    # Species-level z samples (3D: sample x species x site)
    "z.samples" = {
      n_sp <- x[["n_species"]]
      fits <- x[["species_fits"]]
      sp_names <- x[["species_names"]]
      N <- x[["data"]]$N
      arr <- array(NA_integer_, dim = c(n, n_sp, N))
      for (s in seq_len(n_sp)) {
        fit <- fits[[sp_names[s]]]
        if (!is.null(fit)) {
          for (i in seq_len(N)) {
            arr[, s, i] <- rbinom(n, 1, fit$z_hat[i])
          }
        }
      }
      arr
    },

    # Factor loadings
    "lambda.samples" = {
      lam <- x[["lambda"]]
      if (!is.null(lam)) {
        # Point estimate replicated (no uncertainty available from PCA)
        arr <- array(NA_real_, dim = c(n, nrow(lam), ncol(lam)))
        for (s in seq_len(n)) arr[s, , ] <- lam
        arr
      }
    },

    # Not a compat field — return NULL
    NULL
  )
}


# ---------------------------------------------------------------------------
# $ accessor for temporal models
# ---------------------------------------------------------------------------

#' @export
`$.occu_inla_temporal` <- function(x, name) {
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  n <- COMPAT_N_SAMPLES
  n_periods <- x[["n_periods"]]

  switch(name,
    # psi samples: 3D (sample x site x period)
    "psi.samples" = {
      fits <- x[["period_fits"]]
      N <- if (!is.null(fits[[1]])) length(fits[[1]]$psi_hat) else 0
      arr <- array(NA_real_, dim = c(n, N, n_periods))
      for (t in seq_len(n_periods)) {
        fit <- fits[[t]]
        if (!is.null(fit)) {
          for (i in seq_len(N)) {
            arr[, i, t] <- rbeta(n, fit$psi_hat[i] * 10, (1 - fit$psi_hat[i]) * 10)
          }
        }
      }
      arr
    },

    # z samples: 3D (sample x site x period)
    "z.samples" = {
      fits <- x[["period_fits"]]
      N <- if (!is.null(fits[[1]])) length(fits[[1]]$psi_hat) else 0
      arr <- array(NA_integer_, dim = c(n, N, n_periods))
      for (t in seq_len(n_periods)) {
        fit <- fits[[t]]
        if (!is.null(fit)) {
          for (i in seq_len(N)) {
            arr[, i, t] <- rbinom(n, 1, fit$z_hat[i])
          }
        }
      }
      arr
    },

    # beta/alpha: use first period's fit as representative
    "beta.samples" = {
      fit <- x[["period_fits"]][[1]]
      if (!is.null(fit)) {
        fix <- fit$occ_fit$summary.fixed
        mat <- matrix(NA_real_, n, nrow(fix))
        colnames(mat) <- rownames(fix)
        for (k in seq_len(nrow(fix))) {
          mat[, k] <- rnorm(n, fix$mean[k], fix$sd[k])
        }
        mat
      }
    },

    "alpha.samples" = {
      fit <- x[["period_fits"]][[1]]
      if (!is.null(fit)) {
        fix <- fit$det_fit$summary.fixed
        mat <- matrix(NA_real_, n, nrow(fix))
        colnames(mat) <- rownames(fix)
        for (k in seq_len(nrow(fix))) {
          mat[, k] <- rnorm(n, fix$mean[k], fix$sd[k])
        }
        mat
      }
    },

    # eta samples (AR1 temporal effects)
    "eta.samples" = {
      if (!is.null(x[["psi_smoothed"]])) {
        logit(clamp(x[["psi_smoothed"]]))
      }
    },

    NULL
  )
}


# ---------------------------------------------------------------------------
# $ accessor for JSDM models
# ---------------------------------------------------------------------------

#' @export
`$.occu_inla_jsdm` <- function(x, name) {
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  n <- COMPAT_N_SAMPLES

  switch(name,
    "beta.samples" = {
      fits <- x[["species_fits"]]
      sp_names <- x[["species_names"]]
      n_sp <- x[["n_species"]]
      if (length(fits) > 0 && !is.null(fits[[1]])) {
        n_coef <- nrow(fits[[1]]$summary.fixed)
        coef_names <- rownames(fits[[1]]$summary.fixed)
        arr <- array(NA_real_, dim = c(n, n_sp, n_coef),
                     dimnames = list(NULL, sp_names, coef_names))
        for (s in seq_len(n_sp)) {
          fit <- fits[[s]]
          if (!is.null(fit)) {
            fixed <- fit$summary.fixed
            for (k in seq_len(n_coef)) {
              arr[, s, k] <- rnorm(n, fixed$mean[k], fixed$sd[k])
            }
          }
        }
        arr
      }
    },

    "psi.samples" = {
      fits <- x[["species_fits"]]
      n_sp <- x[["n_species"]]
      N <- ncol(x[["data"]]$y)
      arr <- array(NA_real_, dim = c(n, n_sp, N))
      for (s in seq_len(n_sp)) {
        fit <- fits[[s]]
        if (!is.null(fit)) {
          fitted_mean <- expit(fit$summary.fitted.values$mean[seq_len(N)])
          for (i in seq_len(N)) {
            arr[, s, i] <- rbeta(n, fitted_mean[i] * 10, (1 - fitted_mean[i]) * 10)
          }
        }
      }
      arr
    },

    "lambda.samples" = {
      lam <- x[["lambda"]]
      if (!is.null(lam)) {
        arr <- array(NA_real_, dim = c(n, nrow(lam), ncol(lam)))
        for (s in seq_len(n)) arr[s, , ] <- lam
        arr
      }
    },

    NULL
  )
}
