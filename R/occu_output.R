# =============================================================================
# occu_output.R — Summary, print, and plot methods
# =============================================================================

#' Summary method for occu_inla fits
#'
#' @param object fitted occu_inla object
#' @param level character vector of parameter groups to display. Options:
#'   \code{"beta"} (occupancy fixed), \code{"alpha"} (detection fixed),
#'   \code{"sigma.sq.psi"} (occupancy hyperparams), \code{"sigma.sq.p"}
#'   (detection hyperparams). Default shows all available.
#' @param quantiles numeric vector of quantile levels for credible intervals
#'   (default: \code{c(0.025, 0.5, 0.975)}).
#' @param digits number of digits to print (default 4).
#' @param ... additional arguments (ignored)
#' @return Invisibly returns a summary list with occupancy and detection
#'   fixed effects, and estimated probabilities.
#' @export
summary.occu_inla <- function(object, level = NULL,
                               quantiles = c(0.025, 0.5, 0.975),
                               digits = 4L, ...) {
  # Determine which quantile columns exist in INLA output
  q_cols <- paste0(quantiles, "quant")
  available_q <- intersect(q_cols, names(object$occ_fit$summary.fixed))
  show_cols <- c("mean", "sd", available_q)

  show_all <- is.null(level)

  cat("=== Occupancy Model (INLA-Laplace) ===\n\n")
  cat(sprintf("Sites: %d | Max visits: %d\n", object$data$N, object$data$J))
  cat(sprintf("Naive occupancy: %.3f | Naive detection: %.3f\n",
              object$data$naive_occ,
              ifelse(is.na(object$data$naive_det), 0, object$data$naive_det)))
  cat(sprintf("EM iterations: %d | Converged: %s\n\n",
              object$n_iter, object$converged))

  occ_fixed <- object$occ_fit$summary.fixed
  det_fixed <- object$det_fit$summary.fixed

  if (show_all || "beta" %in% level) {
    cat("--- Occupancy (psi) ---\n")
    print(round(occ_fixed[, intersect(show_cols, names(occ_fixed)), drop = FALSE], digits))
    cat("\n")
  }

  if (show_all || "alpha" %in% level) {
    cat("--- Detection (p) ---\n")
    print(round(det_fixed[, intersect(show_cols, names(det_fixed)), drop = FALSE], digits))
    cat("\n")
  }

  if ((show_all || "sigma.sq.psi" %in% level) &&
      length(object$occ_fit$summary.hyperpar) > 0) {
    hyp_occ <- object$occ_fit$summary.hyperpar
    if (nrow(hyp_occ) > 0) {
      cat("--- Hyperparameters (Occupancy) ---\n")
      print(round(hyp_occ[, intersect(show_cols, names(hyp_occ)), drop = FALSE], digits))
      cat("\n")
    }
  }

  if ((show_all || "sigma.sq.p" %in% level) &&
      length(object$det_fit$summary.hyperpar) > 0) {
    hyp_det <- object$det_fit$summary.hyperpar
    if (nrow(hyp_det) > 0) {
      cat("--- Hyperparameters (Detection) ---\n")
      print(round(hyp_det[, intersect(show_cols, names(hyp_det)), drop = FALSE], digits))
      cat("\n")
    }
  }

  if (show_all) {
    cat("--- Model Fit ---\n")
    cat(sprintf("Marginal log-likelihood: %.2f\n",
                tail(object$history, 1)[[1]]$loglik))
    if (!is.null(object$occ_fit$waic))
      cat(sprintf("WAIC (occupancy):  %.2f\n", object$occ_fit$waic$waic))
    if (!is.null(object$det_fit$waic))
      cat(sprintf("WAIC (detection):  %.2f\n", object$det_fit$waic$waic))

    cat(sprintf("\nEstimated occupancy: %.3f (%.3f - %.3f)\n",
                mean(object$psi_hat),
                quantile(object$psi_hat, 0.025),
                quantile(object$psi_hat, 0.975)))
    cat(sprintf("Estimated detection: %.3f (%.3f - %.3f)\n",
                mean(object$p_hat, na.rm = TRUE),
                quantile(object$p_hat, 0.025, na.rm = TRUE),
                quantile(object$p_hat, 0.975, na.rm = TRUE)))
    cat(sprintf("Estimated occupied sites: %.1f / %d\n",
                sum(object$z_hat), object$data$N))
  }

  invisible(list(
    occ_fixed = occ_fixed,
    det_fixed = det_fixed,
    psi_hat   = object$psi_hat,
    p_hat     = object$p_hat,
    z_hat     = object$z_hat
  ))
}


#' Summary for spatial occupancy model
#' @param object fitted occu_inla_spatial object
#' @param ... additional arguments (ignored)
#' @return Invisibly returns the spatial summary via the base summary method.
#' @export
summary.occu_inla_spatial <- function(object, ...) {
  # Call base summary
  NextMethod()

  # Add spatial summary
  if (!is.null(object$spatial)) {
    cat("\n--- Spatial Component (SPDE) ---\n")
    cat(sprintf("Mesh nodes: %d\n", object$spatial$n_mesh))

    # Extract SPDE hyperparameters
    hyp <- object$occ_fit$summary.hyperpar
    sp_rows <- grep("Range|Stdev", rownames(hyp))
    if (length(sp_rows) > 0) {
      print(round(hyp[sp_rows, c("mean", "sd", "0.025quant", "0.975quant")], 4))
    }
  }
}


#' Summary for multi-species model
#' @param object fitted occu_inla_ms object
#' @param ... additional arguments (ignored)
#' @return The multi-species model object, returned invisibly.
#' @export
summary.occu_inla_ms <- function(object, ...) {
  cat("=== Multi-Species Occupancy Model (INLA-Laplace) ===\n\n")
  cat(sprintf("Species: %d | Sites: %d | Max visits: %d\n\n",
              object$n_species, object$data$N, object$data$J))

  # Per-species summary
  cat("--- Per-Species Occupancy ---\n")
  sp_table <- data.frame(
    species     = object$species_names,
    naive_psi   = NA_real_,
    est_psi     = NA_real_,
    est_p       = NA_real_,
    n_detected  = NA_integer_,
    converged   = NA
  )

  for (s in seq_along(object$species_names)) {
    sp <- object$species_names[s]
    fit <- object$species_fits[[sp]]
    dat <- object$data$species_data[[sp]]

    sp_table$naive_psi[s]  <- dat$naive_occ
    sp_table$n_detected[s] <- dat$n_occupied

    if (!is.null(fit)) {
      sp_table$est_psi[s]   <- mean(fit$psi_hat)
      sp_table$est_p[s]     <- mean(fit$p_hat, na.rm = TRUE)
      sp_table$converged[s] <- fit$converged
    }
  }
  print(sp_table, row.names = FALSE)

  # Community-level effects
  if (!is.null(object$community_occ)) {
    cat("\n--- Community-Level Occupancy Effects ---\n")
    print(round(object$community_occ, 4), row.names = FALSE)
  }

  if (!is.null(object$community_det)) {
    cat("\n--- Community-Level Detection Effects ---\n")
    print(round(object$community_det, 4), row.names = FALSE)
  }

  invisible(object)
}


#' Print method for occu_inla
#' @param x fitted occu_inla object
#' @param ... additional arguments (ignored)
#' @return The \code{occu_inla} object \code{x}, returned invisibly.
#' @export
print.occu_inla <- function(x, ...) {
  cat("Occupancy model (INLA-Laplace)\n")
  cat(sprintf("  %d sites, %d max visits\n", x$data$N, x$data$J))
  cat(sprintf("  Occ formula: %s\n", deparse(x$occ_formula)))
  cat(sprintf("  Det formula: %s\n", deparse(x$det_formula)))
  cat(sprintf("  EM: %d iterations, converged = %s\n",
              x$n_iter, x$converged))
  cat(sprintf("  Mean psi = %.3f, Mean p = %.3f\n",
              mean(x$psi_hat), mean(x$p_hat, na.rm = TRUE)))
  invisible(x)
}


#' Plot diagnostics for occu_inla
#'
#' @param x fitted occu_inla object
#' @param which integer vector: which plots to show.
#'   1 = EM convergence, 2 = psi histogram, 3 = p histogram,
#'   4 = psi vs covariates.
#' @param ... additional arguments passed to plot
#' @return The \code{occu_inla} object \code{x}, returned invisibly.
#' @export
plot.occu_inla <- function(x, which = 1:4, ...) {
  if (1 %in% which) {
    # EM convergence trace
    iters <- seq_along(x$history)
    ll <- vapply(x$history, function(h) h$loglik, numeric(1))

    plot(iters, ll, type = "b", pch = 16,
         xlab = "EM Iteration", ylab = "Log-Likelihood",
         main = "EM Convergence", ...)
  }

  if (2 %in% which) {
    # Occupancy probability distribution
    hist(x$psi_hat, breaks = 30, col = "steelblue", border = "white",
         main = "Estimated Occupancy Probabilities",
         xlab = expression(hat(psi)), ...)
    abline(v = mean(x$psi_hat), col = "red", lwd = 2, lty = 2)
  }

  if (3 %in% which) {
    # Detection probability distribution
    p_vec <- as.vector(x$p_hat)
    p_vec <- p_vec[!is.na(p_vec)]
    hist(p_vec, breaks = 30, col = "darkorange", border = "white",
         main = "Estimated Detection Probabilities",
         xlab = expression(hat(p)), ...)
    abline(v = mean(p_vec), col = "red", lwd = 2, lty = 2)
  }

  if (4 %in% which) {
    # z_hat (posterior occupancy) vs naive detection
    d <- x$data
    naive_det_rate <- rowSums(d$y, na.rm = TRUE) / rowSums(!is.na(d$y))

    plot(naive_det_rate, x$z_hat, pch = 16, col = rgb(0, 0, 0, 0.3),
         xlab = "Naive Detection Rate",
         ylab = expression(hat(z)),
         main = "Posterior Occupancy vs Naive Detection", ...)
    abline(h = 0.5, col = "red", lty = 2)
  }

  invisible(x)
}


#' @noRd
plot_spatial <- function(x, what = c("psi", "spatial", "z"),
                         n_grid = 100, ...) {
  what <- match.arg(what)

  if (what == "psi") {
    coords <- x$data$coords
    if (is.null(coords)) stop("No coordinates in model")

    # Scatter plot of estimated psi
    col_ramp <- colorRampPalette(c("white", "steelblue", "darkblue"))
    cols <- col_ramp(100)[cut(x$psi_hat, 100)]

    plot(coords[, 1], coords[, 2], pch = 16, col = cols,
         xlab = "X", ylab = "Y",
         main = expression("Estimated" ~ psi), ...)

  } else if (what == "spatial") {
    if (!inherits(x, "occu_inla_spatial")) {
      stop("Spatial field plot requires occu_inla_spatial object")
    }

    sf <- extract_spatial_field(x$occ_fit, x$spatial)
    grid <- project_spatial_grid(sf, x$spatial, n_grid = n_grid)

    image(grid$x, grid$y, grid$z_mean,
          col = hcl.colors(100, "Viridis"),
          xlab = "X", ylab = "Y",
          main = "Spatial Random Field", ...)
    points(x$data$coords, pch = "+", cex = 0.5)

  } else if (what == "z") {
    coords <- x$data$coords
    if (is.null(coords)) stop("No coordinates in model")

    cols <- ifelse(x$z_hat > 0.5, "steelblue", "grey80")
    plot(coords[, 1], coords[, 2], pch = 16, col = cols,
         xlab = "X", ylab = "Y",
         main = "Estimated Occupancy State (z > 0.5)", ...)
  }
}


# =============================================================================
# Integrated (multi-source) methods
# =============================================================================

#' Summary for integrated occupancy model
#' @param object fitted occu_inla_int object
#' @param ... additional arguments (ignored)
#' @return Invisibly returns a summary list.
#' @export
summary.occu_inla_int <- function(object, ...) {
  N_total <- nrow(object$data$occ.covs)
  cat("=== Integrated Occupancy Model (INLA-Laplace) ===\n\n")
  cat(sprintf("Data sources: %d | Total sites: %d\n", object$n_data, N_total))
  cat(sprintf("EM iterations: %d | Converged: %s\n\n",
              object$n_iter, object$converged))

  for (d in seq_len(object$n_data)) {
    cat(sprintf("  Source %d: %d sites, %d max visits\n",
                d, nrow(object$data$y[[d]]), ncol(object$data$y[[d]])))
  }
  cat("\n")

  # Shared occupancy fixed effects
  cat("--- Occupancy (psi, shared) ---\n")
  occ_fixed <- object$occ_fit$summary.fixed
  print(round(occ_fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
  cat("\n")

  # Per-source detection fixed effects
  for (d in seq_len(object$n_data)) {
    cat(sprintf("--- Detection (p, source %d) ---\n", d))
    if (!is.null(object$det_fits[[d]])) {
      det_fixed <- object$det_fits[[d]]$summary.fixed
      print(round(det_fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
    } else {
      cat("  (not available)\n")
    }
    cat("\n")
  }

  cat(sprintf("Estimated occupancy: %.3f (%.3f - %.3f)\n",
              mean(object$psi_hat),
              quantile(object$psi_hat, 0.025),
              quantile(object$psi_hat, 0.975)))
  cat(sprintf("Estimated occupied sites: %.1f / %d\n",
              sum(object$z_hat), N_total))

  invisible(object)
}


#' Print method for integrated occupancy model
#' @param x fitted occu_inla_int object
#' @param ... additional arguments (ignored)
#' @return The object \code{x}, returned invisibly.
#' @export
print.occu_inla_int <- function(x, ...) {
  N_total <- nrow(x$data$occ.covs)
  cat("Integrated occupancy model (INLA-Laplace)\n")
  cat(sprintf("  %d data sources, %d total sites\n", x$n_data, N_total))
  cat(sprintf("  Occ formula: %s\n", deparse(x$occ_formula)))
  cat(sprintf("  EM: %d iterations, converged = %s\n",
              x$n_iter, x$converged))
  cat(sprintf("  Mean psi = %.3f\n", mean(x$psi_hat)))
  invisible(x)
}


# =============================================================================
# Multi-species methods
# =============================================================================

#' Print method for multi-species occupancy model
#' @param x fitted occu_inla_ms object
#' @param ... additional arguments (ignored)
#' @return The object \code{x}, returned invisibly.
#' @export
print.occu_inla_ms <- function(x, ...) {
  cat("Multi-species occupancy model (INLA-Laplace)\n")
  cat(sprintf("  %d species, %d sites, %d max visits\n",
              x$n_species, x$data$N, x$data$J))
  cat(sprintf("  Occ formula: %s\n", deparse(x$occ.formula)))
  cat(sprintf("  Det formula: %s\n", deparse(x$det.formula)))
  n_ok <- sum(vapply(x$species_fits, function(f) !is.null(f), logical(1)))
  cat(sprintf("  Successfully fitted: %d / %d species\n", n_ok, x$n_species))
  invisible(x)
}


#' Plot diagnostics for multi-species occupancy model
#'
#' @param x fitted occu_inla_ms object
#' @param which integer vector: which plots to show.
#'   1 = per-species occupancy, 2 = per-species detection,
#'   3 = community effects.
#' @param ... additional arguments passed to plot
#' @return The object \code{x}, returned invisibly.
#' @export
plot.occu_inla_ms <- function(x, which = 1:2, ...) {
  psi_means <- vapply(x$species_names, function(sp) {
    fit <- x$species_fits[[sp]]
    if (!is.null(fit)) mean(fit$psi_hat) else NA_real_
  }, numeric(1))

  p_means <- vapply(x$species_names, function(sp) {
    fit <- x$species_fits[[sp]]
    if (!is.null(fit)) mean(fit$p_hat, na.rm = TRUE) else NA_real_
  }, numeric(1))

  if (1 %in% which) {
    ord <- order(psi_means, na.last = TRUE)
    barplot(psi_means[ord], names.arg = x$species_names[ord],
            las = 2, col = "steelblue", border = "white",
            main = "Estimated Occupancy by Species",
            ylab = expression(hat(psi)), ...)
  }

  if (2 %in% which) {
    ord <- order(p_means, na.last = TRUE)
    barplot(p_means[ord], names.arg = x$species_names[ord],
            las = 2, col = "darkorange", border = "white",
            main = "Estimated Detection by Species",
            ylab = expression(hat(p)), ...)
  }

  invisible(x)
}


# =============================================================================
# Temporal (multi-season) methods
# =============================================================================

#' Print method for temporal occupancy model
#' @param x fitted occu_inla_temporal object
#' @param ... additional arguments (ignored)
#' @return The object \code{x}, returned invisibly.
#' @export
print.occu_inla_temporal <- function(x, ...) {
  N <- if (length(x$data_list) > 0) x$data_list[[1]]$N else NA
  cat("Temporal occupancy model (INLA-Laplace)\n")
  cat(sprintf("  %d periods, %d sites\n", x$n_periods, N))
  cat(sprintf("  Occ formula: %s\n", deparse(x$occ.formula)))
  cat(sprintf("  Det formula: %s\n", deparse(x$det.formula)))
  cat(sprintf("  AR(1): %s\n", x$ar1))
  n_ok <- sum(vapply(x$period_fits, function(f) !is.null(f), logical(1)))
  cat(sprintf("  Successfully fitted: %d / %d periods\n", n_ok, x$n_periods))
  invisible(x)
}


#' Summary for temporal occupancy model
#' @param object fitted occu_inla_temporal object
#' @param ... additional arguments (ignored)
#' @return Invisibly returns a per-period summary data.frame.
#' @export
summary.occu_inla_temporal <- function(object, ...) {
  N <- if (length(object$data_list) > 0) object$data_list[[1]]$N else NA
  cat("=== Temporal Occupancy Model (INLA-Laplace) ===\n\n")
  cat(sprintf("Periods: %d | Sites: %d | AR(1): %s\n\n",
              object$n_periods, N, object$ar1))

  period_table <- data.frame(
    period    = seq_len(object$n_periods),
    est_psi   = NA_real_,
    est_p     = NA_real_,
    n_iter    = NA_integer_,
    converged = NA
  )

  for (t in seq_len(object$n_periods)) {
    fit <- object$period_fits[[t]]
    if (!is.null(fit)) {
      period_table$est_psi[t]   <- mean(fit$psi_hat)
      period_table$est_p[t]     <- mean(fit$p_hat, na.rm = TRUE)
      period_table$n_iter[t]    <- fit$n_iter
      period_table$converged[t] <- fit$converged
    }
  }

  cat("--- Per-Period Summary ---\n")
  print(period_table, row.names = FALSE)

  if (!is.null(object$psi_smoothed)) {
    cat(sprintf("\nAR(1)-smoothed occupancy range: %.3f - %.3f\n",
                min(object$psi_smoothed), max(object$psi_smoothed)))
  }

  invisible(period_table)
}


#' Plot diagnostics for temporal occupancy model
#'
#' @param x fitted occu_inla_temporal object
#' @param which integer vector: which plots to show.
#'   1 = occupancy over time, 2 = detection over time.
#' @param ... additional arguments passed to plot
#' @return The object \code{x}, returned invisibly.
#' @export
plot.occu_inla_temporal <- function(x, which = 1:2, ...) {
  periods <- seq_len(x$n_periods)

  psi_means <- vapply(periods, function(t) {
    fit <- x$period_fits[[t]]
    if (!is.null(fit)) mean(fit$psi_hat) else NA_real_
  }, numeric(1))

  p_means <- vapply(periods, function(t) {
    fit <- x$period_fits[[t]]
    if (!is.null(fit)) mean(fit$p_hat, na.rm = TRUE) else NA_real_
  }, numeric(1))

  if (1 %in% which) {
    plot(periods, psi_means, type = "b", pch = 16, col = "steelblue",
         xlab = "Period", ylab = expression(hat(psi)),
         main = "Occupancy Over Time", ylim = c(0, 1), ...)
    if (!is.null(x$psi_smoothed)) {
      smoothed_means <- colMeans(x$psi_smoothed)
      lines(periods, smoothed_means, col = "red", lty = 2, lwd = 2)
    }
  }

  if (2 %in% which) {
    plot(periods, p_means, type = "b", pch = 16, col = "darkorange",
         xlab = "Period", ylab = expression(hat(p)),
         main = "Detection Over Time", ylim = c(0, 1), ...)
  }

  invisible(x)
}


# =============================================================================
# JSDM methods
# =============================================================================

#' Print method for JSDM
#' @param x fitted occu_inla_jsdm object
#' @param ... additional arguments (ignored)
#' @return The object \code{x}, returned invisibly.
#' @export
print.occu_inla_jsdm <- function(x, ...) {
  cat("Joint species distribution model (INLA-Laplace)\n")
  cat(sprintf("  %d species, %d sites, %d latent factors\n",
              x$n_species, nrow(x$data$y[1, , drop = FALSE]), x$n.factors))
  cat(sprintf("  Formula: %s\n", deparse(x$formula)))
  n_ok <- sum(vapply(x$species_fits, function(f) !is.null(f), logical(1)))
  cat(sprintf("  Successfully fitted: %d / %d species\n", n_ok, x$n_species))
  invisible(x)
}


#' Summary for JSDM
#' @param object fitted occu_inla_jsdm object
#' @param ... additional arguments (ignored)
#' @return Invisibly returns a per-species summary data.frame.
#' @export
summary.occu_inla_jsdm <- function(object, ...) {
  cat("=== Joint Species Distribution Model (INLA-Laplace) ===\n\n")
  cat(sprintf("Species: %d | Latent factors: %d\n\n",
              object$n_species, object$n.factors))

  sp_table <- data.frame(
    species  = object$species_names,
    est_psi  = NA_real_
  )

  for (s in seq_along(object$species_names)) {
    fit <- object$species_fits[[s]]
    if (!is.null(fit)) {
      sp_table$est_psi[s] <- mean(expit(fit$summary.fitted.values$mean))
    }
  }

  cat("--- Per-Species Occupancy ---\n")
  print(sp_table, row.names = FALSE)

  cat("\n--- Factor Loadings (lambda) ---\n")
  lam <- object$lambda
  rownames(lam) <- object$species_names
  colnames(lam) <- paste0("F", seq_len(ncol(lam)))
  print(round(lam, 3))

  invisible(sp_table)
}


#' Plot diagnostics for JSDM
#'
#' @param x fitted occu_inla_jsdm object
#' @param which integer vector: which plots to show.
#'   1 = species occupancy, 2 = factor loadings heatmap.
#' @param ... additional arguments passed to plot
#' @return The object \code{x}, returned invisibly.
#' @export
plot.occu_inla_jsdm <- function(x, which = 1:2, ...) {
  if (1 %in% which) {
    psi_means <- vapply(seq_along(x$species_fits), function(s) {
      fit <- x$species_fits[[s]]
      if (!is.null(fit)) mean(expit(fit$summary.fitted.values$mean)) else NA_real_
    }, numeric(1))
    names(psi_means) <- x$species_names
    ord <- order(psi_means, na.last = TRUE)
    barplot(psi_means[ord], names.arg = x$species_names[ord],
            las = 2, col = "steelblue", border = "white",
            main = "Estimated Occupancy by Species",
            ylab = expression(hat(psi)), ...)
  }

  if (2 %in% which) {
    lam <- x$lambda
    image(seq_len(ncol(lam)), seq_len(nrow(lam)), t(lam),
          col = hcl.colors(100, "Blue-Red 3"),
          xlab = "Factor", ylab = "Species",
          main = "Factor Loadings", ...)
  }

  invisible(x)
}


# =============================================================================
# Standard S3 generics: coef, confint, vcov, nobs, update, tidy, ranef
# =============================================================================

#' Extract model coefficients
#'
#' @param object fitted occu_inla object
#' @param process \code{"occupancy"} (default), \code{"detection"}, or
#'   \code{"both"}.
#' @param ... ignored
#'
#' @return Named numeric vector of posterior means. If
#'   \code{process = "both"}, a named list with \code{occ} and \code{det}.
#' @export
coef.occu_inla <- function(object, process = c("both", "occupancy", "detection"),
                           ...) {
  process <- match.arg(process)
  occ <- setNames(object$occ_fit$summary.fixed$mean,
                  rownames(object$occ_fit$summary.fixed))
  det <- setNames(object$det_fit$summary.fixed$mean,
                  rownames(object$det_fit$summary.fixed))
  if (process == "occupancy") return(occ)
  if (process == "detection") return(det)
  list(occ = occ, det = det)
}


#' Compute credible intervals
#'
#' @param object fitted occu_inla object
#' @param parm not used (included for S3 compatibility)
#' @param level credible level (default 0.95)
#' @param process \code{"occupancy"} (default), \code{"detection"}, or
#'   \code{"both"}.
#' @param ... ignored
#'
#' @return Matrix with lower and upper columns, or a list of matrices
#'   if \code{process = "both"}.
#' @export
confint.occu_inla <- function(object, parm, level = 0.95,
                              process = c("both", "occupancy", "detection"),
                              ...) {
  process <- match.arg(process)
  alpha <- (1 - level) / 2
  q_lo <- paste0(alpha, "quant")
  q_hi <- paste0(1 - alpha, "quant")

  make_ci <- function(fixed) {
    lo <- if (q_lo %in% names(fixed)) fixed[[q_lo]] else {
      fixed$mean + qnorm(alpha) * fixed$sd
    }
    hi <- if (q_hi %in% names(fixed)) fixed[[q_hi]] else {
      fixed$mean + qnorm(1 - alpha) * fixed$sd
    }
    ci <- cbind(lo, hi)
    colnames(ci) <- paste0(c(alpha * 100, (1 - alpha) * 100), "%")
    rownames(ci) <- rownames(fixed)
    ci
  }

  occ_ci <- make_ci(object$occ_fit$summary.fixed)
  det_ci <- make_ci(object$det_fit$summary.fixed)

  if (process == "occupancy") return(occ_ci)
  if (process == "detection") return(det_ci)
  list(occ = occ_ci, det = det_ci)
}


#' Approximate variance-covariance matrix
#'
#' Returns the diagonal approximation (independent posteriors) from INLA.
#' For the full posterior covariance, use
#' \code{INLA::inla.posterior.sample()}.
#'
#' @param object fitted occu_inla object
#' @param process \code{"occupancy"} (default) or \code{"detection"}.
#' @param ... ignored
#'
#' @return A diagonal variance-covariance matrix.
#' @export
vcov.occu_inla <- function(object,
                           process = c("occupancy", "detection"), ...) {
  process <- match.arg(process)
  fixed <- if (process == "occupancy") {
    object$occ_fit$summary.fixed
  } else {
    object$det_fit$summary.fixed
  }
  v <- diag(fixed$sd^2, nrow = nrow(fixed))
  dimnames(v) <- list(rownames(fixed), rownames(fixed))
  v
}


#' Number of observations
#'
#' @param object fitted occu_inla object
#' @param ... ignored
#'
#' @return Integer: number of non-NA detection history entries.
#' @export
nobs.occu_inla <- function(object, ...) {
  sum(!is.na(object$data$y))
}


#' Update and refit an occupancy model
#'
#' Modify the formula or arguments of a fitted model and refit.
#'
#' @param object fitted occu_inla object
#' @param occ.formula new occupancy formula (default: keep existing).
#'   Supports \code{update.formula} syntax: \code{. ~ . - term}.
#' @param det.formula new detection formula (default: keep existing).
#' @param data new data (default: keep existing).
#' @param ... additional arguments passed to \code{\link{occu}}.
#' @param evaluate if \code{FALSE}, return the call instead of fitting.
#'
#' @return A new fitted model.
#' @export
update.occu_inla <- function(object, occ.formula = NULL, det.formula = NULL,
                             data = NULL, ..., evaluate = TRUE) {
  # Recover original call
  cl <- object$call
  if (is.null(cl)) stop("Model does not store its call; cannot update")

  # Update formulas

  if (!is.null(occ.formula)) {
    cl$occ.formula <- update.formula(object$occ.formula, occ.formula)
  }
  if (!is.null(det.formula)) {
    cl$det.formula <- update.formula(object$det.formula, det.formula)
  }
  if (!is.null(data)) {
    cl$data <- data
  }

  # Merge extra arguments
  dots <- list(...)
  for (nm in names(dots)) cl[[nm]] <- dots[[nm]]

  if (!evaluate) return(cl)
  eval(cl, parent.frame())
}


#' Tidy model output into a data.frame
#'
#' Returns a tidy data.frame of fixed effect estimates, similar to
#' \code{broom::tidy()}.
#'
#' @param x fitted occu_inla object
#' @param process \code{"occupancy"} (default), \code{"detection"}, or
#'   \code{"both"}.
#' @param conf.level credible level for intervals (default 0.95).
#' @param ... ignored
#'
#' @return A data.frame with columns \code{process}, \code{term},
#'   \code{estimate}, \code{std.error}, \code{conf.low}, \code{conf.high}.
#' @export
tidy.occu_inla <- function(x, process = c("both", "occupancy", "detection"),
                           conf.level = 0.95, ...) {
  process <- match.arg(process)
  alpha <- (1 - conf.level) / 2

  make_tidy <- function(fixed, proc_label) {
    lo <- fixed$mean + qnorm(alpha) * fixed$sd
    hi <- fixed$mean + qnorm(1 - alpha) * fixed$sd
    data.frame(
      process   = proc_label,
      term      = rownames(fixed),
      estimate  = fixed$mean,
      std.error = fixed$sd,
      conf.low  = lo,
      conf.high = hi,
      row.names = NULL
    )
  }

  parts <- list()
  if (process %in% c("both", "occupancy")) {
    parts <- c(parts, list(make_tidy(x$occ_fit$summary.fixed, "occupancy")))
  }
  if (process %in% c("both", "detection")) {
    parts <- c(parts, list(make_tidy(x$det_fit$summary.fixed, "detection")))
  }
  do.call(rbind, parts)
}


#' Extract random effects
#'
#' Returns the posterior summaries of random effect levels, similar to
#' \code{lme4::ranef()}.
#'
#' @param object fitted occu_inla object
#' @param process \code{"occupancy"} (default), \code{"detection"}, or
#'   \code{"both"}.
#' @param ... ignored
#'
#' @return A named list of data.frames (one per random effect group),
#'   each with columns \code{mean}, \code{sd}, \code{0.025quant},
#'   \code{0.975quant}. If \code{process = "both"}, a list with
#'   \code{occ} and \code{det} sub-lists.
#' @export
ranef.occu_inla <- function(object,
                            process = c("both", "occupancy", "detection"),
                            ...) {
  process <- match.arg(process)

  extract_re <- function(fit, prefix) {
    re <- fit$summary.random
    if (is.null(re) || length(re) == 0) return(list())
    # Strip internal prefix from names
    out <- re
    names(out) <- sub(paste0("^", prefix, "_re_"), "", names(out))
    out
  }

  occ_re <- extract_re(object$occ_fit, "occ")
  det_re <- extract_re(object$det_fit, "det")

  if (process == "occupancy") return(occ_re)
  if (process == "detection") return(det_re)
  list(occ = occ_re, det = det_re)
}


#' Extract log-likelihood
#'
#' Returns the observed-data log-likelihood from the final EM iteration.
#' The \code{df} attribute is set to the number of fixed effect
#' coefficients (occupancy + detection) and the \code{nobs} attribute
#' to the number of non-NA detection history entries.
#'
#' @param object fitted occu_inla object
#' @param ... ignored
#'
#' @return An object of class \code{"logLik"}.
#' @export
logLik.occu_inla <- function(object, ...) {
  ll <- NA_real_
  if (!is.null(object$history) && length(object$history) > 0) {
    last <- object$history[[length(object$history)]]
    if (!is.null(last$loglik)) ll <- last$loglik
  }

  # Compute on demand if not cached (e.g. verbose = 0)
  if (is.na(ll) && !is.null(object$psi_hat) && !is.null(object$p_hat)) {
    ll <- occu_loglik(object$data$y, object$psi_hat, object$p_hat)
  }

  n_occ <- nrow(object$occ_fit$summary.fixed)
  n_det <- nrow(object$det_fit$summary.fixed)
  df <- n_occ + n_det

  structure(ll, df = df, nobs = nobs.occu_inla(object),
            class = "logLik")
}


#' Glance at model-level statistics
#'
#' Returns a single-row data.frame of model-level summaries, similar to
#' \code{broom::glance()}.
#'
#' @param x fitted occu_inla object
#' @param ... ignored
#'
#' @return A one-row data.frame with columns \code{nobs}, \code{n_sites},
#'   \code{n_visits}, \code{n_occ_coef}, \code{n_det_coef}, \code{logLik},
#'   \code{WAIC}, \code{n_iter}, \code{converged}.
#' @export
glance.occu_inla <- function(x, ...) {
  ll <- NA_real_
  if (!is.null(x$history) && length(x$history) > 0) {
    last <- x$history[[length(x$history)]]
    if (!is.null(last$loglik)) ll <- last$loglik
  }

  waic_occ <- if (!is.null(x$occ_fit$waic)) x$occ_fit$waic$waic else NA_real_
  waic_det <- if (!is.null(x$det_fit$waic)) x$det_fit$waic$waic else NA_real_
  waic_total <- sum(waic_occ, waic_det, na.rm = TRUE)
  if (is.na(waic_occ) && is.na(waic_det)) waic_total <- NA_real_

  data.frame(
    nobs       = sum(!is.na(x$data$y)),
    n_sites    = x$data$N,
    n_visits   = x$data$J,
    n_occ_coef = nrow(x$occ_fit$summary.fixed),
    n_det_coef = nrow(x$det_fit$summary.fixed),
    logLik     = ll,
    WAIC       = waic_total,
    n_iter     = x$n_iter %||% NA_integer_,
    converged  = x$converged %||% NA,
    row.names  = NULL
  )
}


#' Compare two occupancy models via WAIC
#'
#' @param ... named occu_inla objects to compare
#' @return data.frame with model comparison metrics
#' @export
compare_models <- function(..., criterion = c("waic", "aic", "bic")) {
  criterion <- match.arg(criterion)
  models <- list(...)
  if (is.null(names(models))) {
    names(models) <- paste0("model_", seq_along(models))
  }
  n_models <- length(models)

  comp <- data.frame(
    model     = names(models),
    loglik    = NA_real_,
    df        = NA_integer_,
    AIC       = NA_real_,
    BIC       = NA_real_,
    WAIC      = NA_real_,
    n_iter    = NA_integer_,
    converged = NA
  )

  for (i in seq_len(n_models)) {
    m <- models[[i]]
    if (!inherits(m, "occu_em")) next

    ll <- logLik(m)
    comp$loglik[i] <- as.numeric(ll)
    comp$df[i]     <- attr(ll, "df")
    n_obs          <- attr(ll, "nobs")
    comp$AIC[i]    <- -2 * as.numeric(ll) + 2 * comp$df[i]
    comp$BIC[i]    <- -2 * as.numeric(ll) + log(n_obs) * comp$df[i]

    occ_w <- if (!is.null(m$occ_fit$waic)) m$occ_fit$waic$waic else NA_real_
    det_w <- if (!is.null(m$det_fit$waic)) m$det_fit$waic$waic else NA_real_
    comp$WAIC[i] <- sum(occ_w, det_w, na.rm = TRUE)
    if (is.na(occ_w) && is.na(det_w)) comp$WAIC[i] <- NA_real_

    comp$n_iter[i]    <- m$n_iter
    comp$converged[i] <- m$converged
  }

  # Choose IC for ranking and weights
  ic_col <- switch(criterion, waic = "WAIC", aic = "AIC", bic = "BIC")
  ic_vals <- comp[[ic_col]]

  comp$delta <- ic_vals - min(ic_vals, na.rm = TRUE)
  comp$weight <- exp(-0.5 * comp$delta) / sum(exp(-0.5 * comp$delta),
                                                na.rm = TRUE)

  comp <- comp[order(comp$delta), ]
  rownames(comp) <- NULL
  attr(comp, "criterion") <- criterion
  comp
}


# ---------------------------------------------------------------------------
#  modelAverage  —  multi-model inference
# ---------------------------------------------------------------------------

#' Model-averaged predictions from multiple occupancy models
#'
#' Computes weighted-average occupancy and detection probabilities across
#' a candidate set of models.  Weights are derived from AIC, BIC, or WAIC
#' via \code{\link{compare_models}}.
#'
#' Follows the "full model averaging" approach of Burnham & Anderson (2002):
#' predictions from every model contribute, weighted by information-criterion
#' weights.
#'
#' @param ... named \code{occu_inla} objects (the candidate model set)
#' @param criterion \code{"waic"} (default), \code{"aic"}, or \code{"bic"}
#' @param newdata optional list for out-of-sample prediction (passed to
#'   \code{\link{predict.occu_inla}}).  If \code{NULL}, returns in-sample
#'   averaged psi and p.
#' @param se if \code{TRUE} (default), returns unconditional standard errors
#'   that account for model selection uncertainty
#'
#' @return A list of class \code{"occu_model_avg"} with:
#'   \describe{
#'     \item{psi_hat}{model-averaged occupancy probabilities}
#'     \item{p_hat}{model-averaged detection probabilities (N x J)}
#'     \item{psi_se}{unconditional SE for psi (if \code{se = TRUE})}
#'     \item{weights}{named vector of model weights}
#'     \item{comparison}{data.frame from \code{compare_models()}}
#'     \item{criterion}{IC used for weights}
#'   }
#' @export
modelAverage <- function(..., criterion = c("waic", "aic", "bic"),
                          newdata = NULL, se = TRUE) {
  criterion <- match.arg(criterion)
  models <- list(...)
  if (is.null(names(models))) {
    names(models) <- paste0("model_", seq_along(models))
  }

  comp <- do.call(compare_models, c(models, list(criterion = criterion)))
  w <- setNames(comp$weight, comp$model)
  # Reorder weights to match input model order
  w <- w[names(models)]

  N <- models[[1]]$data$N

  # --- Weighted average of psi and p ---
  psi_avg <- rep(0, N)
  psi_var <- rep(0, N)   # for unconditional SE

  # Detection dimensions can vary across models but J is fixed
  J <- models[[1]]$data$J
  p_avg <- matrix(0, N, J)

  for (nm in names(models)) {
    m <- models[[nm]]
    wi <- w[nm]
    if (is.na(wi) || wi == 0) next

    psi_i <- m$psi_hat
    psi_avg <- psi_avg + wi * psi_i

    p_i <- m$p_hat
    if (!is.null(p_i) && identical(dim(p_i), c(N, J))) {
      p_avg <- p_avg + wi * p_i
    }
  }

  # Unconditional SE (Burnham & Anderson eq. 4.9):
  # var_uncond = sum w_i * (var_i + (theta_i - theta_bar)^2)
  psi_se <- NULL
  if (se) {
    for (nm in names(models)) {
      m <- models[[nm]]
      wi <- w[nm]
      if (is.na(wi) || wi == 0) next

      psi_i <- m$psi_hat
      # Approximate within-model variance from Bernoulli: psi*(1-psi)/nobs
      # Better: use posterior SD if available from INLA
      if (!is.null(m$occ_fit$summary.fitted.values)) {
        sd_i <- m$occ_fit$summary.fitted.values$sd[seq_len(N)]
        if (length(sd_i) == N) {
          var_i <- sd_i^2
        } else {
          var_i <- psi_i * (1 - psi_i)  # fallback
        }
      } else {
        var_i <- psi_i * (1 - psi_i)
      }
      psi_var <- psi_var + wi * (var_i + (psi_i - psi_avg)^2)
    }
    psi_se <- sqrt(psi_var)
  }

  out <- list(
    psi_hat    = psi_avg,
    p_hat      = p_avg,
    psi_se     = psi_se,
    weights    = w,
    comparison = comp,
    criterion  = criterion
  )
  class(out) <- "occu_model_avg"
  out
}


#' @export
print.occu_model_avg <- function(x, ...) {
  cat("Model-averaged occupancy predictions\n")
  cat(sprintf("  Criterion: %s | Models: %d\n", toupper(x$criterion),
              length(x$weights)))
  cat(sprintf("  psi range: [%.3f, %.3f] mean=%.3f\n",
              min(x$psi_hat), max(x$psi_hat), mean(x$psi_hat)))
  if (!is.null(x$psi_se)) {
    cat(sprintf("  Unconditional SE range: [%.3f, %.3f]\n",
                min(x$psi_se), max(x$psi_se)))
  }
  cat("\nModel weights:\n")
  w <- sort(x$weights, decreasing = TRUE)
  for (nm in names(w)) {
    cat(sprintf("  %-20s %.3f\n", nm, w[nm]))
  }
  invisible(x)
}
