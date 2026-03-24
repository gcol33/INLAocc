# =============================================================================
# occu_output.R — Summary, print, and plot methods
# =============================================================================

#' Summary method for occu_inla fits
#'
#' @param object fitted occu_inla object
#' @param ... additional arguments (ignored)
#' @return invisibly returns a summary list
summary.occu_inla <- function(object, ...) {
  cat("=== Occupancy Model (INLA-Laplace) ===\n\n")

  cat(sprintf("Sites: %d | Max visits: %d\n", object$data$N, object$data$J))
  cat(sprintf("Naive occupancy: %.3f | Naive detection: %.3f\n",
              object$data$naive_occ,
              ifelse(is.na(object$data$naive_det), 0, object$data$naive_det)))
  cat(sprintf("EM iterations: %d | Converged: %s\n\n",
              object$n_iter, object$converged))

  # Occupancy fixed effects
  cat("--- Occupancy (psi) ---\n")
  occ_fixed <- object$occ_fit$summary.fixed
  occ_fixed$signif <- ifelse(
    sign(occ_fixed$`0.025quant`) == sign(occ_fixed$`0.975quant`), "*", ""
  )
  print(round(occ_fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
  cat("\n")

  # Detection fixed effects
  cat("--- Detection (p) ---\n")
  det_fixed <- object$det_fit$summary.fixed
  det_fixed$signif <- ifelse(
    sign(det_fixed$`0.025quant`) == sign(det_fixed$`0.975quant`), "*", ""
  )
  print(round(det_fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
  cat("\n")

  # Random effects
  if (length(object$occ_fit$summary.hyperpar) > 0) {
    cat("--- Hyperparameters (Occupancy) ---\n")
    hyp_occ <- object$occ_fit$summary.hyperpar
    if (nrow(hyp_occ) > 0) {
      print(round(hyp_occ[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
    }
    cat("\n")
  }

  if (length(object$det_fit$summary.hyperpar) > 0) {
    cat("--- Hyperparameters (Detection) ---\n")
    hyp_det <- object$det_fit$summary.hyperpar
    if (nrow(hyp_det) > 0) {
      print(round(hyp_det[, c("mean", "sd", "0.025quant", "0.975quant")], 4))
    }
    cat("\n")
  }

  # Model fit
  cat("--- Model Fit ---\n")
  cat(sprintf("Marginal log-likelihood: %.2f\n",
              tail(object$history, 1)[[1]]$loglik))
  if (!is.null(object$occ_fit$waic)) {
    cat(sprintf("WAIC (occupancy):  %.2f\n", object$occ_fit$waic$waic))
  }
  if (!is.null(object$det_fit$waic)) {
    cat(sprintf("WAIC (detection):  %.2f\n", object$det_fit$waic$waic))
  }

  # Derived quantities
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

  invisible(list(
    occ_fixed = occ_fixed,
    det_fixed = det_fixed,
    psi_hat   = object$psi_hat,
    p_hat     = object$p_hat,
    z_hat     = object$z_hat
  ))
}


#' Summary for spatial occupancy model
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
#' @param which integer vector: which plots to show
#'   1 = EM convergence, 2 = psi histogram, 3 = p histogram,
#'   4 = psi vs covariates, 5 = residuals
#' @param ... additional arguments passed to plot
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


#' Plot spatial occupancy map
#'
#' @param x fitted occu_inla_spatial object
#' @param what "psi" (occupancy), "spatial" (spatial field), or "z" (occupancy state)
#' @param n_grid grid resolution for spatial field
#' @param ... additional arguments to image()
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


#' Compare two occupancy models via WAIC
#'
#' @param ... named list of occu_inla objects
#' @return data.frame with model comparison metrics
compare_models <- function(...) {
  models <- list(...)
  if (is.null(names(models))) {
    names(models) <- paste0("model_", seq_along(models))
  }

  comp <- data.frame(
    model      = names(models),
    occ_waic   = NA_real_,
    det_waic   = NA_real_,
    total_waic = NA_real_,
    loglik     = NA_real_,
    n_iter     = NA_integer_,
    converged  = NA
  )

  for (i in seq_along(models)) {
    m <- models[[i]]
    if (inherits(m, "occu_em")) {
      if (!is.null(m$occ_fit$waic)) {
        comp$occ_waic[i] <- m$occ_fit$waic$waic
      }
      if (!is.null(m$det_fit$waic)) {
        comp$det_waic[i] <- m$det_fit$waic$waic
      }
      comp$total_waic[i] <- sum(comp$occ_waic[i], comp$det_waic[i],
                                 na.rm = TRUE)
      comp$loglik[i]   <- tail(m$history, 1)[[1]]$loglik
      comp$n_iter[i]   <- m$n_iter
      comp$converged[i] <- m$converged
    }
  }

  # Delta WAIC
  comp$delta_waic <- comp$total_waic - min(comp$total_waic, na.rm = TRUE)
  comp <- comp[order(comp$delta_waic), ]
  rownames(comp) <- NULL

  comp
}
