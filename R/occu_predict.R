# =============================================================================
# occu_predict.R — Prediction and derived quantities
# =============================================================================

#' Predict from a fitted occupancy model
#'
#' Supports three prediction modes:
#' \enumerate{
#'   \item \strong{terms-based} (like \code{ggpredict}): pass a character vector
#'     of covariate terms to vary, with optional bracket notation for ranges
#'     or specific values. Everything else is held at its mean (continuous)
#'     or mode (factor).
#'   \item \strong{design-matrix}: pass \code{X.0} / \code{X.p.0} directly
#'     (spOccupancy-compatible).
#'   \item \strong{in-sample}: call with no extra arguments.
#' }
#'
#' @param object fitted occu_inla object
#' @param terms character vector of terms to vary, with optional bracket
#'   notation. Examples: \code{"elev"}, \code{"elev [0:100]"},
#'   \code{"elev [0:100 by=5]"}, \code{"elev [1, 5, 10]"},
#'   \code{"habitat [forest, grassland]"}. The first term is the x-axis;
#'   the second (if any) defines groups; the third defines facets.
#' @param process \code{"occupancy"} (default) or \code{"detection"}.
#'   Only used with \code{terms}.
#' @param n_points number of prediction points per continuous term
#'   (default 50). Only used with \code{terms}.
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
#' @return If \code{terms} is provided, an \code{occu_prediction} data.frame
#'   with columns \code{x}, \code{estimate}, \code{lower}, \code{upper},
#'   and optional \code{group}/\code{facet}. Has a \code{plot()} method.
#'   Otherwise, a list with \code{psi.0}, \code{p.0}, \code{z.0}.
#'
#' @examples
#' \donttest{
#' # fit <- occu(~ elev + forest, ~ effort, data)
#' # predict(fit, terms = "elev")
#' # predict(fit, terms = "elev [0:100 by=10]")
#' # predict(fit, terms = c("elev", "forest"))
#' # predict(fit, terms = "effort", process = "detection")
#' }
#'
#' @export
predict.occu_inla <- function(object, X.0 = NULL, X.p.0 = NULL,
                              coords.0 = NULL,
                              ignore.RE = FALSE,
                              type = c("both", "occupancy", "detection"),
                              quantiles = c(0.025, 0.5, 0.975),
                              n.samples = 500,
                              terms = NULL,
                              process = c("occupancy", "detection"),
                              n_points = 50L, ...) {
  # --- terms-based prediction (ggpredict-style) ---
  if (!is.null(terms)) {
    process <- match.arg(process)
    return(predict_terms(object, terms, process, n_points, quantiles))
  }

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

      # Auto-add intercept column if the model has one and X.0 does not
      has_intercept <- "(Intercept)" %in% rownames(object$occ_fit$summary.fixed)
      if (has_intercept && ncol(X.0) == length(beta) - 1L) {
        X.0 <- cbind(1, X.0)
      }

      if (ncol(X.0) != length(beta)) {
        stop(sprintf(
          "X.0 has %d columns but model has %d occupancy coefficients (including intercept).",
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
#' @export
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
#' @export
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
#' @export
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


# =============================================================================
# Predict methods for non-base model classes
# =============================================================================

#' Predict from an integrated occupancy model
#'
#' Returns in-sample predictions. For out-of-sample, use the shared
#' occupancy component via the base \code{predict.occu_inla} method
#' on the \code{occ_fit} sub-object.
#'
#' @param object fitted occu_inla_int object
#' @param ... additional arguments (ignored)
#' @return list with psi.0 (shared occupancy) and per-source p.0
#' @export
predict.occu_inla_int <- function(object, ...) {
  list(
    psi.0 = list(mean = object$psi_hat),
    z.0   = list(mean = object$z_hat),
    p.0   = object$p_hats
  )
}


#' Predict from a multi-species occupancy model
#'
#' Returns per-species predictions by delegating to the base predict method.
#'
#' @param object fitted occu_inla_ms object
#' @param ... additional arguments passed to \code{predict.occu_inla}
#' @return named list of per-species prediction results
#' @export
predict.occu_inla_ms <- function(object, ...) {
  results <- list()
  for (sp in object$species_names) {
    fit <- object$species_fits[[sp]]
    if (!is.null(fit)) {
      results[[sp]] <- predict.occu_inla(fit, ...)
    }
  }
  results
}


#' Predict from a temporal occupancy model
#'
#' Returns per-period predictions. If AR(1) smoothing was applied,
#' includes the smoothed occupancy estimates.
#'
#' @param object fitted occu_inla_temporal object
#' @param ... additional arguments passed to \code{predict.occu_inla}
#' @return list with per-period predictions and optional smoothed psi
#' @export
predict.occu_inla_temporal <- function(object, ...) {
  period_preds <- vector("list", object$n_periods)
  for (t in seq_len(object$n_periods)) {
    fit <- object$period_fits[[t]]
    if (!is.null(fit)) {
      period_preds[[t]] <- predict.occu_inla(fit, ...)
    }
  }

  list(
    period_predictions = period_preds,
    psi_smoothed       = object$psi_smoothed
  )
}


#' Predict from a JSDM
#'
#' Returns per-species occupancy predictions (no detection process).
#'
#' @param object fitted occu_inla_jsdm object
#' @param ... additional arguments (ignored)
#' @return named list of per-species occupancy probability vectors
#' @export
predict.occu_inla_jsdm <- function(object, ...) {
  results <- list()
  for (s in seq_along(object$species_fits)) {
    fit <- object$species_fits[[s]]
    sp <- object$species_names[s]
    if (!is.null(fit)) {
      results[[sp]] <- list(
        psi.0 = list(mean = expit(fit$summary.fitted.values$mean))
      )
    }
  }
  results
}


# =============================================================================
# Terms-based prediction (ggpredict-style)
# =============================================================================

#' @noRd
parse_term <- function(term_str) {
  term_str <- trimws(term_str)
  bracket_pos <- regexpr("\\[", term_str)
  if (bracket_pos == -1L) {
    return(list(name = term_str, spec = NULL))
  }

  name <- trimws(substr(term_str, 1L, bracket_pos - 1L))
  content <- trimws(sub(".*\\[(.*)\\].*", "\\1", term_str))

  # Range: "0:100" or "0:100 by=5" or "-1:1"
  if (grepl("^-?[0-9.]+\\s*:\\s*-?[0-9.]+", content)) {
    by_val <- NULL
    range_str <- content
    by_match <- regexpr("by\\s*=\\s*[0-9.]+", content)
    if (by_match > 0L) {
      by_str <- regmatches(content, by_match)
      by_val <- as.numeric(sub("by\\s*=\\s*", "", by_str))
      range_str <- trimws(sub("by\\s*=\\s*[0-9.]+", "", content))
    }
    parts <- strsplit(range_str, "\\s*:\\s*")[[1]]
    return(list(name = name, spec = list(
      type = "range", from = as.numeric(parts[1]), to = as.numeric(parts[2]),
      by = by_val
    )))
  }

  # Comma-separated values or factor levels
  vals <- trimws(strsplit(content, ",")[[1]])
  nums <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(nums))) {
    return(list(name = name, spec = list(type = "values", values = nums)))
  }
  list(name = name, spec = list(type = "levels", levels = vals))
}

#' @noRd
resolve_term_values <- function(parsed, cov_vec, n_points = 50L) {
  if (is.null(parsed$spec)) {
    if (is.factor(cov_vec) || is.character(cov_vec)) {
      return(if (is.factor(cov_vec)) levels(cov_vec) else sort(unique(cov_vec)))
    }
    return(seq(min(cov_vec, na.rm = TRUE), max(cov_vec, na.rm = TRUE),
               length.out = n_points))
  }
  spec <- parsed$spec
  if (spec$type == "range") {
    if (!is.null(spec$by)) return(seq(spec$from, spec$to, by = spec$by))
    return(seq(spec$from, spec$to, length.out = n_points))
  }
  if (spec$type == "values") return(spec$values)
  spec$levels
}

#' @noRd
predict_terms <- function(object, terms, process, n_points, quantiles) {
  if (process == "occupancy") {
    fit  <- object$occ_fit
    covs <- object$data$occ.covs
    fml  <- object$occ_fixed_formula %||% object$occ.formula
  } else {
    fit  <- object$det_fit
    covs <- if (!is.null(object$det_df)) {
      object$det_df
    } else if (!is.null(object$data$det.covs)) {
      # Flatten det.covs list to data.frame
      as.data.frame(lapply(object$data$det.covs, function(m) {
        if (is.matrix(m)) as.vector(m) else m
      }))
    } else {
      data.frame(.intercept = 1)
    }
    fml <- object$det_fixed_formula %||% object$det.formula
  }

  if (is.null(fit)) stop("No fitted model for process '", process, "'")

  fixed <- fit$summary.fixed
  coef_names <- rownames(fixed)

  # Parse each term
  parsed <- lapply(terms, parse_term)
  term_names <- vapply(parsed, `[[`, character(1), "name")

  # Check all terms exist in covariates
  for (nm in term_names) {
    if (is.null(covs[[nm]])) {
      stop(sprintf("Term '%s' not found in %s covariates. Available: %s",
                   nm, process, paste(names(covs), collapse = ", ")))
    }
  }

  # Resolve values for each term
  term_values <- lapply(seq_along(parsed), function(i) {
    resolve_term_values(parsed[[i]], covs[[term_names[i]]], n_points)
  })
  names(term_values) <- term_names

  # Build prediction grid
  grid <- expand.grid(term_values, stringsAsFactors = TRUE)
  names(grid) <- term_names

  # Fill non-varied covariates at mean (numeric) or mode (factor/character)
  pred_data <- as.data.frame(grid)
  for (nm in names(covs)) {
    if (nm %in% term_names) next
    cv <- covs[[nm]]
    if (is.factor(cv) || is.character(cv)) {
      tab <- table(cv)
      mode_val <- names(tab)[which.max(tab)]
      pred_data[[nm]] <- if (is.factor(cv)) {
        factor(mode_val, levels = levels(cv))
      } else {
        mode_val
      }
    } else {
      pred_data[[nm]] <- mean(cv, na.rm = TRUE)
    }
  }

  # Build design matrix using the fixed-effects formula
  # Strip RE terms if present (parse_re_formula already did this for occ_fixed_formula)
  X <- tryCatch(
    model.matrix(fml, data = pred_data),
    error = function(e) {
      # Fallback: build manually from coefficient names
      X_manual <- matrix(0, nrow(pred_data), length(coef_names))
      colnames(X_manual) <- coef_names
      X_manual[, 1] <- 1  # intercept
      for (nm in term_names) {
        if (nm %in% coef_names) X_manual[, nm] <- pred_data[[nm]]
      }
      X_manual
    }
  )

  if (ncol(X) != length(coef_names)) {
    stop(sprintf(
      "Design matrix has %d columns but model has %d coefficients (%s)",
      ncol(X), length(coef_names), paste(coef_names, collapse = ", ")
    ))
  }

  # Predict on link scale (logit), then transform
  eta <- as.vector(X %*% fixed$mean)
  # Delta-method variance: diag(X V X') where V = diag(sd^2)
  eta_var <- rowSums((X %*% diag(fixed$sd^2, nrow = length(coef_names))) * X)
  eta_sd <- sqrt(pmax(eta_var, 0))

  z_lo <- qnorm(quantiles[1])
  z_hi <- qnorm(quantiles[length(quantiles)])

  result <- data.frame(
    x        = grid[[1]],
    estimate = expit(eta),
    lower    = expit(eta + z_lo * eta_sd),
    upper    = expit(eta + z_hi * eta_sd)
  )

  if (length(term_names) >= 2L) result$group <- grid[[2]]
  if (length(term_names) >= 3L) result$facet <- grid[[3]]

  attr(result, "terms")   <- term_names
  attr(result, "process") <- process
  class(result) <- c("occu_prediction", "data.frame")
  result
}


#' Plot predicted effects from an occupancy model
#'
#' @param x an \code{occu_prediction} object from
#'   \code{predict(fit, terms = ...)}
#' @param ... additional arguments passed to \code{\link[graphics]{plot}}
#'
#' @return The input object, invisibly. Called for side effect of producing
#'   a plot.
#' @export
plot.occu_prediction <- function(x, ...) {
  terms   <- attr(x, "terms")
  process <- attr(x, "process")
  x_label <- terms[1]
  y_label <- if (process == "occupancy") {
    "Occupancy probability"
  } else {
    "Detection probability"
  }

  has_group <- "group" %in% names(x)
  x_is_factor <- is.factor(x$x) || is.character(x$x)

  if (!has_group && !x_is_factor) {
    # Continuous ribbon plot
    plot(x$x, x$estimate, type = "l", lwd = 2,
         xlab = x_label, ylab = y_label, ylim = c(0, 1), ...)
    polygon(c(x$x, rev(x$x)), c(x$lower, rev(x$upper)),
            col = rgb(0, 0, 0, 0.15), border = NA)
    lines(x$x, x$estimate, lwd = 2)
  } else if (!has_group && x_is_factor) {
    # Factor: point + error bar
    lvls <- if (is.factor(x$x)) levels(x$x) else unique(x$x)
    idx <- seq_along(lvls)
    plot(idx, x$estimate, pch = 19, xlim = c(0.5, length(lvls) + 0.5),
         ylim = c(0, 1), xaxt = "n", xlab = x_label, ylab = y_label, ...)
    axis(1, at = idx, labels = lvls)
    arrows(idx, x$lower, idx, x$upper, angle = 90, code = 3, length = 0.05)
  } else {
    # Grouped: one line per group
    groups <- unique(x$group)
    n_groups <- length(groups)
    cols <- hcl.colors(n_groups, palette = "Dark 2")

    plot(range(x$x, na.rm = TRUE), c(0, 1), type = "n",
         xlab = x_label, ylab = y_label, ...)

    for (i in seq_along(groups)) {
      sel <- x$group == groups[i]
      xi <- x$x[sel]
      ord <- order(xi)
      xi <- xi[ord]
      est_i <- x$estimate[sel][ord]
      lo_i  <- x$lower[sel][ord]
      hi_i  <- x$upper[sel][ord]

      polygon(c(xi, rev(xi)), c(lo_i, rev(hi_i)),
              col = adjustcolor(cols[i], alpha.f = 0.15), border = NA)
      lines(xi, est_i, lwd = 2, col = cols[i])
    }

    legend("topright", legend = as.character(groups),
           col = cols, lwd = 2, title = terms[2], bty = "n")
  }

  invisible(x)
}
