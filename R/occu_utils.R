# =============================================================================
# occu_utils.R — Internal utilities for INLA occupancy framework
# =============================================================================

#' @noRd
expit <- function(x) {
  1 / (1 + exp(-x))
}

#' @noRd
logit <- function(p) {
  log(p / (1 - p))
}

#' @noRd
clamp <- function(x, lo = 1e-8, hi = 1 - 1e-8) {
  pmin(pmax(x, lo), hi)
}

#' @noRd
occu_loglik <- function(y, psi, p) {
  N <- nrow(y)
  J <- ncol(y)
  ll <- 0

  for (i in seq_len(N)) {
    visits <- which(!is.na(y[i, ]))
    if (length(visits) == 0) next

    yi <- y[i, visits]
    pi <- p[i, visits]
    detected <- any(yi == 1)

    if (detected) {
      # z_i = 1 certain
      ll_det <- sum(yi * log(clamp(pi)) + (1 - yi) * log(clamp(1 - pi)))
      ll <- ll + log(clamp(psi[i])) + ll_det
    } else {
      # z_i unknown: marginalize
      q_i <- prod(clamp(1 - pi))
      ll <- ll + log(clamp((1 - psi[i]) + psi[i] * q_i))
    }
  }
  ll
}

#' @noRd
compute_weights <- function(y, psi, p) {
  N <- nrow(y)
  w <- rep(1, N)

  for (i in seq_len(N)) {
    visits <- which(!is.na(y[i, ]))
    if (length(visits) == 0) {
      w[i] <- 0.5
      next
    }

    yi <- y[i, visits]
    detected <- any(yi == 1)

    if (!detected) {
      pi <- p[i, visits]
      q_i <- prod(clamp(1 - pi))  # P(all zeros | occupied)
      psi_i <- clamp(psi[i])
      # Bayes: P(z=1 | all zeros) = psi*q / ((1-psi) + psi*q)
      w[i] <- (psi_i * q_i) / ((1 - psi_i) + psi_i * q_i)
    }
  }
  w
}

#' @noRd
build_det_df <- function(y, det_covs = NULL, site_idx = NULL,
                         weights = NULL, visit_idx = NULL) {
  N <- nrow(y)
  J <- ncol(y)

  if (is.null(weights)) weights <- rep(1, N)
  if (is.null(site_idx)) site_idx <- seq_len(N)

  rows <- list()
  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      if (is.na(y[i, j])) next

      row <- list(
        y_det    = y[i, j],
        site     = site_idx[i],
        visit    = j,
        w        = weights[i]
      )

      if (!is.null(det_covs)) {
        for (nm in names(det_covs)) {
          row[[nm]] <- det_covs[[nm]][i, j]
        }
      }

      if (!is.null(visit_idx)) {
        row[["visit_group"]] <- visit_idx[i, j]
      }

      rows[[length(rows) + 1]] <- row
    }
  }

  do.call(rbind, lapply(rows, as.data.frame))
}

#' @noRd
build_occ_df <- function(occ_covs, weights, site_idx = NULL,
                         detected = NULL, M = 1000L) {
  N <- nrow(occ_covs)
  if (is.null(site_idx)) site_idx <- seq_len(N)
  if (is.null(detected)) detected <- weights == 1

  df <- occ_covs
  df$site <- site_idx

  # Encode E-step weights as scaled binomial counts (one row per site)
  df$z       <- ifelse(detected, M, round(weights * M))
  df$Ntrials <- rep(M, N)
  df$w       <- rep(1, N)  # uniform weight — information is in z/Ntrials

  df
}

#' @noRd
parse_re_formula <- function(formula) {
  has_lme4 <- requireNamespace("lme4", quietly = TRUE)

  if (has_lme4) {
    # --- Use lme4's parser (gold standard) ---
    bars <- lme4::findbars(formula)
    fixed <- lme4::nobars(formula)
    if (is.null(fixed) || length(fixed) == 2 && is.null(fixed[[2]])) {
      fixed <- ~ 1
    }

    re_list <- list()
    if (!is.null(bars)) {
      for (bar in bars) {
        # bar is a call like `1 | group` or `x | group`
        grp <- deparse(bar[[3]])  # RHS of |
        lhs <- deparse(bar[[2]])  # LHS of |

        # lme4 parses (1 + x | group) into a single bar with lhs = "1 + x"
        # Split into individual terms
        lhs_terms <- trimws(strsplit(lhs, "\\+")[[1]])

        for (term in lhs_terms) {
          if (term == "1") {
            re_list[[length(re_list) + 1]] <- occu_re("intercept", group = grp)
          } else if (term == "0") {
            next  # (0 + x | group) means no intercept, skip
          } else {
            re_list[[length(re_list) + 1]] <- occu_re("slope", group = grp,
                                                        covariate = term)
          }
        }
      }
    }
  } else {
    # --- Fallback: regex parser ---
    formula_str <- paste(deparse(formula, width.cutoff = 500), collapse = "")

    re_pattern <- "\\(\\s*([^|]+)\\s*\\|\\s*([^)]+)\\s*\\)"
    re_matches <- gregexpr(re_pattern, formula_str, perl = TRUE)
    re_strings <- regmatches(formula_str, re_matches)[[1]]

    re_list <- list()
    if (length(re_strings) > 0) {
      for (rs in re_strings) {
        inner <- sub("^\\(\\s*", "", sub("\\s*\\)$", "", rs))
        parts <- strsplit(inner, "\\s*\\|\\s*")[[1]]
        lhs_terms <- trimws(strsplit(parts[1], "\\+")[[1]])
        group <- trimws(parts[2])

        for (term in lhs_terms) {
          if (term == "1") {
            re_list[[length(re_list) + 1]] <- occu_re("intercept", group = group)
          } else if (term != "0") {
            re_list[[length(re_list) + 1]] <- occu_re("slope", group = group,
                                                        covariate = term)
          }
        }
      }

      # Strip RE terms from formula
      fixed_str <- formula_str
      for (rs in re_strings) {
        fixed_str <- sub(paste0("\\+?\\s*", gsub("([.|()^{}])", "\\\\\\1", rs)),
                         "", fixed_str, perl = TRUE)
      }
      fixed_str <- gsub("~\\s*\\+", "~", fixed_str)
      fixed_str <- gsub("\\+\\s*$", "", fixed_str)
      fixed_str <- trimws(fixed_str)
      if (grepl("~\\s*$", fixed_str)) fixed_str <- paste0(fixed_str, " 1")
    } else {
      fixed_str <- formula_str
    }
    fixed <- as.formula(fixed_str)
  }

  list(
    fixed   = fixed,
    re_list = re_list
  )
}


#' @noRd
parse_occu_formula <- function(occ_formula, det_formula) {
  occ_parsed <- parse_re_formula(occ_formula)
  det_parsed <- parse_re_formula(det_formula)

  list(
    occ       = occ_parsed$fixed,
    det       = det_parsed$fixed,
    occ_terms = attr(terms(occ_parsed$fixed), "term.labels"),
    det_terms = attr(terms(det_parsed$fixed), "term.labels"),
    occ_re    = occ_parsed$re_list,
    det_re    = det_parsed$re_list
  )
}

#' @noRd
inla_fitted_prob <- function(fit, link = "logit") {
  eta <- fit$summary.fitted.values$mean
  if (link == "logit") {
    expit(eta)
  } else if (link == "probit") {
    pnorm(eta)
  } else {
    eta
  }
}

#' Specify priors for INLA occupancy models (spOccupancy-compatible)
#'
#' @param beta.normal list(mean, var) for occupancy fixed effects Normal prior
#' @param alpha.normal list(mean, var) for detection fixed effects Normal prior
#' @param sigma.sq.psi.ig c(shape, scale) for occupancy RE variance Inv-Gamma
#' @param sigma.sq.p.ig c(shape, scale) for detection RE variance Inv-Gamma
#' @param sigma.sq.ig c(shape, scale) for spatial variance (spatial models)
#' @param phi.unif c(lower, upper) for spatial decay Uniform prior
#' @param beta.comm.normal list(mean, var) community occ prior (multi-species)
#' @param alpha.comm.normal list(mean, var) community det prior (multi-species)
#' @param tau.sq.beta.ig list(a, b) species-level occ variance (multi-species)
#' @param tau.sq.alpha.ig list(a, b) species-level det variance (multi-species)
#'
#' @return list of class \code{"occu_priors"} with INLA-compatible prior specs
#' @export
occu_priors <- function(beta.normal = list(mean = 0, var = 2.72),
                        alpha.normal = list(mean = 0, var = 2.72),
                        sigma.sq.psi.ig = c(0.1, 0.1),
                        sigma.sq.p.ig = c(0.1, 0.1),
                        sigma.sq.ig = NULL,
                        phi.unif = NULL,
                        beta.comm.normal = NULL,
                        alpha.comm.normal = NULL,
                        tau.sq.beta.ig = NULL,
                        tau.sq.alpha.ig = NULL) {

  # Convert to INLA prior specifications
  # INLA uses precision, not variance: prec = 1/var
  occ_fixed_prec <- 1 / beta.normal$var
  det_fixed_prec <- 1 / alpha.normal$var

  out <- list(
    # Fixed effects: Gaussian with specified precision
    occ_fixed = list(
      mean     = beta.normal$mean,
      prec     = occ_fixed_prec,
      mean.intercept = beta.normal$mean,
      prec.intercept = occ_fixed_prec
    ),
    det_fixed = list(
      mean     = alpha.normal$mean,
      prec     = det_fixed_prec,
      mean.intercept = alpha.normal$mean,
      prec.intercept = det_fixed_prec
    ),
    # Random effect variance: convert IG(shape, scale) to PC prior for INLA
    # Approximate: PC prior P(sigma > sqrt(scale/shape)) = 0.05
    occ_re_hyper = list(
      prec = list(prior = "pc.prec",
                  param = c(sqrt(sigma.sq.psi.ig[2] / sigma.sq.psi.ig[1]), 0.05))
    ),
    det_re_hyper = list(
      prec = list(prior = "pc.prec",
                  param = c(sqrt(sigma.sq.p.ig[2] / sigma.sq.p.ig[1]), 0.05))
    ),
    # Spatial (if provided)
    spatial = if (!is.null(sigma.sq.ig)) {
      list(
        sigma.sq = sigma.sq.ig,
        phi      = phi.unif
      )
    },
    # Community (multi-species)
    community = if (!is.null(beta.comm.normal)) {
      list(
        beta.comm  = beta.comm.normal,
        alpha.comm = alpha.comm.normal,
        tau.sq.beta = tau.sq.beta.ig,
        tau.sq.alpha = tau.sq.alpha.ig
      )
    }
  )
  class(out) <- "occu_priors"
  out
}


#' @noRd
priors_to_control_fixed <- function(priors, process = c("occ", "det")) {
  process <- match.arg(process)
  if (is.null(priors)) return(list())

  spec <- if (process == "occ") priors$occ_fixed else priors$det_fixed
  list(
    mean            = spec$mean,
    prec            = spec$prec,
    mean.intercept  = spec$mean.intercept,
    prec.intercept  = spec$prec.intercept
  )
}


#' @noRd
kfold_occu <- function(fit_fun, data, k = 5, ...) {
  N <- data$N
  fold_ids <- sample(rep(1:k, length.out = N))

  fold_deviance <- numeric(k)
  fold_results <- list()

  for (f in seq_len(k)) {
    test_sites  <- which(fold_ids == f)
    train_sites <- which(fold_ids != f)

    # Split data
    train_data <- data
    train_data$y         <- data$y[train_sites, , drop = FALSE]
    train_data$occ.covs  <- data$occ.covs[train_sites, , drop = FALSE]
    train_data$N         <- length(train_sites)
    train_data$detected  <- train_data$n_det[train_sites] > 0
    train_data$n_visits  <- data$n_visits[train_sites]
    train_data$n_det     <- data$n_det[train_sites]
    train_data$n_occupied <- sum(train_data$detected)
    train_data$naive_occ <- mean(train_data$detected)
    train_data$naive_det <- if (sum(train_data$n_visits[train_data$detected]) > 0) {
      sum(train_data$n_det) / sum(train_data$n_visits[train_data$detected])
    } else NA_real_

    if (!is.null(data$det.covs)) {
      train_data$det.covs <- lapply(data$det.covs, function(m) {
        m[train_sites, , drop = FALSE]
      })
    }
    if (!is.null(data$coords)) {
      train_data$coords <- data$coords[train_sites, , drop = FALSE]
    }

    # Fit on training data
    fold_fit <- tryCatch({
      fit_fun(data = train_data, verbose = 0, ...)
    }, error = function(e) {
      warning(sprintf("Fold %d failed: %s", f, e$message))
      NULL
    })

    if (!is.null(fold_fit)) {
      # Evaluate on test data
      test_y   <- data$y[test_sites, , drop = FALSE]
      test_psi <- fold_fit$psi_hat[test_sites]
      test_p   <- fold_fit$p_hat[test_sites, , drop = FALSE]

      # Clamp for numerical safety
      test_psi <- clamp(test_psi)
      test_p   <- clamp(test_p)

      fold_deviance[f] <- -2 * occu_loglik(test_y, test_psi, test_p)
    } else {
      fold_deviance[f] <- NA
    }

    fold_results[[f]] <- list(
      test_sites = test_sites,
      deviance   = fold_deviance[f]
    )
  }

  list(
    k.fold.deviance = mean(fold_deviance, na.rm = TRUE),
    fold_deviance   = fold_deviance,
    fold_results    = fold_results,
    k               = k
  )
}


# check_inla() and %||% are defined in zzz.R
