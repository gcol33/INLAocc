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
  p_c <- clamp(p)
  psi_c <- clamp(psi)

  detected <- rowSums(y == 1, na.rm = TRUE) > 0
  valid <- !is.na(y)

  # Detection log-likelihood per cell (0 where NA)
  ll_cell <- ifelse(valid,
                    y * log(p_c) + (1 - y) * log(1 - p_c),
                    0)

  # Detected sites: log(psi) + sum of detection log-lik
  ll_det <- log(psi_c[detected]) + rowSums(ll_cell[detected, , drop = FALSE])

  # Undetected sites: log((1 - psi) + psi * prod(1 - p))
  # prod(1 - p) across visits = exp(sum(log(1-p))) for non-NA cells
  log_q <- rowSums(ifelse(valid, log(1 - p_c), 0))
  q_vec <- exp(log_q)
  ll_undet <- log(clamp((1 - psi_c[!detected]) + psi_c[!detected] * q_vec[!detected]))

  sum(ll_det) + sum(ll_undet)
}

#' @noRd
compute_weights <- function(y, psi, p) {
  N <- nrow(y)
  detected <- rowSums(y == 1, na.rm = TRUE) > 0
  has_visits <- rowSums(!is.na(y)) > 0

  # Default: detected = 1, no visits = 0.5
  w <- rep(0.5, N)
  w[detected] <- 1

  # Undetected with visits: Bayes update
  undet <- !detected & has_visits
  if (any(undet)) {
    p_c <- clamp(p)
    # prod(1 - p) per site via log-sum, only over non-NA cells
    log_q <- rowSums(ifelse(!is.na(y), log(1 - p_c), 0))
    q_vec <- exp(log_q[undet])
    psi_c <- clamp(psi[undet])
    w[undet] <- (psi_c * q_vec) / ((1 - psi_c) + psi_c * q_vec)
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

  # Vectorized: find all non-NA (site, visit) pairs
  valid <- which(!is.na(y), arr.ind = TRUE)
  ri <- valid[, 1]
  ci <- valid[, 2]

  df <- data.frame(
    y_det = y[valid],
    site  = site_idx[ri],
    visit = ci,
    w     = weights[ri]
  )

  if (!is.null(det_covs)) {
    for (nm in names(det_covs)) {
      df[[nm]] <- det_covs[[nm]][valid]
    }
  }

  if (!is.null(visit_idx)) {
    df$visit_group <- visit_idx[valid]
  }

  df
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

#' Build weighted Bernoulli occupancy data for spatial models
#'
#' Uses two rows per undetected site (z=1 with weight w, z=0 with weight 1-w)
#' instead of scaled binomial encoding. This keeps the SPDE prior properly
#' influential relative to the data.
#' @noRd
build_occ_df_weighted <- function(occ_covs, weights, site_idx = NULL,
                                   detected = NULL) {
  N <- nrow(occ_covs)
  if (is.null(site_idx)) site_idx <- seq_len(N)
  if (is.null(detected)) detected <- weights == 1

  n_undet <- sum(!detected)
  n_rows <- N + n_undet

  # Row indices: each site once, plus undetected sites get a second row
  src_rows <- c(seq_len(N), which(!detected))

  df <- occ_covs[src_rows, , drop = FALSE]
  rownames(df) <- NULL
  df$site    <- site_idx[src_rows]
  df$Ntrials <- 1L

  # z and w vectors
  w_clamped <- clamp(weights, 1e-6, 1 - 1e-6)

  z_vec <- integer(n_rows)
  w_vec <- numeric(n_rows)

  # First N rows: one per site
  z_vec[seq_len(N)]  <- 1L           # all z=1 (detected stays 1, undetected gets w)
  w_vec[seq_len(N)]  <- ifelse(detected, 1.0, w_clamped)

  # Extra rows for undetected: z=0 with weight 1-w
  if (n_undet > 0) {
    extra_idx <- (N + 1L):n_rows
    z_vec[extra_idx] <- 0L
    w_vec[extra_idx] <- 1 - w_clamped[!detected]
  }

  df$z <- z_vec
  df$w <- w_vec
  df
}

# =============================================================================
# Formula parser — native AST walker (ported from tulpa)
#
# R formulas are already parse trees. We do structural recursion to find
# random effect terms (| nodes) and separate them from fixed effects.
# No external dependencies (lme4), no regex on deparsed strings.
# =============================================================================

#' Find all bar terms in a formula's parse tree
#'
#' Recursively walks the formula AST and collects all `|` and `||` nodes
#' found inside parentheses.
#' @noRd
findbars_ast <- function(term) {
  if (is.name(term) || !is.language(term)) return(NULL)

  # Parenthesized expression: check if it wraps a bar term
  if (term[[1]] == as.name("(")) {
    inner <- term[[2]]
    if (is.call(inner) &&
        (inner[[1]] == as.name("|") || inner[[1]] == as.name("||"))) {
      return(list(inner))
    }
    return(findbars_ast(inner))
  }

  # Bar term found directly
  if (term[[1]] == as.name("|") || term[[1]] == as.name("||")) {
    return(list(term))
  }

  # Unary operator: recurse into operand
  if (length(term) == 2) return(findbars_ast(term[[2]]))

  # Binary operator: recurse both sides
  c(findbars_ast(term[[2]]), findbars_ast(term[[3]]))
}

#' Remove all bar terms from a formula's parse tree
#'
#' Rewrites the formula AST with any `|`/`||` nodes inside parentheses removed.
#' @noRd
nobars_ast <- function(term) {
  if (is.name(term) || !is.language(term)) return(term)

  # Parenthesized bar term: drop entirely
  if (term[[1]] == as.name("(")) {
    inner <- term[[2]]
    if (is.call(inner) &&
        (inner[[1]] == as.name("|") || inner[[1]] == as.name("||"))) {
      return(NULL)
    }
    nb <- nobars_ast(inner)
    if (is.null(nb)) return(NULL)
    return(call("(", nb))
  }

  # Unary operator
  if (length(term) == 2) {
    nb <- nobars_ast(term[[2]])
    if (is.null(nb)) return(NULL)
    return(call(deparse(term[[1]]), nb))
  }

  # Binary operator: recurse both sides, handle NULL removal
  nb_left  <- nobars_ast(term[[2]])
  nb_right <- nobars_ast(term[[3]])

  if (is.null(nb_left) && is.null(nb_right)) return(NULL)
  if (is.null(nb_left))  return(nb_right)
  if (is.null(nb_right)) return(nb_left)

  call(deparse(term[[1]]), nb_left, nb_right)
}

#' Flatten a + chain into individual terms
#' @noRd
collect_additive_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    c(collect_additive_terms(expr[[2]]), collect_additive_terms(expr[[3]]))
  } else {
    list(expr)
  }
}

#' Parse a single bar term into a list of occu_re objects
#'
#' Takes a `|` or `||` language object and emits one occu_re per effect term.
#' Handles intercept indicators (0, 1, -1), nested grouping (a/b), and
#' uncorrelated syntax (||).
#' @noRd
bar_to_occu_re <- function(bar_term) {
  lhs <- bar_term[[2]]
  rhs <- bar_term[[3]]

  # Resolve grouping variable
  if (is.call(rhs) && identical(rhs[[1]], as.name("/"))) {
    # Nested grouping: (1 | a/b) → group = "a:b"
    group_var <- paste0(as.character(rhs[[2]]), ":", as.character(rhs[[3]]))
  } else {
    group_var <- if (is.name(rhs)) as.character(rhs) else deparse(rhs)
  }

  # Decompose LHS into intercept flag + slope terms
  terms <- collect_additive_terms(lhs)
  has_intercept <- TRUE
  slopes <- list()

  for (term in terms) {
    if (is.numeric(term)) {
      if (term == 0) has_intercept <- FALSE
      # 1 = intercept marker, not a slope
    } else if (is.call(term) && identical(term[[1]], as.name("-")) &&
               length(term) == 2 && is.numeric(term[[2]]) && term[[2]] == 1) {
      # -1 suppresses intercept
      has_intercept <- FALSE
    } else {
      slopes <- c(slopes, list(term))
    }
  }

  # Emit occu_re objects
  re_list <- list()
  if (has_intercept) {
    re_list[[length(re_list) + 1L]] <- occu_re("intercept", group = group_var)
  }
  for (s in slopes) {
    re_list[[length(re_list) + 1L]] <- occu_re("slope", group = group_var,
                                                 covariate = deparse(s))
  }
  re_list
}

#' Parse a mixed-model formula into fixed effects and random effects
#'
#' Walks the formula's abstract syntax tree to find bar terms (`|`, `||`),
#' converts each to occu_re objects, and returns the fixed-effects-only formula.
#' No external dependencies — pure structural recursion on R's own parse tree.
#' @noRd
parse_re_formula <- function(formula) {
  # Find all bar terms via AST recursion
  bars <- findbars_ast(formula)

  # Parse each bar term into occu_re objects
  re_list <- list()
  if (length(bars) > 0) {
    for (bar in bars) {
      re_list <- c(re_list, bar_to_occu_re(bar))
    }
  }

  # Remove bar terms to get fixed-effects formula
  fixed_call <- nobars_ast(formula)
  if (is.null(fixed_call) ||
      (length(fixed_call) == 2 && is.null(fixed_call[[2]]))) {
    fixed <- ~ 1
  } else {
    fixed <- as.formula(fixed_call, env = environment(formula))
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
