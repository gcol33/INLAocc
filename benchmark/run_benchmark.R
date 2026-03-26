# =============================================================================
# INLAocc vs spOccupancy Benchmark
#
# Compares: speed, accuracy (psi recovery), coefficient recovery
# Across: varying N (sites), J (visits), model complexity
# =============================================================================

library(spOccupancy)
devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(2026)
results_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/results.rds"

cat("=== INLAocc vs spOccupancy Benchmark ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# -----------------------------------------------------------------------------
# Scenario definitions
# -----------------------------------------------------------------------------
scenarios <- list(
  # Small: quick validation
  list(name = "small_basic",     N = 100,  J = 4, n_occ = 1, n_det = 1, spatial = FALSE),
  list(name = "small_2cov",      N = 100,  J = 4, n_occ = 2, n_det = 2, spatial = FALSE),

  # Medium: typical field study
  list(name = "medium_basic",    N = 500,  J = 4, n_occ = 1, n_det = 1, spatial = FALSE),
  list(name = "medium_2cov",     N = 500,  J = 4, n_occ = 2, n_det = 2, spatial = FALSE),
  list(name = "medium_6visit",   N = 500,  J = 6, n_occ = 2, n_det = 1, spatial = FALSE),

 # Large: regional study
  list(name = "large_basic",     N = 2000, J = 4, n_occ = 1, n_det = 1, spatial = FALSE),
  list(name = "large_2cov",      N = 2000, J = 4, n_occ = 2, n_det = 2, spatial = FALSE),

  # Very large: pushing limits
  list(name = "xlarge_basic",    N = 5000, J = 4, n_occ = 1, n_det = 1, spatial = FALSE),

  # Spatial models — where INLA should shine
  list(name = "spatial_small",   N = 200,  J = 4, n_occ = 1, n_det = 1, spatial = TRUE),
  list(name = "spatial_medium",  N = 500,  J = 4, n_occ = 1, n_det = 1, spatial = TRUE),
  list(name = "spatial_large",   N = 1000, J = 4, n_occ = 1, n_det = 1, spatial = TRUE),
  list(name = "spatial_xlarge",  N = 2000, J = 4, n_occ = 1, n_det = 1, spatial = TRUE)
)


# -----------------------------------------------------------------------------
# Simulate data for a scenario
# -----------------------------------------------------------------------------
simulate_scenario <- function(sc) {
  N <- sc$N; J <- sc$J
  n_occ <- sc$n_occ; n_det <- sc$n_det

  # True parameters
  beta_occ <- c(0.5, seq(-0.8, 0.8, length.out = n_occ))
  beta_det <- c(0.0, seq(-0.5, 0.5, length.out = n_det))

  # Covariates
  X_occ <- matrix(rnorm(N * n_occ), N, n_occ)
  colnames(X_occ) <- paste0("occ_cov", seq_len(n_occ))
  occ_covs <- as.data.frame(X_occ)

  det_covs <- list()
  for (k in seq_len(n_det)) {
    det_covs[[paste0("det_cov", k)]] <- matrix(rnorm(N * J), N, J)
  }

  coords <- cbind(x = runif(N), y = runif(N))

  # True occupancy
  eta_occ <- beta_occ[1] + X_occ %*% beta_occ[-1]

  # Add spatial effect if needed
  w_true <- rep(0, N)
  if (sc$spatial) {
    D <- as.matrix(dist(coords))
    Sigma <- 1.0 * exp(-D / 0.2)
    w_true <- MASS::mvrnorm(1, rep(0, N), Sigma)
    eta_occ <- eta_occ + w_true
  }

  psi <- plogis(as.vector(eta_occ))
  z <- rbinom(N, 1, psi)

  # Detection
  p_mat <- matrix(NA, N, J)
  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      eta_det <- beta_det[1]
      for (k in seq_len(n_det)) {
        eta_det <- eta_det + beta_det[k + 1] * det_covs[[k]][i, j]
      }
      p_mat[i, j] <- plogis(eta_det)
    }
  }

  y <- matrix(NA, N, J)
  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      y[i, j] <- rbinom(1, 1, z[i] * p_mat[i, j])
    }
  }

  list(
    y = y, occ.covs = occ_covs, det.covs = det_covs, coords = coords,
    truth = list(z = z, psi = psi, beta_occ = beta_occ, beta_det = beta_det,
                 w = w_true)
  )
}


# -----------------------------------------------------------------------------
# Fit with spOccupancy
# -----------------------------------------------------------------------------
fit_spOcc <- function(data, sc) {
  occ_formula <- as.formula(paste("~", paste(names(data$occ.covs), collapse = " + ")))
  det_formula <- as.formula(paste("~", paste(names(data$det.covs), collapse = " + ")))

  # Short MCMC for benchmarking (enough for point estimates)
  n.samples <- 5000
  n.burn <- 2000
  n.thin <- 3

  t0 <- proc.time()

  if (sc$spatial) {
    fit <- tryCatch({
      spPGOcc(
        occ.formula = occ_formula,
        det.formula = det_formula,
        data = list(y = data$y, occ.covs = data$occ.covs,
                    det.covs = data$det.covs, coords = data$coords),
        n.batch = 100, batch.length = 50,
        n.burn = 2000, n.thin = 3,
        cov.model = "exponential",
        NNGP = TRUE, n.neighbors = 8,
        n.omp.threads = 1,
        verbose = FALSE
      )
    }, error = function(e) { message("spOcc spatial failed: ", e$message); NULL })
  } else {
    fit <- tryCatch({
      PGOcc(
        occ.formula = occ_formula,
        det.formula = det_formula,
        data = list(y = data$y, occ.covs = data$occ.covs,
                    det.covs = data$det.covs),
        n.samples = n.samples, n.burn = n.burn, n.thin = n.thin,
        n.omp.threads = 1,
        verbose = FALSE
      )
    }, error = function(e) { message("spOcc failed: ", e$message); NULL })
  }

  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(fit)) return(list(time = elapsed, psi = NULL, beta_occ = NULL, beta_det = NULL))

  # Extract estimates
  psi_est <- apply(fit$psi.samples, 2, mean)
  beta_occ_est <- apply(fit$beta.samples, 2, mean)
  beta_det_est <- apply(fit$alpha.samples, 2, mean)

  list(time = elapsed, psi = psi_est, beta_occ = beta_occ_est, beta_det = beta_det_est)
}


# -----------------------------------------------------------------------------
# Fit with INLAocc
# -----------------------------------------------------------------------------
fit_INLAocc <- function(data, sc) {
  occ_formula <- as.formula(paste("~", paste(names(data$occ.covs), collapse = " + ")))
  det_formula <- as.formula(paste("~", paste(names(data$det.covs), collapse = " + ")))

  t0 <- proc.time()

  fit <- tryCatch({
    if (sc$spatial) {
      occu(occ_formula, det_formula,
           data = list(y = data$y, occ.covs = data$occ.covs,
                       det.covs = data$det.covs, coords = data$coords),
           spatial = data$coords, verbose = 0)
    } else {
      occu(occ_formula, det_formula,
           data = list(y = data$y, occ.covs = data$occ.covs,
                       det.covs = data$det.covs),
           verbose = 0)
    }
  }, error = function(e) { message("INLAocc failed: ", e$message); NULL })

  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(fit)) return(list(time = elapsed, psi = NULL, beta_occ = NULL, beta_det = NULL))

  psi_est <- fit$psi_hat
  beta_occ_est <- fit$occ_fit$summary.fixed$mean
  beta_det_est <- fit$det_fit$summary.fixed$mean

  list(time = elapsed, psi = psi_est, beta_occ = beta_occ_est, beta_det = beta_det_est)
}


# -----------------------------------------------------------------------------
# Run all scenarios
# -----------------------------------------------------------------------------
all_results <- list()

for (i in seq_along(scenarios)) {
  sc <- scenarios[[i]]
  cat(sprintf("[%d/%d] %s (N=%d, J=%d, spatial=%s) ...\n",
              i, length(scenarios), sc$name, sc$N, sc$J, sc$spatial))

  # Simulate
  data <- simulate_scenario(sc)

  # Fit spOccupancy
  cat("  spOccupancy ... ")
  res_spOcc <- fit_spOcc(data, sc)
  cat(sprintf("%.1fs\n", res_spOcc$time))

  # Fit INLAocc
  cat("  INLAocc     ... ")
  res_inla <- fit_INLAocc(data, sc)
  cat(sprintf("%.1fs\n", res_inla$time))

  # Compute accuracy metrics
  metrics <- list(
    scenario   = sc$name,
    N          = sc$N,
    J          = sc$J,
    spatial    = sc$spatial,
    n_occ_covs = sc$n_occ,
    n_det_covs = sc$n_det,

    # Timing
    time_spOcc   = res_spOcc$time,
    time_INLAocc = res_inla$time,
    speedup      = res_spOcc$time / max(res_inla$time, 0.01),

    # Psi recovery (correlation with truth)
    cor_psi_spOcc   = if (!is.null(res_spOcc$psi)) cor(data$truth$psi, res_spOcc$psi) else NA,
    cor_psi_INLAocc = if (!is.null(res_inla$psi))  cor(data$truth$psi, res_inla$psi)  else NA,

    # Psi RMSE
    rmse_psi_spOcc   = if (!is.null(res_spOcc$psi)) sqrt(mean((data$truth$psi - res_spOcc$psi)^2)) else NA,
    rmse_psi_INLAocc = if (!is.null(res_inla$psi))  sqrt(mean((data$truth$psi - res_inla$psi)^2))  else NA,

    # Beta recovery (max absolute error)
    beta_mae_spOcc   = if (!is.null(res_spOcc$beta_occ))
      max(abs(data$truth$beta_occ - res_spOcc$beta_occ)) else NA,
    beta_mae_INLAocc = if (!is.null(res_inla$beta_occ))
      max(abs(data$truth$beta_occ - res_inla$beta_occ)) else NA
  )

  all_results[[i]] <- metrics

  cat(sprintf("  Speedup: %.1fx | Cor(psi): spOcc=%.3f INLAocc=%.3f\n",
              metrics$speedup, metrics$cor_psi_spOcc, metrics$cor_psi_INLAocc))
  cat("\n")

  # Save incrementally
  saveRDS(all_results, results_file)
}


# -----------------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------------
cat("\n=== SUMMARY ===\n\n")

summary_df <- do.call(rbind, lapply(all_results, function(r) {
  data.frame(
    scenario    = r$scenario,
    N           = r$N,
    spatial     = r$spatial,
    t_spOcc     = round(r$time_spOcc, 1),
    t_INLAocc   = round(r$time_INLAocc, 1),
    speedup     = round(r$speedup, 1),
    cor_spOcc   = round(r$cor_psi_spOcc, 3),
    cor_INLAocc = round(r$cor_psi_INLAocc, 3),
    rmse_spOcc  = round(r$rmse_psi_spOcc, 4),
    rmse_INLAocc = round(r$rmse_psi_INLAocc, 4),
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)

# Save final results
saveRDS(all_results, results_file)
write.csv(summary_df, sub("\\.rds$", ".csv", results_file), row.names = FALSE)

cat(sprintf("\nFinished: %s\n", Sys.time()))
cat(sprintf("Results saved to: %s\n", results_file))
