# =============================================================================
# Scaling benchmark: INLAocc vs spOccupancy vs Stan
#
# Measures computation time across increasing dataset sizes.
# Non-spatial: all three methods.  Spatial: INLAocc + spOccupancy.
# =============================================================================

library(cmdstanr)
library(spOccupancy)
devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(42)
out_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results.csv"

cat("=== Scaling Benchmark ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# -----------------------------------------------------------------------------
# 1. Compile Stan model (once)
# -----------------------------------------------------------------------------
cat("Compiling Stan model ... ")
stan_mod <- cmdstan_model(
  "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/occu_model.stan",
  quiet = TRUE
)
cat("done.\n\n")

# -----------------------------------------------------------------------------
# 2. Scenario grid
# -----------------------------------------------------------------------------
# Evenly spaced on log10 scale: 2.0, 2.5, 3.0, 3.5, 4.0
N_grid <- c(100, 300, 1000, 3000, 10000)
J <- 4L
beta_occ <- c(0.5, -0.8)   # intercept + 1 covariate
beta_det <- c(0.0,  0.5)

# Stan is too slow beyond N = 3000 for a single-run benchmark
stan_max_N <- 3000

# Spatial only up to N = 3000 (covariance matrix becomes expensive to simulate)
spatial_max_N <- 3000

# -----------------------------------------------------------------------------
# 3. Simulate data
# -----------------------------------------------------------------------------
simulate_data <- function(N, J, spatial = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Occupancy covariates (with intercept)
  x1 <- rnorm(N)
  X_occ <- cbind(1, x1)
  occ_covs <- data.frame(occ_x1 = x1)

  # Detection covariates (site x visit)
  det_vals <- rnorm(N * J)
  det_mat  <- matrix(det_vals, N, J)
  det_covs <- list(det_x1 = det_mat)
  X_det    <- cbind(1, det_vals)             # stacked: N*J rows

  coords <- cbind(x = runif(N), y = runif(N))

  # True occupancy linear predictor
  eta_occ <- X_occ %*% beta_occ

  # Optional spatial random effect
  if (spatial) {
    D     <- as.matrix(dist(coords))
    Sigma <- 1.0 * exp(-D / 0.2)
    w     <- MASS::mvrnorm(1, rep(0, N), Sigma)
    eta_occ <- eta_occ + w
  }

  psi <- plogis(as.vector(eta_occ))
  z   <- rbinom(N, 1, psi)

  # Detection
  eta_det <- X_det %*% beta_det              # length N*J
  p_vec   <- plogis(as.vector(eta_det))
  z_rep   <- rep(z, each = J)
  y_vec   <- rbinom(N * J, 1, z_rep * p_vec)
  y       <- matrix(y_vec, N, J, byrow = TRUE)

  list(
    y        = y,
    occ.covs = occ_covs,
    det.covs = det_covs,
    coords   = coords,
    psi_true = psi,
    stan_data = list(
      N = N, J = J, y = y,
      K_occ = 2L, K_det = 2L,
      X_occ = X_occ, X_det = X_det
    )
  )
}

# -----------------------------------------------------------------------------
# 4. Fitting wrappers
# -----------------------------------------------------------------------------
time_inlaocc <- function(data, spatial = FALSE) {
  t0 <- proc.time()
  fit <- tryCatch({
    if (spatial) {
      d <- list(y = data$y, occ.covs = data$occ.covs,
                det.covs = data$det.covs, coords = data$coords)
      occu(~ occ_x1, ~ det_x1, data = d, spatial = data$coords, verbose = 0)
    } else {
      d <- list(y = data$y, occ.covs = data$occ.covs, det.covs = data$det.covs)
      occu(~ occ_x1, ~ det_x1, data = d, verbose = 0)
    }
  }, error = function(e) { message("  INLAocc error: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  psi_hat <- if (!is.null(fit)) fit$psi_hat else NULL
  list(time = elapsed, psi = psi_hat)
}

time_spocc <- function(data, spatial = FALSE) {
  d <- if (spatial) {
    list(y = data$y, occ.covs = data$occ.covs,
         det.covs = data$det.covs, coords = data$coords)
  } else {
    list(y = data$y, occ.covs = data$occ.covs, det.covs = data$det.covs)
  }

  t0 <- proc.time()
  fit <- tryCatch({
    if (spatial) {
      spPGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1, data = d,
              n.batch = 100, batch.length = 50,
              n.burn = 2000, n.thin = 3,
              cov.model = "exponential", NNGP = TRUE, n.neighbors = 8,
              n.omp.threads = 1, verbose = FALSE)
    } else {
      PGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1, data = d,
            n.samples = 5000, n.burn = 2000, n.thin = 3,
            n.omp.threads = 1, verbose = FALSE)
    }
  }, error = function(e) { message("  spOcc error: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  psi_hat <- if (!is.null(fit)) apply(fit$psi.samples, 2, mean) else NULL
  list(time = elapsed, psi = psi_hat)
}

time_stan <- function(data, stan_mod) {
  t0 <- proc.time()
  fit <- tryCatch({
    stan_mod$sample(
      data = data$stan_data,
      chains = 1, iter_warmup = 1000, iter_sampling = 1000,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE
    )
  }, error = function(e) { message("  Stan error: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  psi_hat <- if (!is.null(fit)) {
    colMeans(fit$draws("psi", format = "matrix"))
  } else NULL
  list(time = elapsed, psi = psi_hat)
}

# -----------------------------------------------------------------------------
# 5. Run
# -----------------------------------------------------------------------------
results <- data.frame(
  N = integer(), type = character(), method = character(),
  time = numeric(), cor_psi = numeric(),
  stringsAsFactors = FALSE
)

add_row <- function(df, N, type, method, time, psi_hat, psi_true) {
  r <- if (!is.null(psi_hat)) cor(psi_true, psi_hat) else NA_real_
  rbind(df, data.frame(N = N, type = type, method = method,
                       time = time, cor_psi = round(r, 4),
                       stringsAsFactors = FALSE))
}

## Non-spatial
cat("--- Non-spatial ---\n")
for (N in N_grid) {
  cat(sprintf("\nN = %d\n", N))
  data <- simulate_data(N, J, spatial = FALSE, seed = N)

  cat("  INLAocc       ... "); flush.console()
  r <- time_inlaocc(data, FALSE)
  cat(sprintf("%6.1fs\n", r$time))
  results <- add_row(results, N, "Non-spatial", "INLAocc", r$time, r$psi, data$psi_true)

  cat("  spOccupancy   ... "); flush.console()
  r <- time_spocc(data, FALSE)
  cat(sprintf("%6.1fs\n", r$time))
  results <- add_row(results, N, "Non-spatial", "spOccupancy", r$time, r$psi, data$psi_true)

  if (N <= stan_max_N) {
    cat("  Stan          ... "); flush.console()
    r <- time_stan(data, stan_mod)
    cat(sprintf("%6.1fs\n", r$time))
    results <- add_row(results, N, "Non-spatial", "Stan", r$time, r$psi, data$psi_true)
  }

  write.csv(results, out_file, row.names = FALSE)
}

## Spatial
cat("\n--- Spatial ---\n")
for (N in N_grid[N_grid <= spatial_max_N]) {
  cat(sprintf("\nN = %d\n", N))
  data <- simulate_data(N, J, spatial = TRUE, seed = N + 1000)

  cat("  INLAocc       ... "); flush.console()
  r <- time_inlaocc(data, TRUE)
  cat(sprintf("%6.1fs\n", r$time))
  results <- add_row(results, N, "Spatial", "INLAocc", r$time, r$psi, data$psi_true)

  cat("  spOccupancy   ... "); flush.console()
  r <- time_spocc(data, TRUE)
  cat(sprintf("%6.1fs\n", r$time))
  results <- add_row(results, N, "Spatial", "spOccupancy", r$time, r$psi, data$psi_true)

  write.csv(results, out_file, row.names = FALSE)
}

write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\n\nFinished: %s\n", Sys.time()))
cat(sprintf("Results: %s\n", out_file))
print(results)
