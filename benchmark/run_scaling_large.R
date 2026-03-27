# =============================================================================
# Large-N extension: push scaling benchmarks to 30k–100k sites
# Appends results to scaling_results.csv and scaling_results_ms.csv
# =============================================================================

library(spOccupancy)
devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(2026)
root <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark"

# --- Helpers (same as run_scaling.R) ---
J         <- 4L
beta_occ  <- c(0.5, -0.8)
beta_det  <- c(0.0,  0.5)
n_species <- 10L

simulate_data <- function(N, J, spatial = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(N)
  X_occ <- cbind(1, x1)
  occ_covs <- data.frame(occ_x1 = x1)
  det_vals <- rnorm(N * J)
  det_mat  <- matrix(det_vals, N, J)
  det_covs <- list(det_x1 = det_mat)
  coords   <- cbind(x = runif(N), y = runif(N))
  eta_occ  <- X_occ %*% beta_occ

  if (spatial) {
    # Spectral approximation for large N (avoids O(N^2) covariance matrix)
    n_basis <- 100
    w <- rep(0, N)
    for (k in seq_len(n_basis)) {
      freq  <- runif(2, -5, 5) / 0.2
      phase <- runif(1, 0, 2 * pi)
      amp   <- rnorm(1, 0, 1 / sqrt(n_basis))
      w <- w + amp * cos(coords[, 1] * freq[1] + coords[, 2] * freq[2] + phase)
    }
    eta_occ <- eta_occ + w
  }

  psi   <- plogis(as.vector(eta_occ))
  z     <- rbinom(N, 1, psi)
  X_det <- cbind(1, det_vals)
  p_vec <- plogis(as.vector(X_det %*% beta_det))
  y_vec <- rbinom(N * J, 1, rep(z, each = J) * p_vec)
  y     <- matrix(y_vec, N, J, byrow = TRUE)

  list(
    y = y, occ.covs = occ_covs, det.covs = det_covs, coords = coords,
    psi_true = psi
  )
}

time_inlaocc <- function(data, spatial = FALSE) {
  t0 <- proc.time()
  fit <- tryCatch({
    d <- list(y = data$y, occ.covs = data$occ.covs,
              det.covs = data$det.covs, coords = data$coords)
    if (spatial) occu(~ occ_x1, ~ det_x1, data = d, spatial = data$coords, verbose = 0)
    else         occu(~ occ_x1, ~ det_x1, data = d, verbose = 0)
  }, error = function(e) { message("  INLAocc: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  psi <- if (!is.null(fit)) fit$psi_hat else NULL
  list(time = elapsed, psi = psi)
}

time_spocc <- function(data, spatial = FALSE) {
  d <- list(y = data$y, occ.covs = data$occ.covs,
            det.covs = data$det.covs, coords = data$coords)
  t0 <- proc.time()
  fit <- tryCatch({
    if (spatial) {
      spPGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1, data = d,
              n.batch = 100, batch.length = 50, n.burn = 2000, n.thin = 3,
              cov.model = "exponential", NNGP = TRUE, n.neighbors = 8,
              n.omp.threads = 1, verbose = FALSE)
    } else {
      PGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1, data = d,
            n.samples = 5000, n.burn = 2000, n.thin = 3,
            n.omp.threads = 1, verbose = FALSE)
    }
  }, error = function(e) { message("  spOcc: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  psi <- if (!is.null(fit)) apply(fit$psi.samples, 2, mean) else NULL
  list(time = elapsed, psi = psi)
}

add_row <- function(df, N, type, method, time, psi_hat, psi_true) {
  r <- if (!is.null(psi_hat)) cor(psi_true, psi_hat) else NA_real_
  rbind(df, data.frame(N = N, type = type, method = method,
                       time = time, cor_psi = round(r, 4),
                       stringsAsFactors = FALSE))
}

# --- Load existing results ---
ss_file <- file.path(root, "scaling_results.csv")
ms_file <- file.path(root, "scaling_results_ms.csv")
ss <- read.csv(ss_file, stringsAsFactors = FALSE)
ms <- read.csv(ms_file, stringsAsFactors = FALSE)

cat("=== Large-N Extension ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# =========================================================================
# 1. Non-spatial: INLAocc at 30k, 100k; spOcc at 30k
# =========================================================================
cat("--- Non-spatial large N ---\n")

for (N in c(30000, 100000)) {
  cat(sprintf("\nN = %d\n", N))
  data <- simulate_data(N, J, spatial = FALSE, seed = N)

  cat("  INLAocc       ... "); flush.console()
  r <- time_inlaocc(data)
  cat(sprintf("%6.1fs\n", r$time))
  ss <- add_row(ss, N, "Non-spatial", "INLAocc", r$time, r$psi, data$psi_true)

  if (N <= 30000) {
    cat("  spOccupancy   ... "); flush.console()
    r <- time_spocc(data)
    cat(sprintf("%6.1fs\n", r$time))
    ss <- add_row(ss, N, "Non-spatial", "spOccupancy", r$time, r$psi, data$psi_true)
  }

  write.csv(ss, ss_file, row.names = FALSE)
}

# =========================================================================
# 2. Spatial: both at 5k, 10k
# =========================================================================
cat("\n--- Spatial large N ---\n")

for (N in c(5000, 10000)) {
  cat(sprintf("\nN = %d\n", N))
  data <- simulate_data(N, J, spatial = TRUE, seed = N + 1000)

  cat("  INLAocc       ... "); flush.console()
  r <- time_inlaocc(data, spatial = TRUE)
  cat(sprintf("%6.1fs\n", r$time))
  ss <- add_row(ss, N, "Spatial", "INLAocc", r$time, r$psi, data$psi_true)

  cat("  spOccupancy   ... "); flush.console()
  r <- time_spocc(data, spatial = TRUE)
  cat(sprintf("%6.1fs\n", r$time))
  ss <- add_row(ss, N, "Spatial", "spOccupancy", r$time, r$psi, data$psi_true)

  write.csv(ss, ss_file, row.names = FALSE)
}

# =========================================================================
# 3. Multi-species: both at 10k
# =========================================================================
cat("\n--- Multi-species large N ---\n")

N <- 10000
cat(sprintf("\nN = %d, S = %d\n", N, n_species))
sim <- simMsOcc(N = N, J = J, n_species = n_species,
                n_occ_covs = 1, n_det_covs = 1, seed = N + 2000)

cat("  INLAocc       ... "); flush.console()
t0 <- proc.time()
fit <- tryCatch(
  occu(~ occ_x1, ~ det_x1, data = sim$data, multispecies = TRUE, verbose = 0),
  error = function(e) { message("  INLAocc: ", e$message); NULL }
)
elapsed <- (proc.time() - t0)[["elapsed"]]
cor_avg <- NA_real_
if (!is.null(fit) && !is.null(fit$species_fits)) {
  cors <- vapply(sim$data$species_names, function(sp) {
    sp_fit <- fit$species_fits[[sp]]
    if (!is.null(sp_fit)) cor(sim$truth$psi[match(sp, sim$data$species_names), ],
                              sp_fit$psi_hat) else NA_real_
  }, numeric(1))
  cor_avg <- mean(cors, na.rm = TRUE)
}
cat(sprintf("%6.1fs  (cor = %.3f)\n", elapsed, cor_avg))
ms <- rbind(ms, data.frame(N = N, type = "Multi-species", method = "INLAocc",
                           time = elapsed, cor_psi = round(cor_avg, 4),
                           stringsAsFactors = FALSE))
write.csv(ms, ms_file, row.names = FALSE)

cat("  spOccupancy   ... "); flush.console()
y_array <- array(NA, dim = c(n_species, N, J))
for (s in seq_len(n_species)) {
  sp_name <- sim$data$species_names[s]
  y_array[s, , ] <- sim$data$species_data[[sp_name]]$y
}
spocc_data <- list(y = y_array, occ.covs = sim$data$occ.covs, det.covs = sim$data$det.covs)
t0 <- proc.time()
fit <- tryCatch(
  msPGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1, data = spocc_data,
          n.samples = 5000, n.burn = 2000, n.thin = 3,
          n.omp.threads = 1, verbose = FALSE),
  error = function(e) { message("  spOcc: ", e$message); NULL }
)
elapsed <- (proc.time() - t0)[["elapsed"]]
cor_avg <- NA_real_
if (!is.null(fit)) {
  psi_hat <- apply(fit$psi.samples, c(2, 3), mean)
  cors <- vapply(seq_len(n_species), function(s) cor(sim$truth$psi[s, ], psi_hat[s, ]),
                 numeric(1))
  cor_avg <- mean(cors, na.rm = TRUE)
}
cat(sprintf("%6.1fs  (cor = %.3f)\n", elapsed, cor_avg))
ms <- rbind(ms, data.frame(N = N, type = "Multi-species", method = "spOccupancy",
                           time = elapsed, cor_psi = round(cor_avg, 4),
                           stringsAsFactors = FALSE))
write.csv(ms, ms_file, row.names = FALSE)

cat(sprintf("\n\nFinished: %s\n", Sys.time()))
