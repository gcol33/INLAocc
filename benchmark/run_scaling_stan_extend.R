# =============================================================================
# Extend Stan data points to match other methods where feasible
# Non-spatial: add N=10000
# Multi-species: add N=3000
# =============================================================================

library(INLAocc)
library(cmdstanr)

set.seed(42)
beta_occ <- c(0.5, -0.8)
beta_det <- c(0.0,  0.5)
J <- 4L

# --- Compile models ---
cat("Compiling Stan models ... ")
stan_ss <- cmdstan_model(
  "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/occu_model.stan", quiet = TRUE)
stan_ms <- cmdstan_model(
  "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/occu_ms_model.stan", quiet = TRUE)
cat("done.\n\n")

# --- Load existing CSVs ---
ss_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results.csv"
ms_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results_ms.csv"
ss <- read.csv(ss_file, stringsAsFactors = FALSE)
ms <- read.csv(ms_file, stringsAsFactors = FALSE)

# =========================================================================
# 1. Non-spatial Stan at N=10000
# =========================================================================
N <- 10000L
cat(sprintf("Non-spatial Stan, N = %d\n", N))
set.seed(N)
x1 <- rnorm(N)
X_occ <- cbind(1, x1)
det_vals <- rnorm(N * J)
X_det <- cbind(1, det_vals)
eta_occ <- X_occ %*% beta_occ
psi <- plogis(as.vector(eta_occ))
z <- rbinom(N, 1, psi)
p_vec <- plogis(as.vector(X_det %*% beta_det))
y_vec <- rbinom(N * J, 1, rep(z, each = J) * p_vec)
y <- matrix(y_vec, N, J, byrow = TRUE)

stan_data <- list(N = N, J = J, y = y, K_occ = 2L, K_det = 2L,
                  X_occ = X_occ, X_det = X_det)

cat("  Stan          ... "); flush.console()
t0 <- proc.time()
fit <- stan_ss$sample(data = stan_data, chains = 1,
                      iter_warmup = 1000, iter_sampling = 1000,
                      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
elapsed <- (proc.time() - t0)[["elapsed"]]
psi_hat <- colMeans(fit$draws("psi", format = "matrix"))
r <- cor(psi, psi_hat)
cat(sprintf("%6.1fs  (cor = %.3f)\n", elapsed, r))

ss <- rbind(ss, data.frame(N = N, type = "Non-spatial", method = "Stan",
                           time = elapsed, cor_psi = round(r, 4),
                           stringsAsFactors = FALSE))
write.csv(ss, ss_file, row.names = FALSE)

# =========================================================================
# 2. Multi-species Stan at N=3000
# =========================================================================
N <- 3000L; n_species <- 10L
cat(sprintf("\nMulti-species Stan, N = %d, S = %d\n", N, n_species))

sim <- INLAocc::simMsOcc(N = N, J = J, n_species = n_species,
                          n_occ_covs = 1, n_det_covs = 1, seed = N + 2000)

x1 <- sim$data$occ.covs$occ_x1
X_occ <- cbind(1, x1)
det_vals <- as.vector(sim$data$det.covs$det_x1)
X_det <- cbind(1, det_vals)

y_array <- array(NA_integer_, dim = c(n_species, N, J))
for (s in seq_len(n_species)) {
  sp <- sim$data$species_names[s]
  y_array[s, , ] <- sim$data$species_data[[sp]]$y
}

stan_data <- list(S = n_species, N = N, J = J, y = y_array,
                  K_occ = 2L, K_det = 2L, X_occ = X_occ, X_det = X_det)

cat("  Stan          ... "); flush.console()
t0 <- proc.time()
fit <- stan_ms$sample(data = stan_data, chains = 1,
                      iter_warmup = 1000, iter_sampling = 1000,
                      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
elapsed <- (proc.time() - t0)[["elapsed"]]

cors <- vapply(seq_len(n_species), function(s) {
  vars <- paste0("psi[", s, ",", seq_len(N), "]")
  psi_draws <- fit$draws(variables = vars, format = "matrix")
  cor(sim$truth$psi[s, ], colMeans(psi_draws))
}, numeric(1))
cor_avg <- mean(cors, na.rm = TRUE)
cat(sprintf("%6.1fs  (cor = %.3f)\n", elapsed, cor_avg))

ms <- rbind(ms, data.frame(N = N, type = "Multi-species", method = "Stan",
                           time = elapsed, cor_psi = round(cor_avg, 4),
                           stringsAsFactors = FALSE))
write.csv(ms, ms_file, row.names = FALSE)

cat("\nDone.\n")
