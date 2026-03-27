# =============================================================================
# Parameter recovery benchmark: estimated vs true coefficients
# Multi-species at N=1000, all three methods
# Saves benchmark/param_recovery.csv
# =============================================================================

devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")
library(spOccupancy)
library(cmdstanr)

set.seed(42)
out_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/param_recovery.csv"

N <- 1000L; J <- 4L; n_species <- 10L

cat("Simulating multi-species data (N=1000, S=10) ...\n")
sim <- INLAocc::simMsOcc(N = N, J = J, n_species = n_species,
                          n_occ_covs = 1, n_det_covs = 1, seed = 999)

results <- data.frame(
  species = character(), param = character(), true_val = numeric(),
  estimate = numeric(), method = character(), stringsAsFactors = FALSE
)

add_coefs <- function(df, method, species, occ_est, det_est, occ_true, det_true) {
  occ_names <- c("occ_intercept", paste0("occ_beta", seq_len(length(occ_true) - 1)))
  det_names <- c("det_intercept", paste0("det_beta", seq_len(length(det_true) - 1)))
  rbind(df,
    data.frame(species = species, param = occ_names, true_val = occ_true,
               estimate = occ_est, method = method, stringsAsFactors = FALSE),
    data.frame(species = species, param = det_names, true_val = det_true,
               estimate = det_est, method = method, stringsAsFactors = FALSE)
  )
}

# ---- INLAocc ----
cat("INLAocc ... "); flush.console()
fit_inla <- occu(~ occ_x1, ~ det_x1, data = sim$data, multispecies = TRUE, verbose = 0)
for (s in seq_len(n_species)) {
  sp <- sim$data$species_names[s]
  sf <- fit_inla$species_fits[[sp]]
  if (!is.null(sf)) {
    results <- add_coefs(results, "INLAocc", sp,
      sf$occ_fit$summary.fixed$mean, sf$det_fit$summary.fixed$mean,
      sim$truth$betas[s, ], sim$truth$alphas[s, ])
  }
}
cat("done.\n")

# ---- spOccupancy ----
cat("spOccupancy ... "); flush.console()
y_array <- array(NA_integer_, dim = c(n_species, N, J))
for (s in seq_len(n_species)) {
  sp <- sim$data$species_names[s]
  y_array[s, , ] <- sim$data$species_data[[sp]]$y
}
spocc_data <- list(y = y_array, occ.covs = sim$data$occ.covs, det.covs = sim$data$det.covs)

fit_spocc <- msPGOcc(occ.formula = ~ occ_x1, det.formula = ~ det_x1,
                     data = spocc_data,
                     n.samples = 5000, n.burn = 2000, n.thin = 3,
                     n.omp.threads = 1, verbose = FALSE)
# beta.samples is (n_post, n_species * K_occ) mcmc object
# Column order: param1_sp1, param1_sp2, ..., param1_spS, param2_sp1, ...
beta_means <- colMeans(as.matrix(fit_spocc$beta.samples))
alpha_means <- colMeans(as.matrix(fit_spocc$alpha.samples))
K_occ <- ncol(sim$truth$betas)
K_det <- ncol(sim$truth$alphas)
beta_mat <- matrix(beta_means, nrow = n_species, ncol = K_occ, byrow = FALSE)
alpha_mat <- matrix(alpha_means, nrow = n_species, ncol = K_det, byrow = FALSE)
for (s in seq_len(n_species)) {
  results <- add_coefs(results, "spOccupancy", sim$data$species_names[s],
    beta_mat[s, ], alpha_mat[s, ],
    sim$truth$betas[s, ], sim$truth$alphas[s, ])
}
cat("done.\n")

# ---- Stan ----
cat("Stan ... "); flush.console()
x1 <- sim$data$occ.covs$occ_x1
X_occ <- cbind(1, x1)
det_vals <- as.vector(sim$data$det.covs$det_x1)
X_det <- cbind(1, det_vals)

stan_ms <- cmdstan_model(
  "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/occu_ms_model.stan",
  quiet = TRUE
)
stan_data <- list(S = n_species, N = N, J = J, y = y_array,
                  K_occ = 2L, K_det = 2L, X_occ = X_occ, X_det = X_det)

fit_stan <- stan_ms$sample(data = stan_data, chains = 1,
                           iter_warmup = 1000, iter_sampling = 1000,
                           refresh = 0, show_messages = FALSE,
                           show_exceptions = FALSE)
for (s in seq_len(n_species)) {
  beta_vars <- paste0("beta[", s, ",", 1:2, "]")
  alpha_vars <- paste0("alpha[", s, ",", 1:2, "]")
  beta_est <- colMeans(fit_stan$draws(variables = beta_vars, format = "matrix"))
  alpha_est <- colMeans(fit_stan$draws(variables = alpha_vars, format = "matrix"))
  results <- add_coefs(results, "Stan", sim$data$species_names[s],
    beta_est, alpha_est, sim$truth$betas[s, ], sim$truth$alphas[s, ])
}
cat("done.\n")

write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\nSaved %d rows to %s\n", nrow(results), out_file))
