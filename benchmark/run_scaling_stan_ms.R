# =============================================================================
# Add Stan to multi-species benchmark (small N only -- Stan is slow here)
# Appends to scaling_results_ms.csv
# =============================================================================

library(cmdstanr)
library(INLAocc)

set.seed(42)
ms_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results_ms.csv"
ms <- read.csv(ms_file, stringsAsFactors = FALSE)

# Remove any previous Stan MS rows
ms <- ms[!(ms$method == "Stan" & ms$type == "Multi-species"), ]

cat("Compiling multi-species Stan model ... ")
stan_ms <- cmdstan_model(
  "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/occu_ms_model.stan",
  quiet = TRUE
)
cat("done.\n\n")

n_species <- 10L
J <- 4L

for (N in c(100, 300, 1000)) {
  cat(sprintf("N = %d, S = %d\n", N, n_species))

  sim <- INLAocc::simMsOcc(N = N, J = J, n_species = n_species,
                            n_occ_covs = 1, n_det_covs = 1, seed = N + 2000)

  # Build Stan data
  x1 <- sim$data$occ.covs$occ_x1
  X_occ <- cbind(1, x1)
  det_vals <- as.vector(sim$data$det.covs$det_x1)
  X_det <- cbind(1, det_vals)

  y_array <- array(NA_integer_, dim = c(n_species, N, J))
  for (s in seq_len(n_species)) {
    sp <- sim$data$species_names[s]
    y_array[s, , ] <- sim$data$species_data[[sp]]$y
  }

  stan_data <- list(
    S = n_species, N = N, J = J, y = y_array,
    K_occ = 2L, K_det = 2L,
    X_occ = X_occ, X_det = X_det
  )

  cat("  Stan          ... "); flush.console()
  t0 <- proc.time()
  fit <- tryCatch(
    stan_ms$sample(
      data = stan_data,
      chains = 1, iter_warmup = 1000, iter_sampling = 1000,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE
    ),
    error = function(e) { message("  Stan error: ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  cor_avg <- NA_real_
  if (!is.null(fit)) {
    cors <- vapply(seq_len(n_species), function(s) {
      vars <- paste0("psi[", s, ",", seq_len(N), "]")
      psi_draws <- fit$draws(variables = vars, format = "matrix")
      psi_hat <- colMeans(psi_draws)
      cor(sim$truth$psi[s, ], psi_hat)
    }, numeric(1))
    cor_avg <- mean(cors, na.rm = TRUE)
  }
  cat(sprintf("%6.1fs  (cor = %.3f)\n", elapsed, cor_avg))

  ms <- rbind(ms, data.frame(
    N = N, type = "Multi-species", method = "Stan",
    time = elapsed, cor_psi = round(cor_avg, 4),
    stringsAsFactors = FALSE
  ))
  write.csv(ms, ms_file, row.names = FALSE)
}

cat("\nDone.\n")
