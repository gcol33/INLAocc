# =============================================================================
# Multi-species scaling benchmark: INLAocc vs spOccupancy
# Saves to benchmark/scaling_results_ms.csv
# =============================================================================

library(spOccupancy)
devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(123)
out_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results_ms.csv"

cat("=== Multi-species Scaling Benchmark ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
N_grid    <- c(100, 300, 1000, 3000)
n_species <- 10L
J         <- 4L

# -----------------------------------------------------------------------------
# Fitting wrappers
# -----------------------------------------------------------------------------
time_inlaocc_ms <- function(sim) {
  t0 <- proc.time()
  fit <- tryCatch(
    occu(~ occ_x1, ~ det_x1, data = sim$data, multispecies = TRUE, verbose = 0),
    error = function(e) { message("  INLAocc error: ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  # Average correlation across species
  cor_avg <- NA_real_
  if (!is.null(fit) && !is.null(fit$species_fits)) {
    cors <- vapply(sim$data$species_names, function(sp) {
      sp_fit <- fit$species_fits[[sp]]
      if (!is.null(sp_fit) && !is.null(sp_fit$psi_hat)) {
        idx <- match(sp, sim$data$species_names)
        cor(sim$truth$psi[idx, ], sp_fit$psi_hat)
      } else NA_real_
    }, numeric(1))
    cor_avg <- mean(cors, na.rm = TRUE)
  }
  list(time = elapsed, cor_psi = cor_avg)
}

time_spocc_ms <- function(sim) {
  N <- sim$data$N

  # Build 3D array from species_data
  y_array <- array(NA, dim = c(n_species, N, J))
  for (s in seq_len(n_species)) {
    sp_name <- sim$data$species_names[s]
    y_array[s, , ] <- sim$data$species_data[[sp_name]]$y
  }

  spocc_data <- list(
    y        = y_array,
    occ.covs = sim$data$occ.covs,
    det.covs = sim$data$det.covs
  )

  t0 <- proc.time()
  fit <- tryCatch(
    msPGOcc(
      occ.formula = ~ occ_x1, det.formula = ~ det_x1,
      data = spocc_data,
      n.samples = 5000, n.burn = 2000, n.thin = 3,
      n.omp.threads = 1, verbose = FALSE
    ),
    error = function(e) { message("  spOcc error: ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  cor_avg <- NA_real_
  if (!is.null(fit)) {
    psi_hat <- apply(fit$psi.samples, c(2, 3), mean)   # species x N
    cors <- vapply(seq_len(n_species), function(s) {
      cor(sim$truth$psi[s, ], psi_hat[s, ])
    }, numeric(1))
    cor_avg <- mean(cors, na.rm = TRUE)
  }
  list(time = elapsed, cor_psi = cor_avg)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
results <- data.frame(
  N = integer(), type = character(), method = character(),
  time = numeric(), cor_psi = numeric(),
  stringsAsFactors = FALSE
)

for (N in N_grid) {
  cat(sprintf("N = %d, S = %d\n", N, n_species))

  sim <- simMsOcc(N = N, J = J, n_species = n_species,
                  n_occ_covs = 1, n_det_covs = 1, seed = N + 2000)

  cat("  INLAocc       ... "); flush.console()
  r <- time_inlaocc_ms(sim)
  cat(sprintf("%6.1fs  (cor = %.3f)\n", r$time, r$cor_psi))
  results <- rbind(results, data.frame(
    N = N, type = "Multi-species", method = "INLAocc",
    time = r$time, cor_psi = round(r$cor_psi, 4),
    stringsAsFactors = FALSE
  ))

  cat("  spOccupancy   ... "); flush.console()
  r <- time_spocc_ms(sim)
  cat(sprintf("%6.1fs  (cor = %.3f)\n", r$time, r$cor_psi))
  results <- rbind(results, data.frame(
    N = N, type = "Multi-species", method = "spOccupancy",
    time = r$time, cor_psi = round(r$cor_psi, 4),
    stringsAsFactors = FALSE
  ))

  cat("\n")
  write.csv(results, out_file, row.names = FALSE)
}

cat(sprintf("Finished: %s\n", Sys.time()))
cat(sprintf("Results: %s\n", out_file))
print(results)
