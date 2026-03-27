# =============================================================================
# Multi-species parallel benchmark: real wall-clock measurements
# Requires INLAocc to be installed (not just load_all'd)
# =============================================================================

library(INLAocc)
library(spOccupancy)

n_cores <- min(parallel::detectCores(logical = FALSE), 10L)
cat(sprintf("Using %d cores for parallel species fitting\n", n_cores))
options(INLAocc.cores = n_cores)

set.seed(123)
ms_file <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark/scaling_results_ms.csv"

# Load existing, drop any previous parallel rows
ms <- read.csv(ms_file, stringsAsFactors = FALSE)
ms <- ms[ms$type != "Multi-species (parallel)", ]

n_species <- 10L
J <- 4L

for (N in c(1000, 3000, 10000)) {
  cat(sprintf("\nN = %d, S = %d, cores = %d\n", N, n_species, n_cores))

  sim <- INLAocc::simMsOcc(N = N, J = J, n_species = n_species,
                           n_occ_covs = 1, n_det_covs = 1, seed = N + 2000)

  cat("  INLAocc (parallel) ... "); flush.console()
  t0 <- proc.time()
  fit <- tryCatch(
    occu(~ occ_x1, ~ det_x1, data = sim$data, multispecies = TRUE, verbose = 0),
    error = function(e) { message("  error: ", e$message); NULL }
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

  ms <- rbind(ms, data.frame(
    N = N, type = "Multi-species (parallel)",
    method = "INLAocc", time = elapsed,
    cor_psi = round(cor_avg, 4), stringsAsFactors = FALSE
  ))
  write.csv(ms, ms_file, row.names = FALSE)
}

cat("\nDone.\n")
print(ms[ms$type == "Multi-species (parallel)", ])
