# =============================================================================
# Push INLAocc to the 10-minute wall
# Non-spatial: 300k, 1M, 3M sites
# Multi-species: 30k, 100k sites (10 species)
# =============================================================================

devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(2026)
root <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark"
J <- 4L
beta_occ <- c(0.5, -0.8)
beta_det <- c(0.0, 0.5)

simulate_data <- function(N, J, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(N)
  occ_covs <- data.frame(occ_x1 = x1)
  det_mat <- matrix(rnorm(N * J), N, J)
  det_covs <- list(det_x1 = det_mat)
  eta_occ <- cbind(1, x1) %*% beta_occ
  psi <- plogis(as.vector(eta_occ))
  z <- rbinom(N, 1, psi)
  X_det <- cbind(1, as.vector(det_mat))
  p_vec <- plogis(as.vector(X_det %*% beta_det))
  y_vec <- rbinom(N * J, 1, rep(z, each = J) * p_vec)
  y <- matrix(y_vec, N, J, byrow = TRUE)
  list(y = y, occ.covs = occ_covs, det.covs = det_covs, psi_true = psi)
}

ss_file <- file.path(root, "scaling_results.csv")
ms_file <- file.path(root, "scaling_results_ms.csv")
ss <- read.csv(ss_file, stringsAsFactors = FALSE)
ms <- read.csv(ms_file, stringsAsFactors = FALSE)

cat("=== Push INLAocc to 10-min wall ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# ---- Non-spatial: 30k, 100k, 300k, 1M ----
cat("--- Non-spatial (INLAocc only) ---\n")
for (N in c(30000L, 100000L, 300000L, 1000000L)) {
  cat(sprintf("\nN = %s\n", format(N, big.mark = ",")))
  data <- simulate_data(N, J, seed = N)

  cat("  INLAocc       ... "); flush.console()
  t0 <- proc.time()
  fit <- tryCatch(
    occu(~ occ_x1, ~ det_x1,
         data = list(y = data$y, occ.covs = data$occ.covs, det.covs = data$det.covs),
         verbose = 0),
    error = function(e) { message("  error: ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  r <- if (!is.null(fit)) cor(data$psi_true, fit$psi_hat) else NA
  cat(sprintf("%6.1fs  (cor = %.4f)\n", elapsed, r))

  ss <- rbind(ss, data.frame(N = N, type = "Non-spatial", method = "INLAocc",
                             time = elapsed, cor_psi = round(r, 4),
                             stringsAsFactors = FALSE))
  write.csv(ss, ss_file, row.names = FALSE)

  if (elapsed > 600) {
    cat("  Exceeded 10 min — stopping.\n")
    break
  }
}

# ---- Multi-species: 30k, 100k ----
cat("\n--- Multi-species (INLAocc only, 10 species) ---\n")
n_species <- 10L
for (N in c(30000L, 100000L)) {
  cat(sprintf("\nN = %s, S = %d\n", format(N, big.mark = ","), n_species))

  cat("  Simulating... "); flush.console()
  sim <- simMsOcc(N = N, J = J, n_species = n_species,
                  n_occ_covs = 1, n_det_covs = 1, seed = N + 5000)
  cat("done.\n")

  cat("  INLAocc       ... "); flush.console()
  t0 <- proc.time()
  fit <- tryCatch(
    occu(~ occ_x1, ~ det_x1, data = sim$data, multispecies = TRUE, verbose = 0),
    error = function(e) { message("  error: ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  cor_avg <- NA_real_
  if (!is.null(fit) && !is.null(fit$species_fits)) {
    cors <- vapply(sim$data$species_names, function(sp) {
      sf <- fit$species_fits[[sp]]
      if (!is.null(sf)) cor(sim$truth$psi[match(sp, sim$data$species_names), ], sf$psi_hat)
      else NA_real_
    }, numeric(1))
    cor_avg <- mean(cors, na.rm = TRUE)
  }
  cat(sprintf("%6.1fs  (cor = %.4f)\n", elapsed, cor_avg))

  ms <- rbind(ms, data.frame(N = N, type = "Multi-species", method = "INLAocc",
                             time = elapsed, cor_psi = round(cor_avg, 4),
                             stringsAsFactors = FALSE))
  write.csv(ms, ms_file, row.names = FALSE)

  if (elapsed > 600) {
    cat("  Exceeded 10 min — stopping.\n")
    break
  }
}

cat(sprintf("\nFinished: %s\n", Sys.time()))
