# Rerun INLAocc-only rows in scaling_results.csv, keeping competitor data intact
devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")

set.seed(42)
root <- "C:/Users/Gilles Colling/Documents/dev/INLAocc/benchmark"
out_file <- file.path(root, "scaling_results.csv")

df <- read.csv(out_file, stringsAsFactors = FALSE)

# --- simulation function (same as run_scaling.R) ---
beta_occ <- c(0.5, -0.8)
beta_det <- c(0.0,  0.5)
J <- 4L

simulate_data <- function(N, J, spatial = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(N)
  X_occ <- cbind(1, x1)
  occ_covs <- data.frame(occ_x1 = x1)
  det_vals <- rnorm(N * J)
  det_mat  <- matrix(det_vals, N, J)
  det_covs <- list(det_x1 = det_mat)
  coords <- cbind(x = runif(N), y = runif(N))
  eta_occ <- X_occ %*% beta_occ
  if (spatial) {
    D     <- as.matrix(dist(coords))
    Sigma <- 1.0 * exp(-D / 0.2)
    w     <- MASS::mvrnorm(1, rep(0, N), Sigma)
    eta_occ <- eta_occ + w
  }
  psi <- plogis(as.vector(eta_occ))
  z   <- rbinom(N, 1, psi)
  eta_det <- cbind(1, det_vals) %*% beta_det
  p_vec   <- plogis(as.vector(eta_det))
  z_rep   <- rep(z, each = J)
  y_vec   <- rbinom(N * J, 1, z_rep * p_vec)
  y       <- matrix(y_vec, N, J, byrow = TRUE)
  list(y = y, occ.covs = occ_covs, det.covs = det_covs,
       coords = coords, psi_true = psi)
}

# --- rerun INLAocc rows ---
inla_rows <- which(df$method == "INLAocc")
for (i in inla_rows) {
  N <- df$N[i]
  type <- df$type[i]
  spatial <- type == "Spatial"
  seed <- if (spatial) N + 1000 else N

  cat(sprintf("%s N=%d ... ", type, N)); flush.console()
  data <- simulate_data(N, J, spatial = spatial, seed = seed)

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
  }, error = function(e) { message("error: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]

  psi_hat <- if (!is.null(fit)) fit$psi_hat else NULL
  cor_val <- if (!is.null(psi_hat)) round(cor(data$psi_true, psi_hat), 4) else NA
  cat(sprintf("%.1fs (cor=%.4f)\n", elapsed, cor_val))

  df$time[i] <- elapsed
  df$cor_psi[i] <- cor_val
}

write.csv(df, out_file, row.names = FALSE)
cat("\nDone. Updated INLAocc rows in scaling_results.csv\n")
print(df)
