devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")
set.seed(42)

sim <- simulate_occu(N = 100, J = 4, n_occ_covs = 1, n_det_covs = 1, seed = 1)

cat("--- N=100 timing ---\n")
t0 <- system.time({
  fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, verbose = 1)
})
cat(sprintf("Total: %.1f s\n", t0["elapsed"]))
cat(sprintf("Occ intercept: %.3f (true: %.3f)\n",
            fit$occ_fit$summary.fixed$mean[1], sim$truth$beta[1]))
cat(sprintf("Det intercept: %.3f (true: %.3f)\n",
            fit$det_fit$summary.fixed$mean[1], sim$truth$alpha[1]))
