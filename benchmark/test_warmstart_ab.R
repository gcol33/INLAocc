devtools::load_all("C:/Users/Gilles Colling/Documents/dev/INLAocc")
set.seed(42)

for (N in c(100, 300, 1000)) {
  sim <- simulate_occu(N = N, J = 4, n_occ_covs = 1, n_det_covs = 1, seed = 1)

  # Warm-start (current code)
  t1 <- system.time({
    fit1 <- occu(~ occ_x1, ~ det_x1, data = sim$data, verbose = 0)
  })["elapsed"]

  cat(sprintf("N=%5d | %.1f s | occ=[%s] det=[%s]\n",
              N, t1,
              paste(round(fit1$occ_fit$summary.fixed$mean, 3), collapse=", "),
              paste(round(fit1$det_fit$summary.fixed$mean, 3), collapse=", ")))
}
