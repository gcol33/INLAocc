# =============================================================================
# run_demo.R — Full demonstration of the INLA occupancy framework
# =============================================================================
# Exercises all features: basic model, mixed-model formula syntax, random intercepts
# and slopes, spatial SPDE, multi-species, temporal, integrated, latent
# factors, SVC, GOF, WAIC, k-fold CV, marginal effects, prediction.
# =============================================================================

cat("Loading INLA occupancy framework...\n")
base_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else getwd()
  }
)

for (f in c("occu_utils.R", "occu_data.R", "occu_random.R", "occu_spatial.R",
            "occu_engine.R", "occu_predict.R", "occu_output.R",
            "occu_diagnostics.R", "occu_fit.R")) {
  source(file.path(base_dir, f))
}

# ============================================================================
# 1. Simulate + fit basic model
# ============================================================================
cat("\n========================================\n")
cat("1. Basic occupancy model\n")
cat("========================================\n")

sim <- simulate_occu(
  N = 200, J = 4, n_occ_covs = 2, n_det_covs = 1,
  beta_occ = c(0.5, 0.8, -0.6), beta_det = c(0.2, -0.4),
  seed = 42
)
print(sim$data)

fit_basic <- occu_inla(
  occ.formula = ~ occ_x1 + occ_x2,
  det.formula = ~ det_x1,
  data = sim$data
)
summary(fit_basic)

cat("\nTrue vs Estimated:\n")
cat(sprintf("  Occ betas: true = [%.2f, %.2f, %.2f] | est = [%.3f, %.3f, %.3f]\n",
            sim$truth$beta_occ[1], sim$truth$beta_occ[2], sim$truth$beta_occ[3],
            fit_basic$occ_fit$summary.fixed$mean[1],
            fit_basic$occ_fit$summary.fixed$mean[2],
            fit_basic$occ_fit$summary.fixed$mean[3]))

# ============================================================================
# 2. Mixed-model formula with random effects
# ============================================================================
cat("\n========================================\n")
cat("2. Mixed-model formula syntax: (1 | group) and (x | group)\n")
cat("========================================\n")

sim$data$occ.covs$region <- sample(1:5, sim$data$N, replace = TRUE)

# Random intercept via formula (just like spOccupancy)
fit_re <- occu_inla(
  occ.formula = ~ occ_x1 + occ_x2 + (1 | region),
  det.formula = ~ det_x1,
  data = sim$data
)
summary(fit_re)

# Random slope via formula
fit_rs <- occu_inla(
  occ.formula = ~ occ_x1 + occ_x2 + (1 + occ_x1 | region),
  det.formula = ~ det_x1,
  data = sim$data
)
summary(fit_rs)

# ============================================================================
# 3. With priors (spOccupancy-compatible syntax)
# ============================================================================
cat("\n========================================\n")
cat("3. Custom priors\n")
cat("========================================\n")

fit_priors <- occu_inla(
  occ.formula = ~ occ_x1 + occ_x2,
  det.formula = ~ det_x1,
  data = sim$data,
  priors = list(
    beta.normal  = list(mean = 0, var = 2.72),
    alpha.normal = list(mean = 0, var = 2.72)
  )
)
cat("Fitted with custom priors.\n")

# ============================================================================
# 4. Spatial model (SPDE)
# ============================================================================
cat("\n========================================\n")
cat("4. Spatial occupancy model (SPDE)\n")
cat("========================================\n")

sim_sp <- simulate_occu(
  N = 200, J = 4, n_occ_covs = 1, n_det_covs = 1,
  beta_occ = c(0.3, 0.5), beta_det = c(0.0, -0.3),
  spatial_range = 0.2, spatial_var = 1.0, seed = 123
)

fit_spatial <- spatial_occu_inla(
  occ.formula = ~ occ_x1,
  det.formula = ~ det_x1,
  data   = sim_sp$data,
  coords = sim_sp$data$coords,
  spde.args = list(prior.range = c(0.3, 0.5), prior.sigma = c(1, 0.05))
)
summary(fit_spatial)

# ============================================================================
# 5. Multi-species model
# ============================================================================
cat("\n========================================\n")
cat("5. Multi-species model\n")
cat("========================================\n")

sim_ms <- simMsOcc(
  N = 100, J = 3, n_species = 5,
  n_occ_covs = 1, n_det_covs = 1,
  seed = 200
)

fit_ms <- ms_occu_inla(
  occ.formula = ~ occ_x1,
  det.formula = ~ det_x1,
  data = sim_ms$data
)
summary(fit_ms)

rich <- richness(fit_ms)
cat("\nSpecies richness (first 10 sites):\n")
print(head(rich, 10))

# ============================================================================
# 6. Multi-season (temporal) model
# ============================================================================
cat("\n========================================\n")
cat("6. Multi-season model (AR1)\n")
cat("========================================\n")

sim_t <- simTOcc(
  N = 80, J = 3, n_seasons = 4,
  beta_occ = c(0.5, 0.3), beta_det = c(0, -0.3),
  rho = 0.7, sigma_t = 0.5, seed = 300
)

fit_temporal <- temporal_occu_inla(
  occ.formula = ~ occ_x1,
  det.formula = ~ det_x1,
  data = sim_t$data,
  ar1 = TRUE
)

cat("Temporal model fitted.\n")
cat(sprintf("  Seasons: %d\n", fit_temporal$n_periods))
for (t in seq_along(fit_temporal$period_fits)) {
  pf <- fit_temporal$period_fits[[t]]
  if (!is.null(pf)) {
    cat(sprintf("  Season %d: mean psi = %.3f, converged = %s\n",
                t, mean(pf$psi_hat), pf$converged))
  }
}

# ============================================================================
# 7. Integrated (multi-source) model
# ============================================================================
cat("\n========================================\n")
cat("7. Integrated (data fusion) model\n")
cat("========================================\n")

sim_int <- simIntOcc(
  N_total = 150, n_data = 2, J = c(4, 3),
  n_shared = 20,
  beta_occ = c(0.5, 0.3),
  beta_det = list(c(0.2, -0.4), c(-0.1, 0.3)),
  seed = 400
)

fit_int <- intOccu_inla(
  occ.formula = ~ occ_x1,
  det.formula = list(~ det_x1, ~ det_x1),
  data = sim_int$data
)
cat(sprintf("Integrated model: converged = %s, mean psi = %.3f\n",
            fit_int$converged, mean(fit_int$psi_hat)))

# ============================================================================
# 8. Diagnostics: fitted, residuals, GOF, WAIC
# ============================================================================
cat("\n========================================\n")
cat("8. Diagnostics\n")
cat("========================================\n")

# Fitted values
fv <- fitted(fit_basic)
cat(sprintf("Mean fitted y.rep: %.3f\n", mean(fv$y.rep, na.rm = TRUE)))

# Residuals
res <- residuals(fit_basic, type = "deviance")
cat(sprintf("Mean |occ residual|: %.3f\n", mean(abs(res$occ.resids))))
cat(sprintf("Mean |det residual|: %.3f\n", mean(abs(res$det.resids), na.rm = TRUE)))

# Posterior predictive check
gof <- ppcOccu(fit_basic, fit.stat = "freeman-tukey", n.samples = 200)
cat(sprintf("Bayesian p-value: %.3f\n", gof$bayesian.p))

# WAIC
waic <- waicOccu(fit_basic)
cat("WAIC:\n")
print(waic)

# ============================================================================
# 9. K-fold cross-validation
# ============================================================================
cat("\n========================================\n")
cat("9. K-fold CV\n")
cat("========================================\n")

fit_cv <- occu_inla(
  occ.formula = ~ occ_x1 + occ_x2,
  det.formula = ~ det_x1,
  data   = sim$data,
  k.fold = 3,
  verbose = 1
)
cat(sprintf("K-fold deviance: %.2f\n", fit_cv$k.fold$k.fold.deviance))

# ============================================================================
# 10. Prediction (X.0 design matrix)
# ============================================================================
cat("\n========================================\n")
cat("10. Out-of-sample prediction\n")
cat("========================================\n")

# Create prediction design matrix (with intercept)
X_new <- cbind(1, seq(-2, 2, length.out = 50), 0)  # vary occ_x1, fix occ_x2 = 0
pred <- predict(fit_basic, X.0 = X_new, type = "occupancy", n.samples = 100)
cat(sprintf("Prediction range: [%.3f, %.3f]\n",
            min(pred$psi.0$mean), max(pred$psi.0$mean)))

# Marginal effect
me <- marginal_effect(fit_basic, "occ_x1", process = "occupancy")
cat("Marginal effect of occ_x1 (first 5 rows):\n")
print(head(me, 5))

# ============================================================================
# 11. Model comparison
# ============================================================================
cat("\n========================================\n")
cat("11. Model comparison\n")
cat("========================================\n")

comp <- compare_models(
  basic = fit_basic,
  re    = fit_re,
  rs    = fit_rs
)
print(comp)

# ============================================================================
# 12. Plots (if interactive)
# ============================================================================
if (interactive()) {
  cat("\n========================================\n")
  cat("12. Diagnostic plots\n")
  cat("========================================\n")

  par(mfrow = c(2, 2))
  plot(fit_basic)
  par(mfrow = c(1, 1))

  # GOF plot
  plot(gof$fit.y, gof$fit.y.rep, pch = 16, col = rgb(0, 0, 0, 0.3),
       xlab = "Observed fit stat", ylab = "Replicated fit stat",
       main = sprintf("PPC (Bayesian p = %.3f)", gof$bayesian.p))
  abline(0, 1, col = "red", lwd = 2)
}

cat("\nDemo complete.\n")
