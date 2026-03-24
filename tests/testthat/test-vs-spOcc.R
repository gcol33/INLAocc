test_that("INLAocc matches spOccupancy on basic occupancy model", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("spOccupancy")

  set.seed(42)
  sim <- spOccupancy::simOcc(
    J.x = 10, J.y = 10, n.rep = rep(4, 100),
    beta = c(0.5, 0.8), alpha = c(0.2, -0.5)
  )

  data_list <- list(
    y        = sim$y,
    occ.covs = data.frame(occ_cov1 = sim$X[, 2]),
    det.covs = list(det_cov1 = sim$X.p[, , 2])
  )

  # Fit with spOccupancy
  fit_sp <- spOccupancy::PGOcc(
    occ.formula = ~ occ_cov1, det.formula = ~ det_cov1,
    data = data_list, n.samples = 10000, n.burn = 5000,
    n.thin = 5, n.chains = 1, verbose = FALSE
  )

  # Fit with INLAocc
  fit_inla <- occu_inla(
    occ.formula = ~ occ_cov1, det.formula = ~ det_cov1,
    data = data_list, verbose = 0
  )

  psi_true  <- as.vector(sim$psi)
  psi_sp    <- apply(fit_sp$psi.samples, 2, mean)
  psi_inla  <- fit_inla$psi_hat

  # Both should correlate well with truth
  expect_gt(cor(psi_inla, psi_true), 0.95)
  expect_gt(cor(psi_sp, psi_true), 0.95)

  # RMSE should be comparable (within 2x)
  rmse_sp   <- sqrt(mean((psi_sp - psi_true)^2))
  rmse_inla <- sqrt(mean((psi_inla - psi_true)^2))
  expect_lt(rmse_inla / rmse_sp, 2)

  # Fixed effects: both within 0.3 of truth
  beta_sp   <- apply(fit_sp$beta.samples, 2, mean)
  beta_inla <- fit_inla$occ_fit$summary.fixed$mean
  expect_lt(max(abs(beta_inla - c(0.5, 0.8))), 0.3)
  expect_lt(max(abs(beta_sp - c(0.5, 0.8))), 0.3)
})


test_that("INLAocc recovers random effect variance", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("spOccupancy")

  set.seed(123)
  sim <- spOccupancy::simOcc(
    J.x = 10, J.y = 10, n.rep = rep(4, 100),
    beta = c(0.5, 0.8), alpha = c(0.2, -0.5),
    psi.RE = list(levels = 5, sigma.sq.psi = 0.5)
  )

  data_list <- list(
    y        = sim$y,
    occ.covs = data.frame(occ_cov1 = sim$X[, 2], group = sim$X.re[, 1]),
    det.covs = list(det_cov1 = sim$X.p[, , 2])
  )

  fit <- occu_inla(
    occ.formula = ~ occ_cov1 + (1 | group),
    det.formula = ~ det_cov1,
    data = data_list, verbose = 0
  )

  # RE variance should be in reasonable range (true = 0.5)
  prec <- fit$occ_fit$summary.hyperpar$mean[1]
  sigma_sq <- 1 / prec
  expect_gt(sigma_sq, 0.05)
  expect_lt(sigma_sq, 5.0)
})


test_that("data formatting accepts spOccupancy-style lists", {
  y <- matrix(rbinom(200, 1, 0.4), 50, 4)
  covs <- data.frame(elev = rnorm(50))
  det_covs <- list(effort = matrix(runif(200), 50, 4))

  dat <- occu_format(y, covs, det_covs)
  expect_s3_class(dat, "occu_data")
  expect_equal(dat$N, 50)
  expect_equal(dat$J, 4)
})


test_that("3D array accepted for multi-species", {
  y_arr <- array(rbinom(300, 1, 0.3), dim = c(3, 50, 2))
  dimnames(y_arr)[[1]] <- c("sp_a", "sp_b", "sp_c")

  dat <- occu_format_ms(y_arr, data.frame(x = rnorm(50)))
  expect_s3_class(dat, "occu_data_ms")
  expect_equal(dat$n_species, 3)
  expect_equal(dat$species_names, c("sp_a", "sp_b", "sp_c"))
})


test_that("formula parser extracts lme4-style random effects", {
  parsed <- parse_re_formula(~ x1 + x2 + (1 | group) + (x1 | region))

  expect_length(parsed$re_list, 2)
  expect_equal(parsed$re_list[[1]]$type, "intercept")
  expect_equal(parsed$re_list[[1]]$group, "group")
  expect_equal(parsed$re_list[[2]]$type, "slope")
  expect_equal(parsed$re_list[[2]]$group, "region")
  expect_equal(parsed$re_list[[2]]$covariate, "x1")

  # Fixed part should not contain RE terms
  fixed_terms <- attr(terms(parsed$fixed), "term.labels")
  expect_true(all(c("x1", "x2") %in% fixed_terms))
  expect_false(any(grepl("\\|", fixed_terms)))
})


test_that("simulate_occu produces valid data", {
  sim <- simulate_occu(N = 50, J = 3, seed = 1)

  expect_s3_class(sim$data, "occu_data")
  expect_equal(sim$data$N, 50)
  expect_equal(sim$data$J, 3)
  expect_true(all(sim$data$y %in% c(0, 1)))
  expect_equal(length(sim$truth$z), 50)
  expect_true(all(sim$truth$psi >= 0 & sim$truth$psi <= 1))
})
