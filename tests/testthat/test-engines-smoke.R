# =============================================================================
# Smoke tests for all engine types
#
# Each test verifies the engine runs without error on minimal simulated data
# and returns the expected output structure. NOT accuracy tests — just ensures
# no regressions in basic functionality.
# =============================================================================

skip_if_not_installed("INLA")

# --- Helper: minimal simulation sizes ---
N_SMALL <- 30
J_SMALL <- 3
N_SP    <- 3
N_SEAS  <- 2


# ==========================================================================
# Single-species engines
# ==========================================================================

test_that("engine_ss runs without error", {
  sim <- simulate_occu(N = N_SMALL, J = J_SMALL, seed = 1)
  fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, verbose = 0,
              max.iter = 5)
  expect_s3_class(fit, "occu_inla")
  expect_length(fit$psi_hat, N_SMALL)
})

test_that("engine_ss_spatial runs without error", {
  sim <- simulate_occu(N = N_SMALL, J = J_SMALL, seed = 3)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = coords,
              verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_spatial")
})

test_that("engine_temporal runs without error", {
  sim <- simTOcc(N = N_SMALL, J = J_SMALL, n_seasons = N_SEAS, seed = 4)
  fit <- occu(~ 1, ~ 1, data = sim$data, temporal = "ar1",
              verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_temporal")
  expect_equal(fit$n_periods, N_SEAS)
})

test_that("engine_svc runs without error", {
  sim <- simulate_occu(N = N_SMALL, J = J_SMALL, seed = 5)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = coords,
              svc = 1L, verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_svc")
})

test_that("engine_svc_temporal runs without error", {
  sim <- simTOcc(N = N_SMALL, J = J_SMALL, n_seasons = N_SEAS, seed = 10)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ 1, ~ 1, data = sim$data, temporal = "ar1",
              spatial = coords, svc = 1L, verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_svct")
})


# ==========================================================================
# Integrated engines
# ==========================================================================

test_that("engine_int runs without error", {
  sim <- simIntOcc(N_total = 40, n_data = 2, J = c(3, 3), seed = 6)
  fit <- occu(~ 1, ~ 1, data = sim$data, integrated = TRUE,
              verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_int")
})

test_that("engine_int_spatial runs without error", {
  sim <- simIntOcc(N_total = 40, n_data = 2, J = c(3, 3), seed = 11)
  fit <- occu(~ 1, ~ 1, data = sim$data, integrated = TRUE,
              spatial = sim$data$coords, verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_spint")
})

test_that("engine_stint runs without error", {
  # Build spatio-temporal integrated data: 2 periods × 2 sources
  set.seed(17)
  N <- 30
  coords <- cbind(runif(N), runif(N))
  occ_covs <- data.frame(x1 = rnorm(N))
  sites <- list(seq_len(N), seq_len(N))

  y_periods <- lapply(1:2, function(t) {
    lapply(1:2, function(d) matrix(rbinom(N * 3, 1, 0.4), N, 3))
  })
  det_covs <- lapply(1:2, function(t) {
    lapply(1:2, function(d) list(det_x1 = matrix(rnorm(N * 3), N, 3)))
  })

  data <- list(y = y_periods, occ.covs = occ_covs,
               det.covs = det_covs, sites = sites, coords = coords)
  fit <- occu(~ 1, ~ 1, data = data, integrated = TRUE,
              spatial = coords, temporal = "ar1", verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_stint")
})


# ==========================================================================
# Multi-species engines
# ==========================================================================

test_that("engine_ms runs without error", {
  sim <- simMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP, seed = 7)
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_ms")
  expect_equal(fit$n_species, N_SP)
  expect_false(is.null(fit$community_occ))
})

test_that("engine_ms with ensemble pooling runs without error", {
  sim <- simMsOcc(N = N_SMALL, J = J_SMALL, n_species = 5, seed = 8)
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              ensemble = TRUE, verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_ms")
  expect_true("n_subsets" %in% names(fit$community_occ))
})

test_that("engine_ms_lf runs without error", {
  sim <- simMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP, seed = 9)
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              n.factors = 1, verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_lfms")
  expect_false(is.null(fit$lambda))
})

test_that("engine_ms_sf runs without error", {
  sim <- simMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP, seed = 12)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              n.factors = 1, spatial = coords, verbose = 0, max.iter = 5)
  expect_s3_class(fit, "occu_inla_sfms")
})

test_that("engine_ms_temporal runs without error", {
  sim <- simTMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP,
                   n_seasons = N_SEAS, seed = 13)
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              temporal = "ar1", verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_tms")
  expect_equal(fit$n_species, N_SP)
})

test_that("engine_ms_st runs without error", {
  sim <- simTMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP,
                   n_seasons = N_SEAS, seed = 14)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              temporal = "ar1", spatial = coords, verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_stms")
})

test_that("engine_ms_int runs without error", {
  sim <- simIntMsOcc(N_total = 40, n_species = N_SP, n_data = 2,
                     J = c(3, 3), n_shared = 10, seed = 15)
  fit <- occu(~ 1, ~ 1, data = sim$data, multispecies = TRUE,
              integrated = TRUE, verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_intms")
  expect_equal(fit$n_species, N_SP)
})

test_that("engine_ms_svc_t runs without error", {
  sim <- simTMsOcc(N = N_SMALL, J = J_SMALL, n_species = N_SP,
                   n_seasons = N_SEAS, seed = 16)
  coords <- cbind(runif(N_SMALL), runif(N_SMALL))
  fit <- occu(~ occ_x1, ~ 1, data = sim$data, multispecies = TRUE,
              temporal = "ar1", spatial = coords, svc = 1L,
              verbose = 0, max.iter = 3)
  expect_s3_class(fit, "occu_inla_svctms")
})


# ==========================================================================
# JSDM engines
# ==========================================================================

test_that("engine_jsdm_lf runs without error", {
  n_sp <- 5
  N <- 40
  y <- matrix(rbinom(n_sp * N, 1, 0.4), nrow = n_sp)
  rownames(y) <- paste0("sp", seq_len(n_sp))
  data <- list(y = y, occ.covs = data.frame(x = rnorm(N)))
  fit <- occu(~ x, data = data, multispecies = "jsdm",
              n.factors = 2, verbose = 0)
  expect_s3_class(fit, "occu_inla_jsdm")
  expect_equal(fit$n.factors, 2)
})

test_that("engine_jsdm_sf runs without error", {
  n_sp <- 5
  N <- 40
  y <- matrix(rbinom(n_sp * N, 1, 0.4), nrow = n_sp)
  rownames(y) <- paste0("sp", seq_len(n_sp))
  coords <- cbind(runif(N), runif(N))
  data <- list(y = y, occ.covs = data.frame(x = rnorm(N)))
  fit <- occu(~ x, data = data, multispecies = "jsdm",
              n.factors = 2, spatial = coords, verbose = 0)
  expect_s3_class(fit, "occu_inla_sfjsdm")
})


# ==========================================================================
# Formula parser edge cases
# ==========================================================================

test_that("formula parser handles complex RE specifications", {
  # Multiple crossed REs
  parsed <- parse_re_formula(~ x1 + (1 | site) + (1 | observer))
  expect_length(parsed$re_list, 2)
  expect_equal(parsed$re_list[[1]]$group, "site")
  expect_equal(parsed$re_list[[2]]$group, "observer")

  # Nested grouping: (1 | a/b) → interaction group
  parsed2 <- parse_re_formula(~ x1 + (1 | region/site))
  expect_length(parsed2$re_list, 1)
  expect_equal(parsed2$re_list[[1]]$group, "region:site")

  # Suppressed intercept
  parsed3 <- parse_re_formula(~ x1 + (0 + x1 | group))
  expect_length(parsed3$re_list, 1)
  expect_equal(parsed3$re_list[[1]]$type, "slope")

  # Intercept-only (all terms are REs)
  parsed4 <- parse_re_formula(~ (1 | group))
  expect_equal(deparse(parsed4$fixed), "~1")

  # Multiple slopes
  parsed5 <- parse_re_formula(~ x1 + x2 + (1 + x1 + x2 | group))
  expect_length(parsed5$re_list, 3)  # intercept + 2 slopes
})


# ==========================================================================
# New features: areal spatial, ensemble pooling
# ==========================================================================

test_that("occu_areal constructor works", {
  # 4-region grid adjacency
  adj <- matrix(0, 4, 4)
  adj[1, 2] <- adj[2, 1] <- 1
  adj[2, 3] <- adj[3, 2] <- 1
  adj[3, 4] <- adj[4, 3] <- 1
  areal <- occu_areal(adj, model = "bym2")
  expect_s3_class(areal, "occu_areal")
  expect_equal(areal$n_regions, 4)
  expect_equal(areal$model, "bym2")
})

test_that("pool_community_ensemble produces valid output", {
  # Fake beta summaries for 5 species
  beta_list <- lapply(1:5, function(i) {
    data.frame(mean = rnorm(2), sd = abs(rnorm(2, 0, 0.1)),
               row.names = c("(Intercept)", "x1"))
  })
  result <- pool_community_ensemble(beta_list, B = 20, prop = 0.7, seed = 42)
  expect_equal(nrow(result), 2)
  expect_true("n_subsets" %in% names(result))
  expect_equal(result$n_subsets[1], 20)
})

test_that("resolve_group_vals handles interaction groups", {
  df <- data.frame(site = 1:6, year = rep(1:2, 3))
  vals <- resolve_group_vals("site:year", df, NULL)
  expect_length(vals, 6)
  # Each site:year combo should get a unique integer
  expect_equal(length(unique(vals)), length(unique(paste(df$site, df$year))))
})
