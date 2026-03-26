# Quick Start

## Overview

INLAocc fits occupancy models — the kind that estimate species
occurrence while accounting for imperfect detection — using INLA instead
of MCMC. If you’ve used spOccupancy, the data format and model types are
the same. The difference is speed: seconds instead of minutes.

## Fitting a basic model

Simulate some data and fit a single-species occupancy model:

``` r

library(INLAocc)

sim <- simulate_occu(N = 200, J = 4, seed = 42)
fit <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)
summary(fit)
```

The formula interface is `occu(occ.formula, det.formula, data)` —
occupancy covariates on the left, detection covariates on the right.

## Random effects

INLAocc parses mixed-model random effects via its own AST walker:

``` r

# Add a regional random intercept
sim$data$occ.covs$region <- sample(1:5, 200, replace = TRUE)

fit_re <- occu(
  ~ occ_x1 + (1 | region),
  ~ det_x1,
  data = sim$data
)
ranef(fit_re)
```

Supported syntax: `(1 | group)`, `(x | group)`, `(x || group)`
(uncorrelated), `(1 | a/b)` (nested).

## Spatial models

Pass coordinates to enable an SPDE spatial random effect:

``` r

sim_sp <- simulate_occu(N = 200, J = 4, spatial_range = 0.2, seed = 123)

fit_sp <- occu(~ occ_x1, ~ det_x1, data = sim_sp$data,
               spatial = sim_sp$data$coords)
summary(fit_sp)
```

## Model diagnostics

### One-call diagnostic panel

``` r

checkModel(fit_sp)
```

This produces a 2x2 panel: QQ plot of PIT residuals, residuals vs
fitted, dispersion test, and a Moran’s I correlogram (because
coordinates are available).

### Individual tests

Each diagnostic is also available as a standalone function:

``` r

testUniformity(fit)     # KS test: are PIT residuals uniform?
testDispersion(fit)     # over/underdispersion
testOutliers(fit)       # simulation envelope
testZeroInflation(fit)  # excess unoccupied sites
moranI(fit_sp)          # spatial autocorrelation in residuals
```

All return `htest` objects, so `$p.value` and
[`print()`](https://rdrr.io/r/base/print.html) work as expected.

### DHARMa bridge

If you have DHARMa installed, one call gives you the full DHARMa object:

``` r

res <- dharma(fit)
DHARMa::testSpatialAutocorrelation(res,
  x = sim_sp$data$coords[, 1],
  y = sim_sp$data$coords[, 2])
```

## Model comparison and averaging

``` r

m1 <- occu(~ 1, ~ 1, data = sim$data)
m2 <- occu(~ occ_x1, ~ det_x1, data = sim$data)
m3 <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)

# Compare with AIC weights
compare_models(null = m1, elev = m2, full = m3, criterion = "aic")

# Model-averaged predictions (Burnham & Anderson)
avg <- modelAverage(null = m1, elev = m2, full = m3, criterion = "aic")
avg$psi_hat   # weighted-average occupancy
avg$psi_se    # unconditional standard errors
```

[`compare_models()`](https://gillescolling.com/INLAocc/reference/compare_models.md)
supports `"aic"`, `"bic"`, and `"waic"` criteria.

## Prediction

``` r

# ggpredict-style: vary one covariate, hold others at mean
predict(fit, terms = "occ_x1 [-2:2 by=0.5]")

# Design-matrix prediction at new sites
predict(fit, X.0 = data.frame(occ_x1 = 0, occ_x2 = 1))

# Marginal effect plot data
marginal_effect(fit, "occ_x1")
```

## Multi-species models

``` r

sim_ms <- simMsOcc(N = 100, J = 3, n_species = 10, seed = 200)

fit_ms <- occu(~ occ_x1, ~ det_x1, data = sim_ms$data, multispecies = TRUE)
summary(fit_ms)

# Site-level species richness
richness(fit_ms)
```

## Multi-season models

``` r

sim_t <- simTOcc(N = 50, J = 3, n_seasons = 5, seed = 300)

fit_t <- occu(~ 1, ~ 1, data = sim_t$data, temporal = "ar1")
summary(fit_t)
```

## Standard S3 methods

``` r

coef(fit)               # posterior means
confint(fit)            # credible intervals
tidy(fit)               # broom-compatible data.frame
glance(fit)             # model-level summary row
AIC(fit); BIC(fit)      # information criteria
update(fit, ~ . - occ_x2)  # refit without occ_x2
```
