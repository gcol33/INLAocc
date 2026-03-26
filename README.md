<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/INLAocc/graph/badge.svg)](https://app.codecov.io/gh/gcol33/INLAocc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# INLAocc

**Fast occupancy models using INLA — a Laplace-based alternative to MCMC.**

INLAocc fits single-species, multi-species, spatial, temporal, and integrated occupancy models in seconds where spOccupancy takes minutes. Same data format, same model types, same diagnostics — just faster.

## Quick Start

```r
library(INLAocc)

sim <- simulate_occu(N = 200, J = 4, seed = 42)
fit <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)
summary(fit)

# One-call diagnostic panel
checkModel(fit)
```

## Statement of Need

Occupancy models estimate where species occur while accounting for imperfect detection. The standard Bayesian approach (MCMC via spOccupancy) is accurate but slow — minutes per model, hours for multi-species or spatial extensions. INLAocc replaces the MCMC sampler with an EM algorithm that uses INLA at each M-step, giving comparable estimates in a fraction of the time. A novel EM + multiple imputation hybrid corrects for beta attenuation, achieving correlations > 0.99 with MCMC estimates.

## Features

### Model Types

Every model type in spOccupancy has a corresponding INLAocc call through a single unified `occu()` interface:

| spOccupancy | INLAocc | Description |
|---|---|---|
| `PGOcc` | `occu()` | Single-species |
| `spPGOcc` | `occu(..., spatial = coords)` | Spatial SPDE |
| `tPGOcc` | `occu(..., temporal = "ar1")` | Multi-season AR(1) |
| `msPGOcc` | `occu(..., multispecies = TRUE)` | Community |
| `lfMsPGOcc` | `occu(..., multispecies = TRUE, n.factors = k)` | Latent factor |
| `sfMsPGOcc` | `occu(..., multispecies = TRUE, n.factors = k, spatial = coords)` | Spatial factor |
| `svcPGOcc` | `occu(..., spatial = coords, svc = 1)` | Spatially-varying coefficients |
| `intPGOcc` | `occu(..., integrated = TRUE)` | Integrated multi-source |
| `lfJSDM` | `occu(..., multispecies = "jsdm")` | Joint SDM (no detection) |

### Random Effects

Mixed-model formula syntax with native AST parsing:

```r
occu(~ elev + (1 | region), ~ effort, data)              # random intercept
occu(~ elev + (1 + elev | region), ~ effort, data)       # random slope
occu(~ elev + (elev || region), ~ effort, data)           # uncorrelated
occu(~ elev + (1 | site/plot), ~ effort, data)            # nested
```

### Diagnostics

Full diagnostic suite — native implementations with zero external dependencies:

```r
# Simulation-based (DHARMa equivalents)
simulate(fit, nsim = 250)          # posterior predictive simulation
pitResiduals(fit)                  # PIT scaled residuals
testUniformity(fit)                # KS test on PIT residuals
testDispersion(fit)                # over/underdispersion
testOutliers(fit)                  # simulation envelope test
testZeroInflation(fit)             # excess zeros

# Spatial / temporal autocorrelation
moranI(fit)                        # Moran's I (inverse-distance or k-NN)
durbinWatson(fit)                  # Durbin-Watson for temporal models
variogram(fit)                     # empirical semivariogram

# GOF and model comparison
ppcOccu(fit)                       # posterior predictive checks
waicOccu(fit)                      # WAIC
AIC(fit); BIC(fit)                 # information criteria

# One-call diagnostic panel
checkModel(fit)                    # QQ, residuals, dispersion, correlogram
```

If DHARMa is installed, `dharma(fit)` creates a full DHARMa object for access to all DHARMa tests.

### Model Selection & Averaging

```r
comp <- compare_models(null = m1, elev = m2, full = m3, criterion = "aic")
#>   model      AIC delta  weight
#> 1  full  339.793 0.000   0.733
#> 2  elev  342.072 2.278   0.235
#> 3  null  346.010 6.216   0.033

avg <- modelAverage(null = m1, elev = m2, full = m3, criterion = "aic")
avg$psi_hat   # model-averaged occupancy
avg$psi_se    # unconditional SEs (Burnham & Anderson)
```

### Prediction

```r
predict(fit, terms = "elev [0:2000 by=100]")   # ggpredict-style
predict(fit, X.0 = new_covs)                    # design-matrix
predict_spatial(fit, newcoords, newocc.covs)     # spatial interpolation
marginal_effect(fit, "elev")                     # response curves
richness(fit_ms)                                 # multi-species richness
```

### S3 Methods

Full integration with R's generic function system:

```r
coef(fit)              # posterior means
confint(fit)           # credible intervals
vcov(fit)              # variance-covariance
tidy(fit)              # broom-compatible data.frame
glance(fit)            # model-level summary
ranef(fit)             # random effect summaries
update(fit, ~ . - x2)  # refit with modified formula
logLik(fit)            # observed-data log-likelihood
nobs(fit)              # number of observations
```

## Installation

INLAocc requires INLA, which is not on CRAN:

```r
install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

Then install the development version:

```r
# install.packages("pak")
pak::pak("gcol33/INLAocc")
```

## Usage

### Single-species with random effects

```r
sim <- simulate_occu(N = 200, J = 4, seed = 42)
sim$data$occ.covs$region <- sample(1:5, 200, replace = TRUE)

fit <- occu(
  ~ occ_x1 + occ_x2 + (1 | region),
  ~ det_x1,
  data = sim$data
)
summary(fit)
```

### Spatial model

```r
sim <- simulate_occu(N = 200, J = 4, spatial_range = 0.2, seed = 123)

fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = sim$data$coords)
moranI(fit)  # check residual spatial autocorrelation
```

### Multi-species community model

```r
sim_ms <- simMsOcc(N = 100, J = 3, n_species = 10, seed = 200)

fit <- occu(~ occ_x1, ~ det_x1, data = sim_ms$data, multispecies = TRUE)
rich <- richness(fit)
```

### Model averaging

```r
m1 <- occu(~ 1, ~ 1, data = sim$data)
m2 <- occu(~ occ_x1, ~ det_x1, data = sim$data)
m3 <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)

avg <- modelAverage(null = m1, elev = m2, full = m3)
```

## Documentation

- [Quick Start](https://gillescolling.com/INLAocc/articles/quickstart.html)
- [Full Reference](https://gillescolling.com/INLAocc/reference/)

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools should be free and open. I started these projects for my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (see the LICENSE.md file)

## Citation

```bibtex
@software{INLAocc,
  author = {Colling, Gilles},
  title = {INLAocc: Occupancy Models via Integrated Nested Laplace Approximation},
  year = {2026},
  url = {https://github.com/gcol33/INLAocc}
}
```
