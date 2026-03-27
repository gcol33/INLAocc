<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/INLAocc/graph/badge.svg)](https://app.codecov.io/gh/gcol33/INLAocc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# INLAocc

**Occupancy models via INLA — a fast alternative to MCMC.**

Fit single-species, multi-species, spatial, temporal, and integrated occupancy models through one function. INLAocc uses an EM algorithm with INLA at each M-step, giving fast and accurate estimates of species occurrence while accounting for imperfect detection.

## Quick Start

```r
library(INLAocc)

sim <- simulate_occu(N = 200, J = 4, seed = 42)
fit <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)
summary(fit)

# One-call diagnostic panel
checkModel(fit)
```

## Why INLAocc

Occupancy models separate two entangled processes — where a species actually occurs (occupancy) and where you happen to detect it (detection). Getting this right matters for conservation decisions, but fitting these models with MCMC can be slow, especially for spatial or multi-species extensions.

The core problem: at sites where you never detected the species, you don't know if it was absent or just missed. INLAocc handles this by iterating two steps. First, it estimates the probability that each undetected site is truly occupied, using the current model (E-step). Then, treating those probabilities as data, it refits the occupancy and detection models using INLA — a fast Bayesian method (M-step). This back-and-forth repeats until the estimates stabilise. A multiple imputation step then corrects for the bias this approach introduces. The result matches full MCMC estimates (correlation > 0.99) at a fraction of the runtime.

## Features

### Model Types

One `occu()` call handles everything — model type is determined by the arguments you pass:

| Call | Description |
|---|---|
| `occu(~ x, ~ x, data)` | Single-species |
| `occu(..., spatial = coords)` | Spatial (SPDE mesh) |
| `occu(..., temporal = "ar1")` | Multi-season with AR(1) |
| `occu(..., multispecies = TRUE)` | Community model |
| `occu(..., multispecies = TRUE, n.factors = k)` | Latent factor |
| `occu(..., multispecies = TRUE, n.factors = k, spatial = coords)` | Spatial factor |
| `occu(..., spatial = coords, svc = 1)` | Spatially-varying coefficients |
| `occu(..., integrated = TRUE)` | Integrated multi-source |
| `occu(..., multispecies = "jsdm")` | Joint SDM (no detection) |

### Random Effects

Mixed-model formula syntax with native AST parsing:

```r
occu(~ elev + (1 | region), ~ effort, data)              # random intercept
occu(~ elev + (1 + elev | region), ~ effort, data)       # random slope
occu(~ elev + (elev || region), ~ effort, data)           # uncorrelated
occu(~ elev + (1 | site/plot), ~ effort, data)            # nested
```

### Diagnostics

All diagnostics run on model residuals — they check whether the fitted model
captured the structure in the data, or whether patterns remain unexplained.

```r
# Does the model fit? (simulation-based, on residuals)
testUniformity(fit)                # are residuals uniformly distributed?
testDispersion(fit)                # more/less variance in residuals than expected?
testOutliers(fit)                  # sites outside the simulation envelope?
testZeroInflation(fit)             # more unoccupied sites than the model predicts?

# Is there leftover spatial or temporal structure in residuals?
moranI(fit)                        # Moran's I on occupancy residuals
durbinWatson(fit)                  # lag-1 autocorrelation across time periods
variogram(fit)                     # semivariance of residuals vs distance

# One-call panel: QQ plot, residuals vs fitted, dispersion, correlogram
checkModel(fit)

# Before fitting: flag identifiability problems in the data
checkIdentifiability(dat)
checkIdentifiability(fit)          # also works post-fit (boundary estimates, collapsed REs)
```

Also available: `ppcOccu()` (posterior predictive checks), `waicOccu()`,
`AIC()`, `BIC()`, `pitResiduals()`, `simulate()`.
If DHARMa is installed, `dharma(fit)` creates a full DHARMa object.

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

# Auto-scaled mesh from coordinate extent
fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = sim$data$coords)

# Control mesh resolution (units match your CRS — metres for UTM)
fit <- occu(~ occ_x1, ~ det_x1, data = sim$data,
            spatial = sim$data$coords,
            spde.args = list(max.edge = c(5000, 15000)))

# Or build the spatial object explicitly
sp <- occu_spatial(coords,
                   max.edge = c(5000, 15000),    # triangle edges in CRS units
                   cutoff = 1000,                 # min distance between nodes
                   prior.range = c(20000, 0.5),   # P(range < 20km) = 0.5
                   prior.sigma = c(1, 0.5))       # P(sigma > 1) = 0.5
fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = sp)

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

## Coming from spOccupancy?

INLAocc accepts the same `list(y, occ.covs, det.covs, coords)` data format. The main differences: one `occu()` function instead of separate model functions, no MCMC tuning parameters, and random slopes are supported in addition to random intercepts.

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
