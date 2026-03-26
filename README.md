<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/INLAocc/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/INLAocc/graph/badge.svg)](https://app.codecov.io/gh/gcol33/INLAocc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# INLAocc

**Fast occupancy models using INLA — a Laplace-based alternative to MCMC.**

## Quick Start

```r
library(INLAocc)

sim <- simulate_occu(N = 200, J = 4, n_occ_covs = 2, n_det_covs = 1, seed = 42)
fit <- occu(~ occ_x1 + occ_x2, ~ det_x1, data = sim$data)
summary(fit)
```

## Statement of Need

Occupancy models estimate where species occur while accounting for imperfect detection. The standard approach (MCMC via spOccupancy) is accurate but slow, especially for spatial, multi-species, or multi-season extensions. INLAocc replaces the MCMC sampler with an EM algorithm that uses INLA at each M-step, giving comparable estimates in a fraction of the time.

## Features

### Model Types

Every model type in spOccupancy has a corresponding INLAocc call:

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

Native mixed-model formula syntax — no lme4 dependency:

```r
# Random intercept
occu(~ elev + (1 | region), ~ effort, data)

# Random slope
occu(~ elev + (1 + elev | region), ~ effort, data)
```

### Diagnostics

- `ppcOccu()` — posterior predictive checks (Freeman-Tukey, chi-squared)
- `waicOccu()` — WAIC for model comparison
- `fitted()`, `residuals()` — standard S3 methods
- `compare_models()` — side-by-side WAIC/deviance table
- `k.fold` argument — k-fold cross-validation built into `occu()`

### Prediction

- Design-matrix prediction via `X.0` / `coords.0`
- `marginal_effect()` — covariate response curves
- `richness()` — site-level species richness from multi-species models

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
sim <- simulate_occu(N = 200, J = 4, n_occ_covs = 2, n_det_covs = 1, seed = 42)
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
sim <- simulate_occu(N = 200, J = 4, n_occ_covs = 1, n_det_covs = 1,
                     spatial_range = 0.2, spatial_var = 1.0, seed = 123)

fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, spatial = sim$data$coords)
```

### Multi-species community model

```r
sim_ms <- simMsOcc(N = 100, J = 3, n_species = 10, n_occ_covs = 1, n_det_covs = 1, seed = 200)

fit <- occu(~ occ_x1, ~ det_x1, data = sim_ms$data, multispecies = TRUE)
rich <- richness(fit)
```

## Documentation

- [Getting Started](https://gillescolling.com/INLAocc/articles/quickstart.html)
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
