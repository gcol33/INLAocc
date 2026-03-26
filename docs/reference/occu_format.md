# Format data for INLA occupancy models

Accepts data in spOccupancy-compatible format and validates it.

## Usage

``` r
occu_format(y, occ.covs = NULL, det.covs = NULL, coords = NULL, species = NULL)
```

## Arguments

- y:

  N x J detection history matrix (0/1/NA). Rows = sites, cols = visits.

- occ.covs:

  data.frame with N rows of site-level occupancy covariates. Can also be
  a named list of vectors, each length N.

- det.covs:

  named list of detection covariates. Each element is either:

  - a vector of length N (constant across visits)

  - an N x J matrix (varies by visit)

- coords:

  optional N x 2 matrix of site coordinates (for spatial models)

- species:

  optional character or integer species identifier (for multi-species)

## Value

An object of class `"occu_data"` with validated components

## Examples

``` r
y <- matrix(rbinom(200, 1, 0.4), nrow = 50, ncol = 4)
occ_covs <- data.frame(elev = rnorm(50), forest = runif(50))
det_covs <- list(
  effort = matrix(runif(200, 1, 8), 50, 4),
  date   = matrix(rnorm(200), 50, 4)
)
dat <- occu_format(y, occ_covs, det_covs)
```
