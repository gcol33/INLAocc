# Predict occupancy at new spatial locations

Convenience wrapper for spatial models. Constructs X.0 from a covariate
data.frame and handles coordinate projection.

## Usage

``` r
predict_spatial(
  object,
  newcoords,
  newocc.covs = NULL,
  ignore.RE = FALSE,
  n.samples = 500
)
```

## Arguments

- object:

  fitted occu_inla_spatial object

- newcoords:

  n_pred x 2 matrix of prediction coordinates

- newocc.covs:

  data.frame of occupancy covariates at prediction locations

- ignore.RE:

  ignore random effects?

- n.samples:

  number of posterior samples

## Value

list with predicted psi at new locations
