# Compute marginal effects for a covariate

Compute marginal effects for a covariate

## Usage

``` r
marginal_effect(
  object,
  covariate,
  process = c("occupancy", "detection"),
  values = NULL,
  n_points = 100,
  other.means = NULL
)
```

## Arguments

- object:

  fitted occu_inla object

- covariate:

  character: name of covariate

- process:

  "occupancy" or "detection"

- values:

  optional: vector of covariate values

- n_points:

  number of prediction points (if values is NULL)

- other.means:

  named list of values for other covariates (default: means)

## Value

data.frame with covariate values and predicted probability + CI
