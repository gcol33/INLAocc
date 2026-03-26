# Simulate occupancy data for testing

Simulate occupancy data for testing

## Usage

``` r
simulate_occu(
  N = 100,
  J = 4,
  n_occ_covs = 2,
  n_det_covs = 1,
  beta_occ = NULL,
  beta_det = NULL,
  random_occ_sd = 0,
  random_det_sd = 0,
  spatial_range = NULL,
  spatial_var = 1,
  seed = NULL
)
```

## Arguments

- N:

  number of sites

- J:

  number of visits per site

- n_occ_covs:

  number of occupancy covariates

- n_det_covs:

  number of detection covariates

- beta_occ:

  occupancy coefficients (intercept + slopes)

- beta_det:

  detection coefficients (intercept + slopes)

- random_occ_sd:

  SD of site-level random intercept on occupancy (0 = none)

- random_det_sd:

  SD of site-level random intercept on detection (0 = none)

- spatial_range:

  spatial range parameter (NULL = no spatial effect)

- spatial_var:

  spatial variance

- seed:

  random seed

## Value

list with occu_data object and true parameter values
