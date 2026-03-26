# Simulate multi-season occupancy data (cf. simTOcc)

Simulate multi-season occupancy data (cf. simTOcc)

## Usage

``` r
simTOcc(
  N = 100,
  J = 4,
  n_seasons = 5,
  beta_occ = c(0.5, 0.3),
  beta_det = c(0, -0.3),
  ar1 = TRUE,
  rho = 0.7,
  sigma_t = 0.5,
  seed = NULL
)
```

## Arguments

- N:

  number of sites

- J:

  number of visits per season

- n_seasons:

  number of primary periods

- beta_occ:

  occupancy coefficients

- beta_det:

  detection coefficients

- ar1:

  logical: use AR(1) temporal correlation on occupancy

- rho:

  AR(1) correlation parameter

- sigma_t:

  temporal innovation SD

- seed:

  random seed

## Value

list with data suitable for temporal_occu_inla and true values
