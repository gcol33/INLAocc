# Simulate integrated (multi-source) occupancy data (cf. simIntOcc)

Simulate integrated (multi-source) occupancy data (cf. simIntOcc)

## Usage

``` r
simIntOcc(
  N_total = 150,
  n_data = 2,
  J = c(4, 3),
  n_shared = 20,
  beta_occ = c(0.5, 0.3),
  beta_det = list(c(0.2, -0.4), c(-0.1, 0.3)),
  seed = NULL
)
```

## Arguments

- N_total:

  total unique sites

- n_data:

  number of data sources

- J:

  vector of visits per source

- n_shared:

  number of sites shared across sources

- beta_occ:

  occupancy coefficients

- beta_det:

  list of detection coefficients per source

- seed:

  random seed

## Value

list with data for intOccu_inla and true values
