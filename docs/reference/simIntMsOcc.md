# Simulate integrated multi-species occupancy data (cf. simIntMsOcc)

Simulate integrated multi-species occupancy data (cf. simIntMsOcc)

## Usage

``` r
simIntMsOcc(
  N_total = 150,
  n_species = 10,
  n_data = 2,
  J = c(4, 3),
  n_shared = 20,
  beta_comm_mean = NULL,
  beta_comm_sd = NULL,
  beta_det = NULL,
  seed = NULL
)
```

## Arguments

- N_total:

  total unique sites

- n_species:

  number of species

- n_data:

  number of data sources

- J:

  vector of visits per source

- n_shared:

  number of sites shared across sources

- beta_comm_mean:

  community mean for occupancy coefficients

- beta_comm_sd:

  community SD for occupancy coefficients

- beta_det:

  list of detection coefficient vectors (one per source)

- seed:

  random seed

## Value

list with data for integrated multi-species models and true values
