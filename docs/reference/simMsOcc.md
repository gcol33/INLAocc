# Simulate multi-species occupancy data (cf. simMsOcc)

Simulate multi-species occupancy data (cf. simMsOcc)

## Usage

``` r
simMsOcc(
  N = 100,
  J = 4,
  n_species = 10,
  n_occ_covs = 1,
  n_det_covs = 1,
  beta_comm_mean = NULL,
  beta_comm_sd = NULL,
  alpha_comm_mean = NULL,
  alpha_comm_sd = NULL,
  spatial_range = NULL,
  spatial_var = 1,
  seed = NULL
)
```

## Arguments

- N:

  number of sites

- J:

  number of visits

- n_species:

  number of species

- n_occ_covs:

  number of occupancy covariates

- n_det_covs:

  number of detection covariates

- beta_comm_mean:

  community mean for occupancy coefficients

- beta_comm_sd:

  community SD for occupancy coefficients

- alpha_comm_mean:

  community mean for detection coefficients

- alpha_comm_sd:

  community SD for detection coefficients

- spatial_range:

  spatial range (NULL = no spatial effect)

- spatial_var:

  spatial variance

- seed:

  random seed

## Value

list with occu_data_ms object and true parameters
