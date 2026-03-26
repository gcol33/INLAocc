# Simulate temporal multi-species occupancy data (cf. simTMsOcc)

Simulate temporal multi-species occupancy data (cf. simTMsOcc)

## Usage

``` r
simTMsOcc(
  N = 100,
  J = 4,
  n_species = 10,
  n_seasons = 5,
  n_occ_covs = 1,
  n_det_covs = 1,
  beta_comm_mean = NULL,
  beta_comm_sd = NULL,
  alpha_comm_mean = NULL,
  alpha_comm_sd = NULL,
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

- n_species:

  number of species

- n_seasons:

  number of primary periods

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

- ar1:

  logical: use AR(1) temporal correlation

- rho:

  AR(1) correlation parameter

- sigma_t:

  temporal innovation SD

- seed:

  random seed

## Value

list with 4D data and true parameters
