# Create multi-species occupancy data

Create multi-species occupancy data

## Usage

``` r
occu_format_ms(y_list, occ.covs = NULL, det.covs = NULL, coords = NULL)
```

## Arguments

- y_list:

  named list of N x J detection matrices (one per species)

- occ.covs:

  site-level covariates (shared across species)

- det.covs:

  detection covariates (shared or species-specific)

- coords:

  optional coordinates

## Value

object of class `"occu_data_ms"`
