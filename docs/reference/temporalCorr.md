# Extract AR(1) temporal correlation from a fitted temporal model

Returns the posterior summary for the AR(1) autocorrelation parameter
(rho) and the temporal precision/variance.

## Usage

``` r
temporalCorr(object)
```

## Arguments

- object:

  fitted temporal occupancy model (class `occu_inla_temporal`)

## Value

A data.frame with columns `mean`, `sd`, `q025`, `q975`.
