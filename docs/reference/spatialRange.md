# Extract spatial range and standard deviation from a fitted spatial model

Returns the posterior mean, SD, and 95\\ SPDE range and marginal
standard deviation. The range is in the same units as the coordinates
(e.g. metres for UTM).

## Usage

``` r
spatialRange(object)
```

## Arguments

- object:

  fitted spatial occupancy model (class `occu_inla_spatial` or any model
  with a spatial SPDE component)

## Value

A data.frame with columns `mean`, `sd`, `q025`, `q975`, and rows `range`
and `stdev`.
