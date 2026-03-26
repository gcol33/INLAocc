# Compare occupancy models via information criteria

Compare occupancy models via information criteria

## Usage

``` r
compare_models(..., criterion = c("waic", "aic", "bic"))
```

## Arguments

- ...:

  named occu_inla objects to compare

- criterion:

  `"waic"` (default), `"aic"`, or `"bic"` — used for ranking and
  computing model weights

## Value

data.frame with model comparison metrics, delta IC, and weights
