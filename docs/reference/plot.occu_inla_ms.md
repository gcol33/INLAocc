# Plot diagnostics for multi-species occupancy model

Plot diagnostics for multi-species occupancy model

## Usage

``` r
# S3 method for class 'occu_inla_ms'
plot(x, which = 1:2, ...)
```

## Arguments

- x:

  fitted occu_inla_ms object

- which:

  integer vector: which plots to show. 1 = per-species occupancy, 2 =
  per-species detection, 3 = community effects.

- ...:

  additional arguments passed to plot

## Value

The object `x`, returned invisibly.
