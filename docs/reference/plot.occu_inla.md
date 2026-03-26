# Plot diagnostics for occu_inla

Plot diagnostics for occu_inla

## Usage

``` r
# S3 method for class 'occu_inla'
plot(x, which = 1:4, ...)
```

## Arguments

- x:

  fitted occu_inla object

- which:

  integer vector: which plots to show. 1 = EM convergence, 2 = psi
  histogram, 3 = p histogram, 4 = psi vs covariates.

- ...:

  additional arguments passed to plot

## Value

The `occu_inla` object `x`, returned invisibly.
