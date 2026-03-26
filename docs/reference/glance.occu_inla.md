# Glance at model-level statistics

Returns a single-row data.frame of model-level summaries, similar to
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html).

## Usage

``` r
# S3 method for class 'occu_inla'
glance(x, ...)
```

## Arguments

- x:

  fitted occu_inla object

- ...:

  ignored

## Value

A one-row data.frame with columns `nobs`, `n_sites`, `n_visits`,
`n_occ_coef`, `n_det_coef`, `logLik`, `WAIC`, `n_iter`, `converged`.
