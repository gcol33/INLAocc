# Approximate variance-covariance matrix

Returns the diagonal approximation (independent posteriors) from INLA.
For the full posterior covariance, use
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).

## Usage

``` r
# S3 method for class 'occu_inla'
vcov(object, process = c("occupancy", "detection"), ...)
```

## Arguments

- object:

  fitted occu_inla object

- process:

  `"occupancy"` (default) or `"detection"`.

- ...:

  ignored

## Value

A diagonal variance-covariance matrix.
