# Tidy model output into a data.frame

Returns a tidy data.frame of fixed effect estimates, similar to
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html).

## Usage

``` r
# S3 method for class 'occu_inla'
tidy(x, process = c("both", "occupancy", "detection"), conf.level = 0.95, ...)
```

## Arguments

- x:

  fitted occu_inla object

- process:

  `"occupancy"` (default), `"detection"`, or `"both"`.

- conf.level:

  credible level for intervals (default 0.95).

- ...:

  ignored

## Value

A data.frame with columns `process`, `term`, `estimate`, `std.error`,
`conf.low`, `conf.high`.
