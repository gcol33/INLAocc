# Extract model coefficients

Extract model coefficients

## Usage

``` r
# S3 method for class 'occu_inla'
coef(object, process = c("both", "occupancy", "detection"), ...)
```

## Arguments

- object:

  fitted occu_inla object

- process:

  `"occupancy"` (default), `"detection"`, or `"both"`.

- ...:

  ignored

## Value

Named numeric vector of posterior means. If `process = "both"`, a named
list with `occ` and `det`.
