# Extract random effects

Returns the posterior summaries of random effect levels, similar to
[`lme4::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

## Usage

``` r
# S3 method for class 'occu_inla'
ranef(object, process = c("both", "occupancy", "detection"), ...)
```

## Arguments

- object:

  fitted occu_inla object

- process:

  `"occupancy"` (default), `"detection"`, or `"both"`.

- ...:

  ignored

## Value

A named list of data.frames (one per random effect group), each with
columns `mean`, `sd`, `0.025quant`, `0.975quant`. If `process = "both"`,
a list with `occ` and `det` sub-lists.
