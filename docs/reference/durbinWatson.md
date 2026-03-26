# Durbin-Watson test for temporal autocorrelation in residuals

Tests whether occupancy residuals across time periods exhibit
first-order autocorrelation. Designed for temporal occupancy models
where each period yields a site-averaged residual.

## Usage

``` r
durbinWatson(
  object,
  resid.type = "deviance",
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- object:

  fitted `occu_inla_temporal` object, or a numeric vector of
  temporally-ordered residuals

- resid.type:

  residual type (default `"deviance"`)

- alternative:

  `"two.sided"` (default), `"greater"` (positive autocorrelation), or
  `"less"` (negative autocorrelation)

## Value

A list of class `"htest"` with the DW statistic, approximate p-value
(normal approximation), and the lag-1 autocorrelation.
