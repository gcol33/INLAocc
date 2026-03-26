# Empirical semivariogram of occupancy residuals

Computes the empirical semivariogram in distance bins. Useful for
visually assessing whether spatial structure remains in the residuals
after model fitting.

## Usage

``` r
variogram(
  object,
  coords = NULL,
  n.bins = 15L,
  max.dist = NULL,
  resid.type = "deviance"
)
```

## Arguments

- object:

  fitted `occu_inla` object, or a numeric vector of residuals (in which
  case `coords` must be supplied)

- coords:

  optional N x 2 coordinate matrix (defaults to `object$data$coords`)

- n.bins:

  number of distance bins (default 15)

- max.dist:

  maximum distance to consider (default: half the max pairwise distance)

- resid.type:

  residual type (default `"deviance"`)

## Value

A data.frame of class `"occu_variogram"` with columns `dist` (bin
midpoint), `gamma` (semivariance), and `n.pairs` (number of pairs in
bin). Has a [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method.
