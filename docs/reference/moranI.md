# Moran's I test for spatial autocorrelation

Computes Moran's I on occupancy residuals using an inverse-distance
weight matrix (k nearest neighbours). No external dependencies — uses
the normal approximation under the randomisation assumption.

## Usage

``` r
moranI(
  object,
  coords = NULL,
  weights = c("inverse", "knn"),
  k = 10L,
  resid.type = "deviance",
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- object:

  fitted `occu_inla` object, or a numeric vector of residuals (in which
  case `coords` must be supplied)

- coords:

  optional N x 2 coordinate matrix. If `object` is an `occu_inla` fit,
  defaults to `object$data$coords`.

- weights:

  weight scheme: `"inverse"` (default) uses inverse-distance weights on
  all pairs (matches DHARMa), `"knn"` uses k nearest-neighbour binary
  weights

- k:

  number of nearest neighbours (only used when `weights = "knn"`,
  default 10)

- resid.type:

  residual type passed to
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) if `object` is
  a fit (default `"deviance"`)

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`

## Value

A list of class `"htest"` with components `statistic`, `p.value`,
`parameter` (expected I), and `method`.
