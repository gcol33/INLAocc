# Create SPDE mesh and spatial effect for occupancy models

Wraps INLA's mesh construction and SPDE model setup for use in occupancy
models. The spatial effect is added to the occupancy linear predictor.

## Usage

``` r
occu_spatial(
  coords,
  max.edge = NULL,
  cutoff = NULL,
  offset = NULL,
  prior.range = c(0.3, 0.5),
  prior.sigma = c(1, 0.05),
  alpha = 2
)
```

## Arguments

- coords:

  N x 2 matrix of coordinates

- max.edge:

  numeric vector of length 2: max triangle edge length (inner domain,
  outer extension). Default auto-scaled from data extent.

- cutoff:

  minimum distance between mesh nodes. Default: `max.edge[1]/5`.

- offset:

  numeric vector of length 2: inner and outer extension distances.

- prior.range:

  numeric(2): PC prior for range. c(r0, p) means P(range \< r0) = p.

- prior.sigma:

  numeric(2): PC prior for marginal SD. c(s0, p) means P(sigma \> s0) =
  p.

- alpha:

  SPDE smoothness parameter (default 2, Matern nu = 1).

## Value

object of class `"occu_spatial"` with mesh, spde, and A matrix
components
