# Create an areal spatial effect for occupancy models (CAR/BYM2)

For data on grids, administrative units, or any areal structure. Uses
INLA's `besag` or `bym2` model with an adjacency graph.

## Usage

``` r
occu_areal(
  adj,
  model = c("bym2", "besag"),
  prior.sigma = c(1, 0.05),
  scale.model = TRUE
)
```

## Arguments

- adj:

  Adjacency structure: a symmetric matrix, a `nb` object (from `spdep`),
  or a path to an INLA graph file.

- model:

  `"bym2"` (default, recommended) or `"besag"`. BYM2 separates
  structured (spatial) and unstructured (iid) components with a mixing
  parameter.

- prior.sigma:

  numeric(2): PC prior on marginal SD. c(s0, p) means P(sigma \> s0)
  = p. Default c(1, 0.05).

- scale.model:

  Logical: scale the precision matrix so the generalized variance is 1
  (recommended for bym2). Default TRUE.

## Value

object of class `"occu_areal"` with graph and model components
