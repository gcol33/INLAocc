# Specify priors for INLA occupancy models (spOccupancy-compatible)

Specify priors for INLA occupancy models (spOccupancy-compatible)

## Usage

``` r
occu_priors(
  beta.normal = list(mean = 0, var = 2.72),
  alpha.normal = list(mean = 0, var = 2.72),
  sigma.sq.psi.ig = c(0.1, 0.1),
  sigma.sq.p.ig = c(0.1, 0.1),
  sigma.sq.ig = NULL,
  phi.unif = NULL,
  beta.comm.normal = NULL,
  alpha.comm.normal = NULL,
  tau.sq.beta.ig = NULL,
  tau.sq.alpha.ig = NULL
)
```

## Arguments

- beta.normal:

  list(mean, var) for occupancy fixed effects Normal prior

- alpha.normal:

  list(mean, var) for detection fixed effects Normal prior

- sigma.sq.psi.ig:

  c(shape, scale) for occupancy RE variance Inv-Gamma

- sigma.sq.p.ig:

  c(shape, scale) for detection RE variance Inv-Gamma

- sigma.sq.ig:

  c(shape, scale) for spatial variance (spatial models)

- phi.unif:

  c(lower, upper) for spatial decay Uniform prior

- beta.comm.normal:

  list(mean, var) community occ prior (multi-species)

- alpha.comm.normal:

  list(mean, var) community det prior (multi-species)

- tau.sq.beta.ig:

  list(a, b) species-level occ variance (multi-species)

- tau.sq.alpha.ig:

  list(a, b) species-level det variance (multi-species)

## Value

list of class `"occu_priors"` with INLA-compatible prior specs
