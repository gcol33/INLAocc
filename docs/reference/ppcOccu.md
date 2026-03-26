# Posterior predictive checks for occupancy models

Computes a goodness-of-fit statistic for observed and replicated data.
Analogous to spOccupancy::ppcOcc().

## Usage

``` r
ppcOccu(
  object,
  fit.stat = c("freeman-tukey", "chi-squared"),
  group = 1,
  n.samples = 500
)
```

## Arguments

- object:

  fitted occu_inla object

- fit.stat:

  "freeman-tukey" (default) or "chi-squared"

- group:

  1 = aggregate by site, 2 = aggregate by visit

- n.samples:

  number of replicated datasets to generate (default 500)

## Value

list with:

- fit.y:

  vector of fit statistic for observed data across samples

- fit.y.rep:

  vector of fit statistic for replicated data

- bayesian.p:

  Bayesian p-value (proportion fit.y.rep \> fit.y)
