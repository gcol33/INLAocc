# Compute PIT (scaled) residuals for an occupancy model

For each site, computes the quantile of the observed detection count
within the posterior predictive distribution from
[`simulate`](https://rdrr.io/r/stats/simulate.html). If the model is
correct, these residuals are Uniform(0, 1).

## Usage

``` r
pitResiduals(object, nsim = 250L, seed = 123L)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations (default 250)

- seed:

  random seed (default 123)

## Value

A numeric vector of length N with values between 0 and 1.

## Details

For integer-valued responses, a randomisation step avoids discrete
artefacts: the residual is drawn uniformly between P(sim \< obs) and
P(sim \<= obs).
