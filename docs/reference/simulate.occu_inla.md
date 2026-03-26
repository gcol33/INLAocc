# Simulate replicate datasets from a fitted occupancy model

Generates posterior predictive datasets by sampling the latent occupancy
state z from Bernoulli(psi) and observations y from Bernoulli(z \* p) at
each visit. Returns a site-level summary (number of detections per site)
suitable for DHARMa's `createDHARMa()`.

## Usage

``` r
# S3 method for class 'occu_inla'
simulate(object, nsim = 250L, seed = 123L, level = c("site", "obs"), ...)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of replicate datasets (default 250)

- seed:

  random seed for reproducibility (default 123)

- level:

  `"site"` (default) returns an N x nsim matrix of per-site detection
  counts; `"obs"` returns an (N\*J) x nsim matrix of individual visit
  outcomes.

- ...:

  ignored

## Value

A matrix (see `level`).
