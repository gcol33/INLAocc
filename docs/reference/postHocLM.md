# Fit a post-hoc linear model relating covariates to estimated parameters

Useful for exploring drivers of occupancy/detection variation. Supports
both Bayesian (via INLA) and frequentist
([`lm`](https://rdrr.io/r/stats/lm.html)) fitting.

## Usage

``` r
postHocLM(
  formula,
  data,
  weights = NULL,
  method = c("bayes", "freq"),
  n.samples = 1000L
)
```

## Arguments

- formula:

  model formula (e.g., psi_hat ~ trait1 + trait2)

- data:

  data.frame with response and predictors

- weights:

  optional weights (e.g., inverse of psi standard errors)

- method:

  `"bayes"` (default, uses INLA) or `"freq"` (uses lm). Falls back to
  `"freq"` if INLA is not available.

- n.samples:

  number of posterior samples (default 1000, Bayesian only)

## Value

A list of class `"postHocLM"` with:

- beta.samples:

  matrix of posterior coefficient samples (Bayesian only)

- tau.sq.samples:

  posterior residual variance samples (Bayesian only)

- y.hat.samples:

  posterior fitted value samples (Bayesian only)

- bayes.R2:

  posterior samples of Bayesian R-squared (Bayesian only)

- summary:

  data.frame of coefficient summaries

- lm.fit:

  frequentist lm object (always available)

- method:

  character: which method was used
