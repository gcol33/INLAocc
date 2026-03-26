# Specify random effects for occupancy or detection process

Creates a random effects specification that gets translated to INLA f()
terms.

## Usage

``` r
occu_re(
  type = c("intercept", "slope", "iid"),
  group,
  covariate = NULL,
  model = "iid",
  prior = NULL,
  constr = FALSE,
  correlated = TRUE,
  n_groups = NULL
)
```

## Arguments

- type:

  "intercept", "slope", or "iid"

- group:

  character: grouping variable name (column in covariates)

- covariate:

  character: for random slopes, the covariate name

- model:

  INLA model for the random effect ("iid", "ar1", "rw1", "rw2", "besag")

- prior:

  list with prior specification (passed to INLA hyper)

- constr:

  logical: sum-to-zero constraint

- correlated:

  logical: whether this RE is correlated with other REs sharing the same
  grouping factor. Set to `FALSE` when parsed from lme4 `||`
  (double-bar) syntax, indicating uncorrelated random effects.

- n_groups:

  optional: number of groups (auto-detected if NULL)

## Value

object of class `"occu_re"`
