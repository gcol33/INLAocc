# Create a community-level (multi-species) random effect

For multi-species models: species-specific intercepts/slopes drawn from
a community distribution.

## Usage

``` r
occu_community_re(
  type = c("intercept", "slope"),
  covariate = NULL,
  model = "iid",
  prior = NULL
)
```

## Arguments

- type:

  "intercept" or "slope"

- covariate:

  for slopes, the covariate name

- model:

  "iid" for exchangeable species effects, "iid" with group for
  correlated

- prior:

  prior specification

## Value

object of class `"occu_community_re"`
