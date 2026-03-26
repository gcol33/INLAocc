# Predict from a fitted occupancy model

Supports three prediction modes:

1.  **terms-based** (like `ggpredict`): pass a character vector of
    covariate terms to vary, with optional bracket notation for ranges
    or specific values. Everything else is held at its mean (continuous)
    or mode (factor).

2.  **design-matrix**: pass `X.0` / `X.p.0` directly
    (spOccupancy-compatible).

3.  **in-sample**: call with no extra arguments.

## Usage

``` r
# S3 method for class 'occu_inla'
predict(
  object,
  X.0 = NULL,
  X.p.0 = NULL,
  coords.0 = NULL,
  ignore.RE = FALSE,
  type = c("both", "occupancy", "detection"),
  quantiles = c(0.025, 0.5, 0.975),
  n.samples = 500,
  terms = NULL,
  process = c("occupancy", "detection"),
  n_points = 50L,
  ...
)
```

## Arguments

- object:

  fitted occu_inla object

- X.0:

  optional: design matrix for occupancy prediction at new locations.
  Must include intercept column if model has intercept.

- X.p.0:

  optional: design matrix for detection prediction.

- coords.0:

  optional: coordinates for spatial prediction (n_pred x 2).

- ignore.RE:

  logical: set random effects to zero? (default FALSE)

- type:

  "occupancy", "detection", or "both"

- quantiles:

  quantile levels for credible intervals

- n.samples:

  number of posterior samples for uncertainty

- terms:

  character vector of terms to vary, with optional bracket notation.
  Examples: `"elev"`, `"elev [0:100]"`, `"elev [0:100 by=5]"`,
  `"elev [1, 5, 10]"`, `"habitat [forest, grassland]"`. The first term
  is the x-axis; the second (if any) defines groups; the third defines
  facets.

- process:

  `"occupancy"` (default) or `"detection"`. Only used with `terms`.

- n_points:

  number of prediction points per continuous term (default 50). Only
  used with `terms`.

- ...:

  additional arguments

## Value

If `terms` is provided, an `occu_prediction` data.frame with columns
`x`, `estimate`, `lower`, `upper`, and optional `group`/`facet`. Has a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
Otherwise, a list with `psi.0`, `p.0`, `z.0`.

## Examples

``` r
# \donttest{
# fit <- occu(~ elev + forest, ~ effort, data)
# predict(fit, terms = "elev")
# predict(fit, terms = "elev [0:100 by=10]")
# predict(fit, terms = c("elev", "forest"))
# predict(fit, terms = "effort", process = "detection")
# }
```
