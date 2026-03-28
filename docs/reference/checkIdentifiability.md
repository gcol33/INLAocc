# Diagnose identifiability issues in occupancy models

Checks for common sources of non-identifiability in occupancy models,
both before fitting (on data) and after fitting (on the model object).
Each issue is flagged with a severity and a plain-language explanation.

## Usage

``` r
checkIdentifiability(object, occ.re = NULL, det.re = NULL)
```

## Arguments

- object:

  an `occu_data` object (pre-fit) or a fitted `occu_inla` object
  (post-fit)

- occ.re:

  optional list of occupancy random effects (for pre-fit RE checks;
  auto-extracted from fitted models)

- det.re:

  optional list of detection random effects

## Value

A list of class `"occu_identifiability"` with component `issues`, a
data.frame of flagged issues.
