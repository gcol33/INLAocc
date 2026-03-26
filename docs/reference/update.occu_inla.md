# Update and refit an occupancy model

Modify the formula or arguments of a fitted model and refit.

## Usage

``` r
# S3 method for class 'occu_inla'
update(
  object,
  occ.formula = NULL,
  det.formula = NULL,
  data = NULL,
  ...,
  evaluate = TRUE
)
```

## Arguments

- object:

  fitted occu_inla object

- occ.formula:

  new occupancy formula (default: keep existing). Supports
  `update.formula` syntax: `. ~ . - term`.

- det.formula:

  new detection formula (default: keep existing).

- data:

  new data (default: keep existing).

- ...:

  additional arguments passed to
  [`occu`](https://gillescolling.com/INLAocc/reference/occu.md).

- evaluate:

  if `FALSE`, return the call instead of fitting.

## Value

A new fitted model.
