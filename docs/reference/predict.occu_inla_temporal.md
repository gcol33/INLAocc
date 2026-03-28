# Predict from a temporal occupancy model

Returns per-period predictions. If AR(1) smoothing was applied, includes
the smoothed occupancy estimates. If `period` is specified, returns
predictions for that period only.

## Usage

``` r
# S3 method for class 'occu_inla_temporal'
predict(object, period = NULL, ...)
```

## Arguments

- object:

  fitted occu_inla_temporal object

- period:

  optional integer: return predictions for this period only (1-indexed).
  If `NULL` (default), returns all periods.

- ...:

  additional arguments passed to `predict.occu_inla`

## Value

list with per-period predictions and optional smoothed psi. If `period`
is specified, returns predictions for that single period.
