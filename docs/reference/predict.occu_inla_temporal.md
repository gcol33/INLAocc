# Predict from a temporal occupancy model

Returns per-period predictions. If AR(1) smoothing was applied, includes
the smoothed occupancy estimates.

## Usage

``` r
# S3 method for class 'occu_inla_temporal'
predict(object, ...)
```

## Arguments

- object:

  fitted occu_inla_temporal object

- ...:

  additional arguments passed to `predict.occu_inla`

## Value

list with per-period predictions and optional smoothed psi
