# Predict from an integrated occupancy model

Returns in-sample predictions. For out-of-sample, use the shared
occupancy component via the base `predict.occu_inla` method on the
`occ_fit` sub-object.

## Usage

``` r
# S3 method for class 'occu_inla_int'
predict(object, ...)
```

## Arguments

- object:

  fitted occu_inla_int object

- ...:

  additional arguments (ignored)

## Value

list with psi.0 (shared occupancy) and per-source p.0
