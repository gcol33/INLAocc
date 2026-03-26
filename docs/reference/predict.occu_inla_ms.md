# Predict from a multi-species occupancy model

Returns per-species predictions by delegating to the base predict
method.

## Usage

``` r
# S3 method for class 'occu_inla_ms'
predict(object, ...)
```

## Arguments

- object:

  fitted occu_inla_ms object

- ...:

  additional arguments passed to `predict.occu_inla`

## Value

named list of per-species prediction results
