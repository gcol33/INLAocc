# Compute residuals for an integrated occupancy model

Compute residuals for an integrated occupancy model

## Usage

``` r
# S3 method for class 'occu_inla_int'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  fitted occu_inla_int object

- type:

  "deviance" (default), "pearson", or "response"

- ...:

  ignored

## Value

list with per-source detection residuals and shared occupancy residuals
