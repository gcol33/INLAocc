# Compute residuals for a temporal occupancy model

Compute residuals for a temporal occupancy model

## Usage

``` r
# S3 method for class 'occu_inla_temporal'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  fitted occu_inla_temporal object

- type:

  "deviance" (default), "pearson", or "response"

- ...:

  ignored

## Value

list of per-period residual lists
