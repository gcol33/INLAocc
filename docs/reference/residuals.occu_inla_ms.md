# Compute residuals for a multi-species occupancy model

Compute residuals for a multi-species occupancy model

## Usage

``` r
# S3 method for class 'occu_inla_ms'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  fitted occu_inla_ms object

- type:

  "deviance" (default), "pearson", or "response"

- ...:

  ignored

## Value

named list of per-species residual lists
