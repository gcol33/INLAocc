# Compute residuals for an occupancy model

Occupancy residuals (site-level) and detection residuals (visit-level).

## Usage

``` r
# S3 method for class 'occu_inla'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  fitted occu_inla object

- type:

  "deviance" (default), "pearson", or "response"

- ...:

  ignored

## Value

list with:

- occ.resids:

  length-N occupancy residuals

- det.resids:

  N x J detection residuals
