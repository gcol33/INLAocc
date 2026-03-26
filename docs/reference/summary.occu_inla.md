# Summary method for occu_inla fits

Summary method for occu_inla fits

## Usage

``` r
# S3 method for class 'occu_inla'
summary(
  object,
  level = NULL,
  quantiles = c(0.025, 0.5, 0.975),
  digits = 4L,
  ...
)
```

## Arguments

- object:

  fitted occu_inla object

- level:

  character vector of parameter groups to display. Options: `"beta"`
  (occupancy fixed), `"alpha"` (detection fixed), `"sigma.sq.psi"`
  (occupancy hyperparams), `"sigma.sq.p"` (detection hyperparams).
  Default shows all available.

- quantiles:

  numeric vector of quantile levels for credible intervals (default:
  `c(0.025, 0.5, 0.975)`).

- digits:

  number of digits to print (default 4).

- ...:

  additional arguments (ignored)

## Value

Invisibly returns a summary list with occupancy and detection fixed
effects, and estimated probabilities.
