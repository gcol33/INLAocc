# Extract log-likelihood

Returns the observed-data log-likelihood from the final EM iteration.
The `df` attribute is set to the number of fixed effect coefficients
(occupancy + detection) and the `nobs` attribute to the number of non-NA
detection history entries.

## Usage

``` r
# S3 method for class 'occu_inla'
logLik(object, ...)

# S3 method for class 'occu_inla_temporal'
logLik(object, ...)
```

## Arguments

- object:

  fitted occu_inla object

- ...:

  ignored

## Value

An object of class `"logLik"`.
