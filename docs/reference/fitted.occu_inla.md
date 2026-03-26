# Extract fitted values from an occupancy model

Returns fitted detection-level values (y.rep) and detection
probabilities. Analogous to spOccupancy's fitted() method.

## Usage

``` r
# S3 method for class 'occu_inla'
fitted(object, ...)
```

## Arguments

- object:

  fitted occu_inla object

- ...:

  ignored

## Value

list with:

- y.rep:

  N x J matrix of expected detection values (z_hat \* p_hat)

- p:

  N x J matrix of estimated detection probabilities

- psi:

  length-N vector of estimated occupancy probabilities

- z:

  length-N vector of posterior P(z=1)
