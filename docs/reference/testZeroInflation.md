# Test for zero-inflation

Compares the number of all-zero sites (no detections) in the observed
data to the distribution expected under the fitted model. Native
equivalent of
[`DHARMa::testZeroInflation()`](https://rdrr.io/pkg/DHARMa/man/testZeroInflation.html).

## Usage

``` r
testZeroInflation(object, nsim = 250L, seed = 123L)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations (default 250)

- seed:

  random seed (default 123)

## Value

A list of class `"htest"` with the zero-inflation ratio and
simulation-based p-value.
