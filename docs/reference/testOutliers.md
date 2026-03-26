# Test for outliers (simulation envelope)

Counts how many sites have observed detection counts outside the
min-to-max range of all simulations. Under a correct model, the expected
number of such outliers is approximately `2 * N * (1 / (nsim + 1))`. A
binomial test assesses whether more outliers than expected are present.
Native equivalent of
[`DHARMa::testOutliers()`](https://rdrr.io/pkg/DHARMa/man/testOutliers.html).

## Usage

``` r
testOutliers(object, nsim = 250L, seed = 123L)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations (default 250)

- seed:

  random seed (default 123)

## Value

A list of class `"htest"` (binomial test result).
