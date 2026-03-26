# Test uniformity of PIT residuals (KS test)

If the model is correct, PIT residuals should follow Uniform(0, 1). This
test applies a Kolmogorov-Smirnov test against that null. Native
equivalent of
[`DHARMa::testUniformity()`](https://rdrr.io/pkg/DHARMa/man/testUniformity.html).

## Usage

``` r
testUniformity(object, nsim = 250L, seed = 123L, plot = FALSE)
```

## Arguments

- object:

  fitted `occu_inla` object, or a numeric vector of PIT residuals (from
  [`pitResiduals`](https://gillescolling.com/INLAocc/reference/pitResiduals.md))

- nsim:

  number of simulations if `object` is a fit (default 250)

- seed:

  random seed (default 123)

- plot:

  if `TRUE`, draws a QQ plot of residuals vs Uniform

## Value

A list of class `"htest"` (KS test result).
