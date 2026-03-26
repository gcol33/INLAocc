# Test for over- or underdispersion

Compares the variance of observed site-level detection counts to the
variance expected under the fitted model (via simulation). A ratio \> 1
indicates overdispersion; \< 1 indicates underdispersion. Native
equivalent of
[`DHARMa::testDispersion()`](https://rdrr.io/pkg/DHARMa/man/testDispersion.html).

## Usage

``` r
testDispersion(
  object,
  nsim = 250L,
  seed = 123L,
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations (default 250)

- seed:

  random seed (default 123)

- alternative:

  `"two.sided"` (default), `"greater"` (overdispersion), or `"less"`
  (underdispersion)

## Value

A list of class `"htest"` with the dispersion ratio and simulation-based
p-value.
