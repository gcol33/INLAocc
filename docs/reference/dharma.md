# Create a DHARMa residuals object from a fitted occupancy model

Thin convenience wrapper that calls
[`simulate()`](https://rdrr.io/r/stats/simulate.html) on the fitted
model and passes the result to
[`DHARMa::createDHARMa()`](https://rdrr.io/pkg/DHARMa/man/createDHARMa.html).
All DHARMa tests (`testSpatialAutocorrelation`, `testDispersion`, etc.)
then work on the returned object.

## Usage

``` r
dharma(object, nsim = 250L, seed = 123L, ...)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations (default 250)

- seed:

  random seed (default 123)

- ...:

  passed to
  [`DHARMa::createDHARMa()`](https://rdrr.io/pkg/DHARMa/man/createDHARMa.html)

## Value

A `DHARMa` object (see
[`createDHARMa`](https://rdrr.io/pkg/DHARMa/man/createDHARMa.html))

## Details

Requires DHARMa to be installed. For a dependency-free alternative, use
[`moranI`](https://gillescolling.com/INLAocc/reference/moranI.md),
[`durbinWatson`](https://gillescolling.com/INLAocc/reference/durbinWatson.md),
or
[`variogram`](https://gillescolling.com/INLAocc/reference/variogram.md)
directly.
