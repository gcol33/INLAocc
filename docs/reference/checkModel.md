# Visual diagnostic panel for a fitted occupancy model

Produces a 2x2 (or 2x3) panel of diagnostic plots:

1.  QQ plot of PIT residuals vs Uniform(0, 1)

2.  Residuals vs fitted values

3.  Dispersion histogram (observed variance vs simulated)

4.  Moran's I correlogram (if coordinates are available)

## Usage

``` r
checkModel(object, nsim = 250L, seed = 123L)
```

## Arguments

- object:

  fitted `occu_inla` object

- nsim:

  number of simulations for PIT and dispersion (default 250)

- seed:

  random seed (default 123)

## Value

Invisible list of test results.

## Details

Analogous to `performance::check_model()` or
[`DHARMa::plot.DHARMa()`](https://rdrr.io/pkg/DHARMa/man/plot.DHARMa.html).
