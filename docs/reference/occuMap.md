# Spatial occupancy map

Produces a publication-quality map of predicted occupancy probability,
uncertainty, or the spatial random effect from a fitted spatial model.
Survey sites are overlaid with detection status.

## Usage

``` r
occuMap(
  fit,
  type = c("psi", "sd", "spatial", "all"),
  species = NULL,
  n_grid = 100L,
  xlim = NULL,
  ylim = NULL,
  col = NULL,
  sites = TRUE,
  main = NULL,
  ...
)
```

## Arguments

- fit:

  fitted spatial occupancy model (class `occu_inla_spatial` or
  multi-species with spatial component)

- type:

  `"psi"` (default): occupancy probability. `"sd"`: posterior SD.
  `"spatial"`: spatial random effect only. `"all"`: 2-panel (psi + sd
  side by side).

- species:

  for multi-species models, the species name or index

- n_grid:

  grid resolution per axis (default 100)

- xlim, ylim:

  coordinate limits (default: data extent)

- col:

  color palette vector (default chosen per `type`)

- sites:

  logical: overlay survey sites? (default TRUE)

- main:

  plot title (default chosen per `type`)

- ...:

  passed to [`image()`](https://rdrr.io/r/graphics/image.html)

## Value

Invisibly, a list with `x`, `y`, `z` (the gridded values), `type`,
`zlim`, `coords`, `detected`, and `col`.
