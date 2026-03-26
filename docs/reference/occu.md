# Fit occupancy models using INLA

Single entry point for all occupancy model types. The model structure is
controlled by flags: `spatial`, `temporal`, `multispecies`,
`integrated`, `n.factors`, and `svc`. Covers all model types from
spOccupancy (PGOcc, spPGOcc, tPGOcc, msPGOcc, etc.) via flag
combinations.

## Usage

``` r
occu(
  occ.formula,
  det.formula = NULL,
  data,
  spatial = NULL,
  temporal = NULL,
  multispecies = FALSE,
  integrated = FALSE,
  n.factors = NULL,
  svc = NULL,
  spde.args = list(),
  priors = NULL,
  occ.re = NULL,
  det.re = NULL,
  max.iter = 50L,
  tol = 1e-04,
  damping = 0.3,
  ensemble = FALSE,
  num.threads = "1:1",
  k.fold = 0L,
  verbose = 1L
)
```

## Arguments

- occ.formula:

  RHS formula for occupancy (psi). Supports mixed-model random effects
  syntax: `~ elev + forest + (1 | region)`.

- det.formula:

  RHS formula for detection (p). NULL for JSDM models (set
  `multispecies = "jsdm"`).

- data:

  An `occu_data` or `occu_data_ms` object, or a raw list with components
  `y`, `occ.covs`, `det.covs` (spOccupancy-compatible).

- spatial:

  N x 2 coordinate matrix or
  [`occu_spatial()`](https://gillescolling.com/INLAocc/reference/occu_spatial.md)
  object. Enables spatial SPDE models.

- temporal:

  `"ar1"` or `"iid"` for multi-season models. NULL (default) for
  single-season.

- multispecies:

  `FALSE` (default), `TRUE` for community models, or `"jsdm"` for joint
  species distribution models (no detection process).

- integrated:

  `TRUE` for multi-data-source models. Requires `data$y` to be a list of
  detection matrices.

- n.factors:

  Integer number of latent factors. Enables latent factor variants of
  multi-species models.

- svc:

  Integer vector of occupancy design matrix columns that get
  spatially-varying coefficients (1 = intercept). Requires `spatial`.

- spde.args:

  List of arguments passed to
  [`occu_spatial()`](https://gillescolling.com/INLAocc/reference/occu_spatial.md)
  when `spatial` is a coordinate matrix.

- priors:

  [`occu_priors()`](https://gillescolling.com/INLAocc/reference/occu_priors.md)
  object or spOccupancy-compatible named list.

- occ.re:

  Explicit list of
  [`occu_re()`](https://gillescolling.com/INLAocc/reference/occu_re.md)
  specs for occupancy.

- det.re:

  Explicit list of
  [`occu_re()`](https://gillescolling.com/INLAocc/reference/occu_re.md)
  specs for detection.

- max.iter:

  Maximum EM iterations (default 50).

- tol:

  Convergence tolerance (default 1e-4).

- damping:

  EM damping factor 0–1 (default 0.3).

- ensemble:

  Logical; if `TRUE`, use ensemble averaging across multiple imputation
  chains (default `FALSE`).

- num.threads:

  Thread specification for INLA, as `"A:B"` where A is outer and B is
  inner OpenMP threads (default `"1:1"`).

- k.fold:

  Number of cross-validation folds (default 0, no CV).

- verbose:

  0 = silent, 1 = iteration summaries, 2 = full INLA output.

## Value

An S3 object whose class depends on the model type (e.g., `"occu_inla"`,
`"occu_inla_spatial"`, `"occu_inla_ms"`).

## Examples

``` r
# \donttest{
# Single-species (cf. PGOcc)
# fit <- occu(~ elev, ~ effort, data)

# Spatial (cf. spPGOcc)
# fit <- occu(~ elev, ~ effort, data, spatial = coords)

# Multi-species (cf. msPGOcc)
# fit <- occu(~ elev, ~ effort, ms_data, multispecies = TRUE)

# Temporal (cf. tPGOcc)
# fit <- occu(~ elev, ~ effort, data, temporal = "ar1")
# }
```
