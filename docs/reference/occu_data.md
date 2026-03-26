# Convert a data.frame to occupancy data

Converts a long-format data.frame (one row per site-visit) into the
structured list format required by
[`occu`](https://gillescolling.com/INLAocc/reference/occu.md).
Pipe-friendly: the data.frame is the first argument.

## Usage

``` r
occu_data(df, y, site, visit, occ.covs = NULL, det.covs = NULL, coords = NULL)
```

## Arguments

- df:

  a data.frame in long format (one row per site-visit).

- y:

  character: name of the detection column (0/1/NA).

- site:

  character: name of the site identifier column.

- visit:

  character: name of the visit/replicate column (integer or factor
  indicating visit number within each site).

- occ.covs:

  character vector of column names for site-level occupancy covariates.
  If `NULL` (default), all columns that are constant within each site
  (excluding `y`, `site`, `visit`, and `det.covs`) are used.

- det.covs:

  character vector of column names for visit-level detection covariates.
  If `NULL` (default), all columns that vary within at least one site
  are used.

- coords:

  character vector of length 2 giving coordinate column names, or
  `NULL`.

## Value

An object of class `"occu_data"` ready for
[`occu`](https://gillescolling.com/INLAocc/reference/occu.md).

## Examples

``` r
# Long-format data
df <- data.frame(
  site = rep(1:50, each = 4),
  visit = rep(1:4, 50),
  detected = rbinom(200, 1, 0.4),
  elev = rep(rnorm(50), each = 4),
  effort = runif(200, 1, 8)
)
dat <- occu_data(df, y = "detected", site = "site", visit = "visit")

# With pipe
# dat <- df |> occu_data(y = "detected", site = "site", visit = "visit")
```
