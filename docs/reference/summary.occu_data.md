# Summary statistics for occupancy detection histories

Computes descriptive statistics for the detection matrix before model
fitting: detection frequencies, visit completeness, per-visit detection
rates, and spatial patterns (if coordinates are available).

## Usage

``` r
# S3 method for class 'occu_data'
summary(object, ...)
```

## Arguments

- object:

  an `occu_data` object or a raw data list with component `y`

- ...:

  ignored

## Value

A list of class `"occu_data_summary"` with:

- N:

  number of sites

- J:

  max visits per site

- n_obs:

  total non-NA observations

- n_missing:

  total NAs in detection matrix

- naive_psi:

  proportion of sites with at least one detection

- naive_p:

  detection rate conditional on occupied sites

- det_freq:

  table of detection counts per site (0, 1, 2, ...)

- visit_rate:

  proportion of non-NA cells per visit column

- det_per_visit:

  detection rate per visit column

- has_coords:

  logical
