# Compute WAIC for occupancy models

Returns WAIC for the occupancy and detection components separately and
combined. Analogous to spOccupancy::waicOcc().

## Usage

``` r
waicOccu(object, by.sp = FALSE)
```

## Arguments

- object:

  fitted occu_inla object

- by.sp:

  logical: if TRUE and object is multi-species, return per-species WAIC

## Value

data.frame with elpd, pD, and WAIC columns
