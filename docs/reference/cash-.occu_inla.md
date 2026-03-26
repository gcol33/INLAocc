# Access spOccupancy-compatible fields from INLAocc fits

Allows accessing spOccupancy-style posterior sample fields (e.g.,
`$beta.samples`, `$psi.samples`) on INLAocc model objects. Samples are
generated lazily via
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
on first access.

## Usage

``` r
# S3 method for class 'occu_inla'
x$name
```

## Arguments

- x:

  an INLAocc model object

- name:

  field name to access

## Value

The requested field value
