# Compute credible intervals

Compute credible intervals

## Usage

``` r
# S3 method for class 'occu_inla'
confint(
  object,
  parm,
  level = 0.95,
  process = c("both", "occupancy", "detection"),
  ...
)
```

## Arguments

- object:

  fitted occu_inla object

- parm:

  not used (included for S3 compatibility)

- level:

  credible level (default 0.95)

- process:

  `"occupancy"` (default), `"detection"`, or `"both"`.

- ...:

  ignored

## Value

Matrix with lower and upper columns, or a list of matrices if
`process = "both"`.
