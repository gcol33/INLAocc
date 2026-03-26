# Model-averaged predictions from multiple occupancy models

Computes weighted-average occupancy and detection probabilities across a
candidate set of models. Weights are derived from AIC, BIC, or WAIC via
[`compare_models`](https://gillescolling.com/INLAocc/reference/compare_models.md).

## Usage

``` r
modelAverage(
  ...,
  criterion = c("waic", "aic", "bic"),
  newdata = NULL,
  se = TRUE
)
```

## Arguments

- ...:

  named `occu_inla` objects (the candidate model set)

- criterion:

  `"waic"` (default), `"aic"`, or `"bic"`

- newdata:

  optional list for out-of-sample prediction (passed to
  [`predict.occu_inla`](https://gillescolling.com/INLAocc/reference/predict.occu_inla.md)).
  If `NULL`, returns in-sample averaged psi and p.

- se:

  if `TRUE` (default), returns unconditional standard errors that
  account for model selection uncertainty

## Value

A list of class `"occu_model_avg"` with:

- psi_hat:

  model-averaged occupancy probabilities

- p_hat:

  model-averaged detection probabilities (N x J)

- psi_se:

  unconditional SE for psi (if `se = TRUE`)

- weights:

  named vector of model weights

- comparison:

  data.frame from
  [`compare_models()`](https://gillescolling.com/INLAocc/reference/compare_models.md)

- criterion:

  IC used for weights

## Details

Follows the "full model averaging" approach of Burnham & Anderson
(2002): predictions from every model contribute, weighted by
information-criterion weights.
