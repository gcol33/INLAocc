# Changelog

## INLAocc 0.1.0

Initial release.

### Core Engine

- EM-INLA algorithm with scaled-binomial occupancy M-step.
- Novel EM + multiple imputation hybrid for beta debiasing (cor \> 0.99
  with MCMC).
- Adaptive damping and GLM warm start for convergence stability.
- Unified
  [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) API
  dispatching to all model types.

### Model Types

- Single-species:
  [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md),
  [`occu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Spatial SPDE: `occu(..., spatial = coords)`,
  [`spatial_occu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Multi-season AR(1): `occu(..., temporal = "ar1")`,
  [`temporal_occu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Multi-species community: `occu(..., multispecies = TRUE)`,
  [`ms_occu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Latent factor multi-species:
  [`lfMsOccu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md),
  [`sfMsOccu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Spatially-varying coefficients:
  [`svcOccu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Integrated multi-source:
  [`intOccu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md),
  [`spIntOccu_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).
- Joint SDM (no detection):
  [`lfJSDM_inla()`](https://gillescolling.com/INLAocc/reference/occu.md).

### Random Effects

- Native AST-based formula parser for `(1 | group)`, `(x | group)`,
  `(x || group)` syntax.
- Supports nested grouping `(1 | a/b)`, operator expansion
  `(a*b | group)`.
- Pure AST walker — no regex on deparsed code.

### Diagnostics

- [`ppcOccu()`](https://gillescolling.com/INLAocc/reference/ppcOccu.md):
  posterior predictive checks (Freeman-Tukey, chi-squared).
- [`waicOccu()`](https://gillescolling.com/INLAocc/reference/waicOccu.md):
  WAIC for model comparison (per-species with `by.sp`).
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html): S3 methods
  for all model types.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html): posterior
  predictive simulation (site-level and observation-level).
- [`pitResiduals()`](https://gillescolling.com/INLAocc/reference/pitResiduals.md):
  PIT scaled residuals with integer randomisation.
- [`testUniformity()`](https://gillescolling.com/INLAocc/reference/testUniformity.md):
  KS test on PIT residuals.
- [`testDispersion()`](https://gillescolling.com/INLAocc/reference/testDispersion.md):
  simulation-based over/underdispersion test.
- [`testOutliers()`](https://gillescolling.com/INLAocc/reference/testOutliers.md):
  simulation envelope exceedance test.
- [`testZeroInflation()`](https://gillescolling.com/INLAocc/reference/testZeroInflation.md):
  observed vs expected zero-detection sites.
- [`moranI()`](https://gillescolling.com/INLAocc/reference/moranI.md):
  Moran’s I with inverse-distance (default) or k-NN weights.
- [`durbinWatson()`](https://gillescolling.com/INLAocc/reference/durbinWatson.md):
  Durbin-Watson test for temporal autocorrelation.
- [`variogram()`](https://gillescolling.com/INLAocc/reference/variogram.md):
  empirical semivariogram with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
- [`checkModel()`](https://gillescolling.com/INLAocc/reference/checkModel.md):
  one-call 2x2 diagnostic panel (QQ, residuals, dispersion,
  correlogram).
- [`dharma()`](https://gillescolling.com/INLAocc/reference/dharma.md):
  optional DHARMa convenience bridge.

### Model Selection & Averaging

- [`compare_models()`](https://gillescolling.com/INLAocc/reference/compare_models.md):
  AIC, BIC, and WAIC with delta and Akaike/Schwarz weights.
- [`modelAverage()`](https://gillescolling.com/INLAocc/reference/modelAverage.md):
  Burnham & Anderson full model averaging with unconditional SEs.
- [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) via
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) with correct df and
  nobs attributes.

### S3 Generics

- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html).
- [`update()`](https://rdrr.io/r/stats/update.html) with formula
  modification support.
- [`tidy()`](https://gillescolling.com/INLAocc/reference/tidy.md),
  [`glance()`](https://gillescolling.com/INLAocc/reference/glance.md):
  broom-compatible output.
- [`ranef()`](https://gillescolling.com/INLAocc/reference/ranef.md):
  random effect posterior summaries.
- [`summary()`](https://rdrr.io/r/base/summary.html),
  [`print()`](https://rdrr.io/r/base/print.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for all model
  types.

### Prediction

- [`predict()`](https://rdrr.io/r/stats/predict.html) with
  ggpredict-style `terms` argument.
- Design-matrix prediction via `X.0`, `X.p.0`, `coords.0`.
- [`predict_spatial()`](https://gillescolling.com/INLAocc/reference/predict_spatial.md):
  convenience wrapper for spatial interpolation.
- [`marginal_effect()`](https://gillescolling.com/INLAocc/reference/marginal_effect.md):
  single-covariate response curves.
- [`richness()`](https://gillescolling.com/INLAocc/reference/richness.md):
  site-level species richness from multi-species models.

### Data Formatting

- [`occu_format()`](https://gillescolling.com/INLAocc/reference/occu_format.md):
  spOccupancy-compatible data validation.
- [`occu_data()`](https://gillescolling.com/INLAocc/reference/occu_data.md):
  pipe-friendly data construction from long-format data.frames.
- [`occu_format_ms()`](https://gillescolling.com/INLAocc/reference/occu_format_ms.md):
  multi-species data formatting.

### Post-Hoc Analysis

- [`postHocLM()`](https://gillescolling.com/INLAocc/reference/postHocLM.md):
  Bayesian or frequentist linear model on estimated parameters.
- [`getSVCSamples()`](https://gillescolling.com/INLAocc/reference/getSVCSamples.md):
  extract spatially-varying coefficient posteriors.

### Simulation

- [`simulate_occu()`](https://gillescolling.com/INLAocc/reference/simulate_occu.md):
  single-species with optional spatial correlation.
- [`simMsOcc()`](https://gillescolling.com/INLAocc/reference/simMsOcc.md):
  multi-species community data.
- [`simTOcc()`](https://gillescolling.com/INLAocc/reference/simTOcc.md):
  multi-season data with AR(1).
- [`simIntOcc()`](https://gillescolling.com/INLAocc/reference/simIntOcc.md),
  [`simIntMsOcc()`](https://gillescolling.com/INLAocc/reference/simIntMsOcc.md):
  integrated multi-source data.
- [`simTMsOcc()`](https://gillescolling.com/INLAocc/reference/simTMsOcc.md):
  temporal multi-species data.
