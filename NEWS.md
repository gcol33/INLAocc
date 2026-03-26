# INLAocc 0.1.0

Initial release.

## Core Engine
* EM-INLA algorithm with scaled-binomial occupancy M-step.
* Novel EM + multiple imputation hybrid for beta debiasing (cor > 0.99 with MCMC).
* Adaptive damping and GLM warm start for convergence stability.
* Unified `occu()` API dispatching to all model types.

## Model Types
* Single-species: `occu()`, `occu_inla()`.
* Spatial SPDE: `occu(..., spatial = coords)`, `spatial_occu_inla()`.
* Multi-season AR(1): `occu(..., temporal = "ar1")`, `temporal_occu_inla()`.
* Multi-species community: `occu(..., multispecies = TRUE)`, `ms_occu_inla()`.
* Latent factor multi-species: `lfMsOccu_inla()`, `sfMsOccu_inla()`.
* Spatially-varying coefficients: `svcOccu_inla()`.
* Integrated multi-source: `intOccu_inla()`, `spIntOccu_inla()`.
* Joint SDM (no detection): `lfJSDM_inla()`.

## Random Effects
* Native AST-based formula parser for `(1 | group)`, `(x | group)`, `(x || group)` syntax.
* Supports nested grouping `(1 | a/b)`, operator expansion `(a*b | group)`.
* Pure AST walker — no regex on deparsed code.

## Diagnostics
* `ppcOccu()`: posterior predictive checks (Freeman-Tukey, chi-squared).
* `waicOccu()`: WAIC for model comparison (per-species with `by.sp`).
* `fitted()`, `residuals()`: S3 methods for all model types.
* `simulate()`: posterior predictive simulation (site-level and observation-level).
* `pitResiduals()`: PIT scaled residuals with integer randomisation.
* `testUniformity()`: KS test on PIT residuals.
* `testDispersion()`: simulation-based over/underdispersion test.
* `testOutliers()`: simulation envelope exceedance test.
* `testZeroInflation()`: observed vs expected zero-detection sites.
* `moranI()`: Moran's I with inverse-distance (default) or k-NN weights.
* `durbinWatson()`: Durbin-Watson test for temporal autocorrelation.
* `variogram()`: empirical semivariogram with `plot()` method.
* `checkModel()`: one-call 2x2 diagnostic panel (QQ, residuals, dispersion, correlogram).
* `dharma()`: optional DHARMa convenience bridge.

## Model Selection & Averaging
* `compare_models()`: AIC, BIC, and WAIC with delta and Akaike/Schwarz weights.
* `modelAverage()`: Burnham & Anderson full model averaging with unconditional SEs.
* `AIC()`, `BIC()` via `logLik()` with correct df and nobs attributes.

## S3 Generics
* `coef()`, `confint()`, `vcov()`, `nobs()`, `logLik()`.
* `update()` with formula modification support.
* `tidy()`, `glance()`: broom-compatible output.
* `ranef()`: random effect posterior summaries.
* `summary()`, `print()`, `plot()` for all model types.

## Prediction
* `predict()` with ggpredict-style `terms` argument.
* Design-matrix prediction via `X.0`, `X.p.0`, `coords.0`.
* `predict_spatial()`: convenience wrapper for spatial interpolation.
* `marginal_effect()`: single-covariate response curves.
* `richness()`: site-level species richness from multi-species models.

## Data Formatting
* `occu_format()`: spOccupancy-compatible data validation.
* `occu_data()`: pipe-friendly data construction from long-format data.frames.
* `occu_format_ms()`: multi-species data formatting.

## Post-Hoc Analysis
* `postHocLM()`: Bayesian or frequentist linear model on estimated parameters.
* `getSVCSamples()`: extract spatially-varying coefficient posteriors.

## Simulation
* `simulate_occu()`: single-species with optional spatial correlation.
* `simMsOcc()`: multi-species community data.
* `simTOcc()`: multi-season data with AR(1).
* `simIntOcc()`, `simIntMsOcc()`: integrated multi-source data.
* `simTMsOcc()`: temporal multi-species data.
