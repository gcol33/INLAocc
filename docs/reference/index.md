# Package index

## Unified Interface

One function for all model types.

- [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) : Fit
  occupancy models using INLA

## Single-Species Models

Non-spatial and spatial occupancy models for one species.

- [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) : Fit
  occupancy models using INLA

## Multi-Species Models

Community occupancy models with species correlations.

- [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) : Fit
  occupancy models using INLA

## Temporal & Integrated Models

Multi-season and multi-data-source models.

- [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) : Fit
  occupancy models using INLA

## Data Formatting

Prepare detection histories and covariates.

- [`occu_format()`](https://gillescolling.com/INLAocc/reference/occu_format.md)
  : Format data for INLA occupancy models
- [`occu_format_ms()`](https://gillescolling.com/INLAocc/reference/occu_format_ms.md)
  : Create multi-species occupancy data
- [`occu_data()`](https://gillescolling.com/INLAocc/reference/occu_data.md)
  : Convert a data.frame to occupancy data
- [`occu_areal()`](https://gillescolling.com/INLAocc/reference/occu_areal.md)
  : Create an areal spatial effect for occupancy models (CAR/BYM2)
- [`print(`*`<occu_data>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_data.md)
  : Print method for occu_data
- [`summary(`*`<occu_data>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_data.md)
  : Summary statistics for occupancy detection histories
- [`plot(`*`<occu_data>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_data.md)
  : Plot detection history patterns
- [`print(`*`<occu_spatial>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_spatial.md)
  : Print method for occu_spatial

## Model Specification

Random effects, priors, and spatial components.

- [`occu_re()`](https://gillescolling.com/INLAocc/reference/occu_re.md)
  : Specify random effects for occupancy or detection process
- [`occu_community_re()`](https://gillescolling.com/INLAocc/reference/occu_community_re.md)
  : Create a community-level (multi-species) random effect
- [`occu_priors()`](https://gillescolling.com/INLAocc/reference/occu_priors.md)
  : Specify priors for INLA occupancy models (spOccupancy-compatible)
- [`occu_spatial()`](https://gillescolling.com/INLAocc/reference/occu_spatial.md)
  : Create SPDE mesh and spatial effect for occupancy models

## Prediction & Derived Quantities

Predict occupancy at new locations, marginal effects, richness.

- [`predict(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/predict.occu_inla.md)
  : Predict from a fitted occupancy model
- [`predict(`*`<occu_inla_int>`*`)`](https://gillescolling.com/INLAocc/reference/predict.occu_inla_int.md)
  : Predict from an integrated occupancy model
- [`predict(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/predict.occu_inla_ms.md)
  : Predict from a multi-species occupancy model
- [`predict(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/predict.occu_inla_temporal.md)
  : Predict from a temporal occupancy model
- [`predict(`*`<occu_inla_jsdm>`*`)`](https://gillescolling.com/INLAocc/reference/predict.occu_inla_jsdm.md)
  : Predict from a JSDM
- [`predict_spatial()`](https://gillescolling.com/INLAocc/reference/predict_spatial.md)
  : Predict occupancy at new spatial locations
- [`marginal_effect()`](https://gillescolling.com/INLAocc/reference/marginal_effect.md)
  : Compute marginal effects for a covariate
- [`richness()`](https://gillescolling.com/INLAocc/reference/richness.md)
  : Compute species richness from multi-species model
- [`plot(`*`<occu_prediction>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_prediction.md)
  : Plot predicted effects from an occupancy model

## S3 Methods

Standard R generics for occupancy model objects.

- [`coef(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/coef.occu_inla.md)
  : Extract model coefficients
- [`confint(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/confint.occu_inla.md)
  : Compute credible intervals
- [`vcov(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/vcov.occu_inla.md)
  : Approximate variance-covariance matrix
- [`logLik(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/logLik.occu_inla.md)
  [`logLik(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/logLik.occu_inla.md)
  : Extract log-likelihood
- [`nobs(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/nobs.occu_inla.md)
  : Number of observations
- [`tidy(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/tidy.occu_inla.md)
  : Tidy model output into a data.frame
- [`glance(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/glance.occu_inla.md)
  : Glance at model-level statistics
- [`ranef(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/ranef.occu_inla.md)
  : Extract random effects
- [`update(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/update.occu_inla.md)
  : Update and refit an occupancy model
- [`tidy()`](https://gillescolling.com/INLAocc/reference/tidy.md) : Tidy
  model output into a data.frame
- [`ranef()`](https://gillescolling.com/INLAocc/reference/ranef.md) :
  Extract random effects
- [`glance()`](https://gillescolling.com/INLAocc/reference/glance.md) :
  Glance at model-level statistics
- [`` `$`( ``*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/cash-.occu_inla.md)
  : Access spOccupancy-compatible fields from INLAocc fits

## Summary, Print & Plot

Model output display methods.

- [`summary(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla.md)
  : Summary method for occu_inla fits
- [`summary(`*`<occu_inla_spatial>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla_spatial.md)
  : Summary for spatial occupancy model
- [`summary(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla_ms.md)
  : Summary for multi-species model
- [`summary(`*`<occu_inla_int>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla_int.md)
  : Summary for integrated occupancy model
- [`summary(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla_temporal.md)
  : Summary for temporal occupancy model
- [`summary(`*`<occu_inla_jsdm>`*`)`](https://gillescolling.com/INLAocc/reference/summary.occu_inla_jsdm.md)
  : Summary for JSDM
- [`print(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_inla.md)
  : Print method for occu_inla
- [`print(`*`<occu_inla_int>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_inla_int.md)
  : Print method for integrated occupancy model
- [`print(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_inla_ms.md)
  : Print method for multi-species occupancy model
- [`print(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_inla_temporal.md)
  : Print method for temporal occupancy model
- [`print(`*`<occu_inla_jsdm>`*`)`](https://gillescolling.com/INLAocc/reference/print.occu_inla_jsdm.md)
  : Print method for JSDM
- [`plot(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_inla.md)
  : Plot diagnostics for occu_inla
- [`plot(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_inla_ms.md)
  : Plot diagnostics for multi-species occupancy model
- [`plot(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_inla_temporal.md)
  : Plot diagnostics for temporal occupancy model
- [`plot(`*`<occu_inla_jsdm>`*`)`](https://gillescolling.com/INLAocc/reference/plot.occu_inla_jsdm.md)
  : Plot diagnostics for JSDM

## Fitted Values & Residuals

Extract fitted values and residuals for all model types.

- [`fitted(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/fitted.occu_inla.md)
  : Extract fitted values from an occupancy model
- [`fitted(`*`<occu_inla_int>`*`)`](https://gillescolling.com/INLAocc/reference/fitted.occu_inla_int.md)
  : Extract fitted values from an integrated occupancy model
- [`fitted(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/fitted.occu_inla_ms.md)
  : Extract fitted values from a multi-species occupancy model
- [`fitted(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/fitted.occu_inla_temporal.md)
  : Extract fitted values from a temporal occupancy model
- [`fitted(`*`<occu_inla_jsdm>`*`)`](https://gillescolling.com/INLAocc/reference/fitted.occu_inla_jsdm.md)
  : Extract fitted values from a JSDM
- [`residuals(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/residuals.occu_inla.md)
  : Compute residuals for an occupancy model
- [`residuals(`*`<occu_inla_int>`*`)`](https://gillescolling.com/INLAocc/reference/residuals.occu_inla_int.md)
  : Compute residuals for an integrated occupancy model
- [`residuals(`*`<occu_inla_ms>`*`)`](https://gillescolling.com/INLAocc/reference/residuals.occu_inla_ms.md)
  : Compute residuals for a multi-species occupancy model
- [`residuals(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/residuals.occu_inla_temporal.md)
  : Compute residuals for a temporal occupancy model

## Residual Diagnostics

Simulation-based diagnostics.

- [`simulate(`*`<occu_inla>`*`)`](https://gillescolling.com/INLAocc/reference/simulate.occu_inla.md)
  [`simulate(`*`<occu_inla_temporal>`*`)`](https://gillescolling.com/INLAocc/reference/simulate.occu_inla.md)
  : Simulate replicate datasets from a fitted occupancy model
- [`pitResiduals()`](https://gillescolling.com/INLAocc/reference/pitResiduals.md)
  : Compute PIT (scaled) residuals for an occupancy model
- [`testUniformity()`](https://gillescolling.com/INLAocc/reference/testUniformity.md)
  : Test uniformity of PIT residuals (KS test)
- [`testDispersion()`](https://gillescolling.com/INLAocc/reference/testDispersion.md)
  : Test for over- or underdispersion
- [`testOutliers()`](https://gillescolling.com/INLAocc/reference/testOutliers.md)
  : Test for outliers (simulation envelope)
- [`testZeroInflation()`](https://gillescolling.com/INLAocc/reference/testZeroInflation.md)
  : Test for zero-inflation
- [`checkModel()`](https://gillescolling.com/INLAocc/reference/checkModel.md)
  : Visual diagnostic panel for a fitted occupancy model
- [`dharma()`](https://gillescolling.com/INLAocc/reference/dharma.md) :
  Create a DHARMa residuals object from a fitted occupancy model

## Spatial & Temporal Diagnostics

Autocorrelation testing and semivariograms.

- [`moranI()`](https://gillescolling.com/INLAocc/reference/moranI.md) :
  Moran's I test for spatial autocorrelation
- [`durbinWatson()`](https://gillescolling.com/INLAocc/reference/durbinWatson.md)
  : Durbin-Watson test for temporal autocorrelation in residuals
- [`variogram()`](https://gillescolling.com/INLAocc/reference/variogram.md)
  : Empirical semivariogram of occupancy residuals
- [`spatialRange()`](https://gillescolling.com/INLAocc/reference/spatialRange.md)
  : Extract spatial range and standard deviation from a fitted spatial
  model
- [`temporalCorr()`](https://gillescolling.com/INLAocc/reference/temporalCorr.md)
  : Extract AR(1) temporal correlation from a fitted temporal model
- [`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
  : Diagnose identifiability issues in occupancy models
- [`occuMap()`](https://gillescolling.com/INLAocc/reference/occuMap.md)
  : Spatial occupancy map

## Model Comparison & Averaging

Information criteria, weights, and multi-model inference.

- [`compare_models()`](https://gillescolling.com/INLAocc/reference/compare_models.md)
  : Compare occupancy models via information criteria
- [`modelAverage()`](https://gillescolling.com/INLAocc/reference/modelAverage.md)
  : Model-averaged predictions from multiple occupancy models
- [`waicOccu()`](https://gillescolling.com/INLAocc/reference/waicOccu.md)
  : Returns WAIC for the occupancy and detection components separately
  and combined. Analogous to spOccupancy::waicOcc().
- [`ppcOccu()`](https://gillescolling.com/INLAocc/reference/ppcOccu.md)
  : Posterior predictive checks for occupancy models

## Post-Hoc Analysis

Secondary analyses on fitted model outputs.

- [`postHocLM()`](https://gillescolling.com/INLAocc/reference/postHocLM.md)
  : Fit a post-hoc linear model relating covariates to estimated
  parameters
- [`getSVCSamples()`](https://gillescolling.com/INLAocc/reference/getSVCSamples.md)
  : Extract spatially-varying coefficient summaries

## Simulation

Generate synthetic occupancy data for testing.

- [`simulate_occu()`](https://gillescolling.com/INLAocc/reference/simulate_occu.md)
  : Simulate occupancy data for testing
- [`simMsOcc()`](https://gillescolling.com/INLAocc/reference/simMsOcc.md)
  : Simulate multi-species occupancy data (cf. simMsOcc)
- [`simTOcc()`](https://gillescolling.com/INLAocc/reference/simTOcc.md)
  : Simulate multi-season occupancy data (cf. simTOcc)
- [`simTMsOcc()`](https://gillescolling.com/INLAocc/reference/simTMsOcc.md)
  : Simulate temporal multi-species occupancy data (cf. simTMsOcc)
- [`simIntOcc()`](https://gillescolling.com/INLAocc/reference/simIntOcc.md)
  : Simulate integrated (multi-source) occupancy data (cf. simIntOcc)
- [`simIntMsOcc()`](https://gillescolling.com/INLAocc/reference/simIntMsOcc.md)
  : Simulate integrated multi-species occupancy data (cf. simIntMsOcc)
