# spOccupancy → INLAocc Compatibility Matrix

Status: ✅ = covered \| ⚠️ = partial/different \| ❌ = missing \| N/A =
not applicable (MCMC-specific)

## 1. Output Fields (`$` accessors)

### Core Parameter Samples

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$beta.samples` | Occ coefficient MCMC samples (coda) | `$occ_fit$summary.fixed` (summary only) | ⚠️ no samples |
| `$alpha.samples` | Det coefficient MCMC samples (coda) | `$det_fit$summary.fixed` (summary only) | ⚠️ no samples |
| `$z.samples` | Latent occupancy state samples | `$z_hat` (point estimate only) | ⚠️ no samples |
| `$psi.samples` | Occupancy probability samples | `$psi_hat` (point estimate only) | ⚠️ no samples |
| `$p.samples` | Detection probability samples | `$p_hat` (point estimate only) | ⚠️ no samples |
| `$sigma.sq.psi.samples` | Occ RE variance samples | `$occ_fit$summary.hyperpar` | ⚠️ no samples |
| `$sigma.sq.p.samples` | Det RE variance samples | `$det_fit$summary.hyperpar` | ⚠️ no samples |
| `$beta.star.samples` | Occ random effect realizations | `$occ_fit$summary.random` | ⚠️ no samples |
| `$alpha.star.samples` | Det random effect realizations | `$det_fit$summary.random` | ⚠️ no samples |

### Spatial Fields

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$theta.samples` | Spatial covariance params (coda) | `$occ_fit$summary.hyperpar` (Range/Stdev rows) | ⚠️ no samples |
| `$w.samples` | Spatial random effects (coda) | `$occ_fit$summary.random$spatial` | ⚠️ no samples |

### Community (Multi-Species) Fields

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$beta.comm.samples` | Community occ coefficients | `$community_occ` (data.frame summary) | ⚠️ no samples |
| `$alpha.comm.samples` | Community det coefficients | `$community_det` (data.frame summary) | ⚠️ no samples |
| `$tau.sq.beta.samples` | Community occ variance | `$community_occ$species_sd` | ⚠️ no samples |
| `$tau.sq.alpha.samples` | Community det variance | `$community_det$species_sd` | ⚠️ no samples |

### Latent Factor Fields

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$lambda.samples` | Factor loadings (3D: sample×sp×factor) | `$lambda` (matrix, point estimate) | ⚠️ no samples |
| `$w.samples` (factor) | Latent factor scores | `$factors` (matrix, point estimate) | ⚠️ no samples |

### Temporal Fields

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$eta.samples` | AR(1) temporal effects | `$psi_smoothed` (smoothed estimates) | ⚠️ different form |

### Diagnostics & Metadata

| spOccupancy field | Description | INLAocc equivalent | Status |
|----|----|----|----|
| `$rhat` | Gelman-Rubin convergence | — | N/A (no MCMC) |
| `$ESS` | Effective sample size | — | N/A (no MCMC) |
| `$like.samples` | Likelihood values for WAIC | `$occ_fit$waic`, `$det_fit$waic` | ⚠️ computed internally |
| `$run.time` | Computation time | — | ❌ missing |
| `$k.fold.deviance` | CV deviance | `$k.fold$k.fold.deviance` | ✅ (nested deeper) |

------------------------------------------------------------------------

## 2. S3 Methods

### summary()

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| Posterior quantiles (2.5%, 50%, 97.5%) | Fixed effect mean, sd, quantiles from INLA | ✅ |
| Rhat column | — | N/A |
| ESS column | — | N/A |
| `level` parameter (select which params) | Not supported | ❌ missing |
| `quantiles` parameter | Hardcoded 2.5%/97.5% | ❌ not configurable |

### print()

| spOccupancy feature   | INLAocc | Status |
|-----------------------|---------|--------|
| Concise model summary | ✅      | ✅     |

### plot()

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| Trace plots of MCMC chains | Histogram/convergence plots | ⚠️ different (no MCMC) |
| Parameter selection | `which` argument (plot type, not param) | ⚠️ different interface |

### fitted()

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| `$y.rep.samples` (3D: sample×site×visit) | `$y.rep` (matrix, point estimate) | ⚠️ no samples |
| `$p.samples` (3D: sample×site×visit) | `$p` (matrix, point estimate) | ⚠️ no samples |

### predict()

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| `$psi.0.samples` (coda) | `$psi.0$samples` (matrix) | ✅ (via inla.posterior.sample) |
| `$z.0.samples` (coda) | `$z.0$mean` (point only) | ⚠️ no samples |
| `$p.0.samples` (coda) | `$p.0$mean` (point only) | ⚠️ no samples |
| `$w.0.samples` (spatial RE at new locs) | Not returned separately | ❌ missing |
| `X.0` argument | ✅ | ✅ |
| `coords.0` argument | ✅ | ✅ |
| `ignore.RE` argument | ✅ (present but not fully used) | ⚠️ stub |
| `type` argument | ✅ (“occupancy”, “detection”, “both”) | ✅ |
| `n.omp.threads` argument | — | N/A |
| `n.report` argument | — | N/A |

### residuals()

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| `$occ.resids` (sample×site matrix) | `$occ.resids` (vector, point) | ⚠️ no samples |
| `$det.resids` (3D: sample×site×visit) | `$det.resids` (matrix, point) | ⚠️ no samples |
| `type = "marginal"` | `type = "deviance"/"pearson"/"response"` | ⚠️ different types |

------------------------------------------------------------------------

## 3. Diagnostic Functions

### ppcOcc / ppcOccu

| spOccupancy feature       | INLAocc                      | Status     |
|---------------------------|------------------------------|------------|
| `$fit.y`                  | `$fit.y`                     | ✅         |
| `$fit.y.rep`              | `$fit.y.rep`                 | ✅         |
| `$fit.y.group.quants`     | —                            | ❌ missing |
| `$fit.y.rep.group.quants` | —                            | ❌ missing |
| `type = "chi-square"`     | `fit.stat = "chi-squared"`   | ✅         |
| `type = "freeman-tukey"`  | `fit.stat = "freeman-tukey"` | ✅         |
| `$bayesian.p`             | `$bayesian.p`                | ✅         |

### waicOcc / waicOccu

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| Returns `c(elpd, pD, WAIC)` vector | Returns data.frame with component column | ⚠️ different format |
| `by.sp = TRUE` per-species | ✅ | ✅ |

### getSVCSamples

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| Returns list of coda::mcmc objects (sample×site) | Returns list of data.frames (mean, sd, quantiles) | ⚠️ no samples |

### postHocLM

| spOccupancy feature | INLAocc | Status |
|----|----|----|
| Returns Bayesian model with `$beta.samples`, `$bayes.R2` | Returns `lm` object | ⚠️ frequentist, not Bayesian |

### updateMCMC

| spOccupancy feature | INLAocc | Status        |
|---------------------|---------|---------------|
| Extend MCMC chains  | —       | N/A (no MCMC) |

------------------------------------------------------------------------

## 4. Simulation Functions

| spOccupancy | INLAocc | Status |
|----|----|----|
| `simOcc()` | [`simulate_occu()`](https://gillescolling.com/INLAocc/reference/simulate_occu.md) | ✅ |
| [`simMsOcc()`](https://gillescolling.com/INLAocc/reference/simMsOcc.md) | [`simMsOcc()`](https://gillescolling.com/INLAocc/reference/simMsOcc.md) | ✅ |
| [`simTOcc()`](https://gillescolling.com/INLAocc/reference/simTOcc.md) | [`simTOcc()`](https://gillescolling.com/INLAocc/reference/simTOcc.md) | ✅ |
| [`simIntOcc()`](https://gillescolling.com/INLAocc/reference/simIntOcc.md) | [`simIntOcc()`](https://gillescolling.com/INLAocc/reference/simIntOcc.md) | ✅ |
| [`simTMsOcc()`](https://gillescolling.com/INLAocc/reference/simTMsOcc.md) | — | ❌ missing |
| [`simIntMsOcc()`](https://gillescolling.com/INLAocc/reference/simIntMsOcc.md) | — | ❌ missing |
| `simBinom()` | — | ❌ (non-occupancy, low priority) |
| `simTBinom()` | — | ❌ (non-occupancy, low priority) |

------------------------------------------------------------------------

## 5. Datasets

| spOccupancy  | INLAocc | Status     |
|--------------|---------|------------|
| `hbef2015`   | —       | ❌ missing |
| `hbefTrends` | —       | ❌ missing |
| `hbefElev`   | —       | ❌ missing |
| `neon2015`   | —       | ❌ missing |

------------------------------------------------------------------------

## 6. Summary of Gaps

### Systemic gap: No posterior samples

spOccupancy returns full MCMC sample matrices for every parameter.
INLAocc returns point estimates and INLA summary statistics. This
affects: - Every `$*.samples` field (20+ fields across model types) -
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) return values -
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return values -
[`getSVCSamples()`](https://gillescolling.com/INLAocc/reference/getSVCSamples.md)
return values - [`predict()`](https://rdrr.io/r/stats/predict.html)
partially (psi.0 has samples, but z.0 and p.0 don’t)

**Fix:** Add a `$.occu_inla` accessor that lazily generates samples via
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
when `$beta.samples`, `$alpha.samples`, etc. are accessed. This provides
backwards compatibility without pre-computing unused samples.

### Individual gaps

| Gap | Priority | Fix |
|----|----|----|
| No `$run.time` field | Low | Add [`proc.time()`](https://rdrr.io/r/base/proc.time.html) tracking in [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) |
| [`summary()`](https://rdrr.io/r/base/summary.html) not configurable (level, quantiles) | Medium | Add `level` and `quantiles` args |
| [`ppcOccu()`](https://gillescolling.com/INLAocc/reference/ppcOccu.md) missing `$fit.y.group.quants` | Low | Add quantile computation |
| [`waicOccu()`](https://gillescolling.com/INLAocc/reference/waicOccu.md) returns data.frame not vector | Medium | Match spOcc format or document difference |
| [`postHocLM()`](https://gillescolling.com/INLAocc/reference/postHocLM.md) is frequentist not Bayesian | Medium | Rewrite with INLA or document as intentional |
| [`predict()`](https://rdrr.io/r/stats/predict.html) missing `$w.0.samples` | Low | Extract from spatial field |
| [`residuals()`](https://rdrr.io/r/stats/residuals.html) missing `type = "marginal"` | Low | Add marginal type |
| Missing [`simTMsOcc()`](https://gillescolling.com/INLAocc/reference/simTMsOcc.md), [`simIntMsOcc()`](https://gillescolling.com/INLAocc/reference/simIntMsOcc.md) | Medium | Write simulation functions |
| Missing example datasets | Medium | Include or create equivalent datasets |
| `ignore.RE` in predict is a stub | Low | Implement properly |
