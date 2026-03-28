# Migrating from spOccupancy

## Introduction

spOccupancy (Doser et al. 2022) is the standard R package for Bayesian
occupancy models. It covers single-species, multi-species, spatial,
temporal, and integrated designs, all fitted via Polya-Gamma MCMC. If
you have used spOccupancy for any of these model types, you already know
the data format and the modelling workflow that INLAocc builds on.

INLAocc replaces the MCMC backend with INLA (integrated nested Laplace
approximation), an analytic approximation to the posterior that avoids
sampling entirely. The user-facing consequence is speed: models that
take minutes under MCMC finish in seconds under INLA, with no burn-in,
no thinning, and no convergence diagnostics. The statistical consequence
is that posteriors are approximate rather than exact — but for
well-identified models with moderate to large data, the difference is
negligible (benchmark correlations between INLA and MCMC posteriors
exceed 0.99).

This vignette covers everything you need to move an existing spOccupancy
workflow to INLAocc: data compatibility, function mapping, side-by-side
code examples, posterior sample access, diagnostics equivalences, and
honest caveats about when MCMC remains the better choice. For
model-specific details, see
[`vignette("spatial-models")`](https://gillescolling.com/INLAocc/articles/spatial-models.md),
[`vignette("multi-species")`](https://gillescolling.com/INLAocc/articles/multi-species.md),
and
[`vignette("diagnostics")`](https://gillescolling.com/INLAocc/articles/diagnostics.md).

## Why switch?

**Speed.** INLAocc is 5–20x faster than spOccupancy across model types.
On 300 sites with 4 visits, a spatial occupancy model runs in 15–30
seconds instead of 3–5 minutes. On 1000 sites, 1–2 minutes instead of
20–30 minutes. The gap widens with dataset size because INLA scales
sublinearly with the number of sites.

**Deterministic results.** The EM-INLA algorithm converges to the same
point estimate and posterior approximation every run. No chain-to-chain
variability, no need to set random seeds for reproducibility.

**No MCMC tuning.** spOccupancy requires setting `n.batch`,
`batch.length`, `n.burn`, `n.thin`, `n.chains`, `n.neighbors`, and
`cov.model` — each a tuning decision that affects results and runtime.
INLAocc has no equivalent parameters. Fewer researcher degrees of
freedom means fewer opportunities for accidental p-hacking of posterior
summaries.

**Unified interface.** One function
([`occu()`](https://gillescolling.com/INLAocc/reference/occu.md))
handles all 17 model types via keyword arguments. No need to remember
whether the spatial multi-species version is
[`spMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/spMsPGOcc.html)
or
[`stMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/stMsPGOcc.html).

**When NOT to switch.** Three situations favour staying with
spOccupancy:

1.  You need exact posterior samples for custom derived quantities that
    go beyond what [`predict()`](https://rdrr.io/r/stats/predict.html)
    provides. INLAocc generates posterior samples via INLA’s Gaussian
    approximation, which is accurate for marginals but does not capture
    the full joint posterior.
2.  You need informative inverse-Gamma priors on specific variance
    components. INLAocc uses PC priors (penalizing complexity), which
    are conceptually different and not easily mapped to inverse-Gamma
    hyperparameters.
3.  The INLA approximation quality matters for your problem: very small
    datasets (\< 30 sites), weakly identified models (many covariates,
    low detection probability), or models where posterior skewness is
    scientifically important.

## Data format compatibility

INLAocc accepts spOccupancy-formatted data lists directly, with no
conversion or reformatting:

``` r

library(INLAocc)

# spOccupancy data format --- used by both packages
data <- list(
  y = y_matrix,                          # N x J detection history
  occ.covs = data.frame(elev, forest),   # site-level covariates
  det.covs = list(effort = effort_mat),  # visit-level covariates (N x J)
  coords = coords_matrix                 # N x 2 coordinates (if spatial)
)

# Works directly in INLAocc --- no conversion needed
fit <- occu(~ elev + forest, ~ effort, data = data, verbose = 0)
```

**Multi-species data.** The detection array `y` is a 3D array with
dimensions `[species, sites, visits]`, or a named list of site-by-visit
matrices (one per species). Both formats work in INLAocc without
modification.

**Multi-season (temporal) data.** The detection array is a 3D array with
dimensions `[sites, seasons, visits]`. Same format as spOccupancy.

**Structured data objects.** INLAocc also provides
[`occu_data()`](https://gillescolling.com/INLAocc/reference/occu_data.md)
and
[`occu_format()`](https://gillescolling.com/INLAocc/reference/occu_format.md)
for constructing validated data objects with built-in checks. These are
optional — raw spOccupancy lists work without conversion. But if you are
starting a new analysis (rather than migrating an existing one), the
structured constructors catch formatting errors earlier.

## Function mapping

Every spOccupancy model function maps to a single
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) call
with different keyword arguments. The table below covers all 17 model
types:

| spOccupancy function | INLAocc equivalent | Notes |
|:---|:---|:---|
| [`PGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/PGOcc.html) | `occu(~ occ, ~ det, data)` | Basic single-species |
| [`spPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/spPGOcc.html) | `occu(..., spatial = coords)` | Spatial (SPDE mesh) |
| [`tPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/tPGOcc.html) | `occu(..., temporal = "ar1")` | Multi-season |
| [`stPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/stPGOcc.html) | `occu(..., spatial = coords, temporal = "ar1")` | Spatio-temporal |
| [`svcPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/svcPGOcc.html) | `occu(..., spatial = coords, svc = k)` | Spatially varying coefficients |
| [`svcTPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/svcTPGOcc.html) | `occu(..., spatial = coords, svc = k, temporal = "ar1")` | SVC + temporal |
| [`msPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/msPGOcc.html) | `occu(..., multispecies = TRUE)` | Community (multi-species) |
| [`spMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/spMsPGOcc.html) | `occu(..., multispecies = TRUE, spatial = coords)` | Spatial community |
| [`tMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/tMsPGOcc.html) | `occu(..., multispecies = TRUE, temporal = "ar1")` | Temporal community |
| [`stMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/stMsPGOcc.html) | `occu(..., multispecies = TRUE, spatial = coords, temporal = "ar1")` | Spatio-temporal community |
| [`svcTMsPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/svcTMsPGOcc.html) | `occu(..., multispecies = TRUE, spatial = coords, svc = k, temporal = "ar1")` | MS + SVC + temporal |
| `lfMsOcc()` | `occu(..., multispecies = TRUE, n.factors = k)` | Latent factor community |
| `sfMsOcc()` | `occu(..., multispecies = TRUE, n.factors = k, spatial = coords)` | Spatial latent factor |
| [`lfJSDM()`](https://www.doserlab.com/files/spoccupancy-web/reference/lfJSDM.html) | `occu(~ occ, data, multispecies = "jsdm", n.factors = k)` | Joint species distribution |
| [`sfJSDM()`](https://www.doserlab.com/files/spoccupancy-web/reference/sfJSDM.html) | `occu(~ occ, data, multispecies = "jsdm", n.factors = k, spatial = coords)` | Spatial JSDM |
| [`intPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/intPGOcc.html) | `occu(..., integrated = TRUE)` | Integrated data sources |
| [`spIntPGOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/spIntPGOcc.html) | `occu(..., spatial = coords, integrated = TRUE)` | Spatial integrated |

The pattern is consistent: `spatial = coords` adds a spatial field,
`temporal = "ar1"` adds temporal structure, `multispecies = TRUE`
enables community modelling, and `n.factors = k` activates latent factor
dimension reduction.

## Side-by-side example: single-species spatial model

The clearest way to see the difference is to fit the same model in both
packages. We use
[`spOccupancy::simOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/simOcc.html)
to generate data with known parameters, then compare the code and
output.

### spOccupancy version

``` r

library(spOccupancy)

# Simulate spatially structured occupancy data
dat <- simOcc(J.x = 15, J.y = 15, n.rep = 4,
              beta = c(0.5, -0.8), alpha = c(0, -0.5),
              sp = TRUE, sigma.sq = 1, phi = 5)

# Fit spatial occupancy model via Polya-Gamma MCMC
out <- spPGOcc(
  occ.formula = ~ occ.cov.1,
  det.formula  = ~ det.cov.1,
  data         = dat,
  n.batch      = 400,
  batch.length = 25,
  n.burn       = 5000,
  n.thin       = 5,
  n.chains     = 3,
  NNGP         = TRUE,
  n.neighbors  = 15,
  cov.model    = "exponential"
)
summary(out)
```

Note the seven tuning parameters (`n.batch`, `batch.length`, `n.burn`,
`n.thin`, `n.chains`, `NNGP`, `n.neighbors`) plus the covariance model
choice. Each affects runtime and posterior quality.

### INLAocc version

``` r

library(INLAocc)

# Same data, no conversion needed
fit <- occu(~ occ.cov.1, ~ det.cov.1, data = dat,
            spatial = dat$coords, verbose = 0)
summary(fit)
```

No tuning parameters. The EM algorithm converges deterministically.
Coefficient estimates should match the MCMC posterior means closely —
the true values are \\\beta = (0.5, -0.8)\\ and \\\alpha = (0, -0.5)\\,
and both packages should recover them.

## Formula syntax differences

The main syntactic difference is how formulas are passed:

``` r

# spOccupancy: named formula arguments
PGOcc(occ.formula = ~ elev + (1 | region),
      det.formula  = ~ effort,
      data = data, ...)

# INLAocc: positional formulas (occupancy first, detection second)
occu(~ elev + (1 | region),
     ~ effort,
     data = data, verbose = 0)
```

Both packages support lme4-style random effects in either formula:

- `(1 | group)` — random intercept

- `(x | group)` — correlated random intercept and slope

- `(x || group)` — uncorrelated random intercept and slope

The formula parsing is otherwise identical: standard R formula syntax
with `+`, `:`, `*`, [`I()`](https://rdrr.io/r/base/AsIs.html),
[`poly()`](https://rdrr.io/r/stats/poly.html), and `ns()` all working as
expected. The only restriction is that interaction terms between
occupancy and detection covariates are not supported (this is a
limitation of occupancy models generally, not of either package).

## Posterior sample compatibility

INLAocc provides a `$` accessor that lazily generates posterior samples
in the same format as spOccupancy:

``` r

fit <- occu(~ elev, ~ effort, data = data, verbose = 0)

# Access samples exactly like spOccupancy
fit$beta.samples       # occupancy coefficients (n_samples x n_occ_covs)
fit$alpha.samples      # detection coefficients (n_samples x n_det_covs)
fit$z.samples          # latent occupancy states (n_samples x N)
fit$psi.samples        # occupancy probabilities (n_samples x N)
```

Samples are generated on first access via
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
and cached in the fit object. This means downstream code that operates
on spOccupancy posterior matrices — computing derived quantities,
feeding into custom plotting functions, or passing to other packages —
works with INLAocc fits without modification.

However, INLAocc also provides native S3 methods that are often more
convenient than working with raw sample matrices:

``` r

# Tidier alternatives to manual posterior summaries
coef(fit)              # posterior means (vs. colMeans(out$beta.samples))
confint(fit)           # credible intervals
tidy(fit)              # broom-compatible data.frame with term, estimate, CI
```

If you are writing new analysis code (rather than migrating existing
scripts), the S3 methods are the recommended interface. They handle
multi-species indexing, temporal indexing, and spatial predictions
correctly without manual bookkeeping.

## Diagnostics comparison

The table below maps common diagnostic tasks between the two packages:

| Task | spOccupancy | INLAocc |
|:---|:---|:---|
| Posterior predictive check | `ppcOcc(out, "freeman-tukey", 1)` | `ppcOccu(fit, "freeman-tukey", 1)` |
| WAIC | `waicOcc(out)` | `waicOccu(fit)` |
| K-fold CV | `out$k.fold.deviance` | `fit$k.fold` |
| Summary | `summary(out)` | `summary(fit)` |
| Coefficients | `out$beta.samples` | `coef(fit)` or `fit$beta.samples` |
| Predictions | `predict(out, X.0)` | `predict(fit, X.0 = X.0)` |
| Traceplots | `plot(out$beta.samples)` | Not needed (deterministic) |
| Rhat / ESS | `out$rhat`, `out$ESS` | Not needed (deterministic) |

### Diagnostics you can drop

Because INLA is deterministic, the entire family of MCMC convergence
diagnostics becomes irrelevant:

- **Traceplots** — these diagnose chain mixing. No chains, no mixing to
  diagnose.
- **Rhat** (potential scale reduction factor) — this detects
  between-chain disagreement. A single deterministic run cannot disagree
  with itself.
- **Effective sample size (ESS)** — this measures how many independent
  samples the chain produced. INLA does not sample.
- **Multiple chain comparison** — running three chains and checking
  agreement is standard MCMC practice. With INLA, a single run is
  sufficient.

### Diagnostics you gain

INLAocc provides simulation-based diagnostics that go beyond what
spOccupancy offers out of the box:

``` r

# One-call diagnostic panel (residuals, QQ, spatial, temporal)
checkModel(fit)

# Individual simulation-based tests
testUniformity(fit)       # PIT uniformity (randomized quantile residuals)
testDispersion(fit)       # over/underdispersion
testOutliers(fit)         # outlier detection
testZeroInflation(fit)    # excess zeros beyond what the model predicts

# Spatial residual diagnostics
moranI(fit)               # Moran's I on residuals
variogram(fit)            # empirical variogram of residuals

# Temporal residual diagnostics
durbinWatson(fit)         # Durbin-Watson test for temporal autocorrelation

# DHARMa integration
dharma(fit)               # returns a DHARMa object for custom diagnostics

# Parameter identifiability
checkIdentifiability(fit) # flags weakly identified parameters
```

These tests are particularly useful for model criticism — detecting
whether the fitted model captures the data-generating process — which is
orthogonal to whether the estimation algorithm converged. See
[`vignette("diagnostics")`](https://gillescolling.com/INLAocc/articles/diagnostics.md)
for worked examples.

## What you gain — detailed

**Speed.** The table below gives rough wall-clock times for common model
types on a standard laptop (4-core, 16 GB RAM). Exact times depend on
data dimensions and model complexity.

| Model type            | Sites       | spOccupancy | INLAocc    | Speedup |
|:----------------------|:------------|:------------|:-----------|:--------|
| Basic (`PGOcc`)       | 300         | ~30 sec     | ~3 sec     | 10x     |
| Spatial (`spPGOcc`)   | 300         | ~3–5 min    | ~15–30 sec | 8–12x   |
| Spatial (`spPGOcc`)   | 1000        | ~20–30 min  | ~1–2 min   | 15–20x  |
| Multi-species spatial | 300, 20 spp | ~30–60 min  | ~3–5 min   | 10x     |

Speed advantages compound when fitting many models (e.g., model
selection across covariate sets, cross-validation folds, or
multi-species analyses where each species is a separate fit).

**No MCMC tuning.** Removing `n.batch`, `batch.length`, `n.burn`,
`n.thin`, `n.chains`, `n.neighbors`, and `cov.model` from the user’s
decision space eliminates a common source of analytical variability. Two
researchers fitting the same model to the same data with different MCMC
settings can get different posteriors; two researchers using INLAocc
cannot.

**Deterministic convergence.** Run the same
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) call
twice and get identical results. This simplifies debugging, makes
analyses reproducible without seed management, and removes “did my
chains converge?” from the review checklist.

**Richer diagnostics.** The simulation-based tests (`testUniformity`,
`testDispersion`, etc.) and spatial residual checks (`moranI`,
`variogram`) provide model criticism tools that spOccupancy does not
include. These complement posterior predictive checks by testing
specific aspects of model fit rather than a single omnibus statistic.

## What to watch for — honest caveats

### Posterior approximation quality

INLA approximates the posterior; MCMC (given convergence) samples from
it exactly. For most occupancy models — moderate number of sites,
reasonable detection probability, a handful of covariates — the
approximation is excellent. Benchmarks on simulated data show
correlations above 0.99 between INLA and MCMC posterior means, and
credible interval coverage within 1–2% of nominal.

Where the approximation can degrade:

- **Very small datasets** (\< 30 sites). The Gaussian approximation to
  the posterior becomes less accurate when the likelihood is weak.
- **Weakly identified models.** Many covariates relative to sites, or
  very low detection probability, can produce posteriors that are
  substantially non-Gaussian. MCMC handles this naturally; INLA may
  underestimate uncertainty.
- **Heavy posterior skewness.** If the scientific question depends on
  the tail shape of a posterior (e.g., probability that an effect
  exceeds a threshold), MCMC’s exact samples may be more reliable than
  INLA’s Gaussian approximation.

### Gibbs data augmentation

INLAocc’s post-EM bias correction uses a short Gibbs data-augmentation
chain (default 100 iterations). This corrects the attenuation bias that
arises from treating the EM point estimates of latent occupancy states
as fixed during INLA fitting. The correction is fast (negligible runtime
overhead) and accurate for fixed-effect coefficients.

For variance components estimated from few groups (e.g., a random
intercept with 5 levels), the short Gibbs chain may not fully propagate
uncertainty. If variance component inference is the primary goal and the
number of groups is small (\< 10), full MCMC via spOccupancy is the
safer choice.

### Prior specification

The two packages use fundamentally different prior frameworks:

- **spOccupancy** uses conjugate priors: normal priors on regression
  coefficients, inverse-Gamma priors on variance components. These can
  be made informative or vague by tuning hyperparameters.
- **INLAocc** uses PC priors (penalizing complexity) for spatial range,
  spatial variance, and random effect standard deviations. PC priors
  shrink toward the base model (no spatial effect, no random variation),
  which is statistically conservative.

For moderate-to-large datasets, the prior choice has negligible impact
on posteriors — the likelihood dominates. For small datasets, priors
matter, and the two frameworks can give meaningfully different results.
If you have strong prior information expressed as inverse-Gamma
hyperparameters, translating to PC prior parameters requires care.

### Custom derived quantities

spOccupancy gives you full MCMC samples from the joint posterior. You
can compute any derived quantity (species richness at a site,
probability of co-occurrence, etc.) by applying a function to each
posterior draw.

INLAocc provides posterior samples via the `$` accessors, but these come
from INLA’s Gaussian approximation to the marginals, not from exact
joint posterior sampling. For simple derived quantities (posterior
means, credible intervals, posterior probabilities of direction), the
difference is negligible. For complex nonlinear functions of multiple
parameters, MCMC samples from the exact joint posterior may be more
appropriate.

## Migration checklist

A step-by-step recipe for converting an existing spOccupancy script:

1.  **Install INLAocc and INLA.** INLAocc requires the INLA package,
    which is not on CRAN. See
    [`vignette("quickstart")`](https://gillescolling.com/INLAocc/articles/quickstart.md)
    for installation instructions.

2.  **Load your existing data.** spOccupancy data lists (with `y`,
    `occ.covs`, `det.covs`, `coords`) work directly. No reformatting.

3.  **Replace the model call.** Swap the spOccupancy function and its
    MCMC arguments for a single
    [`occu()`](https://gillescolling.com/INLAocc/reference/occu.md)
    call:

    ``` r

    # Before
    out <- spPGOcc(occ.formula = ~ elev, det.formula = ~ effort,
                   data = data, n.batch = 400, batch.length = 25,
                   n.burn = 5000, n.thin = 5, n.chains = 3,
                   NNGP = TRUE, n.neighbors = 15,
                   cov.model = "exponential")

    # After
    fit <- occu(~ elev, ~ effort, data = data,
                spatial = data$coords, verbose = 0)
    ```

4.  **Replace `summary(out)` with `summary(fit)`.** Same information
    (coefficient estimates, credible intervals), different layout.

5.  **Replace diagnostic functions.**
    [`ppcOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/ppcOcc.html)
    becomes
    [`ppcOccu()`](https://gillescolling.com/INLAocc/reference/ppcOccu.md),
    [`waicOcc()`](https://www.doserlab.com/files/spoccupancy-web/reference/waicOcc.html)
    becomes
    [`waicOccu()`](https://gillescolling.com/INLAocc/reference/waicOccu.md).

6.  **Drop MCMC diagnostics.** Remove traceplot inspection, Rhat checks,
    and ESS checks from your script. These are not needed.

7.  **Add INLAocc diagnostics.** Replace MCMC convergence checks with
    model criticism:

    ``` r

    checkModel(fit)
    testUniformity(fit)
    ```

8.  **Check downstream code.** If your script uses `out$beta.samples` or
    `out$alpha.samples` for custom calculations, replace `out` with
    `fit` — the `$` accessor returns matrices in the same format.
