# Model Identifiability

## Introduction

Occupancy models have a fundamental identifiability challenge: a site
with no detections could be truly unoccupied, or it could be occupied
with low detection probability. Separating occupancy from detection
requires enough visits per site and enough variation in
detection/non-detection patterns across sites and visits. Without these,
the model cannot determine which process — occupancy or detection — is
responsible for the observed data, and parameter estimates become
unreliable.

When models are not identifiable, parameters have flat posteriors (the
data cannot distinguish between different parameter values), estimates
sit on the boundary of the parameter space, or small changes in the data
produce large changes in estimates. These symptoms are not always
obvious. A model can converge, return point estimates, and produce
confidence intervals that look reasonable — while the estimates are
driven almost entirely by the prior or by numerical regularization
rather than by the data.

INLAocc provides
[`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
for automated screening of common identifiability issues, both before
and after fitting a model. This vignette explains the conceptual
framework behind occupancy identifiability, demonstrates the built-in
checks, and presents a simulation-based stress-testing approach for
deeper assessment. For model fitting basics, see
[`vignette("quickstart")`](https://gillescolling.com/INLAocc/articles/quickstart.md);
for post-fit diagnostics, see
[`vignette("diagnostics")`](https://gillescolling.com/INLAocc/articles/diagnostics.md).

## The identifiability problem

### The occupancy likelihood

Consider a site \\i\\ surveyed on \\J\\ visits. The probability of the
observed detection history depends on two unknowns: the occupancy
probability \\\psi_i\\ and the detection probabilities \\p\_{ij}\\. For
a site where the species was never detected, the likelihood contribution
is:

\\L_i = (1 - \psi_i) + \psi_i \prod\_{j=1}^{J} (1 - p\_{ij})\\

The first term is the probability that the site is truly unoccupied. The
second is the probability that the site is occupied but the species went
undetected on every visit. The data cannot distinguish between these two
contributions — both produce the same observation (all zeros).

### Why one visit is not enough

With \\J = 1\\ visit, the likelihood for an undetected site simplifies
to:

\\L_i = (1 - \psi_i) + \psi_i (1 - p_i) = 1 - \psi_i p_i\\

Only the product \\\psi_i p_i\\ is identifiable — not \\\psi_i\\ and
\\p_i\\ individually. Any combination of occupancy and detection that
produces the same product fits the data equally well. A species with 80%
occupancy and 50% detection looks identical to one with 40% occupancy
and 100% detection.

### Repeated visits break the degeneracy

With \\J \geq 2\\ independent visits, the pattern of detections and
non-detections across visits provides information to separate the two
processes. Sites where the species is detected on some visits but not
others are particularly informative: they demonstrate that detection is
imperfect, which allows the model to estimate how much non-detection is
due to absence versus missed detections.

As a rule of thumb: \\J \geq 3\\ visits per site, with detection
probability \\p \> 0.15\\, gives reasonable identifiability for simple
models. More complex models (with many covariates or random effects)
require more data.

## Using checkIdentifiability()

INLAocc’s
[`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
function screens for common identifiability problems. It works in two
modes: pre-fit (data only) and post-fit (fitted model), catching
progressively more issues.

### Pre-fit checks (data only)

Pre-fit checks examine the detection history matrix and covariates
without fitting any model. They are fast and should be run before any
analysis.

``` r

library(INLAocc)

sim <- simulate_occu(N = 200, J = 4,
                     beta_occ = c(0.5, -0.8),
                     beta_det = c(0, -0.5),
                     n_occ_covs = 1, n_det_covs = 1,
                     seed = 42)

issues <- checkIdentifiability(sim$data)
issues
#> No identifiability issues detected.
```

For well-identified data, the output should confirm no issues. The
summary statistics provide a quick sanity check: naive occupancy
(proportion of sites with at least one detection), naive detection
(proportion of visits with detections at occupied sites), visits per
site, and the number of unique detection histories.

### Problematic data

Now simulate data with known identifiability issues — few sites, few
visits, and very low detection:

``` r

sim_bad <- simulate_occu(N = 50, J = 2,
                         beta_occ = c(0.5, -0.8),
                         beta_det = c(-2, -0.5),  # very low detection
                         n_occ_covs = 1, n_det_covs = 1,
                         seed = 42)

issues <- checkIdentifiability(sim_bad$data)
issues
#> Identifiability check: 4 issue(s) found (3 HIGH, 1 MEDIUM)
#> 
#> [!!] Few visits (J <= 2)
#>     J = 2. With 1-2 visits per site, occupancy and detection are confounded.
#>     -> Increase survey effort to >= 3 visits per site.
#> 
#> [!!] Very low naive detection rate
#>     Naive p = 0.090. Below 10%, the model struggles to separate absence from non-detection.
#>     -> Consider pooling visits, using a simpler detection model, or increasing survey effort.
#> 
#> [!!] Most sites have no detections
#>     42 of 50 sites (84%) have zero detections.
#>     -> Very sparse data. Consider whether the species is too rare for occupancy modelling at this scale.
#> 
#> [!] Low detection history diversity
#>     Only 4 unique detection histories across 50 sites.
#>     -> Many sites have identical patterns. The model may not have enough variation to estimate parameters.
```

Each issue has a concrete interpretation:

- **Visits \\\leq\\ 2**: With only two visits, the detection/occupancy
  trade-off is severely underdetermined. The model has very little
  information to separate the two processes, especially when covariates
  are included.
- **Low detection**: When species are rarely detected (\\p \< 0.10\\),
  almost all sites have all-zero histories. The model cannot reliably
  distinguish low occupancy from low detection because both produce the
  same data pattern.
- **Over-parameterization**: Five parameters for 50 sites is aggressive.
  The rule of thumb is at least 10 sites per parameter (20 is better).
  With too many parameters relative to data, individual coefficients are
  poorly determined.

### Post-fit checks

When passed a fitted model,
[`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
runs additional checks that require posterior estimates:

``` r

fit <- occu(~ occ_x1, ~ det_x1, data = sim_bad$data, verbose = 0)
issues_post <- checkIdentifiability(fit)
issues_post
#> Identifiability check: 5 issue(s) found (4 HIGH, 1 MEDIUM)
#> 
#> [!!] Few visits (J <= 2)
#>     J = 2. With 1-2 visits per site, occupancy and detection are confounded.
#>     -> Increase survey effort to >= 3 visits per site.
#> 
#> [!!] Very low naive detection rate
#>     Naive p = 0.090. Below 10%, the model struggles to separate absence from non-detection.
#>     -> Consider pooling visits, using a simpler detection model, or increasing survey effort.
#> 
#> [!!] Most sites have no detections
#>     42 of 50 sites (84%) have zero detections.
#>     -> Very sparse data. Consider whether the species is too rare for occupancy modelling at this scale.
#> 
#> [!!] EM did not converge
#>     Ran 50 iterations without converging.
#>     -> Try increasing max.iter, adjusting damping, or simplifying the model.
#> 
#> [!] Low detection history diversity
#>     Only 4 unique detection histories across 50 sites.
#>     -> Many sites have identical patterns. The model may not have enough variation to estimate parameters.
```

Post-fit checks also flag:

- **Non-convergence**: The EM algorithm did not reach the convergence
  criterion within the maximum number of iterations.
- **Large posterior SDs**: A posterior standard deviation more than 10
  times the absolute posterior mean indicates the data provides almost
  no information about that parameter.
- **Boundary estimates**: Occupancy or detection estimates at 0 or 1 for
  a large fraction of sites suggest the model is hitting the boundary of
  the parameter space.
- **High correlation between processes**: Strong posterior correlation
  between occupancy and detection parameters indicates the model cannot
  separate the two processes.

## Covariate confounding

The most common identifiability issue in applied work is including the
same covariate in both the occupancy and detection submodels. This is
tempting when a variable (e.g., habitat quality, elevation) plausibly
affects both processes, but it creates a fundamental identifiability
problem.

``` r

# Share a covariate between both processes
sim$data$det.covs$occ_x1 <- matrix(rep(sim$data$occ.covs$occ_x1, 4), nrow = 200, ncol = 4)

# Same variable in both submodels
fit_confound <- occu(~ occ_x1, ~ occ_x1, data = sim$data, verbose = 0)
checkIdentifiability(fit_confound)
#> Identifiability check: 1 issue(s) found (1 HIGH)
#> 
#> [!!] Same covariate in occupancy and detection
#>     Shared covariates: occ_x1
#>     -> Using the same variable for both processes causes confounding. Remove from one process or use different transformations.
```

Why this is problematic: if a covariate affects both occupancy and
detection, the model cannot determine which process it operates through.
The estimated effect gets split arbitrarily between the two submodels,
and the split depends on the data realization, the starting values, and
the prior — not on the underlying biology.

**Fixes:**

- Include the covariate in the process where it has the strongest
  theoretical justification. Elevation affects where a species lives
  (occupancy) more directly than whether a surveyor detects it
  (detection).
- If there is genuine reason to include a variable in both, use a
  different functional form in each: linear in one process, quadratic in
  the other, or a categorical version in one and continuous in the
  other.
- Random effects with the same grouping variable in both processes cause
  the same problem.
  [`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
  flags this as well.

## Simulation-based stress testing

Automated checks catch common problems, but they cannot cover every
scenario. Simulation-based stress testing provides a deeper assessment:
simulate data from known parameters, fit the model, and check whether
the true values are recovered. If recovery fails under conditions
matching your real data, the model is not identifiable for your
application.

### Basic parameter recovery

The simplest stress test simulates multiple datasets from the same
parameters and checks whether the estimates are centered on the truth:

``` r

n_sims <- 20
results <- data.frame()

for (i in 1:n_sims) {
  sim_i <- simulate_occu(N = 200, J = 4,
                         beta_occ = c(0.5, -0.8),
                         beta_det = c(0, -0.5),
                         n_occ_covs = 1, n_det_covs = 1,
                         seed = i)

  fit_i <- occu(~ occ_x1, ~ det_x1, data = sim_i$data, verbose = 0)
  cf <- coef(fit_i)

  results <- rbind(results, data.frame(
    sim = i,
    beta0_occ = cf$occ$mean[1],
    beta1_occ = cf$occ$mean[2],
    alpha0_det = cf$det$mean[1],
    alpha1_det = cf$det$mean[2]
  ))
}

# Estimates should be centered on truth
colMeans(results[, -1])
```

Good recovery: the means across simulations should be close to the true
values (0.5, -0.8, 0, -0.5). Substantial bias would indicate an
identifiability problem under these conditions.

You can also check coverage — the proportion of simulations where the
95% credible interval contains the true value. Nominal coverage should
be close to 0.95.

### Varying sample size

A natural question is: how many sites do I need? Simulate across a range
of sample sizes to find the threshold where estimation breaks down:

``` r

for (N in c(50, 100, 200, 500)) {
  sim_i <- simulate_occu(N = N, J = 4,
                         beta_occ = c(0.5, -0.8),
                         beta_det = c(0, -0.5),
                         n_occ_covs = 1, n_det_covs = 1,
                         seed = 42)

  fit_i <- occu(~ occ_x1, ~ det_x1, data = sim_i$data, verbose = 0)
  est <- coef(fit_i)$occ$mean[2]
  se <- coef(fit_i)$occ$sd[2]

  cat(sprintf("N = %3d: beta1 = %6.2f (SE = %.2f, true = -0.80)\n",
              N, est, se))
}
```

With few sites the estimate is biased and the standard error is large.
By 200 sites, recovery is good. This tells you the minimum sample size
for reliable estimation under your specific model structure.

### Varying detection probability

Detection probability is the single most important factor for occupancy
model identifiability. Test how your model behaves across a range of
detection levels:

``` r

for (det_int in c(-2, -1, 0, 1)) {
  sim_i <- simulate_occu(N = 200, J = 4,
                         beta_occ = c(0.5, -0.8),
                         beta_det = c(det_int, -0.5),
                         n_occ_covs = 1, n_det_covs = 1,
                         seed = 42)

  fit_i <- occu(~ occ_x1, ~ det_x1, data = sim_i$data, verbose = 0)
  naive_p <- mean(sim_i$data$y, na.rm = TRUE)

  cat(sprintf("det_int = %2d (naive p = %.2f): beta1 = %6.2f (true = -0.80)\n",
              det_int, naive_p, coef(fit_i)$occ$mean[2]))
}
```

The pattern should be clear: low detection produces substantial bias in
the occupancy coefficient. Above \\p \approx 0.20\\, recovery is
adequate. This reveals the detection threshold below which your model
breaks down.

### Varying number of visits

The same approach works for the number of visits:

``` r

for (J in c(2, 3, 4, 6, 10)) {
  sim_i <- simulate_occu(N = 200, J = J,
                         beta_occ = c(0.5, -0.8),
                         beta_det = c(0, -0.5),
                         n_occ_covs = 1, n_det_covs = 1,
                         seed = 42)

  fit_i <- occu(~ occ_x1, ~ det_x1, data = sim_i$data, verbose = 0)
  est <- coef(fit_i)$occ$mean[2]
  se <- coef(fit_i)$occ$sd[2]

  cat(sprintf("J = %2d: beta1 = %6.2f (SE = %.2f)\n", J, est, se))
}
```

Diminishing returns set in quickly: going from 2 to 4 visits is far more
valuable than going from 4 to 10.

## When identifiability fails: symptoms and fixes

The table below summarizes common symptoms of identifiability failure,
their likely causes, and practical fixes.

| Symptom | Likely cause | Fix |
|:---|:---|:---|
| Posterior SD \>\> posterior mean | Parameter unidentifiable | Reduce model complexity |
| Coefficients flip sign across runs | Flat likelihood surface | More data or fewer parameters |
| EM fails to converge | Likelihood too flat | Increase visits, simplify model |
| \\\hat{\psi} \approx 1\\ for all sites | Detection absorbing all variation | Check detection model specification |
| \\\hat{p} \approx 0\\ for all sites | Occupancy absorbing all variation | Ensure \\J \geq 3\\ |
| Same covariate in both processes | Confounding | Move to one process |

A useful diagnostic when you suspect identifiability problems is to fit
progressively simpler models: drop covariates from detection, then from
occupancy, until the remaining parameters are stable. The point where
parameters start shifting substantially as you add complexity marks the
boundary of what your data can support.

## Design recommendations

Based on identifiability theory and practical experience:

- **Minimum 3 visits per site**, ideally 4–5. With \\J = 2\\, only the
  simplest models (intercept-only or single covariate) are identifiable.
  More visits always help, but with diminishing returns beyond \\J =
  5\\.
- **Detection probability above 0.15.** Below this threshold, occupancy
  and detection are hard to separate without informative priors. If your
  study species has very low detection, consider increasing survey
  effort per visit (longer transects, more traps) rather than adding
  more visits.
- **At least 10 sites per parameter** in the occupancy model; 20 is
  better. Count intercepts and random effect variance components as
  parameters.
- **Avoid the same covariate in both processes** unless you have strong
  theoretical justification and adequate data. If you must include a
  shared variable, use different functional forms (linear vs. quadratic)
  or check that the simulation recovery test above succeeds with your
  specific design.
- **Scale covariates.** Unscaled covariates with wildly different
  magnitudes create numerical identifiability issues even when the model
  is theoretically identifiable. Center and scale continuous covariates
  to unit variance.
- **Start simple.** Fit the simplest defensible model first
  (intercept-only for detection, single covariate for occupancy). Add
  complexity only when diagnostics or WAIC comparisons justify it. Each
  added parameter requires the data to support one more degree of
  freedom.

## Interaction with priors

INLA uses penalized complexity (PC) priors that provide implicit
regularization. This regularization can mask identifiability issues in a
way that is easy to miss.

A weakly identified parameter with a strong prior will produce a
posterior that looks well-behaved — finite mean, reasonable standard
deviation, no boundary issues. But the posterior is driven by the prior,
not by the data. The model appears to work, but the estimates reflect
your assumptions rather than your observations.

To check whether the data is actually informing a parameter: compare the
posterior to the prior. If they are nearly identical in location and
spread, the data provides little information about that parameter.
INLAocc’s
[`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
post-fit mode flags cases where the posterior standard deviation is
close to the prior standard deviation.

This concern is less acute with INLAocc’s default weakly informative
priors, which are intentionally vague. But it becomes important when
using informative priors to “fix” convergence issues. If a model only
converges with a tight prior on a parameter, that parameter is not
identified by the data — the prior is doing the work. In that case, the
correct response is to simplify the model, not to tighten the prior.

A practical test: fit the model twice with different priors (e.g., prior
SD = 1 vs. prior SD = 5 on the same coefficient). If the posterior
changes substantially, the data is not constraining that parameter.

## Summary

Identifiability is not a binary property — it exists on a continuum from
well-identified (data strongly constrains parameters) to completely
unidentifiable (data provides no information). Most real datasets fall
somewhere in between, and the goal is to match model complexity to the
information content of the data.

The recommended workflow:

1.  **Before fitting**: run
    [`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
    on the data to catch obvious issues (too few visits, low detection,
    covariate confounding).
2.  **After fitting**: run
    [`checkIdentifiability()`](https://gillescolling.com/INLAocc/reference/checkIdentifiability.md)
    on the fitted model to check for large posterior SDs, boundary
    estimates, and non-convergence.
3.  **For critical analyses**: run the simulation-based stress test with
    parameters and sample sizes matching your real data to verify that
    your specific model structure is identifiable under your specific
    conditions.

If identifiability checks reveal problems, simplify the model rather
than adding stronger priors. A simple model that is well-identified
produces more reliable inference than a complex model that is not.
