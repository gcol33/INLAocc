# Integrated Occupancy Models

## Introduction

Many species are monitored simultaneously by multiple programs that
differ in protocol, effort, and spatial coverage. A state wildlife
agency may run structured point-count surveys at a fixed set of sites,
while a network of camera traps operates year-round at a partially
overlapping set of locations, and citizen-science platforms like eBird
accumulate opportunistic checklists across the broader landscape. Each
data source captures the same underlying quantity — whether the species
occupies a site — but filters it through a different detection process.
Point counts detect vocalisations during brief windows; cameras detect
physical presence over extended periods; eBird checklists reflect
variable observer skill and effort. Analysing any single source in
isolation discards the information contained in the others.

Integrated occupancy models solve this problem by combining all sources
into a single inference. The true occupancy state \\z_i\\ at site \\i\\
is shared across every data source. What differs is how likely each
source is to detect the species when it is present. By modelling each
source’s detection process separately while linking them through a
common occupancy layer, the model borrows strength across programs.
Sites surveyed by multiple sources are particularly informative, because
detections (or non-detections) from one source constrain what the other
source’s data imply about occupancy.

INLAocc fits integrated occupancy models through the `integrated = TRUE`
argument in
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md). The EM
algorithm iterates over a shared occupancy posterior while fitting
source-specific detection models via INLA, exactly paralleling the
single-source workflow. This vignette covers the model formulation, the
data structure required, simulation, model fitting, and extensions to
spatial and multi-species integrated models. It assumes familiarity with
single-species models; see
[`vignette("quickstart")`](https://gillescolling.com/INLAocc/articles/quickstart.md)
for an introduction to
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) and the
INLAocc data format.

## The integrated occupancy model

For \\D\\ data sources observing the same set of sites, the integrated
occupancy model has a shared occupancy layer and source-specific
detection layers.

**Occupancy.** At each of \\N\\ unique sites, the latent occupancy state
is:

\\z_i \sim \text{Bernoulli}(\psi_i), \quad \text{logit}(\psi_i) =
\mathbf{x}\_i^\top \boldsymbol{\beta}\\

where \\\mathbf{x}\_i\\ is a vector of site-level covariates and
\\\boldsymbol{\beta}\\ the occupancy coefficients. This is identical to
the standard single-source model — every source shares the same
\\\psi_i\\.

**Detection.** Each source \\d = 1, \ldots, D\\ has its own detection
process. Source \\d\\ surveys a subset of sites (not necessarily all
\\N\\), and at each surveyed site \\i\\, it makes \\J_d\\ visits. The
observation model is:

\\y\_{ijd} \mid z_i \sim \text{Bernoulli}(z_i \cdot p\_{ijd}), \quad
\text{logit}(p\_{ijd}) = \mathbf{w}\_{ijd}^\top \boldsymbol{\alpha}\_d\\

where \\\mathbf{w}\_{ijd}\\ are visit-level covariates specific to
source \\d\\, and \\\boldsymbol{\alpha}\_d\\ is a source-specific
detection coefficient vector. Each source has its own intercept, its own
covariates, and its own regression coefficients — reflecting the fact
that, say, a camera trap and a point count have fundamentally different
detection mechanisms.

**E-step.** The EM algorithm updates the posterior probability that each
site is occupied, conditional on all observed data. For a site surveyed
by sources \\\mathcal{D}\_i \subseteq \\1, \ldots, D\\\\:

\\P(z_i = 1 \mid \mathbf{y}) = \frac{\psi_i \prod\_{d \in
\mathcal{D}\_i} \prod\_{j=1}^{J_d} p\_{ijd}^{y\_{ijd}} (1 -
p\_{ijd})^{1 - y\_{ijd}}}{\psi_i \prod\_{d \in \mathcal{D}\_i}
\prod\_{j=1}^{J_d} p\_{ijd}^{y\_{ijd}} (1 - p\_{ijd})^{1 - y\_{ijd}} +
(1 - \psi_i) \\ \mathbb{I}\[\mathbf{y}\_{i \cdot} = \mathbf{0}\]}\\

The product runs over all sources that surveyed site \\i\\ and all
visits within each source. For sites with at least one detection from
any source, the posterior is 1 (the species is present). For sites with
zero detections across all sources, the posterior depends on \\\psi_i\\
and the detection probabilities from every source. Each additional
source of non-detection further lowers the posterior occupancy
probability, giving more precise estimates of true absence than any
single source alone.

## Data structure

Integrated models require data in a specific list format. The key
difference from standard models is that detection data (`y`) and
detection covariates (`det.covs`) are lists of length \\D\\, one element
per source, and a `sites` mapping tells the model which rows of the
shared occupancy covariates correspond to each source’s sites.

``` r

data <- list(
  y = list(
    point_counts = matrix(..., nrow = 80, ncol = 4),    # 80 sites, 4 visits
    cameras      = matrix(..., nrow = 60, ncol = 10)     # 60 sites, 10 occasions
  ),
  occ.covs = data.frame(
    elevation = rnorm(110),
    forest    = runif(110)
  ),
  det.covs = list(
    list(date = matrix(rnorm(80 * 4), 80, 4)),           # source 1 covariates
    list(effort = matrix(rnorm(60 * 10), 60, 10))        # source 2 covariates
  ),
  sites = list(1:80, 51:110),
  coords = matrix(runif(110 * 2), ncol = 2)              # all unique sites
)
```

Each component plays a specific role:

- **`y`**: a list of \\D\\ detection matrices. Source \\d\\ has an \\N_d
  \times J_d\\ matrix where rows are sites and columns are visits.
  Sources can have different numbers of sites and visits — the only
  requirement is that each matrix has rows matching the corresponding
  entry in `sites`.

- **`occ.covs`**: a data.frame with one row per unique site (\\N = 110\\
  in this example). All occupancy covariates live here, shared across
  sources.

- **`det.covs`**: a list of \\D\\ lists. Each inner list contains named
  matrices of visit-level detection covariates for that source. Source 1
  has a `date` matrix with 80 rows and 4 columns; source 2 has an
  `effort` matrix with 60 rows and 10 columns.

- **`sites`**: a list of \\D\\ integer vectors. Each vector maps the
  rows of the corresponding `y` matrix to rows of `occ.covs`. In this
  example, source 1 covers `occ.covs` rows 1–80 and source 2 covers rows
  51–110. Sites 51–80 appear in both vectors — these 30 overlap sites
  are surveyed by both programs.

- **`coords`**: an \\N \times 2\\ matrix of coordinates for all unique
  sites. Required only for spatial models.

The overlap sites are particularly valuable. They are the only sites
where the model directly observes the same occupancy state through two
different detection lenses, which is what allows it to calibrate the
detection parameters across sources.

## Simulating integrated data

INLAocc provides
[`simIntOcc()`](https://gillescolling.com/INLAocc/reference/simIntOcc.md)
to generate integrated occupancy datasets with known parameters, which
is essential for validating model performance and understanding the data
structure.

``` r

library(INLAocc)

sim <- simIntOcc(
  N_total  = 150,
  n_data   = 2,
  J        = c(4, 3),
  n_shared = 20,
  beta_occ = c(0.5, -0.8),
  beta_det = list(c(0, -0.5), c(-0.3, 0.2)),
  seed     = 42
)

str(sim$data)
#> List of 5
#>  $ y       :List of 2
#>   ..$ : int [1:85, 1:4] 0 1 1 1 0 1 0 1 0 0 ...
#>   ..$ : int [1:85, 1:3] 0 1 1 1 0 0 0 0 1 0 ...
#>  $ occ.covs:'data.frame':    150 obs. of  1 variable:
#>   ..$ occ_x1: num [1:150] -0.0407 -1.5515 1.1672 -0.2736 -0.4678 ...
#>  $ det.covs:List of 2
#>   ..$ :List of 1
#>   .. ..$ det_x1: num [1:85, 1:4] 0.931 1.3349 -0.8693 0.0555 0.0491 ...
#>   ..$ :List of 1
#>   .. ..$ det_x1: num [1:85, 1:3] 0.3703 0.2791 0.2019 -0.013 -0.0909 ...
#>  $ sites   :List of 2
#>   ..$ : int [1:85] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ : int [1:85] 7 8 11 17 20 32 37 39 46 54 ...
#>  $ coords  : num [1:150, 1:2] 0.915 0.937 0.286 0.83 0.642 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "x" "y"
```

The simulator creates a dataset with 150 unique sites. Source 1 surveys
85 sites with 4 visits each, using occupancy coefficients
\\\boldsymbol{\beta} = (0.5, -0.8)\\ and detection coefficients
\\\boldsymbol{\alpha}\_1 = (0, -0.5)\\. Source 2 surveys a different set
of 85 sites with 3 visits, with detection coefficients
\\\boldsymbol{\alpha}\_2 = (-0.3, 0.2)\\. Twenty sites appear in both
sources, providing the overlap needed to calibrate detection across
programs. The occupancy formula is shared — both sources inform the same
\\\psi_i\\.

## Fitting the integrated model

Fitting an integrated model uses the same
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md)
interface, with `integrated = TRUE`:

``` r

fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, integrated = TRUE, verbose = 0)
summary(fit)
#> === Integrated Occupancy Model (INLA-Laplace) ===
#> 
#> Data sources: 2 | Total sites: 150
#> EM iterations: 15 | Converged: TRUE
#> 
#>   Source 1: 85 sites, 4 max visits
#>   Source 2: 85 sites, 3 max visits
#> 
#> --- Occupancy (psi, shared) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept)  1.5932 0.0076     1.5783     1.6082
#> occ_x1      -0.9398 0.0079    -0.9554    -0.9242
#> 
#> --- Detection (p, source 1) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept) -0.6409 0.1167    -0.8696    -0.4121
#> det_x1      -0.3830 0.1196    -0.6174    -0.1486
#> 
#> --- Detection (p, source 2) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept) -0.8902 0.1381    -1.1609    -0.6195
#> det_x1       0.1105 0.1291    -0.1425     0.3635
#> 
#> Estimated occupancy: 0.799 (0.420 - 0.963)
#> Estimated occupied sites: 119.9 / 150
```

Several things to note in the output:

- The occupancy estimates integrate information from both sources. The
  estimated coefficients should be close to the true values (0.5 and
  -0.8), reflecting the combined power of 150 sites observed through two
  detection processes.

- Source 1 and source 2 should have visibly different detection
  characteristics, matching the simulation parameters.

- The 20 overlap sites are doing the heavy lifting for detection
  calibration. At those sites, the model sees the same \\z_i\\ filtered
  through two different detection lenses, which pins down the relative
  detection rates.

## Source-specific detection formulas

When a single detection formula is provided, it is applied to all
sources. But sources often have different covariates — point counts have
observer ID and wind speed, while cameras have trigger sensitivity and
battery voltage. In INLAocc, you can pass a list of formulas, one per
source:

``` r

fit2 <- occu(
  ~ occ_x1,
  list(~ det_x1, ~ det_x1),
  data = sim$data,
  integrated = TRUE,
  verbose = 0
)
summary(fit2)
```

The first formula applies to source 1, the second to source 2. The
occupancy formula is always shared — there is no mechanism for
source-specific occupancy effects, because occupancy is a property of
the site, not of the observation method.

This is the recommended approach when sources genuinely differ in what
drives detection. Start with a shared detection formula to establish a
baseline, then try source-specific formulas and compare via WAIC (see
below). If the source-specific model has lower WAIC despite the
additional parameters, the detection processes are sufficiently
different to warrant separate treatment.

## Spatial integrated models

Integrated models extend naturally to spatial settings. Adding a spatial
random field to the occupancy layer accounts for unmeasured
environmental variation that is spatially structured:

``` r

fit_sp <- occu(
  ~ occ_x1, ~ det_x1,
  data    = sim$data,
  integrated = TRUE,
  spatial = sim$data$coords,
  verbose = 0
)
summary(fit_sp)
#> === Integrated Occupancy Model (INLA-Laplace) ===
#> 
#> Data sources: 2 | Total sites: 150
#> EM iterations: 15 | Converged: TRUE
#> 
#>   Source 1: 85 sites, 4 max visits
#>   Source 2: 85 sites, 3 max visits
#> 
#> --- Occupancy (psi, shared) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept)  1.5932 0.0076     1.5783     1.6082
#> occ_x1      -0.9398 0.0079    -0.9554    -0.9242
#> 
#> --- Detection (p, source 1) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept) -0.6409 0.1167    -0.8696    -0.4121
#> det_x1      -0.3830 0.1196    -0.6174    -0.1486
#> 
#> --- Detection (p, source 2) ---
#>                mean     sd 0.025quant 0.975quant
#> (Intercept) -0.8902 0.1381    -1.1609    -0.6195
#> det_x1       0.1105 0.1291    -0.1425     0.3635
#> 
#> Estimated occupancy: 0.799 (0.420 - 0.963)
#> Estimated occupied sites: 119.9 / 150
```

The spatial field is shared across all data sources — it captures
spatial variation in occupancy that is common regardless of how the
species is detected. This is conceptually correct: spatial
autocorrelation in occupancy arises from environmental processes
(habitat connectivity, dispersal limitation, unmeasured gradients) that
operate independently of the observation method.

The spatial integrated model works best when sources have complementary
spatial coverage. Structured surveys may cover a dense grid in a core
area, while citizen-science data fill in the periphery. The shared
spatial field allows information to flow between regions dominated by
different sources, smoothing the occupancy surface across the entire
study area.

For details on mesh construction, PC priors for the spatial range and
variance, and occupancy mapping, see
[`vignette("spatial-models")`](https://gillescolling.com/INLAocc/articles/spatial-models.md).

## Multi-species integrated models

Community-level integrated models combine multiple species and multiple
data sources under a single hierarchical framework. Species-level
occupancy and detection parameters are drawn from community
distributions, and each source retains its own detection process.

``` r

sim_ms <- simIntMsOcc(
  N_total   = 150,
  n_data    = 2,
  J         = c(4, 3),
  n_species = 8,
  n_shared  = 20,
  seed      = 100
)

fit_ms <- occu(
  ~ occ_x1, ~ det_x1,
  data         = sim_ms$data,
  multispecies = TRUE,
  integrated   = TRUE,
  verbose      = 0
)
summary(fit_ms)
```

This fits 8 species simultaneously across 2 data sources. Each species
has its own occupancy and detection parameters, but all species share
the community-level hyperparameters (mean and variance of the occupancy
and detection coefficient distributions). Rare species that appear in
only one data source borrow strength from the community mean and from
the overlap sites observed across sources.

Multi-species integrated models are particularly useful for biodiversity
monitoring programs that combine structured surveys (few species
observed reliably) with broad-coverage citizen science (many species
observed inconsistently). The hierarchical structure stabilises
estimates for the long tail of rarely detected species.

## WAIC for integrated models

Model comparison via WAIC works as expected:

``` r

w <- waicOccu(fit)
w
#>   component      elpd       pD    WAIC
#> 1 occupancy -757521.4 752147.7 1515043
#> 2     total -757521.4 752147.7 1515043
```

The WAIC is decomposed by component: one entry for the occupancy layer
and one for each detection source. This decomposition reveals which part
of the model drives complexity. Typically, the occupancy component
contributes the most to both the effective number of parameters
(\\p_D\\) and the overall WAIC, while the source with fewer visits and
sites is the simplest component.

Use the total WAIC to compare competing integrated models — for
instance, a model with shared versus source-specific detection formulas,
or a non-spatial versus spatial occupancy layer. The per-component
breakdown helps diagnose where a model is over- or under-fitting.

## When to use integrated models

Integrated models are appropriate when multiple survey programs cover
the same species in overlapping regions with different observation
protocols:

- **Structured + citizen science.** Structured surveys provide high
  detection reliability at a small number of sites, while
  citizen-science platforms like eBird or iNaturalist offer broad
  spatial coverage with variable and often lower detection probability.
  Combining them yields better occupancy estimates than either source
  alone.

- **Multiple sensor types.** Camera traps, acoustic recorders, and
  visual surveys detect different behaviours (movement, vocalisation,
  presence) and operate on different temporal scales. An integrated
  model acknowledges these differences explicitly.

- **Historical + modern data.** Older surveys using different methods
  (e.g., mist netting vs. point counts) can be combined with
  contemporary data, provided the occupancy state can be treated as
  shared within a season.

- **Multi-scale monitoring.** Regional monitoring programs with sparse
  but standardised coverage can be integrated with local intensive
  studies to improve estimates at both scales.

### When NOT to use integrated models

- **No spatial overlap.** If two sources cover completely disjoint
  regions with no shared sites, the model cannot calibrate detection
  across sources. Fit them separately and combine predictions post hoc.

- **One source has no detections.** A source with zero detections across
  all sites and visits contributes no information about detection and
  only adds parameters. Drop it.

- **Identical detection processes.** If two sources use the same
  protocol, observers, and equipment, they have the same detection
  process and should simply be pooled into a single detection matrix
  with more visits. Fitting them as separate sources wastes degrees of
  freedom on estimating differences that do not exist.

## Practical guidance

**Overlap is critical.** Sites surveyed by multiple sources are the
anchor points that let the model learn how detection differs across
programs. With zero overlap, the model must extrapolate detection
calibration entirely from the occupancy structure, which is fragile. Aim
for at least 10–20 shared sites. More overlap is always better, but even
a modest number of shared sites substantially improves identifiability.

**Unequal contributions.** Sources with more visits or higher per-visit
detection probability contribute more to the occupancy estimates. This
is desirable — it means the model automatically up-weights the more
informative source rather than treating all sources equally.

**Start simple.** Begin with a shared detection formula across sources.
This baseline model has fewer parameters and is easier to diagnose. Then
try source-specific detection formulas and compare via WAIC. Only add
source-specific terms when the data support them.

**Scaling considerations.** The total number of site-visit combinations
grows as \\\sum_d N_d \times J_d\\, which can become large with many
sources. The EM algorithm handles this efficiently because each
detection source is fit as a separate INLA call, but memory requirements
scale with the total data volume. For very large datasets, consider
reducing the number of visits (e.g., collapsing weekly camera occasions
into monthly summaries) before fitting.

**Convergence.** Integrated models typically require more EM iterations
than single-source models because the E-step must reconcile information
from multiple detection processes. Monitor convergence using the EM
trace (available via `plot(fit, type = "convergence")`) and increase
`max_iter` if needed. Warm-starting from a simpler model (e.g., the
non-spatial version) can accelerate convergence for spatial integrated
models.
