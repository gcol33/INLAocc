# Multi-Season Occupancy Models

## Introduction

Many ecological studies span multiple years or seasons. A species that
occupies a site in one year may go locally extinct the next, and
previously unoccupied sites may be colonized. Single-season occupancy
models treat each site as a snapshot, discarding the temporal
information embedded in repeated surveys across years. When the data
contain temporal structure, ignoring it wastes information at best and
produces biased estimates at worst — occupancy trends are confounded
with detection variation, and standard errors are too narrow because the
model treats non-independent observations as independent.

Multi-season (dynamic) occupancy models address this by decomposing
temporal variation into a structured process. In the classical
formulation (MacKenzie et al. 2003), colonization and extinction rates
are modelled explicitly. INLAocc takes a different approach: it models
the temporal dynamics as an AR(1) random effect on occupancy, which
captures the same persistence structure (high occupancy tends to follow
high occupancy) without requiring separate colonization/extinction
parameters. This is more parsimonious when the goal is inference on
occupancy probability and its drivers rather than the transition rates
themselves.

INLAocc fits multi-season models by stacking all site-period
combinations into a joint occupancy model with an AR(1) temporal random
effect. Detection is fit per-period, so each season gets its own
detection parameters — survey conditions in July 2019 need not resemble
those in July 2023. The temporal random effect \\\eta_t\\ captures
shared year-to-year fluctuation: good years (favourable weather, high
resource availability) raise occupancy across all sites, and bad years
lower it. Fixed effects on occupancy and detection are estimated
simultaneously, with the AR(1) structure absorbing residual temporal
autocorrelation that would otherwise inflate standard errors.

## The multi-season occupancy model

For site \\i\\ at time period \\t\\, the occupancy process is:

\\z\_{it} \sim \text{Bernoulli}(\psi\_{it}), \quad
\text{logit}(\psi\_{it}) = \mathbf{x}\_i^\top \boldsymbol{\beta} +
\eta_t\\

where \\z\_{it}\\ is the true (latent) occupancy state and
\\\psi\_{it}\\ is the occupancy probability. The temporal random effect
follows an AR(1) process:

\\\eta_t = \rho \cdot \eta\_{t-1} + \epsilon_t, \quad \epsilon_t \sim
\mathcal{N}(0, \sigma^2\_\eta)\\

where \\\rho \in (-1, 1)\\ is the temporal autocorrelation parameter and
\\\sigma^2\_\eta\\ is the innovation variance.

The observation process (within each period) is:

\\y\_{ijt} \mid z\_{it} \sim \text{Bernoulli}(z\_{it} \cdot p\_{ijt}),
\quad \text{logit}(p\_{ijt}) = \mathbf{w}\_{ijt}^\top
\boldsymbol{\alpha}\\

where \\y\_{ijt}\\ is the detection/non-detection outcome at site \\i\\,
visit \\j\\, period \\t\\, and \\p\_{ijt}\\ is the detection
probability.

Three features of this formulation are worth highlighting:

- **\\\eta_t\\ captures shared temporal trends.** All sites experience
  the same temporal random effect in a given period. This represents
  “good years vs. bad years” for the species — shared environmental
  conditions that raise or lower occupancy across the landscape.

- **\\\rho\\ controls persistence.** A high \\\rho\\ (e.g., 0.8) means
  occupancy patterns persist across periods: if this year is good, next
  year is likely good too. A low \\\rho\\ means each year is nearly
  independent of the last.

- **Detection varies across periods.** The detection covariates
  \\\mathbf{w}\_{ijt}\\ can include period-specific variables (e.g.,
  survey date, observer effort, weather during surveys). Each period’s
  detection model is estimated from that period’s data, so changing
  survey protocols across years are handled naturally.

What does \\\rho = 0.7\\ mean in ecological terms? It means that if
occupancy is above average this year, about 70% of that excess carries
over to the next year. A site network where most occupied sites remain
occupied and most empty sites remain empty produces high \\\rho\\. This
is typical for species with strong site fidelity or slow population
dynamics — think territorial birds, long-lived amphibians, or perennial
plants. A low \\\rho\\ (say 0.2) indicates that knowing this year’s
occupancy tells you little about next year: the species may be nomadic,
or local populations may be driven by stochastic events like flooding or
fire.

The innovation standard deviation \\\sigma\_\eta\\ captures the
magnitude of year-to-year environmental stochasticity after accounting
for temporal persistence. A large \\\sigma\_\eta\\ with high \\\rho\\
means the system has strong shocks that then persist (a bad year stays
bad for several years). A large \\\sigma\_\eta\\ with low \\\rho\\ means
large but transient fluctuations. The classical colonization-extinction
parameterization (MacKenzie et al. 2003) models the same persistence
through explicit transition probabilities: colonization \\\gamma\\ and
extinction \\\epsilon\\. Both formulations capture temporal dependence,
but the AR(1) is more parsimonious when colonization and extinction
rates are not themselves the inferential target. If you need to estimate
how quickly empty sites are recolonized, the dynamic parameterization is
more natural. If you need to estimate occupancy trends and their
drivers, the AR(1) is typically sufficient and avoids identifiability
issues that arise when both \\\gamma\\ and \\\epsilon\\ interact with
covariates.

## Data structure

Multi-season data in INLAocc uses a 3D array with dimensions sites
\\\times\\ periods \\\times\\ visits:

    y[site, period, visit]
    dim: N × T × J

Each entry is 1 (detected), 0 (not detected), or `NA` (visit not
conducted). Here is an example for a single site surveyed over 5 years
with 3 visits per year:

``` r

# 100 sites, surveyed over 5 years, 3 visits per year
dim(y)
#> [1] 100   5   3

y[1, , ]  # site 1, all years:
#>      [,1] [,2] [,3]
#> [1,]    1    0    1    # year 1: detected twice
#> [2,]    0    0    0    # year 2: never detected
#> [3,]    1    1    0    # year 3: detected twice
#> [4,]    0    0   NA    # year 4: detected zero, one visit missing
#> [5,]    1    1    1    # year 5: detected all visits
```

Covariates follow the same dimensional logic:

- **`occ.covs`**: Site-level covariates can be either a length-\\N\\
  vector (static across periods, e.g., elevation) or an \\N \times T\\
  matrix (time-varying, e.g., yearly temperature anomaly). Static
  covariates are recycled across all periods internally.

- **`det.covs`**: Visit-level covariates are \\N \times T \times J\\
  arrays. For example, Julian survey date varies by site, year, and
  visit.

- **`coords`**: An \\N \times 2\\ matrix of site coordinates, constant
  across periods. Used only if a spatial component is added (see Section
  8).

Missing visits (`NA` in the detection array) are handled automatically.
Entire missing periods for a site (all visits `NA`) effectively remove
that site-period combination from the likelihood.

Common temporal data collection designs map onto this structure in
different ways. Annual surveys (e.g., breeding bird surveys repeated
yearly) produce one period per year, with visits as replicate surveys
within each breeding season. Seasonal surveys (spring and autumn counts
within each year) can be treated as separate periods if occupancy is
expected to change between seasons, or as additional visits within a
single period if the species’ presence is assumed stable within a year.
Before-after-control-impact (BACI) designs fit naturally into the 3D
array, with a period index that spans the pre- and post-intervention
years. The key question is whether the species’ true occupancy state can
change between two time points: if yes, they should be separate periods;
if no, they are replicate visits within the same period.

## Simulating temporal data

The
[`simTOcc()`](https://gillescolling.com/INLAocc/reference/simTOcc.md)
function generates multi-season data with known parameters, which is
useful for testing model recovery and understanding the data structure:

``` r

library(INLAocc)

sim <- simTOcc(
  N = 100, J = 3, n_seasons = 5,
  beta_occ = c(0.5, 0.3),
  beta_det = c(0, -0.3),
  rho = 0.7, sigma_t = 0.5,
  seed = 42
)

str(sim$data)
#> List of 4
#>  $ y       : int [1:100, 1:5, 1:3] 1 0 0 0 0 1 0 0 1 0 ...
#>  $ occ.covs:List of 1
#>   ..$ occ_x1: num [1:100] 1.201 1.045 -1.003 1.848 -0.667 ...
#>  $ det.covs:List of 1
#>   ..$ det_x1: num [1:100, 1:5, 1:3] -2.001 0.334 1.171 2.06 -1.377 ...
#>  $ coords  : num [1:100, 1:2] 0.915 0.937 0.286 0.83 0.642 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "x" "y"
```

The returned object also contains the true parameter values and latent
states, so you can verify that the model recovers them:

``` r

sim$truth$temporal_re
#>               [,1]        [,2]         [,3]        [,4]         [,5]
#>   [1,]  1.01394398  0.66683296 -0.320101831  0.42216680  0.337750636
#>   [2,]  0.67706762  1.44350383  1.550389003  0.90847554  0.568471315
#>   [3,] -0.69113094 -0.89464573 -0.381486552 -0.41667932  0.037451928
#>   [4,]  0.58280326  0.58969391  0.213462908 -0.59954096 -0.362364584
#>   [5,] -0.72661697 -0.18564244  0.256472990  0.33753817 -0.282779692
#>   [6,] -0.52870185 -0.61584345 -0.357477544  0.92801724  0.189277685
#>   [7,] -1.51009550 -1.11726687 -1.201708976 -1.92741624 -1.786345354
#>   [8,] -0.26155786 -0.73418936 -0.156102621  0.34454832  0.078781256
#>   [9,] -1.03731764 -0.42903678 -0.727763087 -0.92439293 -0.572521643
#>  [10,] -0.68343893 -0.62574946 -0.823512076 -0.38031660 -0.689374572
#>  [11,] -0.59075011 -0.89446930 -0.620603684  0.17305190 -0.323302925
#>  [12,]  0.15270143 -0.22116750 -0.699075175  0.05621283 -0.245264682
#>  [13,]  1.13479407  0.45086163  0.644101492  1.00621715  1.774712371
#>  [14,] -0.58179664 -0.36905472  0.209573920 -0.10413161  0.224534607
#>  [15,] -0.30151613 -0.77773082 -1.662997351 -0.90677836 -0.805698409
#>  [16,]  1.20275478  0.28373939 -0.045622238  0.37173058  0.129707035
#>  [17,]  0.34171029  0.91558911  1.729320176  1.23621768  0.446304652
#>  [18,] -0.75938240 -0.80483406 -0.776638697 -1.14914587 -0.505505607
#>  [19,]  0.99692646  0.47619921 -0.186827288 -0.76772764 -0.152333542
#>  [20,]  0.02207025 -0.93697820 -0.717762329 -0.50226944 -0.773477782
#>  [21,] -0.85510440 -0.82506607 -0.928706218 -0.98728455 -1.842610777
#>  [22,] -0.10455395  0.40115703  0.858991153  1.63847005  1.279574020
#>  [23,] -0.73655001 -0.26441389  0.364832792  1.19863537  1.244671509
#>  [24,]  0.99296928  0.91466345  0.001513608 -0.42962708  0.642753984
#>  [25,]  0.49295972  0.01121154 -0.345009092 -0.70186565 -0.417069312
#>  [26,] -0.17504717 -0.59341749  0.191707264  0.14347830  0.479350661
#>  [27,] -0.26609711 -0.09879362  0.447606118  0.65118628  0.926471018
#>  [28,] -1.27889686 -0.80662561 -0.660541857 -0.78747865 -0.831583229
#>  [29,]  1.17327073  1.03655681  0.653286459  0.92681338  0.376020989
#>  [30,]  0.09721470  0.79148145  0.387722969 -0.12817206 -0.363866813
#>  [31,] -0.73570450 -0.38866190 -0.108132406  0.35214352  1.147381597
#>  [32,]  0.98386536  0.06330986  1.043871326  0.80034917  0.366909869
#>  [33,]  1.29920485  1.12387109  1.652104831  1.13471154  0.428755177
#>  [34,] -1.37645017 -0.40008246 -0.724099918 -0.07119031  0.010113565
#>  [35,]  0.33077062  0.86485939  0.362743384  0.02315358 -0.099671666
#>  [36,]  0.74239613  0.28228444  0.109843969  0.47505303  0.721909777
#>  [37,]  0.55787338  0.41503715  0.671832475  0.61327927  0.112887968
#>  [38,]  1.26332183 -0.04718021  0.389738045  0.44365047  0.321944134
#>  [39,] -0.42806561 -0.42990932  0.601898636  0.37282941  0.068737238
#>  [40,] -0.61497651 -0.67827282 -0.516937198 -1.04071707 -0.157153042
#>  [41,]  0.06849263  0.78045535  1.041536673  1.43626905  0.337439448
#>  [42,]  0.40578971  0.10314219  0.048773310 -0.26559742 -0.078365947
#>  [43,]  0.06555515 -1.01402860 -0.283598476 -0.12133165 -0.570862973
#>  [44,] -0.77115373 -0.88493006 -1.411268500 -0.91900492  0.107868536
#>  [45,]  0.13504412  0.35810825 -0.490927368 -1.11324952 -1.227491107
#>  [46,]  0.50484997 -0.08885600  0.217100066  0.04260129  0.221632647
#>  [47,] -0.72752960 -0.98376684 -0.740717267 -0.62285764 -1.330365444
#>  [48,]  0.37603431 -0.71162402 -0.705754044 -0.46646300 -0.677523956
#>  [49,] -0.98288392 -0.88743071 -0.676412396 -0.93552790 -1.333491032
#>  [50,] -0.83639043 -0.84971719 -0.821115478  0.13496176 -0.548367339
#>  [51,]  0.26685763 -0.02193396  0.031076941 -0.31003841  0.184386808
#>  [52,] -1.19239635 -0.75772391 -0.517368989  0.16577340  0.003807821
#>  [53,]  0.34768702  0.44319789 -0.013217853 -0.05113388 -0.102604316
#>  [54,]  0.78387949  0.66832613  0.613129399  1.60099129  1.328842958
#>  [55,]  0.25378752  0.44349257  1.084465921  1.23987030  1.037127296
#>  [56,] -0.42775266 -0.12090084 -0.013876247 -0.92785904 -1.420016854
#>  [57,] -1.25631929 -0.17216605  0.188652710 -1.04317080 -0.496703765
#>  [58,]  0.04469130 -0.09684631 -0.440053439 -0.94255005 -0.228471028
#>  [59,] -0.40393507 -0.34888822  0.049192296  0.30193475  0.014649354
#>  [60,] -0.37850230 -0.13522615 -0.683797687 -0.32558914 -0.141171785
#>  [61,]  0.17543974 -0.01615421 -0.873675812 -1.61492554 -1.776352044
#>  [62,]  0.25613799  0.10319497 -0.294810568 -0.59735343 -0.142363881
#>  [63,] -1.12781427 -0.41538178 -0.066119237 -0.48131473 -0.385794274
#>  [64,]  0.62529522  1.11316420  1.523700139  1.25212370  0.769249971
#>  [65,]  0.31891000  0.09418145 -0.435534864  0.11719734  0.429367524
#>  [66,]  0.88287922  1.45129510  1.751835537  1.62094273  1.116432774
#>  [67,] -0.06391108 -0.25116922 -0.566349056 -0.59149073 -0.107124701
#>  [68,]  0.26932779 -0.41429888  0.308916397  0.86181211  1.344804279
#>  [69,] -0.85668802 -0.80626893 -0.769178685 -0.76400744 -0.817154683
#>  [70,] -0.39208566 -0.75973899 -0.813229070 -0.73527857 -0.417869944
#>  [71,]  0.42403481  0.77576117  0.142623600  0.12111085  1.252628024
#>  [72,] -0.32337190  0.14385462 -0.248474094  0.35875258  0.523885220
#>  [73,] -0.64287511  0.06258309  0.606677367  0.86607435 -0.777787261
#>  [74,]  0.10069875  0.83043069  0.495643441  0.01787791 -0.064632765
#>  [75,] -1.66846125 -1.73707142 -1.259577759 -0.46404107 -0.218773316
#>  [76,]  0.44903223 -0.20172793  0.601675779  1.01641117  1.020101184
#>  [77,] -0.12108105  0.34603349 -0.067819520 -0.48658019 -1.454903705
#>  [78,]  0.13531132 -1.06161043 -0.755294199 -0.30803564  0.179163426
#>  [79,]  0.71541684  0.80842237  0.136106728  0.54064920  0.410091812
#>  [80,] -1.30552387 -0.97043513 -1.570997782 -1.62644367 -1.448187906
#>  [81,] -0.41093920  0.17583043  0.090378376 -0.40529102 -0.600633428
#>  [82,] -0.40273756 -0.18780577  0.138363950 -0.11358727 -0.249244369
#>  [83,]  0.36272882 -0.18078363 -0.206884000 -0.21612719 -0.533894519
#>  [84,] -0.46195280 -0.47706968 -0.774912501 -0.59011133 -1.051062618
#>  [85,]  0.26296209 -0.05054173  0.364087394  0.46800528  0.071537686
#>  [86,] -0.24330954  0.13437969  0.461797944  0.56161284  0.557271872
#>  [87,] -1.24713668 -0.81329747 -1.278455914 -0.85690195 -0.930247987
#>  [88,] -0.07386901 -0.64175125 -0.531405468  0.38773656  1.183998131
#>  [89,] -1.74239641 -2.04683861 -2.143214669 -2.13839855 -0.457328586
#>  [90,]  0.71101743  0.44290579  0.800793497  1.69122905  1.502405298
#>  [91,]  1.24312423  0.29643619 -0.332155966  0.30146575  0.656902623
#>  [92,]  0.60918941  0.66931382  0.232972826  0.63789597  1.370390978
#>  [93,]  0.86845891 -0.08248587  0.179598328  0.13189201  0.135874375
#>  [94,] -0.12768327 -0.72221343 -0.339907441  0.43792911  1.124833243
#>  [95,]  0.33601977 -0.13294412 -0.384569191 -0.57491442 -0.100765075
#>  [96,] -0.16506735  1.33656090  0.801446989  0.63505932  0.341144030
#>  [97,] -0.44988208  0.10886845 -0.187578145 -0.13011932 -0.232023643
#>  [98,] -1.98032836 -1.06076470 -0.663655701 -0.89127750  0.057655792
#>  [99,] -0.88420418 -1.50791620 -1.399001161 -1.57146857 -0.705189351
#> [100,]  0.06018668 -0.14290005 -0.095361849 -0.54228309 -0.288739307
sim$truth$rho
#> [1] 0.7
sim$truth$sigma_t
#> [1] 0.5
```

The temporal random effects show the year-to-year fluctuation. Years
with the highest values had elevated occupancy across all sites, while
years with negative values were below-average.

## Fitting the temporal model

Pass `temporal = "ar1"` to
[`occu()`](https://gillescolling.com/INLAocc/reference/occu.md) to fit a
multi-season model. The function detects the 3D array and sets up the
stacked likelihood with the AR(1) temporal random effect automatically:

``` r

fit <- occu(~ occ_x1, ~ det_x1, data = sim$data, temporal = "ar1", verbose = 0)
summary(fit)
```

The output has three blocks. The occupancy and detection fixed effects
are interpreted exactly as in a single-season model. The temporal
component block is new:

- **Fixed effects should be recovered well.** The true occupancy
  intercept was 0.5 and the true slope was 0.3. Both 95% credible
  intervals should contain the true values.

- **AR(1) correlation \\\rho\\** should be close to the true value of
  0.7, indicating strong temporal persistence — occupancy conditions in
  one year are a good predictor of the next.

- **Innovation SD \\\sigma\\** should be close to the true value of 0.5,
  the standard deviation of the year-to-year shocks after accounting for
  the autocorrelation.

- **Occupancy varies across periods.** The temporal random effect
  creates meaningful variation in landscape-level occupancy across
  years.

The estimated \\\rho\\ has direct management implications. A high value
(above 0.6) means occupancy patterns are predictable: sites that are
occupied this year are likely occupied next year, and population
monitoring can rely on less frequent surveys because the system changes
slowly. Conservation planning can target currently occupied sites with
reasonable confidence that they will remain occupied in the near term. A
low \\\rho\\ (below 0.3) signals high environmental stochasticity — the
species’ distribution reshuffles substantially from year to year, and
snapshot surveys from a single year may not represent the long-term
pattern. In this regime, protecting a fixed set of sites is less
effective than maintaining landscape connectivity so the species can
track shifting conditions.

The range of period-specific occupancy estimates (the minimum and
maximum across years) gives a sense of overall population dynamics. If
occupancy ranges from 0.3 to 0.7 across five years, the species
undergoes large swings that may reflect boom-bust dynamics, cyclic prey
availability, or irregular disturbance events. A narrow range (0.45 to
0.55) suggests a stable population with the AR(1) absorbing only minor
fluctuations. Compare these estimates against any known events in the
study area — a sharp drop in a particular year may coincide with a
drought, a harsh winter, or a habitat disturbance.

### Extracting temporal correlation

The
[`temporalCorr()`](https://gillescolling.com/INLAocc/reference/temporalCorr.md)
function extracts the AR(1) parameters with full posterior summaries:

``` r

temporalCorr(fit)
```

The period-specific random effects should track the true values
reasonably well, with wider credible intervals for the endpoints
(periods 1 and 5 have fewer temporal neighbours).

## Checking temporal residuals

The Durbin-Watson test checks whether residual temporal autocorrelation
remains after fitting the AR(1). If the model captures the temporal
structure adequately, residuals should show no remaining autocorrelation
(DW statistic near 2):

``` r

durbinWatson(fit)
```

A DW statistic close to 2 with a non-significant p-value confirms the
AR(1) structure is capturing the temporal dependence adequately — no
residual autocorrelation remains.

Values substantially below 2 (e.g., DW \< 1.5 with a small p-value)
would indicate remaining positive autocorrelation, suggesting that the
AR(1) model is insufficient. In that case, consider adding a linear time
trend as a fixed effect alongside the AR(1), or investigate whether a
specific year is an outlier.

Residual temporal autocorrelation means the model’s temporal structure
is not flexible enough. Several ecological scenarios produce this.
Regime shifts — where the system flips between two stable states (e.g.,
forested vs. open) — create autocorrelation patterns that an AR(1)
cannot capture because the process is non-stationary. Non-stationary
trends (accelerating decline, exponential growth during invasion) also
leave residual autocorrelation because the AR(1) assumes constant
\\\rho\\ and \\\sigma\\. Seasonal cycles within a year (if periods are
months rather than years) produce periodic autocorrelation at the
seasonal lag. In these cases, adding a fixed time trend, a second-order
AR term, or restructuring the period definition may be needed.

The general diagnostic suite also works with temporal models:

``` r

checkModel(fit)
```

## Visualizing temporal trends

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
temporal models shows occupancy trajectories across periods:

``` r

plot(fit)
```

This produces a ribbon plot with the mean occupancy per period and 95%
credible intervals. The x-axis is the period index (1 to \\T\\) and the
y-axis is the estimated occupancy probability.

To extract period-specific occupancy estimates programmatically:

``` r

psi_by_period <- fitted(fit)$psi  # N x T matrix
colMeans(psi_by_period)
```

These values represent the average estimated occupancy across all sites
for each period. The trajectory should show variation across years,
consistent with the temporal random effects estimated for each period.

When reading temporal trajectories, look for three broad patterns. A
monotonic decline (occupancy dropping steadily across periods) suggests
ongoing habitat loss or increasing pressure from an invasive competitor;
this pattern often warrants adding a linear time trend as a fixed effect
to separate the directional signal from stochastic fluctuation. Cyclical
variation (occupancy oscillating with a regular period) may reflect
predator-prey dynamics, resource pulses, or climatic oscillations like
ENSO. Step changes — a sudden drop or jump between two consecutive
periods — often coincide with discrete management interventions (e.g.,
habitat restoration, species reintroduction) or catastrophic events
(wildfire, flooding). Overlaying known events on the temporal trajectory
plot helps distinguish real ecological signals from sampling noise.

You can also compute site-specific trajectories:

``` r

# Occupancy trajectory for site 1
psi_by_period[1, ]

# Sites with the most variable occupancy across periods
site_var <- apply(psi_by_period, 1, var)
head(order(site_var, decreasing = TRUE))
```

## Prediction at new sites and times

To predict occupancy at new covariate values across all periods:

``` r

new_sites <- data.frame(occ_x1 = c(-1, 0, 1))
preds <- predict(fit, X.0 = new_sites, type = "occupancy")
head(preds)
```

Predictions are period-specific: the temporal random effect shifts all
predictions up or down in each period, so a site’s predicted occupancy
depends on both its covariate value and which year you are predicting
for.

To predict at a specific period only (e.g., the most recent year):

``` r

predict(fit, period = 5)
```

## K-fold cross-validation

For out-of-sample validation of the temporal model:

``` r

fit_cv <- occu(~ occ_x1, ~ det_x1, data = sim$data,
               temporal = "ar1", k.fold = 5, verbose = 0)
fit_cv$k.fold
```

K-fold cross-validation for temporal models holds out entire sites (not
individual site-periods), so the held-out sites contribute no
information during fitting. This tests whether the model generalises to
new locations, not whether it interpolates between observed time points.

## Space-time models

Temporal and spatial structure often coexist. A species may be more
common in some regions (spatial structure) *and* more common in some
years (temporal structure). INLAocc fits joint space-time models by
combining the SPDE spatial random field with the AR(1) temporal random
effect:

``` r

fit_st <- occu(
  ~ occ_x1, ~ det_x1, data = sim$data,
  spatial = sim$data$coords,
  temporal = "ar1", verbose = 0
)
summary(fit_st)
```

The spatial field is shared across all periods (the same spatial pattern
applies every year), and the temporal field is shared across all sites
(all sites experience the same temporal fluctuation). The two components
are additive on the logit scale:

\\\text{logit}(\psi\_{it}) = \mathbf{x}\_i^\top \boldsymbol{\beta} +
w(\mathbf{s}\_i) + \eta_t\\

Since this simulated dataset has no spatial structure (coordinates were
not used in simulation), the spatial variance is estimated near zero
(0.18), and the temporal parameters are essentially unchanged. In real
data where spatial structure exists, including it alongside the temporal
component typically sharpens fixed-effect estimates by absorbing an
additional source of confounding variation.

For details on mesh construction, spatial priors, and mapped
predictions, see
[`vignette("spatial-models")`](https://gillescolling.com/INLAocc/articles/spatial-models.md).

## Multi-species temporal models

Multi-species temporal models combine the multi-species framework with
the AR(1) temporal structure. Each species gets its own AR(1) process,
capturing the fact that good years for one species may not be good years
for another:

``` r

sim_mt <- simTMsOcc(
  N = 50, J = 3, n_species = 5, n_seasons = 4,
  seed = 100
)

fit_mt <- occu(
  ~ 1, ~ 1, data = sim_mt$data,
  multispecies = TRUE, temporal = "ar1", verbose = 0
)
summary(fit_mt)
```

Species with the highest temporal autocorrelation have the most
persistent occupancy patterns — once they colonize a landscape, they
tend to stay. Species with low autocorrelation indicate more volatile
dynamics. The community-level summary averages across species, giving an
overall sense of how persistent occupancy is in this assemblage.

For details on the multi-species framework, see
[`vignette("multi-species")`](https://gillescolling.com/INLAocc/articles/multi-species.md).

## Comparing temporal vs. static models

When is the temporal model worth the extra complexity? Model comparison
via WAIC provides a principled answer:

``` r

# Null temporal model: intercept only
fit_null <- occu(~ 1, ~ 1, data = sim$data, temporal = "ar1", verbose = 0)

# Covariate temporal model
fit_cov <- occu(~ occ_x1, ~ det_x1, data = sim$data, temporal = "ar1", verbose = 0)

compare_models(null = fit_null, covariate = fit_cov, criterion = "waic")
#>       model loglik df AIC BIC WAIC n_iter converged delta weight
#> 1      null     NA  3  NA  NA   NA     26      TRUE    NA     NA
#> 2 covariate     NA  4  NA  NA   NA     39      TRUE    NA     NA
```

The temporal model should have a substantially lower WAIC and receive
the majority of the model weight. This is expected, since the data were
generated with \\\rho = 0.7\\ — strong temporal structure exists and
ignoring it degrades model fit.

When would the static model be preferred? If \\\rho \approx 0\\ (no
temporal autocorrelation) and the temporal random effects are
negligible, the static model will have a similar or lower WAIC because
it avoids estimating the unnecessary AR(1) parameters. As a rule of
thumb: if the credible interval for \\\rho\\ includes zero and the WAIC
difference is small (\< 2), the simpler static model is sufficient.

## Practical guidance

**Minimum number of seasons.** At least 3–4 seasons are needed for
stable AR(1) estimation. With only 2 seasons, \\\rho\\ is very weakly
identified because there is only one transition to inform the
autocorrelation. Five or more seasons produce substantially tighter
credible intervals on the temporal parameters.

**Study design considerations.** For trend estimation, the number of
periods matters more than the number of visits per period, up to a
point. Five periods with 3 visits each gives a reasonable AR(1)
estimate; 10 or more periods are needed for nonlinear trends or for
detecting changes in \\\rho\\ over time. Within each period, the timing
of visits affects detection parameters — surveys conducted early and
late in a season may yield different detection probabilities if the
species’ activity varies seasonally. If visit timing is inconsistent
across years, include Julian date as a detection covariate to absorb
this variation rather than letting it leak into the occupancy estimates.

**When to use the static model instead.** If \\\rho\\ is estimated near
zero and the temporal model does not improve WAIC over the static model,
use the static model. It is simpler, faster, and the temporal random
effects are adding noise rather than signal.

**Time-varying covariates.** Occupancy covariates that change across
periods (e.g., yearly mean temperature, annual land-use change) should
be supplied as \\N \times T\\ matrices in `occ.covs`. The AR(1) temporal
effect then captures residual temporal variation *after* accounting for
the measured covariates. This is analogous to including both fixed and
random effects for time in a mixed model.

**Detection varies by season.** Detection conditions often change across
seasons — different observers, different phenological timing, different
weather. INLAocc estimates separate detection parameters for each period
by default when it detects a 3D detection array. This is usually the
right behaviour. If you want to constrain detection to be constant
across periods (e.g., standardized protocol), collapse the detection
array to 2D before fitting.

**Non-stationary temporal patterns.** The AR(1) model assumes
stationarity — the same autocorrelation and innovation variance apply
throughout the time series. If you suspect a monotonic trend (e.g.,
consistent decline due to habitat loss), the AR(1) will partially
capture it through a sequence of same-sign innovations, but it is more
efficient to add a linear time trend as a fixed effect:

``` r

# Add a linear time index as a covariate
sim$data$occ.covs$time_index <- matrix(
  rep(1:5, each = 100), nrow = 100, ncol = 5
)

fit_trend <- occu(
  ~ occ_x1 + time_index, ~ det_x1,
  data = sim$data, temporal = "ar1", verbose = 0
)
```

The fixed time trend captures the directional component, and the AR(1)
captures the residual year-to-year fluctuation around the trend.

**Computational considerations.** Space-time models (spatial + temporal)
stack \\N \times T\\ site-period combinations, so the effective sample
size grows quickly. For large \\N\\ (\> 500 sites) and long time series
(\\T\\ \> 10), consider:

- Fitting spatial-only and temporal-only models first to assess which
  component matters more.

- Using a coarser mesh for the spatial SPDE to reduce the number of mesh
  nodes.

- Subsetting to a region of interest if the full dataset is too large.

The temporal component alone (without spatial) is fast even for large
\\N\\, because the AR(1) precision matrix is tridiagonal and INLA
exploits its sparsity.
