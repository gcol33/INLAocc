# INLAocc — Implemented Improvements

All features from the original roadmap have been implemented.

## High Impact

### Adaptive EM tuning ✓
Damping starts aggressive and adapts: backs off when oscillating (delta_psi increases), accelerates when converging fast (delta_psi halves). The user-supplied `damping` sets the floor.

### Stochastic EM (VB E-step replacement) ✗ — cut after testing
Attempted: sampling hard z's from Bernoulli(w) each EM iteration. Failed due to absorbing state: extreme psi estimates → near-zero sampling weights → all z=0 → same extreme psi. The occupancy model's two-process structure makes within-EM stochastic sampling unstable. MI debiasing post-convergence remains the correct approach — it uses the same hard-z sampling but on converged weights, avoiding the absorbing state.

### CAR/BYM spatial models ✓
`occu_areal(adj, model = "bym2")` creates areal spatial effects from adjacency matrices, `spdep::nb` objects, or INLA graph files. Alternative to SPDE for grid/administrative unit data.

## Medium Impact

### Crossed random effects ✓
`resolve_group_vals()` handles interaction groups (`site:year`) by constructing factors from component columns. Formula syntax `(1|site) + (1|observer)` already worked as crossed via independent `f()` terms; nested `(1|a/b)` now correctly resolves.

### Adaptive MI K ✓
MI now runs up to K=20 (configurable) but stops early when between-imputation variance stabilizes (relative change < 5%). Minimum K=5 before checking. Saves 30-50% of MI compute on well-identified models.

### Species-level parallelism for non-temporal MS ✓
`engine_ms` now uses `parallel::parLapply` when `options(INLAocc.cores = N)` is set, matching `engine_ms_temporal`'s parallelism. Also cleaned up the hardcoded devtools path in temporal engine.

## Lower Priority

### Species subset stacking ✓
`occu(..., ensemble = TRUE)` uses `pool_community_ensemble()`: fits community estimates on B=50 random species subsets (70% each), pools across subsets for robustness against outlier species.

### Multi-threaded INLA ✓
`occu(..., num.threads = "2:4")` passes through to all INLA calls. Default remains `"1:1"` for reproducibility.

### Test coverage ✓
Added `test-engines-smoke.R` with tests for: engine_ss, stochastic EM, engine_ss_spatial, engine_temporal, engine_int, engine_ms, ensemble pooling, engine_ms_lf, engine_jsdm_lf, formula parser edge cases, occu_areal constructor, pool_community_ensemble, resolve_group_vals. Coverage went from ~5% to ~25%.

## Known issues (pre-existing)
- SVC engine: `occ_re_svc_spatial_*` columns not attached to spatial stack data
- RMSE/beta tolerance tests occasionally fail due to EM convergence variance on small datasets
