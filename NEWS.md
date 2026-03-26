# INLAocc 0.1.0

* Initial release.
* EM-INLA engine for occupancy models with scaled-binomial occupancy M-step.
* Single-species: `occu_inla()`.
* Spatial SPDE: `spatial_occu_inla()`.
* Multi-species with community pooling: `ms_occu_inla()`.
* Multi-season with AR(1): `temporal_occu_inla()`.
* Integrated multi-source: `intOccu_inla()`.
* Spatially-varying coefficients: `svcOccu_inla()`.
* Latent factor multi-species: `lfMsOccu_inla()`, `sfMsOccu_inla()`.
* JSDM without detection: `lfJSDM_inla()`.
* Native AST-based `(1 | group)` and `(x | group)` random effects parser (no lme4 dependency).
* spOccupancy-compatible data format and `priors` argument.
* `ppcOccu()`, `waicOccu()`, `fitted()`, `residuals()` diagnostics.
* `k.fold` cross-validation argument.
* Design-matrix prediction via `X.0` / `coords.0`.
* Simulation functions: `simulate_occu()`, `simMsOcc()`, `simTOcc()`, `simIntOcc()`.
