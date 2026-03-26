# Contributing to INLAocc

Thank you for your interest in contributing to INLAocc! This document
outlines how to get started.

Please note that this project has a [Code of
Conduct](https://gillescolling.com/INLAocc/CODE_OF_CONDUCT.md). By
participating, you agree to abide by its terms.

## Installation from Source

``` bash
git clone https://github.com/gcol33/INLAocc.git
cd INLAocc
```

Then in R:

``` r

install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
devtools::install_deps(dependencies = TRUE)
devtools::load_all()
```

## Testing

``` r

devtools::test()
devtools::check()
```

## Documentation

Rebuild documentation after editing roxygen comments:

``` r

devtools::document()
```

## Project Organization

    INLAocc/
    ├── R/                  # Source files
    │   ├── occu.R          # Unified API entry point
    │   ├── occu_dispatch.R # Flag-to-engine dispatch
    │   ├── occu_engines.R  # Engine functions per model type
    │   ├── occu_engine.R   # Core EM-INLA algorithm
    │   ├── occu_fit.R      # Legacy direct-call wrappers
    │   ├── occu_data.R     # Data formatting
    │   ├── occu_random.R   # Random effects parsing
    │   ├── occu_spatial.R  # SPDE mesh construction
    │   ├── occu_predict.R  # Prediction methods
    │   ├── occu_output.R   # print/summary/plot methods
    │   ├── occu_diagnostics.R # ppcOccu, waicOccu, fitted, residuals
    │   ├── occu_utils.R    # Shared helpers
    │   └── occu_compat.R   # spOccupancy compatibility
    ├── man/                # Auto-generated (roxygen2)
    ├── tests/testthat/     # Test files
    └── inst/examples/      # Runnable demos

## Contributing Workflow

1.  Fork the repository and create a feature branch
2.  Make your changes, keeping commits focused
3.  Run
    [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
    and
    [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
4.  Submit a pull request with a clear description

## Style Guidelines

- R functions: camelCase for exports (e.g., `occu_inla`, `ppcOccu`)
- Arguments: dot.case following spOccupancy conventions (e.g.,
  `occ.formula`, `det.covs`)
- Internal helpers: snake_case with `@noRd`
- Use roxygen2 with markdown enabled

## Reporting Bugs

Open an issue at <https://github.com/gcol33/INLAocc/issues> with:

- A minimal reproducible example
- Output of [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
- Expected vs. actual behavior

## License

By contributing, you agree that your contributions will be licensed
under the MIT License.
