# R implementation

The R workflow is a complete, independent implementation of the analyses. It reads only the curated aggregate inputs in `data/` and writes generated files to `results/R/`.

## Run

From the repository root:

```bash
Rscript R/install.R
Rscript R/run.R
```

Use `Rscript R/run.R --quick` to retain only the one-year and four-year deletion-block summaries used by the figures.

## Scripts

- `01_main.R`: national primary models and direct road-user comparisons.
- `02_sensitivity.R`: reference, calendar, dynamic, conditional, curvature, Alaska-inclusion and common-target analyses.
- `03_registry.R`: fixed-state and Mississippi county analyses.
- `04_traffic.R`: matched FHWA traffic summaries.
- `05_outputs.R`: manuscript tables and ten figures.
- `functions/`: reusable estimation, inference, registry and plotting functions.
- `tests/`: numerical and end-to-end checks.

The allocation estimator uses grouped binomial counts with policy-year CR1 covariance. The fixed-margin sensitivity uses a separately implemented finite-support conditional likelihood. The year-specific-slope sensitivity conditions on the three-center sufficient statistics. None of these calculations call the Python implementation.
