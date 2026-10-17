# Validation

The public `v1.0.0` tagged version was validated from the curated aggregate inputs and published on 17 October 2026.

## Independent workflows

- The complete R workflow generated 33 estimate tables, two manuscript tables and ten figures.
- The complete Python workflow generated the same analytical output set.
- Each implementation ran directly from `data/`; neither implementation calls the other.
- The five-center curvature, exact policy-year multiplier, seasonal, policy-era, leave-one-year, Arizona-support and Alaska-inclusion diagnostics were independently implemented in both languages.

## Automated checks

- R: all `testthat` checks passed; 23 R source files parsed and ran in the complete workflow.
- Python: all 25 `pytest` checks passed.
- Python formatting and static checks passed with Ruff.
- Each workflow generated ten PDF figures, ten PNG figures and ten plot-data CSV files.

## Cross-language agreement

The validator compared 33 keyed output files containing 2,222 rows, including the Alaska inclusion sensitivity.

The independently reproduced main results were:

| Analysis | Evening RCR (95% CI) | Morning RCR (95% CI) |
|---|---:|---:|
| Fatal pedestrian crashes | 0.155 (0.084–0.286) | 0.204 (0.085–0.489) |
| Pedestrian versus other fatal crashes | 0.212 (0.116–0.389) | 0.214 (0.083–0.553) |
| Five-center quadratic pedestrian model | 0.144 (0.076–0.274) | 0.191 (0.080–0.456) |
| Main pedestrian model with Alaska included | 0.154 (0.083–0.285) | 0.206 (0.086–0.493) |

Exact policy-year Rademacher calibration gave pedestrian intervals of 0.082–0.291 and 0.083–0.499, and direct-comparison intervals of 0.114–0.395 and 0.081–0.566. The workflow enumerates all `2^23` unique sign patterns for the 24 policy-year influence contributions.

The executable validation commands are documented in `README.md`. Generated outputs remain outside version control.
