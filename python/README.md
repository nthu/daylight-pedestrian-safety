# Python implementation

The Python workflow is a complete, independent implementation of the analyses. It reads only the curated aggregate inputs in `data/` and writes generated files to `results/python/`.

## Run

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r python/requirements.txt
python python/run.py
```

Use `python python/run.py --quick` to retain only the one-year and four-year deletion-block summaries used by the figures.

## Scripts

- `01_main.py`: national primary models and direct road-user comparisons.
- `02_sensitivity.py`: reference, calendar, dynamic, conditional, curvature, Alaska-inclusion and common-target analyses.
- `03_registry.py`: fixed-state and Mississippi county analyses.
- `04_traffic.py`: matched FHWA traffic summaries.
- `05_outputs.py`: manuscript tables and ten figures.
- `functions/`: reusable estimation, inference, registry and plotting functions.
- `tests/`: numerical and end-to-end checks.

The workflow does not call R. Cross-language agreement is checked separately by `validate_cross_language.py` after both workflows have completed.
