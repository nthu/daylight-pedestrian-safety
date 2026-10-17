"""Run the complete analysis in manuscript order."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parent
STAGES = ["01_main.py", "02_sensitivity.py", "03_registry.py", "04_traffic.py", "05_outputs.py"]


def _load(path: Path) -> ModuleType:
    module_name = path.stem.replace("-", "_")
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise ImportError(f"Could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def run(*, quick: bool = False) -> None:
    for filename in STAGES:
        print(f"\n{'=' * 72}\n{filename}\n{'=' * 72}")
        module = _load(ROOT / filename)
        if filename == "02_sensitivity.py":
            module.run(full=not quick)
        else:
            module.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Fit only the one-year and four-year deletion summaries needed for the figures.",
    )
    arguments = parser.parse_args()
    run(quick=arguments.quick)
