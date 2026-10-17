"""Load the curated aggregate inputs used by the study.

The public repository starts from analysis-level tables. Provider download
archives and event-level working databases are intentionally outside this
package; see ``data/README.md`` for the release boundary.
"""

from __future__ import annotations

from pathlib import Path
import os

import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[2]
DATA = REPOSITORY / "data"
RESULTS = Path(os.environ.get("DAYLIGHT_RESULTS_DIR", REPOSITORY / "results" / "python"))


def read_csv(name: str, *, folder: Path | None = None) -> pd.DataFrame:
    """Read a CSV or gzip-compressed CSV by repository-relative name."""

    base = DATA if folder is None else folder
    path = base / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def national_counts() -> pd.DataFrame:
    return read_csv("national.csv.gz")


def state_counts() -> pd.DataFrame:
    return read_csv("states.csv.gz")


def mississippi_counts() -> pd.DataFrame:
    return read_csv("mississippi.csv.gz")


def reference_counts() -> pd.DataFrame:
    return read_csv("reference_holiday.csv.gz")


def traffic_counts() -> pd.DataFrame:
    return read_csv("traffic.csv")


def primary_tables() -> pd.DataFrame:
    return read_csv("primary_tables.csv")


def policy_calendar() -> pd.DataFrame:
    return read_csv("calendar.csv")


def sensitivity_table(name: str) -> pd.DataFrame:
    return read_csv(name, folder=DATA / "sensitivity")


def support_table(name: str) -> pd.DataFrame:
    return read_csv(name, folder=DATA / "support")


def ensure_results() -> dict[str, Path]:
    """Create and return the standard result directories."""

    paths = {
        "root": RESULTS,
        "estimates": RESULTS / "estimates",
        "tables": RESULTS / "tables",
        "figures": RESULTS / "figures",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    """Write a deterministic CSV with stable column order and no index."""

    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")
