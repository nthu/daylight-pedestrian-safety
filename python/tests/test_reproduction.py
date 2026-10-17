"""End-to-end checks of the numbered analysis workflow."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PYTHON_ROOT = Path(__file__).resolve().parents[1]
REFERENCE = PYTHON_ROOT / "reference"


@pytest.fixture(scope="module")
def quick_results(tmp_path_factory: pytest.TempPathFactory) -> Path:
    results = tmp_path_factory.mktemp("python-results")
    environment = os.environ.copy()
    environment["DAYLIGHT_RESULTS_DIR"] = str(results)
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[variable] = "1"

    completed = subprocess.run(
        [sys.executable, "-B", "run.py", "--quick"],
        cwd=PYTHON_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stderr[-3000:]
    return results


def _compare(frame: pd.DataFrame, expected: pd.DataFrame, keys: list[str], fields: list[str], tolerance=3e-6):
    left = frame.set_index(keys).sort_index()
    right = expected.set_index(keys).sort_index()
    assert left.index.equals(right.index)
    assert np.allclose(left[fields], right[fields], rtol=tolerance, atol=1e-9, equal_nan=True)


def test_quick_pipeline_reproduces_key_outputs(quick_results: Path):
    national = pd.read_csv(quick_results / "estimates" / "national_models.csv")
    expected_national = pd.read_csv(REFERENCE / "national_models.csv").rename(columns={"label": "specification"})
    _compare(
        national,
        expected_national,
        ["specification"],
        ["ratio", "low", "high", "se", "df"],
    )

    direct = pd.read_csv(quick_results / "estimates" / "direct_contrasts.csv")
    expected_direct = pd.read_csv(REFERENCE / "direct_contrasts.csv")
    _compare(direct, expected_direct, ["contrast"], ["ratio", "low", "high", "se", "df"])

    state = pd.read_csv(quick_results / "estimates" / "state_effects.csv")
    expected_state = pd.read_csv(REFERENCE / "state_stratified.csv")
    _compare(state, expected_state, ["band", "endpoint"], ["ratio", "low", "high"], tolerance=2e-6)

    mississippi = pd.read_csv(quick_results / "estimates" / "mississippi_effects.csv")
    expected_mississippi = pd.read_csv(REFERENCE / "mississippi_stratified.csv")
    expected_mississippi = expected_mississippi[
        expected_mississippi.band.isin(["evening", "morning"])
        & expected_mississippi.endpoint.isin(mississippi.endpoint.unique())
    ]
    _compare(
        mississippi,
        expected_mississippi,
        ["endpoint", "band"],
        ["ratio", "low", "high"],
        tolerance=2e-6,
    )

    traffic = pd.read_csv(quick_results / "estimates" / "traffic_support.csv")
    expected_traffic = pd.read_csv(REFERENCE / "traffic_support.csv").rename(
        columns={
            "min_sites": "minimum_sites",
            "states": "jurisdictions",
            "ratio_geo_mean_state": "ratio_equal_jurisdiction",
            "ci_low": "low",
            "ci_high": "high",
        }
    )
    _compare(
        traffic,
        expected_traffic,
        ["metric", "band", "minimum_sites"],
        ["jurisdictions", "sites", "ratio_equal_jurisdiction", "low", "high"],
    )

    alaska = pd.read_csv(quick_results / "estimates" / "alaska_inclusion.csv")
    expected_alaska = pd.read_csv(REFERENCE / "alaska_inclusion.csv")
    _compare(
        alaska,
        expected_alaska,
        ["geography", "endpoint", "band"],
        ["ratio", "low", "high", "se", "df"],
    )


def test_all_ten_figures_are_bound_to_calculated_data(quick_results: Path):
    manifest_path = quick_results / "figures" / "FIGURE_DATA_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    assert len(manifest) == 10
    assert set(manifest) == {
        "Figure_1_design",
        "Figure_2_aligned_dynamics",
        "Figure_3_primary_contrasts",
        "Figure_S1_estimator_sensitivity",
        "Figure_S2_traffic_support",
        "Figure_S3_sample_flow",
        "Figure_S4_seasonal_counts",
        "Figure_S5_original_event_weeks",
        "Figure_S6_calendar",
        "Figure_S7_policy_year_influence",
    }
    for name, item in manifest.items():
        assert item["rows"] > 0
        for suffix in (".csv", ".pdf", ".png"):
            path = quick_results / "figures" / f"{name}{suffix}"
            minimum = 50 if suffix == ".csv" else 5_000
            assert path.is_file() and path.stat().st_size > minimum


def test_analysis_code_does_not_read_reference_outputs():
    paths = [*PYTHON_ROOT.glob("0*.py"), PYTHON_ROOT / "run.py", *PYTHON_ROOT.joinpath("functions").glob("*.py")]
    for path in paths:
        text = path.read_text()
        assert "reference/" not in text, path
        assert "expected/" not in text, path
