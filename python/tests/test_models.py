"""Numerical tests for the core and targeted estimators."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from functions.allocation import contrast, estimate
from functions.data import national_counts, primary_tables
from functions.exact import jackknife_interval, make_model, scipy_probability_check
from functions.targeted import (
    ENDPOINTS,
    build_triplets,
    fit_triplets,
    interval,
    standardized_estimate,
    standardized_inputs,
)

REFERENCE = Path(__file__).resolve().parents[1] / "reference"


def _assert_close(actual: dict, expected: pd.Series, fields=("ratio", "low", "high"), tolerance=3e-6):
    for field in fields:
        assert actual[field] == pytest.approx(expected[field], rel=tolerance, abs=1e-9), field


@pytest.mark.parametrize("band", ["evening", "morning"])
def test_primary_estimates_match_reference_results(band):
    counts = national_counts()
    actual = estimate(counts, band=band)
    expected = pd.read_csv(REFERENCE / "national_models.csv").set_index("label").loc[f"primary_{band}"]
    _assert_close(actual, expected, fields=("ratio", "low", "high", "se", "df"))


@pytest.mark.parametrize("band", ["evening", "morning"])
def test_direct_road_user_contrasts_match_reference_results(band):
    counts = national_counts()
    pedestrian = estimate(counts, band=band, endpoint="fatal_pedestrian_crashes")
    other = estimate(counts, band=band, endpoint="other_fatal_crashes")
    actual = contrast(pedestrian, other, f"pedestrian_vs_other_{band}")
    expected = pd.read_csv(REFERENCE / "direct_contrasts.csv").set_index("contrast").loc[f"pedestrian_vs_other_{band}"]
    _assert_close(actual, expected, fields=("ratio", "low", "high", "se", "df"))


@pytest.mark.parametrize("band", ["evening", "morning"])
def test_fixed_margin_conditional_likelihood(band):
    tables = primary_tables()
    model, _ = make_model(tables, band)
    full, beta = model.fit()
    assert scipy_probability_check(model, beta) < 1e-9

    deleted = []
    for year in range(2001, 2025):
        reduced, _ = make_model(tables, band, years=[candidate for candidate in range(2001, 2025) if candidate != year])
        deleted.append(reduced.fit()[0]["log_ratio"])
    actual = jackknife_interval(full["log_ratio"], deleted)
    expected = (
        pd.read_csv(REFERENCE / "conditional_blocks.csv").query("band == @band and block_size_years == 1").iloc[0]
    )
    _assert_close(actual, expected, fields=("ratio", "low", "high", "se", "df"))


@pytest.mark.parametrize("band", ["evening", "morning"])
def test_year_specific_slope_conditional_model(band):
    counts = national_counts()
    expected = pd.read_csv(REFERENCE / "annual_slope_conditional.csv")

    fitted = {}
    deleted = {}
    for endpoint in ENDPOINTS:
        triplets = build_triplets(counts, band, endpoint)
        full = fit_triplets(triplets)
        omissions = np.asarray(
            [
                fit_triplets([triplet for triplet in triplets if triplet.year != year])["log_ratio"]
                for year in range(2001, 2025)
            ]
        )
        fitted[endpoint] = full
        deleted[endpoint] = omissions
        actual = {**full, **interval(full["log_ratio"], omissions)}
        row = expected.query("band == @band and endpoint == @endpoint").iloc[0]
        _assert_close(actual, row)

    actual_direct = interval(
        fitted[ENDPOINTS[0]]["log_ratio"] - fitted[ENDPOINTS[1]]["log_ratio"],
        deleted[ENDPOINTS[0]] - deleted[ENDPOINTS[1]],
    )
    expected_direct = expected.query("band == @band and endpoint == 'direct'").iloc[0]
    _assert_close(actual_direct, expected_direct)


@pytest.mark.parametrize("band", ["evening", "morning"])
def test_common_time_standardization(band):
    counts = national_counts()
    frame, years, excluded = standardized_inputs(counts, band)
    assert years == list(range(2001, 2025))
    assert excluded == []

    full, parts = standardized_estimate(frame, years)
    deleted = [standardized_estimate(frame, [year for year in years if year != omitted])[0] for omitted in years]
    expected = pd.read_csv(REFERENCE / "common_time_standardized.csv")

    for endpoint in (*ENDPOINTS, "direct"):
        actual = interval(full[endpoint], [result[endpoint] for result in deleted])
        row = expected.query("band == @band and endpoint == @endpoint").iloc[0]
        _assert_close(actual, row)

    for endpoint in ENDPOINTS:
        reconstructed = sum(part["contribution"] for part in parts if part["endpoint"] == endpoint)
        assert reconstructed == pytest.approx(full[endpoint], abs=1e-12)


def test_five_center_quadratic_models_match_references():
    from functions.curvature import five_center_analysis

    actual = five_center_analysis(full=False)["five_center_curvature"]
    expected = pd.read_csv(REFERENCE / "five_center_curvature.csv")
    expected = expected[expected.model.eq("five_center_quadratic")]
    for band in ("evening", "morning"):
        for endpoint in (*ENDPOINTS, "direct"):
            observed = actual.query(
                "model == 'five_center_quadratic' and band == @band and endpoint == @endpoint"
            ).iloc[0]
            reference = expected.query("band == @band and endpoint == @endpoint").iloc[0]
            _assert_close(observed.to_dict(), reference, fields=("ratio", "low", "high", "se", "df"))


def test_exact_wild_cluster_multiplier_matches_references():
    from functions.curvature import wild_cluster_calibration

    actual = wild_cluster_calibration().set_index("comparison")
    expected = pd.read_csv(REFERENCE / "wild_cluster_calibration.csv").set_index("comparison")
    fields = (
        "ratio",
        "wild_p_value",
        "wild_critical_95",
        "wild_low",
        "wild_high",
        "clusters",
        "unique_sign_patterns",
    )
    for comparison in expected.index:
        _assert_close(actual.loc[comparison].to_dict(), expected.loc[comparison], fields=fields, tolerance=1e-6)


def test_season_and_policy_era_diagnostics_match_references():
    from functions.review_diagnostics import policy_era_interaction, seasonal_heterogeneity

    season = seasonal_heterogeneity()
    season_expected = pd.read_csv(REFERENCE / "season_heterogeneity.csv")
    for _, reference in season_expected.iterrows():
        observed = season.query(
            "band == @reference.band and result_type == @reference.result_type and label == @reference.label"
        ).iloc[0]
        _assert_close(observed.to_dict(), reference, fields=("ratio", "low", "high", "se", "p_value", "df"))

    era = policy_era_interaction()
    era_expected = pd.read_csv(REFERENCE / "policy_era_interaction.csv")
    for _, reference in era_expected.iterrows():
        observed = era.query(
            "band == @reference.band and result_type == @reference.result_type and label == @reference.label"
        ).iloc[0]
        _assert_close(observed.to_dict(), reference, fields=("ratio", "low", "high", "se", "p_value", "df"))
