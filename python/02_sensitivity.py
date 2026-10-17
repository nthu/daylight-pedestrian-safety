"""Run the principal calendar, estimator and target-definition sensitivities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from functions.allocation import contrast, estimate
from functions.curvature import five_center_analysis, wild_cluster_calibration
from functions.data import (
    ensure_results,
    national_counts,
    primary_tables,
    reference_counts,
    sensitivity_table,
    write_csv,
)
from functions.exact import jackknife_interval, make_model
from functions.review_diagnostics import review_diagnostics
from functions.specifications import BANDS, CALENDAR_SPECS, DYNAMIC_WEEKS, OTHER_FATAL, PEDESTRIAN, SEASONS
from functions.targeted import (
    ENDPOINTS,
    build_triplets,
    fit_triplets,
    interval,
    standardized_estimate,
    standardized_inputs,
)
from scipy import stats


def _calendar_sensitivity() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for specification, filename in CALENDAR_SPECS.items():
        frame = sensitivity_table(filename)
        geography = frame.geo.iloc[0]
        for band in BANDS:
            fitted = estimate(frame, band=band, geo=geography)
            rows.append(
                {
                    "specification": specification,
                    "band": band,
                    **{key: fitted[key] for key in ("ratio", "low", "high", "se", "df")},
                }
            )
    return pd.DataFrame(rows)


def _dynamic_sensitivity() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for week in DYNAMIC_WEEKS:
        frame = sensitivity_table(f"dynamic_counts_week_{week:+d}.csv")
        geography = frame.geo.iloc[0]
        for band in BANDS:
            for season in SEASONS:
                fitted = estimate(frame, band=band, seasons=[season], geo=geography)
                rows.append(
                    {
                        "week": week,
                        "band": band,
                        "season": season,
                        **{key: fitted[key] for key in ("ratio", "low", "high", "se", "df")},
                    }
                )
    return pd.DataFrame(rows)


def _reference_sensitivity() -> pd.DataFrame:
    frame = reference_counts()
    specifications = [
        ("evening_ref_early_evening", {"band": "evening", "reference": "early_evening"}),
        ("evening_ref_late_evening", {"band": "evening", "reference": "late_evening"}),
        ("evening_DST_holidaypaired", {"band": "evening", "geo": "DST_holidaypaired"}),
        ("evening_DST_weekday", {"band": "evening", "geo": "DST_weekday"}),
        ("evening_DST_weekend", {"band": "evening", "geo": "DST_weekend"}),
        ("morning_ref_early_morning", {"band": "morning", "reference": "early_morning"}),
        ("morning_ref_late_morning", {"band": "morning", "reference": "late_morning"}),
        ("morning_DST_holidaypaired", {"band": "morning", "geo": "DST_holidaypaired"}),
        ("morning_DST_weekday", {"band": "morning", "geo": "DST_weekday"}),
        ("control_night_vs_noon", {"band": "night", "reference": "noon"}),
        (
            "control_late_vs_early_evening",
            {"band": "late_evening", "reference": "early_evening"},
        ),
    ]
    rows = []
    for name, kwargs in specifications:
        fitted = estimate(frame, **kwargs)
        rows.append(
            {
                "specification": name,
                **{key: fitted[key] for key in ("ratio", "low", "high", "p_value", "se", "df")},
            }
        )
    return pd.DataFrame(rows)


def _conditional_likelihood(full: bool = True) -> pd.DataFrame:
    tables = primary_tables()
    rows: list[dict[str, object]] = []
    for band in BANDS:
        model, _ = make_model(tables, band)
        fitted, _ = model.fit()
        for block_size in (1, 2, 3, 4):
            if not full and block_size not in (1, 4):
                continue
            deleted: list[float] = []
            for start in range(2001, 2025, block_size):
                omitted = set(range(start, start + block_size))
                years = [year for year in range(2001, 2025) if year not in omitted]
                reduced, _ = make_model(tables, band, years=years)
                deleted.append(reduced.fit()[0]["log_ratio"])
            result = jackknife_interval(fitted["log_ratio"], deleted)
            rows.append({"band": band, "block_size_years": block_size, **result})
    return pd.DataFrame(rows)


def _whole_day(counts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    endpoints = (PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")
    for endpoint in endpoints:
        subset = counts[
            counts.geo.eq("DST_main")
            & counts.width_minutes.eq(40)
            & counts.band.eq("all_day")
            & counts.endpoint.eq(endpoint)
        ]
        tables = subset.pivot(index=["year", "season"], columns=["sham", "period"], values="events")
        for season in SEASONS:
            annual = tables.xs(season, level="season")
            if not (annual > 0).all().all():
                raise ValueError(f"Whole-day counts contain zero cells for {endpoint}, {season}")
            log_changes = {center: np.log(annual[center, "post"] / annual[center, "pre"]) for center in (-35, 0, 35)}
            methods = {
                "raw_actual": log_changes[0],
                "symmetric_sham_adjusted": log_changes[0] - 0.5 * (log_changes[-35] + log_changes[35]),
            }
            for method, values in methods.items():
                mean = float(values.mean())
                standard_error = float(values.std(ddof=1) / np.sqrt(len(values)))
                critical = stats.t.ppf(0.975, len(values) - 1)
                rows.append(
                    {
                        "endpoint": endpoint,
                        "season": season,
                        "method": method,
                        "ratio": float(np.exp(mean)),
                        "low": float(np.exp(mean - critical * standard_error)),
                        "high": float(np.exp(mean + critical * standard_error)),
                    }
                )
    return pd.DataFrame(rows)


def _targeted_models(counts: pd.DataFrame) -> dict[str, pd.DataFrame]:
    years = list(range(2001, 2025))

    pre_change: list[dict[str, object]] = []
    for band in BANDS:
        fitted: dict[str, dict[str, object]] = {}
        for endpoint in ENDPOINTS:
            result = estimate(counts, band=band, endpoint=endpoint, years=range(2001, 2021))
            fitted[endpoint] = result
            pre_change.append(
                {
                    "band": band,
                    "endpoint": endpoint,
                    **{key: result[key] for key in ("ratio", "low", "high", "p_value", "se", "df")},
                }
            )
        direct = contrast(fitted[PEDESTRIAN], fitted[OTHER_FATAL], f"pre2022_{band}")
        pre_change.append({"band": band, "endpoint": "direct", **direct})

    conditional: list[dict[str, object]] = []
    support: list[dict[str, object]] = []
    stored: dict[tuple[str, str], tuple[float, np.ndarray]] = {}
    for band in BANDS:
        for endpoint in ENDPOINTS:
            triplets = build_triplets(counts, band, endpoint)
            fitted = fit_triplets(triplets)
            deleted = np.asarray(
                [fit_triplets([triplet for triplet in triplets if triplet.year != year])["log_ratio"] for year in years]
            )
            conditional.append(
                {
                    "band": band,
                    "endpoint": endpoint,
                    **fitted,
                    **interval(fitted["log_ratio"], deleted),
                }
            )
            stored[band, endpoint] = fitted["log_ratio"], deleted
            for triplet in triplets:
                support.append(
                    {
                        "band": band,
                        "endpoint": endpoint,
                        "year": triplet.year,
                        "season": triplet.season,
                        "support_size": len(triplet.support),
                        "informative": len(triplet.support) > 1,
                    }
                )
        ped, ped_deleted = stored[band, PEDESTRIAN]
        other, other_deleted = stored[band, OTHER_FATAL]
        conditional.append(
            {
                "band": band,
                "endpoint": "direct",
                **interval(ped - other, ped_deleted - other_deleted),
            }
        )

    standardized: list[dict[str, object]] = []
    decomposition: list[dict[str, object]] = []
    cell_support: list[pd.DataFrame] = []
    for band in BANDS:
        frame, supported_years, excluded_years = standardized_inputs(counts, band)
        full, parts = standardized_estimate(frame, supported_years)
        deleted = [
            standardized_estimate(frame, [year for year in supported_years if year != omitted])[0]
            for omitted in supported_years
        ]
        for endpoint in (*ENDPOINTS, "direct"):
            values = [result[endpoint] for result in deleted]
            standardized.append(
                {
                    "band": band,
                    "endpoint": endpoint,
                    **interval(full[endpoint], values),
                    "year_seasons": len(supported_years) * 2,
                    "excluded_years": ";".join(map(str, excluded_years)),
                }
            )
        decomposition.extend({"band": band, **part} for part in parts)
        export = frame[["year", "season", "endpoint", "sham", "period", "total", "fraction", "gain_side"]].copy()
        export["band"] = band
        export["weight"] = 1 / (2 * len(supported_years))
        cell_support.append(export)

    return {
        "pre2022": pd.DataFrame(pre_change),
        "annual_slope_conditional": pd.DataFrame(conditional),
        "triplet_support": pd.DataFrame(support),
        "common_time": pd.DataFrame(standardized),
        "sham_decomposition": pd.DataFrame(decomposition),
        "standardization_support": pd.concat(cell_support, ignore_index=True),
    }


def _clock_shift() -> pd.DataFrame:
    frame = sensitivity_table("clock_shift.csv.gz")
    rows: list[dict[str, object]] = []
    for minutes in sorted(frame.clock_shift_minutes.unique()):
        shifted = frame[frame.clock_shift_minutes.eq(minutes)]
        for band in BANDS:
            fitted = {endpoint: estimate(shifted, band=band, endpoint=endpoint) for endpoint in ENDPOINTS}
            for endpoint, result in fitted.items():
                rows.append(
                    {
                        "clock_shift_minutes": int(minutes),
                        "band": band,
                        "endpoint": endpoint,
                        **{key: result[key] for key in ("ratio", "low", "high", "se", "df")},
                    }
                )
            direct = contrast(fitted[PEDESTRIAN], fitted[OTHER_FATAL], "clock_shift")
            rows.append({"clock_shift_minutes": int(minutes), "band": band, "endpoint": "direct", **direct})
    return pd.DataFrame(rows)


def run(full: bool = True) -> dict[str, pd.DataFrame]:
    counts = national_counts()
    paths = ensure_results()

    curvature = five_center_analysis(full=full)
    diagnostics = review_diagnostics(counts)
    outputs = {
        "reference_sensitivity": _reference_sensitivity(),
        "calendar_sensitivity": _calendar_sensitivity(),
        "dynamic_estimates": _dynamic_sensitivity(),
        "conditional_blocks": _conditional_likelihood(full=full),
        "whole_day": _whole_day(counts),
        "clock_shift": _clock_shift(),
        "wild_cluster_calibration": wild_cluster_calibration(counts),
        **curvature,
        **diagnostics,
        **_targeted_models(counts),
    }

    for name, frame in outputs.items():
        write_csv(frame, paths["estimates"] / f"{name}.csv")

    print("\nKey sensitivity estimates")
    print(outputs["calendar_sensitivity"].to_string(index=False))
    print("\nCommon-time direct comparisons")
    print(outputs["common_time"].query("endpoint == 'direct'")[["band", "ratio", "low", "high"]].to_string(index=False))
    return outputs


if __name__ == "__main__":
    run(full=True)
